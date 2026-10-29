################################################################################
# cancel.jl
#
# Cancellation tokens and interruptible waits.
#
# A `CancelToken` stands for "the work this task is doing has been cancelled".
# It is created lazily per task, may be replaced or cleared with
# `setcancaltoken!`, and may be cancelled from any task with `cancel!`. Every
# blocking wait in Air that accepts a token registers a small `Waiter` with the
# token, so a cancellation wakes it; `take` registers the same waiter with the
# promise it is waiting on.
#
# Lock order (the soundness argument for the whole scheme): a task holds its
# waiter's lock only while it is parked or checking whether it was woken, and it
# never acquires a promise mutex, a stream mutex or a token lock while holding
# it. Promise mutexes, stream mutexes and token locks may therefore wake waiters
# while they hold their own state lock, and no cycle is possible.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

"""
    Cancelled(token, reason)

An exception thrown by a blocking Air primitive when the `CancelToken` it was
given has been cancelled. `token` is the token and `reason` is the reason given
to [`cancel!`](@ref), which may be `nothing`.

See also: [`CancelToken`](@ref), [`cancel!`](@ref), [`checkcancel`](@ref)
"""
struct Cancelled <: Exception
    token::Any
    reason::Any
end
export Cancelled
function Base.showerror(io::IO, e::Cancelled)
    print(io, "cancelled")
    reason = getfield(e, :reason)
    if reason === nothing
        print(io, " (no reason given)")
    else
        print(io, ": ")
        showerror(io, reason isa Exception ? reason : ErrorException(string(reason)))
    end
    return nothing
end

"""
    Waiter

A one-task wake-up channel used by Air's blocking primitives. A waiter has its
own lock and condition, and two flags: `done` (it was woken) and `timedout` (its
timer expired). Waiting tasks park on the condition, and every path that should
wake them — a promise settlement, a cancellation, or a timeout — calls `_wake!`
or `_timeout!`, which take the waiter's lock and notify.

A waiter is created by the blocking call, so it has exactly one waiter task.

This struct is part of Air's private/internal implementation details.
"""
mutable struct Waiter
    lock::ReentrantLock
    cond::Threads.Condition
    done::Bool
    timedout::Bool
end
Waiter() = begin
    l = ReentrantLock()
    Waiter(l, Threads.Condition(l), false, false)
end
function _wake!(w::Waiter)
    lock(getfield(w, :lock))
    try
        setfield!(w, :done, true)
        notify(getfield(w, :cond); all=true)
    finally
        unlock(getfield(w, :lock))
    end
    return nothing
end
function _timeout!(w::Waiter)
    lock(getfield(w, :lock))
    try
        setfield!(w, :timedout, true)
        setfield!(w, :done, true)
        notify(getfield(w, :cond); all=true)
    finally
        unlock(getfield(w, :lock))
    end
    return nothing
end
# Park on a waiter until it is woken, and yield whether its timer expired. A
# wake that arrived before the park is consumed immediately rather than missed,
# and the waiter's lock is released on return, so the caller never holds it
# while checking state or taking another lock (see the lock-order note at the
# top of this file). This is the slow path of every blocking primitive.
function _park!(w::Waiter)
    timedout = false
    lock(getfield(w, :lock))
    try
        if getfield(w, :timedout)
            timedout = true
        else
            getfield(w, :done) || wait(getfield(w, :cond))
            timedout = getfield(w, :timedout)
        end
        setfield!(w, :done, false)
    finally
        unlock(getfield(w, :lock))
    end
    return timedout
end

"""
    CancelToken()

A cancellation token: a flag plus the waiters parked under it. Blocking Air
primitives accept one as a `token` keyword argument and throw
[`Cancelled`](@ref) when it has been cancelled. A token is per task by default;
[`canceltoken`](@ref) yields the calling task's token, [`setcancaltoken!`](@ref)
replaces it, and [`cancel!`](@ref) cancels it from any task.

Cancellation is one-way: a token that has been cancelled stays cancelled.

See also: [`canceltoken`](@ref), [`cancel!`](@ref), [`iscancelled`](@ref),
[`checkcancel`](@ref), [`Cancelled`](@ref)
"""
mutable struct CancelToken
    lock::ReentrantLock
    cancelled::Bool
    reason::Any
    waiters::Base.IdSet{Waiter}
end
CancelToken() = CancelToken(ReentrantLock(), false, nothing, Base.IdSet{Waiter}())
export CancelToken

# The task-local cancellation-token slot. The stored value is a `CancelToken`
# or `nothing`; `canceltoken` creates a token when the slot holds neither.
const _CANCELTOKEN_KEY = :AirCancelToken

"""
    canceltoken()

Yields the calling task's `CancelToken`, creating it if the task does not have
one yet. This is the default token of every blocking Air call, so a task needs
to do nothing to be cancellable; to cancel it, another task needs the token
itself (typically recorded when the task is created).

See also: [`setcancaltoken!`](@ref), [`cancel!`](@ref)
"""
function canceltoken()
    d = task_local_storage()
    tok = get(d, _CANCELTOKEN_KEY, nothing)
    (tok isa CancelToken) && return tok
    tok = CancelToken()
    task_local_storage(_CANCELTOKEN_KEY, tok)
    return tok
end
export canceltoken

"""
    setcancaltoken!(token)

Sets the calling task's cancellation token to `token`, which may be `nothing`
to give the task a fresh token on its next `canceltoken()` call, and yields it.

See also: [`canceltoken`](@ref), [`CancelToken`](@ref)
"""
function setcancaltoken!(token::Union{Nothing,CancelToken})
    task_local_storage(_CANCELTOKEN_KEY, token)
    return token
end
export setcancaltoken!

# Park `w` under `tok`. The token's lock is the only lock taken; see the lock
# order note at the top of this file.
function _register_waiter!(tok::CancelToken, w::Waiter)
    lock(getfield(tok, :lock))
    try
        push!(getfield(tok, :waiters), w)
    finally
        unlock(getfield(tok, :lock))
    end
    return nothing
end
function _unregister_waiter!(tok::CancelToken, w::Waiter)
    lock(getfield(tok, :lock))
    try
        delete!(getfield(tok, :waiters), w)
    finally
        unlock(getfield(tok, :lock))
    end
    return nothing
end

"""
    cancel!(token, reason = nothing)

Cancels the given `CancelToken` and wakes every task waiting on it; the reason
is available as the `reason` field of the [`Cancelled`](@ref) exception those
tasks throw. Cancelling a token that is already cancelled changes nothing and
is not an error; the first reason given is kept. Yields the token.

See also: [`CancelToken`](@ref), [`iscancelled`](@ref), [`checkcancel`](@ref)
"""
function cancel!(tok::CancelToken, reason = nothing)
    # Collect the waiters under the lock, but wake them without it: waking takes
    # each waiter's lock, and although the lock order permits holding this one
    # while doing so, releasing first keeps the critical section short.
    waiters = Waiter[]
    lock(getfield(tok, :lock))
    try
        getfield(tok, :cancelled) && return tok
        setfield!(tok, :cancelled, true)
        setfield!(tok, :reason, reason)
        for w in getfield(tok, :waiters)
            push!(waiters, w)
        end
        empty!(getfield(tok, :waiters))
    finally
        unlock(getfield(tok, :lock))
    end
    for w in waiters
        _wake!(w)
    end
    return tok
end
export cancel!

"""
    iscancelled(token)

Yields `true` if the given `CancelToken` has been cancelled. This is the
non-throwing form of [`checkcancel`](@ref).

See also: [`CancelToken`](@ref), [`cancel!`](@ref), [`checkcancel`](@ref)
"""
function iscancelled(tok::CancelToken)
    lock(getfield(tok, :lock))
    try
        return getfield(tok, :cancelled)
    finally
        unlock(getfield(tok, :lock))
    end
end
export iscancelled

"""
    checkcancel(token)

Throws [`Cancelled`](@ref) if the given `CancelToken` has been cancelled and
yields `nothing` otherwise. A `nothing` token is always allowed and does
nothing, so `checkcancel(token)` may be called unconditionally when a token is
optional. This is the explicit, non-blocking way to observe cancellation; every
blocking Air call checks its token itself.

See also: [`CancelToken`](@ref), [`iscancelled`](@ref), [`cancel!`](@ref)
"""
checkcancel(::Nothing) = nothing
function checkcancel(tok::CancelToken)
    lock(getfield(tok, :lock))
    try
        getfield(tok, :cancelled) || return nothing
        throw(Cancelled(tok, getfield(tok, :reason)))
    finally
        unlock(getfield(tok, :lock))
    end
end
export checkcancel

"""
    sleep_interruptible(seconds; token = canceltoken())

Sleeps for the given number of seconds, but throws [`Cancelled`](@ref) as soon
as the given token is cancelled. Unlike `Base.sleep`, this cannot leave a task
uninterruptible. Cancellation wins if it races the end of the sleep.

This is part of Air's internal interface: it is not exported, because exporting
`sleep` would clash with `Base.sleep` for anyone who does `using Air`. AirLang
calls it qualified as `Air.sleep_interruptible`.

See also: [`CancelToken`](@ref), [`canceltoken`](@ref)
"""
function sleep_interruptible(seconds::Real; token::Union{Nothing,CancelToken} = canceltoken())
    (token === nothing) || checkcancel(token)
    w = Waiter()
    (token === nothing) || _register_waiter!(token, w)
    timer = Timer(_ -> _timeout!(w), max(seconds, 0))
    try
        lock(getfield(w, :lock))
        try
            # Re-check after registering: a cancellation between the check above
            # and the registration wakes this waiter.
            getfield(w, :done) || wait(getfield(w, :cond))
        finally
            unlock(getfield(w, :lock))
        end
        # Cancellation wins over the timer, including one that arrives between
        # the wake and here.
        (token === nothing) || checkcancel(token)
        return nothing
    finally
        (token === nothing) || _unregister_waiter!(token, w)
        close(timer)
    end
end
