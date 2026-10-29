################################################################################
# stream.jl
#
# Streams: reusable, optionally bounded FIFOs that block their readers and
# writers and that take part in transactions.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

# A taken item, wrapped so that an item that is itself `nothing` is distinct
# from "the stream was empty"; cf. `_Settled` in util.jl.
struct _StreamItem{T}
    value::T
end

"""
    StreamClosed()

An exception thrown when a closed [`Stream`](@ref) is written to, or read from
once it is empty. A closed stream can still be drained by `take` until no items
remain.

See also: [`Stream`](@ref), `close`
"""
struct StreamClosed <: Exception end
export StreamClosed
Base.showerror(io::IO, ::StreamClosed) = print(io, "stream is closed")

"""
    StreamDeadlock()

An exception thrown when an operation would block on a [`Stream`](@ref) that
was created in the current transaction. No other task can have a reference to
such a stream until the transaction commits, so waiting for one is a guaranteed
deadlock; the stream itself is left unchanged.

See also: [`Stream`](@ref), [`tx`](@ref)
"""
struct StreamDeadlock <: Exception end
export StreamDeadlock
function Base.showerror(io::IO, ::StreamDeadlock)
    print(io, "cannot block on a stream created in this transaction: ")
    return print(io, "no other task can reach it until the transaction commits")
end

"""
    Stream{T}(; maxsize = nothing)
    Stream(; maxsize = nothing)

A `Stream{T}` is a FIFO of values of type `T` that may be written and read
repeatedly, unlike the one-shot [`Promise`](@ref). `put!` appends a value and
[`take`](@ref) removes the oldest one. Both block when they cannot proceed — a
read when the stream is empty, a write when it is full — and both accept the
`timeout` and `token` keywords of `take(::Promise)`.

`maxsize` bounds the stream: with `maxsize = 3`, three values may be waiting,
and a fourth `put!` waits for a `take`. The default, `nothing`, leaves the
stream unbounded, so `put!` never waits. `Stream()` is a `Stream{Any}`.

Outside a transaction, `put!` and `take` happen immediately, under the stream's
lock. Inside one they are staged, and happen only if the transaction commits;
the transaction that staged them sees its own writes and reads, and a concurrent
change to the stream invalidates the transaction, which retries. Blocking
inside a transaction is allowed, but waiting for a change that only the
transaction's own body could make is a deadlock; a stream created in the current
transaction therefore raises [`StreamDeadlock`](@ref) instead of blocking. See
the [Streams](streams.md) page for the full rules.

`close` ends the stream: a later `put!` throws
[`StreamClosed`](@ref), and a `take` drains the values that remain and then
throws it as well. Streams compare by identity, and `isopen`, `isempty`,
`length` and `isready` inspect one without blocking.

See also: [`take`](@ref), [`StreamClosed`](@ref), [`StreamDeadlock`](@ref),
[`maxsize`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> s = Stream{Int}();

julia> put!(s, 1); put!(s, 2);

julia> take(s)
1

julia> take(s)
2

julia> isempty(s)
true
```

A bounded stream makes `put!` wait, and `timeout` bounds that wait:

```jldoctest
julia> s = Stream{Int}(maxsize = 1);

julia> put!(s, 1);

julia> try
           put!(s, 2; timeout = 0.01)
       catch e
           sprint(showerror, e)
       end
"reader timed out after 0.01 seconds"

julia> take(s)
1
```

`close` lets a reader drain what is left and then stop:

```jldoctest
julia> s = Stream{Int}();

julia> put!(s, 7); close(s);

julia> take(s)
7

julia> try
           take(s)
       catch e
           e
       end
StreamClosed()
```
"""
mutable struct Stream{T}
    mutex::ReentrantLock
    # The live items are `items[head:end]`; the dead prefix is compacted away
    # when it is at least as large as the live suffix.
    items::Vector{T}
    head::Int
    # `nothing` for an unbounded stream, otherwise the maximum logical length.
    maxsize::Union{Nothing,Int}
    # Bumped once per mutation: on every change outside a transaction, and once
    # per committed transaction that changes the stream. A transaction validates
    # it at commit and uses it to detect changes while it is running.
    version::UInt64
    closed::Bool
    # Created on the first blocking operation and dropped when the last one
    # leaves. Guarded by `mutex`.
    waiters::Union{Nothing,Base.IdSet{Waiter}}
    function Stream{T}(maxsize::Union{Nothing,Integer}) where {T}
        (maxsize === nothing || maxsize >= 1) ||
            throw(ArgumentError("maxsize must be nothing or at least 1"))
        s = new{T}(ReentrantLock(), T[], 1,
                   maxsize === nothing ? nothing : Int(maxsize),
                   UInt64(0), false, nothing)
        _stream_register_born!(s)
        return s
    end
end
export Stream
Stream{T}(; maxsize::Union{Nothing,Integer} = nothing) where {T} = Stream{T}(maxsize)
Stream(maxsize::Union{Nothing,Integer}) = Stream{Any}(maxsize)
Stream(; maxsize::Union{Nothing,Integer} = nothing) = Stream{Any}(maxsize)

"""
    maxsize(stream)

Yields the maximum number of values a [`Stream`](@ref) may hold, or `nothing`
if it is unbounded.

See also: [`Stream`](@ref)
"""
maxsize(s::Stream) = getfield(s, :maxsize)
export maxsize

# The number of live items; the caller holds the stream's mutex.
@inline _live_count(s::Stream) = length(getfield(s, :items)) - getfield(s, :head) + 1
# The number of staged puts still to be read from a transaction's staging data.
@inline function _staged_count(d)
    return length(getfield(d, :puts)) - getfield(d, :puts_head) + 1
end

# Drop the consumed prefix of the live items. The caller holds the mutex, and
# the compaction never changes the logical contents, so it never bumps the
# version: a transaction's `ntaken` counts items, not indices.
function _stream_compact!(s::Stream)
    h = getfield(s, :head)
    (h == 1) && return nothing
    items = getfield(s, :items)
    live = length(items) - h + 1
    if live <= 0
        empty!(items)
        setfield!(s, :head, 1)
    elseif h - 1 >= live
        deleteat!(items, 1:(h - 1))
        setfield!(s, :head, 1)
    end
    return nothing
end

# Record that this stream was created in the current transaction attempt. The
# marker lives in the transaction (and is dropped by `tx_clear!`), so the
# stream does not keep a reference to the transaction, and an attempt that
# retries does not see a stale marker.
function _stream_register_born!(s::Stream{T}) where {T}
    t = currtx()
    (t === nothing) && return nothing
    getfield(t, :streams)[s] =
        StreamTxData{T}(getfield(s, :version), 0, T[], 1, false, true, false)
    return nothing
end

# Register `w` so that any mutation wakes it. The caller holds the mutex; unlike
# a promise's waiters, a stream's registration survives a wake, because a woken
# waiter may lose the race for the item (or the space) and have to wait again.
function _stream_register_locked!(s::Stream, w::Waiter)
    ws = getfield(s, :waiters)
    if ws === nothing
        ws = Base.IdSet{Waiter}()
        setfield!(s, :waiters, ws)
    end
    push!(ws, w)
    return nothing
end
function _stream_unregister!(s::Stream, w::Waiter)
    mux = getfield(s, :mutex)
    lock(mux)
    try
        ws = getfield(s, :waiters)
        if ws !== nothing
            delete!(ws, w)
            isempty(ws) && setfield!(s, :waiters, nothing)
        end
    finally
        unlock(mux)
    end
    return nothing
end
# Wake every registered waiter and keep the registrations. The caller holds the
# mutex. See `_wake_waiters!(::Promise)` in util.jl for the one-shot variant,
# which drops the registry instead.
function _wake_waiters!(s::Stream)
    ws = getfield(s, :waiters)
    (ws === nothing) && return nothing
    for w in ws
        _wake!(w)
    end
    return nothing
end

# Append `x` if the stream is open and has room; yield whether it did. The
# caller holds the mutex, so the wake-up sees the change.
function _stream_push!(s::Stream{T}, x::T) where {T}
    getfield(s, :closed) && throw(StreamClosed())
    m = getfield(s, :maxsize)
    if m !== nothing
        _live_count(s) < m || return false
    end
    push!(getfield(s, :items), x)
    setfield!(s, :version, getfield(s, :version) + UInt64(1))
    _wake_waiters!(s)
    return true
end
# Remove and yield the oldest item, `nothing` if the stream is empty, and throw
# if it is empty and closed. The caller holds the mutex.
function _stream_pop!(s::Stream{T}) where {T}
    if _live_count(s) > 0
        head = getfield(s, :head)
        x = getfield(s, :items)[head]
        setfield!(s, :head, head + 1)
        _stream_compact!(s)
        setfield!(s, :version, getfield(s, :version) + UInt64(1))
        _wake_waiters!(s)
        return _StreamItem{T}(x)
    elseif getfield(s, :closed)
        throw(StreamClosed())
    else
        return nothing
    end
end

"""
    take(stream; timeout = nothing, token = canceltoken())

Removes and yields the oldest value of the given [`Stream`](@ref), suspending
the current task until a value is available if necessary. A closed stream that
still holds values yields them, and then throws [`StreamClosed`](@ref).

`timeout` bounds the wait in seconds: if no value is available when it expires,
`take` throws [`ReaderTimeout`](@ref). `token` is the [`CancelToken`](@ref) to
observe, as for `take(::Promise)`; cancellation wins even over an available
value. Inside a transaction, `take` is staged and the value is removed only if
the transaction commits; the transaction also sees its own staged writes. A
`take` that has to wait inside a transaction waits normally, but a stream
created in that transaction raises [`StreamDeadlock`](@ref) instead.

See also: [`Stream`](@ref), `put!`, `close`
"""
function take(s::Stream{T};
              timeout::Union{Nothing,Real} = nothing,
              token::Union{Nothing,CancelToken} = canceltoken()) where {T}
    # Cancellation wins over availability: a task whose work was cancelled
    # should not go on to use a value it was waiting for. `token = nothing`
    # opts out.
    (token === nothing) || checkcancel(token)
    t = currtx()
    (t === nothing) || return _stream_take_tx(s, t, timeout, token)
    mux = getfield(s, :mutex)
    # The common case, and the only one that allocates nothing: a value waiting
    # at the front of the stream.
    lock(mux)
    r = try
        _stream_pop!(s)
    finally
        unlock(mux)
    end
    (r === nothing) || return r.value
    (timeout === nothing) || (timeout >= 0) ||
        throw(ArgumentError("timeout must be non-negative"))
    w = Waiter()
    timer = nothing
    try
        # Register before the deciding check, so a mutation that happens after
        # it cannot be missed.
        lock(mux)
        try
            _stream_register_locked!(s, w)
        finally
            unlock(mux)
        end
        (token === nothing) || _register_waiter!(token, w)
        (timeout === nothing) || (timer = Timer(_ -> _timeout!(w), timeout))
        while true
            (token === nothing) || checkcancel(token)
            lock(mux)
            r = try
                _stream_pop!(s)
            finally
                unlock(mux)
            end
            (r === nothing) || return r.value
            if _park!(w)
                # The timer may have raced a value or a close; one last check
                # lets them win.
                lock(mux)
                r = try
                    _stream_pop!(s)
                finally
                    unlock(mux)
                end
                (r === nothing) || return r.value
                throw(ReaderTimeout(timeout))
            end
        end
    finally
        _stream_unregister!(s, w)
        (token === nothing) || _unregister_waiter!(token, w)
        (timer === nothing) || close(timer)
    end
end

"""
    put!(stream, value; timeout = nothing, token = canceltoken())

Appends `value` to the given [`Stream`](@ref) and yields it, converted to the
stream's element type. If the stream has a `maxsize` and is full, `put!`
suspends the current task until a `take` makes room. A closed stream throws
[`StreamClosed`](@ref).

`timeout` bounds the wait in seconds and throws [`ReaderTimeout`](@ref) when it
expires; `token` is the [`CancelToken`](@ref) to observe, as for `take`. Inside
a transaction, the append is staged and happens only if the transaction
commits, and a `put!` after a staged `close` throws [`StreamClosed`](@ref). A
`put!` that has to wait on a stream created in that transaction raises
[`StreamDeadlock`](@ref) instead.

See also: [`Stream`](@ref), [`take`](@ref), `close`
"""
function Base.put!(s::Stream{T}, x;
                   timeout::Union{Nothing,Real} = nothing,
                   token::Union{Nothing,CancelToken} = canceltoken()) where {T}
    (token === nothing) || checkcancel(token)
    xT = convert(T, x)
    t = currtx()
    (t === nothing) || return _stream_put_tx(s, xT, t, timeout, token)
    mux = getfield(s, :mutex)
    lock(mux)
    ok = try
        _stream_push!(s, xT)
    finally
        unlock(mux)
    end
    ok && return xT
    (timeout === nothing) || (timeout >= 0) ||
        throw(ArgumentError("timeout must be non-negative"))
    w = Waiter()
    timer = nothing
    try
        lock(mux)
        try
            _stream_register_locked!(s, w)
        finally
            unlock(mux)
        end
        (token === nothing) || _register_waiter!(token, w)
        (timeout === nothing) || (timer = Timer(_ -> _timeout!(w), timeout))
        while true
            (token === nothing) || checkcancel(token)
            lock(mux)
            ok = try
                _stream_push!(s, xT)
            finally
                unlock(mux)
            end
            ok && return xT
            if _park!(w)
                # The timer may have raced a take or a close; one last check
                # lets them win.
                lock(mux)
                ok = try
                    _stream_push!(s, xT)
                finally
                    unlock(mux)
                end
                ok && return xT
                throw(ReaderTimeout(timeout))
            end
        end
    finally
        _stream_unregister!(s, w)
        (token === nothing) || _unregister_waiter!(token, w)
        (timer === nothing) || close(timer)
    end
end

"""
    close(stream)

Closes the given [`Stream`](@ref), waking every task waiting on it. Closing is
idempotent. Afterwards a `put!` throws [`StreamClosed`](@ref), and a `take`
yields the values that remain and then throws it too. Inside a transaction the
close is staged and happens only if the transaction commits; `isopen` sees the
staged close in the transaction that staged it.

See also: [`Stream`](@ref), [`StreamClosed`](@ref), `isopen`
"""
function Base.close(s::Stream{T}) where {T}
    t = currtx()
    mux = getfield(s, :mutex)
    lock(mux)
    try
        if t === nothing
            getfield(s, :closed) && return nothing
            setfield!(s, :closed, true)
            setfield!(s, :version, getfield(s, :version) + UInt64(1))
            _wake_waiters!(s)
            return nothing
        end
        d = get(getfield(t, :streams), s, nothing)
        (d !== nothing && getfield(d, :closed)) && return nothing
        getfield(s, :closed) && return nothing
        (d === nothing) && (d = _stream_txdata_locked!(s, t))
        _stream_tx_recheck!(s, d)
        setfield!(d, :closed, true)
        setfield!(d, :touched, true)
        return nothing
    finally
        unlock(mux)
    end
end

"""
    isopen(stream)

Yields `false` if the given [`Stream`](@ref) has been closed and `true`
otherwise. Inside a transaction, a close staged by that transaction is visible,
so this is `false` there as well.

See also: [`Stream`](@ref), `close`
"""
function Base.isopen(s::Stream)
    mux = getfield(s, :mutex)
    lock(mux)
    try
        getfield(s, :closed) && return false
        t = currtx()
        (t === nothing) && return true
        d = get(getfield(t, :streams), s, nothing)
        return !(d !== nothing && getfield(d, :closed))
    finally
        unlock(mux)
    end
end

# The logical number of values in the stream. Inside a transaction this
# includes the transaction's own staged puts and excludes its staged takes, so
# it is a read-your-writes view; unlike a `Volatile` read it is an advisory
# snapshot rather than a validated read (see the Streams page).
function Base.length(s::Stream)
    mux = getfield(s, :mutex)
    lock(mux)
    try
        live = _live_count(s)
        t = currtx()
        (t === nothing) && return live
        d = get(getfield(t, :streams), s, nothing)
        (d === nothing) && return live
        return max(0, live - getfield(d, :ntaken)) + _staged_count(d)
    finally
        unlock(mux)
    end
end
Base.isempty(s::Stream) = length(s) == 0
# Ready means a `take` would return or throw immediately rather than block; a
# closed stream is ready even when it is empty, as a failed promise is.
Base.isready(s::Stream) = length(s) > 0 || !isopen(s)
Base.eltype(::Type{Stream{T}}) where {T} = T

"""
    snapshot(stream) -> Vector{T}

Yields a copy of the values the given [`Stream`](@ref) is holding, oldest
first, without blocking and without removing them. The copy is independent of
the stream, so it may be kept, mutated or serialized; whether the stream was
open and what its bound is are reported separately by [`isopen`](@ref) and
[`maxsize`](@ref).

Inside a transaction the snapshot is the transaction's own view: it includes
the values the body has staged with `put!` and excludes those it has staged
with `take`, so it agrees with `length`, `isempty` and `isready`. Like those,
it is an advisory snapshot rather than a read that is validated at commit.

See also: [`Stream`](@ref), [`take`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> s = Stream{Int}(maxsize = 4); put!(s, 1); put!(s, 2);

julia> snapshot(s)
2-element Vector{Int64}:
 1
 2

julia> take(s)          # a snapshot consumes nothing
1

julia> snapshot(s)
1-element Vector{Int64}:
 2
```
"""
function snapshot(s::Stream{T}) where {T}
    mux = getfield(s, :mutex)
    lock(mux)
    try
        items = getfield(s, :items)
        # Indexing a `Vector` with a range copies, which is the snapshot.
        out = items[getfield(s, :head):end]
        t = currtx()
        (t === nothing) && return out
        d = get(getfield(t, :streams), s, nothing)
        (d === nothing) && return out
        # Staged takes consume the front of the live items, then the staged
        # puts; `ntaken` may exceed the live count if another task has taken
        # some since, in which case nothing live remains.
        skip = min(getfield(d, :ntaken), length(out))
        (skip > 0) && deleteat!(out, 1:skip)
        ph = getfield(d, :puts_head)
        ps = getfield(d, :puts)
        (ph <= length(ps)) && append!(out, @view ps[ph:end])
        return out
    finally
        unlock(mux)
    end
end
export snapshot

# Streams are handles, so equality and hashing are by identity, as for a
# promise; a stream is a perfectly good dictionary key.
Base.isequal(a::Stream, b::Stream) = a === b
Base.hash(s::Stream, h::UInt) = hash(objectid(s), hash(0xb3d1f4a07c2e5911, h))

function Base.show(io::IO, s::Stream)
    mux = getfield(s, :mutex)
    lock(mux)
    try
        print(io, "Stream{", eltype(s), "}(", _live_count(s), " items")
        m = getfield(s, :maxsize)
        (m === nothing) || print(io, ", maxsize ", m)
        print(io, getfield(s, :closed) ? ", closed)" : ", open)")
    finally
        unlock(mux)
    end
    return nothing
end
Base.show(io::IO, ::MIME"text/plain", s::Stream) = show(io, s)

# #Transactions
#
# A transaction's staged view of one stream. The logical contents are the live
# items after `ntaken` of them, followed by the staged puts from `puts_head`
# onward. `touched` records that this attempt has actually staged a take, a put
# or a close: while it is false the attempt has no stream-derived state, so a
# blocking operation may adopt a newer version instead of restarting.

"""
    StreamTxData{T}

A transaction's staged view of a [`Stream`](@ref): the version it first
observed, how many live items it has taken, its staged puts, a staged close,
and whether it created the stream. This is part of `Air`'s private/internal
implementation details.
"""
mutable struct StreamTxData{T}
    start_version::UInt64
    ntaken::Int
    puts::Vector{T}
    puts_head::Int
    closed::Bool
    born::Bool
    touched::Bool
end

# The transaction's staging data for `s`, created on first use. The caller holds
# the stream's mutex. While the attempt has staged nothing for this stream, a
# newer version is adopted: the attempt has no earlier stream-derived state to
# invalidate, so waiting on the stream is a legitimate observation of the new
# state.
function _stream_txdata_locked!(s::Stream{T}, t) where {T}
    d = get(getfield(t, :streams), s, nothing)
    if d === nothing
        d = StreamTxData{T}(getfield(s, :version), 0, T[], 1, false, false, false)
        getfield(t, :streams)[s] = d
    else
        d = d::StreamTxData{T}
        getfield(d, :touched) ||
            setfield!(d, :start_version, getfield(s, :version))
    end
    return d
end
# A staged attempt whose stream has changed since its first staged operation
# must restart: its staged counts no longer describe the live items, and the
# body may already have branched on what it saw.
@inline function _stream_tx_recheck!(s::Stream, d::StreamTxData)
    (getfield(d, :touched) &&
     getfield(d, :start_version) != getfield(s, :version)) &&
        throw(TxRetryException())
    return nothing
end

# Take the next logical item, staging the take. The caller holds the mutex.
function _try_tx_take_locked!(s::Stream{T}, d::StreamTxData{T}) where {T}
    live = _live_count(s)
    nt = getfield(d, :ntaken)
    if live > nt
        x = getfield(s, :items)[getfield(s, :head) + nt]
        setfield!(d, :ntaken, nt + 1)
    else
        ph = getfield(d, :puts_head)
        ps = getfield(d, :puts)
        if ph > length(ps)
            (getfield(s, :closed) || getfield(d, :closed)) &&
                throw(StreamClosed())
            return nothing
        end
        x = ps[ph]
        setfield!(d, :puts_head, ph + 1)
        if ph + 1 > length(ps)
            empty!(ps)
            setfield!(d, :puts_head, 1)
        end
    end
    if !getfield(d, :touched)
        setfield!(d, :start_version, getfield(s, :version))
    end
    setfield!(d, :touched, true)
    return _StreamItem{T}(x)
end
# Stage a put if the stream is open and the logical queue has room; the caller
# holds the mutex. Yields whether it did, and throws if the stream is closed.
function _try_tx_push_locked!(s::Stream{T}, d::StreamTxData{T}, x::T) where {T}
    (getfield(s, :closed) || getfield(d, :closed)) && throw(StreamClosed())
    m = getfield(s, :maxsize)
    if m !== nothing
        logical = max(0, _live_count(s) - getfield(d, :ntaken)) + _staged_count(d)
        logical < m || return false
    end
    push!(getfield(d, :puts), x)
    if !getfield(d, :touched)
        setfield!(d, :start_version, getfield(s, :version))
    end
    setfield!(d, :touched, true)
    return true
end

# The born marker is installed by the constructor and never changes, and the
# staging dictionary is task-local, so reading it outside the stream's mutex is
# safe.
@inline function _stream_born(t, s)
    d = get(getfield(t, :streams), s, nothing)
    return d !== nothing && getfield(d, :born)
end

function _stream_take_tx(s::Stream{T}, t, timeout,
                         token::Union{Nothing,CancelToken}) where {T}
    mux = getfield(s, :mutex)
    w = nothing
    timer = nothing
    try
        while true
            (token === nothing) || checkcancel(token)
            lock(mux)
            r = try
                d = _stream_txdata_locked!(s, t)
                _stream_tx_recheck!(s, d)
                _try_tx_take_locked!(s, d)
            finally
                unlock(mux)
            end
            (r === nothing) || return r.value
            _stream_born(t, s) && throw(StreamDeadlock())
            (timeout === nothing) || (timeout >= 0) ||
                throw(ArgumentError("timeout must be non-negative"))
            if w === nothing
                w = Waiter()
                lock(mux)
                try
                    _stream_register_locked!(s, w)
                finally
                    unlock(mux)
                end
                (token === nothing) || _register_waiter!(token, w)
                (timeout === nothing) || (timer = Timer(_ -> _timeout!(w), timeout))
            end
            if _park!(w)
                lock(mux)
                r = try
                    d = _stream_txdata_locked!(s, t)
                    _stream_tx_recheck!(s, d)
                    _try_tx_take_locked!(s, d)
                finally
                    unlock(mux)
                end
                (r === nothing) || return r.value
                throw(ReaderTimeout(timeout))
            end
        end
    finally
        (w === nothing) || _stream_unregister!(s, w)
        ((w === nothing) || (token === nothing)) ||
            _unregister_waiter!(token, w)
        (timer === nothing) || close(timer)
    end
end

function _stream_put_tx(s::Stream{T}, xT::T, t, timeout,
                        token::Union{Nothing,CancelToken}) where {T}
    mux = getfield(s, :mutex)
    w = nothing
    timer = nothing
    try
        while true
            (token === nothing) || checkcancel(token)
            lock(mux)
            ok = try
                d = _stream_txdata_locked!(s, t)
                _stream_tx_recheck!(s, d)
                _try_tx_push_locked!(s, d, xT)
            finally
                unlock(mux)
            end
            ok && return xT
            _stream_born(t, s) && throw(StreamDeadlock())
            (timeout === nothing) || (timeout >= 0) ||
                throw(ArgumentError("timeout must be non-negative"))
            if w === nothing
                w = Waiter()
                lock(mux)
                try
                    _stream_register_locked!(s, w)
                finally
                    unlock(mux)
                end
                (token === nothing) || _register_waiter!(token, w)
                (timeout === nothing) || (timer = Timer(_ -> _timeout!(w), timeout))
            end
            if _park!(w)
                lock(mux)
                ok = try
                    d = _stream_txdata_locked!(s, t)
                    _stream_tx_recheck!(s, d)
                    _try_tx_push_locked!(s, d, xT)
                finally
                    unlock(mux)
                end
                ok && return xT
                throw(ReaderTimeout(timeout))
            end
        end
    finally
        (w === nothing) || _stream_unregister!(s, w)
        ((w === nothing) || (token === nothing)) ||
            _unregister_waiter!(token, w)
        (timer === nothing) || close(timer)
    end
end

# Apply a committed transaction's staged operations to `s`. The caller holds the
# mutex and has validated `start_version`, so the live prefix it counts is still
# intact. Removals happen before appends, so the FIFO order is preserved.
function _stream_commit!(s::Stream{T}, d::StreamTxData{T}) where {T}
    changed = false
    n = getfield(d, :ntaken)
    if n > 0
        setfield!(s, :head, getfield(s, :head) + n)
        _stream_compact!(s)
        changed = true
    end
    ph = getfield(d, :puts_head)
    ps = getfield(d, :puts)
    if ph <= length(ps)
        append!(getfield(s, :items), @view ps[ph:end])
        changed = true
    end
    if getfield(d, :closed) && !getfield(s, :closed)
        setfield!(s, :closed, true)
        changed = true
    end
    if changed
        setfield!(s, :version, getfield(s, :version) + UInt64(1))
        _wake_waiters!(s)
    end
    return nothing
end
