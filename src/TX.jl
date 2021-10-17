################################################################################
# tx.jl
#
# Functional multi-threading tools built around atomic transactions.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

import DataStructures

# The beginning of this file contains the various types and simple functions
# that depend only on those types (like show). These types end with the
# Transaction type.
# The remainder of this file contains methods that depend on transaction
# details. These include all read and write operations for Actors and Volatiles.


# #Actor #######################################################################
# The (private) type of a message that is queued for an actor to eventually run.
"""
    ActorMsg

A message, queued for an `Actor` via the `send()` function.

The `ActorMsg` struct is considered part of `Air`'s private/internal
implementation details.
"""
struct ActorMsg
    fn::Function
end
const ActorMsgQueue = DataStructures.Queue{ActorMsg}
"""
    ActorException{T}(error, value, argno, args)

An ActorException object is thrown whenever one attempts to obtain the value
of or send a function to an actor that is in an error state. An error state
occurs when an unhandled exception is raised while an actor is processing a
sent function. In such a case the `geterror()` and `restart()` functions may
be used. The `geterror()` function yields an `ActorException` object in which
the exception that was raised is stored as `error`, the value of the actor when
the send was run is stored in `value`, and the `args` and `argno` give the
full contex of the `send` call (i.e., the fucntion followed by the arguments is
`args` and the argument number that corresponds to the actor is `argno`).

See also: [`Actor`](@ref), [`geterror`](@ref), [`reset`](@ref), [`send`](@ref)
"""
struct ActorException{T} <: Exception
    error::Any
    value::T
    fn::Function
    queue::Tuple
end
export ActorException
Base.show(io::IO, ::MIME"text/plain", a::ActorException{T}) where {T} = begin
    print(io, "ActorException{$T}($(typeof(a.error)), $(a.value), ...)")
end
Base.show(io::IO, a::ActorException{T}) where {T} = begin
    print(io, "ActorException{$T}($(typeof(a.error)))")
end
const ActorValue{T} = Union{Some{T}, ActorException{T}} where {T}
"""
    Actor{T}

An actor is an object that represents a worker-thread to which tasks can be
scheduled. Any scheduled function is guaranteed to be evaluated at some point in
the future in some other thread, and the return value of that function will
become the new value held by the actor. Each function, when it is run, is passed
the actor's value as one of its argument.

Like with `Ref`s, you can access an actor's current value using `getindex`, 
e.g., `actor[]`. Within a transaction, this will be guaranteed to remain fixed
for the duration of the transaction; outside of a transaction, this may change
at arbitrary times.

`Actor`s may be initialized with post-processing functions. These functions are
called immediately after executing any scheduled function and are given as their
first argument the original actor value and as their second argument the return
value of the scheduled function. Instead of saving this return value, the
return value of the post-processing function is saved in the actor. This is
useful for actors whose job is, for example, to serialize output to a log file
being written to across many threads. If we run something like 
`send(println, log_actor, "Some log message")` the `println` function will
return `nothing`, which we would like to convert back into the log stream so
that subsequent `send` calls can continue to use the `println` fuction.

Actors may additionally be given an error handler. Whenever an exception occurs
during a scheduled function, the error handler will be called with the arguments
of (1) the actor object, (2) the current actor value, and (3) the exception that
was caught. When an error has occurred, any attempt to schedule a function to
the actor or to read from the actor will cause an exception to be raised. The
error may be examined using `geterror(actor)` and restarted using
`reset(actor)`.

All fields in an actor object should be considered strictly private. These
fields are likely to change between releases, and changing the values will break
your code.

See also: [`tx`](@ref), [`@tx`](@ref), [`send`](@ref), [`geterror`](@ref),
[`reset`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Actor{Symbol}\\(@[0-9a-zA-Z]+: :start_sym\\)"
julia> a = Actor{Symbol}(:start_sym)
Actor{Symbol}(@JIm7aUS2sOl: :start_sym)

julia> send(a) do val; Symbol("new_\$(val)") end
Actor{Symbol}(@JIm7aUS2sOl: :start_sym)

julia> sleep(1); a[]
:new_start_sym
```
"""
mutable struct Actor{T} <: TransactionalRef{T}
    mutex::ReentrantLock
    queue::ActorMsgQueue
    task::Union{Nothing, Task}
    value::ActorValue{T}
    function Actor{T}(initval) where {T}
        return new{T}(ReentrantLock(), ActorMsgQueue(), nothing, Some{T}(initval))
    end
end
export Actor
Actor(t::T) where {T} = Actor{T}(t)
Base.propertynames(::Actor) = (:state,:value,:error)
Base.getproperty(a::Actor{T}, k::Symbol) where {T} = begin
    if k == :state
        return (geterror(a) === nothing ? :okay : :error)
    elseif k == :value
        return a[]
    elseif k == :error
        return geterror(a)
    else
        error("type $(typeof(a)) has no field $k")
    end
end
Base.show(io::IO, ::MIME"text/plain", a::Actor{T}) where {T} = begin
    print(io, "Actor{$T}(@")
    print(io, string(objectid(a), base=62))
    print(io, ": ")
    s = getfield(a, :value)
    if isa(s, ActorException{T})
        print(io, "--error--")
    else
        show(io, something(s))
    end
    print(io, ")")
end
Base.show(io::IO, a::Actor{T}) where {T} = begin
    print(io, "Actor{$T}(@")
    print(io, string(objectid(a), base=62))
    print(io, ")")
end
Base.setindex!(a::Actor{T}, args...) where {T} = error(
    "Actor objects cannot be assigned---they can only be sent functions")
# The private main-loop that processes an actor's messages.
"""
    actor_main(actor)

The actor main function is what actually manages everything behind the scenes in
the actor tasks. Note that this function and the task/thread it is running in is
the only place that the actor's value can be legally changed. Changing the value
anywhere else can result in undefined behavior. The actor's condition should
never change.

This function is designed to process all the messages in an actor's queue then
to return. When a new thread wants to enqueue a message to an actor that has an
empty queue, that thread must spin up the `actor_main` function to process the
items in the queue (this is done by the `send` function automatically).

This function is part of `Air`'s internal/private implementation details.
"""
actor_main(a::Actor{T}) where {T} = begin
    mux = getfield(a, :mutex)
    msg = nothing
    args = nothing
    rval = nothing
    # For starters, lock the condition!
    lock(mux)
    try
        queue = getfield(a, :queue)
        value = getfield(a, :value)
        isa(value, ActorException{T}) && return nothing
        # This function is a loop! We process everything in the queue until it
        # is empty, then we exit. If the queue is empty right now, then there's
        # nothing to process (which would be weird).
        while !isempty(queue)
            val = value.value # (value is a Some{T})
            # get the first item in the queue
            msg = DataStructures.dequeue!(queue)
            # Having grabbed that message, we can go ahead and unlock.
            unlock(mux)
            # Okay, we've grabbed the message we want to process, and we've
            # unlocked so that we can do the processing without halting other
            # tasks from queuing up new messages.
            # Note that even though we don't hold the condition, this task is
            # the only one that's allowed to update the value, so the value is
            # effectively task-local.
            f = msg.fn
            try
                value = Some{T}(f(val))
            catch e
                # There was an exception, so we need to lock the actor and set
                # the error state.
                lock(mux)
                e = ActorException{T}(e, val, f, (queue...,))
                setfield!(a, :value, e)
                empty!(queue)
                return nothing
            end
            # Relock the condition to update the value.
            lock(mux)
            setfield!(a, :value, value)
            # At this point we have updated the actor's value and we still
            # hold the lock. We can go ahead and enter the next loop.
        end
        # If we reach this point, then the queue is empty, so we can exit.
    finally
        setfield!(a, :task, nothing)
        unlock(mux)
    end
    return nothing
end
"""
    actor_start!(actor)

If the given `Actor` object is not in an error state, is not running, and does
not have an empty queue, then `actor_start!(actor)` will start `actor`'s
processing task. This should only be called by `send!` when the first message to
an empty queue is finalized.

If either the actor is already running or the task is successfully started, then
that task is yielded. If the actor has an empty queue, then nothing is returned,
and if the actor is in an error state, then an exception is thrown.

Note: It is required that the mutx for the given actor be held when this
function is called; otherwise race conditions could be inadvertantely created.
The [`actor_start`](@ref) method locks the actor's mutex then calls this
function. Both methods are private to `Air` and should not generally be used
outside of it.

This function is considered part of `Air`'s internal/private implementation
details.
"""
actor_start!(a::Actor{T}) where {T} = begin
    task = getfield(a, :task)
    (task === nothing) || return task
    queue = getfield(a, :queue)
    isempty(queue) && return nothing
    isa(getfield(a, :value), ActorException{T}) && return nothing
    # It needs a task! Start one up!
    task = Threads.@spawn actor_main(a)
    setfield!(a, :task, task)
    return task
end
"""
    actor_start(actor)

If the given `Actor` object is not in an error state, is not running, and does
not have an empty queue, `actor_start(actor)` will start `actor`'s processing
task. This should only be called by `send!` when the first message to an empty
queue is finalized.

If either the actor is already running or the task is successfully started, then
that task is yielded. If the actor has an empty queue, then nothing is returned,
and if the actor is in an error state, then an exception is thrown.

This function is considered part of `Air`'s internal/private implementation
details.
"""
actor_start(a::Actor{T}) where {T} = begin
    mux = getfield(a, :mutex)
    lock(mux)
    try
        return actor_start!(a)
    finally
        unlock(mux)
    end
end
"""
    actor_send!(actor, message)

Enqueues the given message to the given actor and starts the actor's task, if
necessary. This function requires that the actor's condition be held. Yields the
actor.

This function is considered part of `Air`'s internal/private implementation
details.
"""
actor_send!(a::Actor{T}, msg::ActorMsg) where {T} = begin
    val = getfield(a, :value)
    isa(val, ActorException{T}) && throw(val)
    DataStructures.enqueue!(getfield(a, :queue), msg)
    (getfield(a, :task) === nothing) && actor_start!(a)
    return a
end
actor_send!(a::Actor{T}, qq::ActorMsgQueue) where {T} = begin
    val = getfield(a, :value)
    isa(val, ActorException{T}) && throw(val)
    if !isempty(qq)
        for msg in qq
            DataStructures.enqueue!(getfield(a, :queue), msg)
        end
        (getfield(a, :task) === nothing) && actor_start!(a)
    end
    return a
end
"""
    actor_reset!(actor, val)

Resets the given actor to have the given value. If the actor is not in an error
state, yields `nothing` and does nothing. Otherwise, yields the `ActorExcption`
object that was cleared. It is required that the actor's mutex be held at ths
time this function is called.

This function is considered part of `Air`'s internal/private implementation
details.
"""
actor_reset!(a::Actor{T}, val) where {T} = begin
    val0 = getfield(a, :value)
    isa(val0, ActorException{T}) || return nothing
    setfield!(a, :value, Some{T}(val))
    empty!(getfield(a, :queue))
    # Since the queue is, as of yet, empty, there's no need to start a task.
    return val0
end


# #Volatile ####################################################################
# Volatile's, like Actor's, have a single immutable state structure.
struct VolatileData{T}
    value::T
    filter::Union{Nothing, Function}
    finalize::Union{Nothing, Function}
    function VolatileData{T}(v::S, flt::Function, fin::Function) where {T,S}
        return new{T}(flt(v), flt, fin)
    end
    function VolatileData{T}(v::S, flt::Function, ::Nothing) where {T,S}
        return new{T}(flt(v), flt, nothing)
    end
    function VolatileData{T}(v::S, ::Nothing, fin::Function) where {T,S}
        return new{T}(v, nothing, fin)
    end
    function VolatileData{T}(v::S, ::Nothing, ::Nothing) where {T,S}
        return new{T}(v, nothing, nothing)
    end
    function VolatileData{T}(v::S, flt::Function) where {T,S}
        return new{T}(flt(v), flt, nothing)
    end
    function VolatileData{T}(v::S, ::Nothing) where {T,S}
        return new{T}(v, nothing, nothing)
    end
    function VolatileData{T}(v::S) where {T,S}
        return new{T}(v, nothing, nothing)
    end
end
"""
    voldata_finalize(voldata)

Calls the finalize function of the given volatile state object (of type
`VolatileData`) if necessary. Yields the new `VolatileData` object, which may be
`voldata` if there is no finalizer or nothing has changed.

This function is a private/internal implementation detail of the `Air` library.
"""
voldata_finalize(v::VolatileData{T}) where {T} = begin
    (v.finalize === nothing) && return v
    u = v.finalize(v.value)
    return (u === v.value ? v : VolatileData{T}(u, v.filter, v.finalize))
end
"""
    Volatile{T}(value)
    Volatile{T}(value, filter_fn)
    Volatile{T}(value, filter_fn, finalize_fn)

`Volatile` objects are `Ref` objects that must be used in conjunction with
transaction blocks (see [`@tx`](@ref) and [`tx`](@ref)). The value of a
`Volatile` can be accessed at any time, and there is no particular guarantee
that a `Volatile`'s value won't be changed by another thread outside of a
transaction, thus reading them outside of a transaction can potentially create
race conditions. Critically, `Volatile` objects can only be set inside of a
transaction, and, within a transaction, a `Volatile`'s value is guaranteed to be
constant. In other words, `Volatile` objects are `Ref` objects that behave as if
they are locked for the executing thread whenever that thread is inside a
transaction block.

The arguments `filter_fn` and `finnalize_fn` are functions for making sure that
the value of a volatile conforms to some standard. The two functions are similar
but are intended for slightly different use cases:
 * `filter_fn` is called every time the value of a volatile is set. The value
   saved in the volatile is `filter_fn(value)` instead of `value` itself. Values
   are only filtered when the volatile is set directly (`vol[] = value`)---the
   filter is not rerun when the `filter_fn` is changed (via `setfilter!`) nor
   when the `finalize_fn` returns a value.
 * `finalize_fn` is called with the volatile's value immediately after the body
   of the transaction has completed, but before a commit is attempted. The
   `finalize_fn` function acts like the `filter_fn` in that `finalize_fn(value)`
   replaces `value` in the volatile.
Both `filter_fn` and `finalize_fn` can throw exceptions, which will abort the
running transaction.

For a `Volatile` `v`, one may set `v`'s stored value via `v[] = value`. This must
be done inside a transaction (via `tx` or `@tx`). To change the filter or the
finalize functions, use the `setfilter!` and `setfinalize!` functions, both
of which also must be run in transactions as well.

See also: [`tx`](@ref), [`@tx`](@ref), [`setfilter!`](@ref),
[`setfinalize!`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Volatile{Symbol}\\(@[0-9a-zA-Z]+: :startval\\)"
julia> v = Volatile{Symbol}(:startval)
Volatile{Symbol}(@JIm7aUS2sOl: :startval)

julia> @tx v[] = :newval
:newval

julia> v[]
:newval
```
"""
mutable struct Volatile{T} <: TransactionalRef{T}
    mutex::ReentrantLock
    value::VolatileData{T}
    function Volatile{T}(val::S) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, nothing, nothing))
    end
    function Volatile{T}(val::S, f::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, f, nothing))
    end
    function Volatile{T}(val::S, ::Nothing) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, nothing, nothing))
    end
    function Volatile{T}(val::S, f::Function, ::Nothing) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, f, nothing))
    end
    function Volatile{T}(val::S, f::Function, ff::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, f, ff))
    end
    function Volatile{T}(val::S, ::Nothing, ff::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, nothing, ff))
    end
    function Volatile{T}(val::S, ::Nothing, ::Nothing) where {T,S}
        return new{T}(ReentrantLock(), VolatileData{T}(val, nothing, nothing))
    end
end
export Volatile
Volatile(t::T) where {T} = Volatile{T}(t)
Volatile(t::T, filterfn::Union{Nothing,Function},
         finfn::Union{Nothing,Function}=nothing) where {T} = Volatile{T}(t, filterfn, finfn)
Base.show(io::IO, ::MIME"text/plain", v::Volatile) = begin
    show(io, typeof(v))
    print(io, "(@")
    print(io, string(objectid(v), base=62))
    print(io, ": ")
    show(io, v[])
    print(io, ")")
end
Base.show(io::IO, v::Volatile) = begin
    show(io, typeof(v))
    print(io, "(@")
    print(io, string(objectid(v), base=62))
    print(io, ")")
end
Base.propertynames(::Volatile) = (:value,)
Base.getproperty(v::Volatile, p) = begin
    if p == :value
        return v[]
    else
        error("type $(typeof(v)) has no field $p")
    end
end


# #Transaction #################################################################
"""
    ActorTxData{T}

The data tracked for an `Actor{T}` during a transaction by the transaction's 
`Transaction` struct.

The `ActorTxData` struct is considered part of the internal/private code
of Air and generally should not be used outside of the Air library.
"""
mutable struct ActorTxData{T}
    start_value::ActorValue{T}
    tx_value::ActorValue{T}
    reset_value::Union{Nothing, Some{T}}
    msgs::ActorMsgQueue
    function ActorTxData{T}(w::S) where {T,S}
        return new{T}(w, w, nothing, ActorMsgQueue())
    end
end
"""
    TxRetryException

A `TxRetryException` is throws when a transaction needs to be retried but an
error hasn't necessarily been generated otherwise. This can be thrown during a
transaction to restart the transaction; though doing this can easily lead to
infinite loops.
"""
struct TxRetryException <: Exception end
export TxRetryException
"""
    Transaction

A transaction object keeps track of what is going on during a particular
transaction. These are generally low-level objects that shouldn't be
touched directly.

Transactions have the following propertiies:
* `state` is either `:running`, `:validating`, or `:error`;
* `rvolatiles` is the set of all `Volatile` objects that have been read during
  the transaction;
* `wvolatiles` is the set of all `Volatile` objects that have been changed;
* `actors` is the set of all `Actor` objects to which messagees have been sent
  during the transaction; and
* `sources` is the set of all `Source` objects from which items have been
  popped.
"""
mutable struct Transaction
    state::Symbol
    reads::IdDict{Volatile, VolatileData}
    writes::IdDict{Volatile, NTuple{2,VolatileData}}
    actors::IdDict{Actor, ActorTxData}
    Transaction() = new(:running,
                        IdDict{Volatile, VolatileData}(),
                        IdDict{Volatile, NTuple{2,VolatileData}}(),
                        IdDict{Actor, ActorTxData}())
end
export Transaction
Base.propertynames(::Transaction) =
    (:state, :rvolatiles, :wvolatiles, :actors)
Base.getproperty(t::Transaction, p) = begin
    if p == :state
        return getfield(t, :state)
    elseif p == :rvolatiles
        return Base.IdSet(keys(getfield(p, :reads)))
    elseif p == :wvolatiles
        return Base.IdSet(keys(getfield(p, :writes)))
    elseif p == :actors
        return Base.IdSet(keys(getfield(p, :actors)))
    else
        error("type $(typeof(t)) has no field $p")
    end
end
tx_clear!(t::Transaction) = begin
    setfield!(t, :state, :running)
    empty!(getfield(t, :reads))
    empty!(getfield(t, :writes))
    empty!(getfield(t, :actors))
    return nothing
end

# The current running transaction; see currtx() below.
"""
    current_tx

The `current_tx` constant is a `Var{T}` that stores the current task's running
`Transaction`, or `nothing` if there is no transaction running.

See also: [`currtx`](@ref)
"""
@var current_tx = nothing::Union{Nothing, Transaction}

"""
    currtx()

Yields the current `Transaction` object, or `nothing` if there is no transaction
running in the current task.

See also: [`Transaction`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> currtx() === nothing
true

julia> @tx (currtx() === nothing)
false
```
"""
currtx() = current_tx[]
export currtx

"""
    TX_MAX_ATTEMPTS

The maximum number of retry attempts that a transaction will make before
aborting the transaction.
"""
const TX_MAX_ATTEMPTS = 2^14

"""
    tx(fn)

Runs the given function in a transaction and yields the result. The function
fn is called as `fn()` without arguments.
"""
function tx(fn::F) where {F <: Function}
    success = false
    res = nothing
    the_tx = current_tx[]
    # If there's already a transaction running, we needn't make a new one---
    # This transaction function will just get rolled up into the current one.
    (the_tx === nothing) || return fn()
    the_tx = Transaction()
    # We need to make a new transaction for this task. We might have to do
    # this a few times if the transaction fails.
    for attempt in 1:TX_MAX_ATTEMPTS
        # Clear the transaction
        tx_clear!(the_tx)
        try
            res = withvars(fn, current_tx => the_tx)
        catch e
            if isa(e, TxRetryException)
                continue
            else
                rethrow()
            end
        end
        # At this point the transaction is over and we need to ensure that the
        # current volatile values are consistent with the transaction we just
        # executed. Before we do that, we should finalize the volatiles.
        writes = getfield(the_tx, :writes)
        setfield!(the_tx, :state, :finalizing)
        for (v,(x0,x1)) in writes
            x = voldata_finalize(x1)
            if x !== x1
                setfield!(the_tx, :state, :running)
                volatile_setindex!(v, x, the_tx)
                setfield!(the_tx, :state, :finalizing)
            end
        end
        # We're now ready to attempt a commit. To do this we must lock the vols
        # and actors in a sorted order. We do volatiles first, since their
        # values must all be tested, then actors.
        reads = getfield(the_tx, :reads)
        actors = getfield(the_tx, :actors)
        setfield!(the_tx, :state, :locking)
        nw = length(writes)
        nr = length(reads)
        n  = nw + nr
        m  = length(actors)
        vols = Vector{Volatile}(undef, n)
        vols[1:nr] .= keys(reads)
        vols[nr+1:end] .= keys(writes)
        sort!(vols, by=objectid)
        acts = Vector{Actor}(undef, m)
        acts[1:m] .= keys(actors)
        sort!(acts, by=objectid)
        locked = 0
        checked = 0
        success = false
        try
            while !success
                # First, go through and lock all the volatiles in order.
                for vol in vols
                    lock(getfield(vol, :mutex))
                    locked += 1
                end
                # Now check their values. First do the read volatiles.
                for (vol,x0) in reads
                    (x0 === getfield(vol, :value)) || break
                    checked += 1
                end
                (checked == nr) || break
                # Next, check the written volatiles.
                for (vol,(x0,x)) in writes
                    (x0 === getfield(vol, :value)) || break
                    checked += 1
                end
                (checked == n) || break
                # Next we need to lock and check the actors that were read from.
                for act in acts
                    lock(getfield(act, :mutex))
                    locked += 1
                end
                for (act,dat) in actors
                    val = getfield(act, :value)
                    (dat.start_value === val) || break
                    checked += 1
                end
                (checked == n + m) || break
                # If we reach this point, then everything checks out and we can
                # commit all of the changes and complete the transaction.
                setfield!(the_tx, :state, :committing)
                for (vol,(x0,x)) in writes
                    setfield!(vol, :value, x)
                end
                # We need to commit the actors as well.
                for (act,dat) in actors
                    q = dat.msgs
                    r = dat.reset_value
                    (r === nothing) || actor_reset!(act, r.value)
                    isempty(q) || actor_send!(act, q)
                end
                # That's all we need to do, aside from unlock.
                success = true
            end
        finally
            for mux in 1:locked
                (mux > n) && break
                unlock(getfield(vols[mux], :mutex))
            end
            locked -= n
            if locked > 0
                for mux in 1:locked
                    #(mux > m) && break
                    unlock(getfield(acts[mux], :mutex))
                end
            end
        end
        success && break
    end
    # It's possible we got here because we were suceessful, but it might be that
    # we failed too many times.
    if success
        return res
    else
        error("transaction aborted after failing $TX_MAX_ATTEMPTS times")
    end
end
export tx

"""
    @tx expr

The macro @tx should be followed by an expression; that expression is run in an
atomic transaction.
"""
macro tx(expr)
    return :(tx(() -> $(esc(expr))))
end
export @tx

# Now that we have definend transactions, we can define the volatile access
# methods.
_volatile_getindex(v::Volatile{T}, t::Transaction) where {T} = begin
    w = get(getfield(t, :writes), v, nothing)
    (w === nothing) || return w[2]
    w = get(getfield(t, :reads), v, nothing)
    (w === nothing) || return w
    w = getfield(v, :value)
    getfield(t, :reads)[v] = w
    return w
end
_volatile_getindex(v::Volatile{T}, ::Nothing) where {T} = getfield(v, :value)
Base.getindex(v::Volatile{T}) where {T} = _volatile_getindex(v, currtx()).value
volatile_setindex!(v::Volatile{T}, x::VolatileData{T}, t::Transaction) where {T} = begin
    (getfield(t, :state) === :running) || error(
        "volatiles can only be set before finalizing the transaction")
    # See if it this ref has already been written to in this transaction.
    w = get(getfield(t, :writes), v, nothing)
    if w !== nothing
        (x0,x1) = w
        # If this doesn't change anything, don't do aything.
        (x === x1) && return x
        # If this changes the volatile back to its inintial value, we convert
        # this to a read operation. Otherwise we record the new value.
        if x === x0
            delete!(getfield(t, :writes), v)
            getfield(t, :reads)[v] = x0
        else
            getfield(t, :writes)[v] = (x0, x)
        end
        return x
    end
    # See if this volatile has already been read.
    w = get(getfield(t, :reads), v, nothing)
    if w !== nothing
        # There is a previously-read value already. If it hasn't changed, it's
        # still just a read value.
        (w === x) && return x 
        # Otherwise, we are going to promote it to a written value.
        delete!(getfield(t, :reads), v)
        getfield(t, :writes)[v] = (w, x)
    else
        # There is no previously-read value; grab the current value.
        w = getfield(v, :value)
        if w === x
            getfield(t, :reads)[v] = w
        else
            getfield(t, :writes)[v] = (w, x)
        end
    end
    return x
end
volatile_setindex!(v::Volatile{T}, x::S, ::Nothing) where {T,S} = error(
    "cannot set volatile outside of transaction")
Base.setindex!(v::Volatile{T}, x::S) where {T,S} = begin
    value = getfield(v, :value)
    # Run the filter on the new value.
    (value.filter === nothing) || (x = value.filter(x))
    newdat = VolatileData{T}(x, value.filter, value.finalize)
    return volatile_setindex!(v, newdat, currtx())
end
"""
    getfilter(v)

Yields the filter-function for the volatile v.
"""
getfilter(v::Volatile{T}) where {T} = _volatile_getindex(v, currtx()).filter
export getfilter

"""
    getfinalize(v)

Yields the finalize-function for the volatile v.
"""
getfinalize(v::Volatile{T}) where {T} = _volatile_getindex(v, currtx()).finalize
export getfinalize
"""
    setfilter!(vol, fn)

Sets the filter-function associated with the Volatile object vol. Any time that
the vol is set (`vol[] = x`) the filter-function is called and the value saved
in `vol` is instead `fn(x)`. This must be called within a transaction.
"""
setfilter!(vol::Volatile{T}, f::Function) where {T} = begin
    value = getfield(v, :value)
    newdat = VolatileData{T}(value.value, f, value.finalize)
    return volatile_setindex!(v, newdat, currtx())
end
export setfilter!
"""
    setfinalize!(vol, fn)

Sets the finalize-function associated with the Volatile object vol. Any time
that a transaction contains a change to `vol`, immediately prior to making an
attempt at committing the transaction, the finalize function is called and
the value committed to `vol` is instead `fn(x)` where `x` is the value set to
`vol` in the transaction. This funvtion must also be called within a
transaction.
"""
setfinalize!(vol::Volatile{T}, f::Function) where {T} = begin
    value = getfield(v, :value)
    newdat = VolatileData{T}(value.value, value.filter, f)
    return volatile_setindex!(v, newdat, currtx())
end
export setfinalize!

# Same for actors: we can put the actor code here.
# First, we definne some private functions for handling things on the back
# end.
"""
    tx_actordata(a)

Yields the in-transaction data for the `Actor` object `a`. If `a` is in an
error-state, then throws an `ActorException`. If there is no transaction
currently running, yields `nothing`.

This function is considered part of the internal/private interface of `Air` and
shouldn't generally be called outside of it.
"""
tx_actordata(a::Actor{T}) where {T} = tx_actordata(a, currtx())
tx_actordata(a::Actor{T}, t::Transaction) where {T} = begin
    # See if we've already grabbed it in the transaction.
    actors = getfield(t, :actors)
    w = get(actors, a, nothing)
    if w === nothing
        val = getfield(a, :value)
        w = ActorTxData{T}(val)
        actors[a] = w
    end
    return w
end
tx_actordata(a::Actor{T}, ::Nothing) where {T} = nothing

"""
    actor_value(a)

Yields the in-transaction value for the given `Actor` object `a` if there is a
running transaction; otherwise yields the current out-of-transaction value of
`a`.

This function is considered part of the internal/private interface of `Air` and
shouldn't generally be called outside of it.
"""
actor_value(a::Actor{T}) where {T} = actor_value(a, currtx())
actor_value(a::Actor{T}, t::Transaction) where {T} = tx_actordata(a, t).tx_value
actor_value(a::Actor{T}, ::Nothing) where {T} = getfield(a, :value)
Base.getindex(a::Actor{T}) where {T} = begin
    x = actor_value(a, currtx())
    isa(x, ActorException{T}) && throw(x)
    return x.value
end

"""
    actor_reset(a, s)

If the given `Actor` object `a` is in an error state, this function resets the
actor to have the new value `s` then yields the `ActorException` object that was
just cleared. If `a` is not in an error state, this function simply yields
`nothing`.

This function is considered part of the internal/private interface of `Air` and
shouldn't generally be called outside of it.
"""
actor_reset(a::Actor{T}, x::S) where {T,S} = actor_reset(a, x, currtx())
actor_reset(a::Actor{T}, newval::S, t::Transaction) where {T,S} = begin
    (getfield(t, :state) === :running) || error(
        "actors can only be reset before finalizing the transaction")
    w = tx_actordata(a, t)
    val = w.tx_value
    isa(val, ActorException{T}) || return nothing
    newval = Some{T}(newval)
    w.reset_value = newval
    w.tx_value = newval
    return val
end
actor_reset(a::Actor{T}, newval::S, ::Nothing) where {T,S} = begin
    # We're not in a transaction, so we just lock and queue.
    mux = getfield(a, :mutex)
    lock(mux)
    try
        val = getfield(a, :value)
        isa(val, ActorException{T}) || return false
        return actor_reset!(a, newval)
    finally
        unlock(mux)
    end
end
"""
    reset(actor, x)

If the given actor is in an error state, this (1) resets it, meaning it will
start handling sent messages again, (2) gives it the new initial value `x`, and
(3) yields the `ActorException` object that was just clared.  If the actor is
not in an error state, this just yields `nothing`.
"""
reset(a::Actor{T}, s::S) where {T,S} = actor_reset(a, s, currtx())
export reset

"""
    actor_send(a, msg)

Adds the given `ActorMsg` message object `msg` to the queue of the given `Actor`
object `a` such that it will be processed in a separate thread. Yields
`nothing`.

This function is considered part of the internal/private interface of `Air` and
shouldn't generally be called outside of it.
"""
actor_send(a::Actor{T}, msg::ActorMsg) where {T} = actor_send(a, msg, currtx())
actor_send(a::Actor{T}, msg::ActorMsg, t::Transaction) where {T} = begin
    (getfield(t, :state) === :running) || error(
        "actors can only be sent functions before finalizing the transaction")
    w = tx_actordata(a, t)
    isa(w.tx_value, ActorException{T}) && throw(w.tx_value)
    DataStructures.enqueue!(w.msgs, msg)
    return nothing
end
actor_send(a::Actor{T}, msg::ActorMsg, ::Nothing) where {T} = begin
    # We're not in a transaction, so we just lock and queue.
    mux = getfield(a, :mutex)
    lock(mux)
    try
        # There could be an error on the actor; if so we want to throw an
        # exception ourself.
        val = getfield(a, :value)
        isa(val, ActorException{T}) && throw(val)
        # The actor isn't in an error state; we can schedule the function.
        return actor_send!(a, msg)
    finally
        unlock(mux)
    end
end
"""
   send(fn, actor)

The `send()` method can be used to send a function to an actor. The function
`fn` is put in the actor's queue, which is processed sequentially by a separate
thread.  In the actual function call that is evaluated in this separate thread,
the function is passed the stored value of the actor, and the return value of
the function becomes the new stored value of the actor.

Note that the value of an actor when `send()` is called does not make any
guarantee about its value when the function that is sent gets evaluated---other
functions may be processed in the interim.

If `send()` is called inside of a transaction, then the function is not
immediately queued but rather is held until the transaction successfully
completes. Once this happens, all sends are dispatched simultaneously such that
two sends to the same actor during the same transaction will always run
back-to-back in the actor's processing thread.

**Important**. Due to the design pattern described in the previous paragraph,
one should never wait on the result of a function sent to an actor during a
transaction. This mistake might be made, for example, by sendng an actor, during
a transactoin, a function that delivers some result to a [`Promise`](@ref)
object. If the sending thread then waits on that promise while still in the same
transaction then it will deadlock because the sent function is not put in the
actor's queue until the enclosing transaction is finished.  Thus the promise
will never recieve the result from the actor.

TL;DR---When you touch an actor in a transaction, it's value freezes inside the
transaction. So if the transaction somehow waits for the actor to update or to
do something, that transaction will deadlock.
"""
send(f::Function, a::Actor{T}) where {T} = actor_send(a, ActorMsg(f), currtx())
export send

"""
    geterror(actor)

If the given `Actor` object is currently in an error state, then yields the
`ActorException` object that describes the error. Otherwise, yields `nothing`.

See also: [`Actor`](@ref), [`reset`](@ref), [`send`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Actor{Symbol}\\(@[0-9a-zA-Z]+: (:start_sym|--error--)\\)"
julia> a = Actor{Symbol}(:start_sym)
Actor{Symbol}(@JIm7aUS2sOl: :start_sym)

julia> geterror(a) === nothing
true

julia> send(a) do val; error("test") end
Actor{Symbol}(@JIm7aUS2sOl: :start_sym)

julia> sleep(1); a
Actor{Symbol}(@JIm7aUS2sOl: --error--)

julia> geterror(a) isa ActorException
true
```
"""
geterror(a::Actor{T}) where {T} = begin
    w = actor_value(a, currtx())
    return isa(w, ActorException{T}) ? w : nothing
end
export geterror

"""
    TxIO

The `TxIO` type is a derivative of the `IO` class and is intended for handling
output to streams across multiple threads in a thread-safe way. `TxIO` objects
can be constructed from `IO` objects and can be printed/written to like normal
`IO` objects. However, all output is performed in a separate asynchronous
thread, and writes that occur during a transaction always occur only if the
transaction succeeds.

See also: [`Actor`](@ref), [`airout`](@ref)
"""
struct TxIO <: IO
    io::IO
    actor::Actor{IO}
end
TxIO(io::IO) = begin
    iswritable(io) || error("TxIO works only with writable IO objects")
    return TxIO(io, Actor{IO}(io))
end
export TxIO
Base.isreadable(io::TxIO) = false
Base.iswritable(io::TxIO) = true
Base.eof(io::TxIO) = eof(io.io)
Base.write(io::TxIO, x::UInt8) = send(io.actor) do io;
    write(io, x)
    return io
end
Base.write(io::TxIO, x::Char) = send(io.actor) do io;
    write(io, x)
    return io
end
Base.write(io::TxIO, x::Union{String,SubString{String}}) = send(io.actor) do io;
    write(io, x)
    return io
end

"""
    AirOut

The `AirOut` type is a derivative of the `IO` class and is intended for handling
output to `stdout` across multiple threads in a thread-safe way. It is similar
`TxIO` except that it only writes to `stdout` and is only really intended for
use with the `airout` object.

See also: [`Actor`](@ref), [`airout`](@ref), [`TxIO`](@ref)
"""
struct AirOut <: IO
    actor::Actor{Nothing}
    function AirOut()
        return new(Actor{Nothing}(nothing))
    end
end
Base.isreadable(io::AirOut) = false
Base.iswritable(io::AirOut) = true
Base.eof(io::AirOut) = eof(stdout)
Base.write(io::AirOut, x::UInt8) = let so = stdout
    return send(io.actor) do io;
        write(so, x)
        return io
    end
end
Base.write(io::AirOut, x::Char) = let so = stdout
    return send(io.actor) do io;
        write(so, x)
        return io
    end
end
Base.write(io::AirOut, x::Union{String,SubString{String}}) = let so = stdout
    return send(io.actor) do io;
        write(so, x)
        return io
    end
end

"""
    airout

The `airout` object is a `AirOut` object that wraps the `stdout` object.

See also: [`AirOut`](@ref)
"""
const airout = AirOut()
export airout
