################################################################################
# tx.jl
#
# Functional multi-threading tools built around atomic transactions.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

using DataStructures: DataStructures

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
    # Deliberately unparameterized. An actor accepts arbitrary closures, so its
    # queue holds messages of many different function types and the call the
    # actor makes is dynamic however this struct is declared; parameterizing it
    # would not change that. Measured, sending is ~52 ns with about one
    # allocation, so there is nothing here to win.
    fn::Function
end
const ActorMsgQueue = DataStructures.Queue{ActorMsg}
"""
    ActorException{T}(error, value, fn, queue)

An `ActorException` object is thrown whenever one attempts to obtain the value
of or send a function to an actor that is in an error state. An error state
occurs when an exception is raised while an actor is processing a sent function
or its post-processing function. In such a case the `geterror` and `reset`
functions may be used. `geterror` yields an `ActorException` object in which the
exception that was raised is stored as `error`, the actor's value when the
failing function was scheduled is stored as `value`, the scheduled function
itself is `fn`, and `queue` is the tuple of messages that were still pending
when the error occurred.

See also: [`Actor`](@ref), [`geterror`](@ref), [`reset`](@ref), [`send`](@ref)
"""
struct ActorException{T} <: Exception
    error::Any
    value::T
    fn::Function
    queue::Tuple
end
export ActorException
function Base.show(io::IO, ::MIME"text/plain", a::ActorException{T}) where {T}
    return print(io, "ActorException{$T}($(typeof(a.error)), $(a.value), ...)")
end
function Base.show(io::IO, a::ActorException{T}) where {T}
    return print(io, "ActorException{$T}($(typeof(a.error)))")
end
const ActorValue{T} = Union{Some{T},ActorException{T}} where {T}
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

`Actor`s may be given two hooks, either at construction (with the `postprocess`
and `errorhandler` keywords) or with [`setpostprocess!`](@ref) and
[`seterrorhandler!`](@ref).

A *post-processing* function is called immediately after each scheduled function
has run, and is given the actor's value before the call as its first argument and
the scheduled function's return value as its second. The value it returns,
rather than the scheduled function's own return value, is stored in the actor
(and converted to `T`); the post-processing function is not called if the
scheduled function throws. This is useful for actors whose job is, for example,
to serialize output to a log file written to across many threads: if the actor's
value is the `IO` stream, a scheduled function such as
`io -> println(io, "message")` returns `nothing`, and a post-processor of
`(io, _) -> io` puts the stream back in the actor so that subsequent `send`
calls can keep writing to it.

An *error handler* is called when a scheduled function---or the post-processing
function---throws. It is given (1) the actor, (2) the actor's value when the
failing function was scheduled, and (3) the exception. The actor still enters
its error state: the handler is for observing or logging the failure, not for
recovering the actor, and its return value is ignored. An exception thrown by
the error handler itself is not caught and ends the actor's processing task.

Both hooks run on the actor's own task. Setting a hook affects the next function
dequeued, never one that is already running. When an error has occurred, any
attempt to schedule a function to the actor or to read from the actor will
cause an exception to be raised. The error may be examined using
`geterror(actor)` and cleared using `reset(actor)`.

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

julia> send(a) do val; sleep(0.1); Symbol("new_\$(val)") end
Actor{Symbol}(@JIm7aUS2sOl: :start_sym)

julia> sleep(1); a[]
:new_start_sym
```
"""
mutable struct Actor{T} <: TransactionalRef{T}
    mutex::ReentrantLock
    queue::ActorMsgQueue
    task::Union{Nothing,Task}
    value::ActorValue{T}
    postprocess::Union{Nothing,Function}
    errorhandler::Union{Nothing,Function}
    function Actor{T}(
        initval;
        postprocess::Union{Nothing,Function} = nothing,
        errorhandler::Union{Nothing,Function} = nothing,
    ) where {T}
        return new{T}(
            ReentrantLock(),
            ActorMsgQueue(),
            nothing,
            Some{T}(initval),
            postprocess,
            errorhandler,
        )
    end
end
export Actor
function Actor(t::T; postprocess = nothing, errorhandler = nothing) where {T}
    return Actor{T}(t; postprocess=postprocess, errorhandler=errorhandler)
end
Base.propertynames(::Actor) = (:state, :value, :error)
function Base.getproperty(a::Actor{T}, k::Symbol) where {T}
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
function Base.show(io::IO, ::MIME"text/plain", a::Actor{T}) where {T}
    print(io, "Actor{$T}(@")
    print(io, string(objectid(a); base=62))
    print(io, ": ")
    s = getfield(a, :value)
    if isa(s, ActorException{T})
        print(io, "--error--")
    else
        show(io, something(s))
    end
    return print(io, ")")
end
Base.show(io::IO, a::Actor{T}) where {T} = begin
    print(io, "Actor{$T}(@")
    print(io, string(objectid(a); base=62))
    print(io, ")")
end
function Base.setindex!(a::Actor{T}, args...) where {T}
    return error("Actor objects cannot be assigned---they can only be sent functions")
end
# As with `Delay`: an `Actor` is a `Ref`, so Base's 0-dimensional `a[i] = x`
# form (with a `CartesianIndex{0}`) is otherwise ambiguous with the method above.
Base.setindex!(a::Actor{T}, x, ::CartesianIndex{0}) where {T} = setindex!(a, x)
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
function actor_main(a::Actor{T}) where {T}
    mux = getfield(a, :mutex)
    err_args = nothing
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
            # Snapshot the hooks before releasing the lock, so that a hook set
            # while this message is in flight affects the next message and not
            # this one.
            pp = getfield(a, :postprocess)
            eh = getfield(a, :errorhandler)
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
            rval = nothing
            err = nothing
            # The three user functions are called with `invokelatest` rather
            # than directly. An actor's task outlives the world it was created
            # in: a message sent while the queue is still being drained may be a
            # closure defined afterwards---at the REPL, for instance, or in any
            # top-level code---and a direct call would raise a `MethodError`
            # ("the applicable method may be too new") inside the actor. The
            # cost is small here because the call is dynamic either way, the
            # message type being unknown: measured on an M1 Max, +2 ns on an
            # actor message that already costs ~165 ns, and no extra allocation.
            try
                rval = Base.invokelatest(f, val)
                (pp === nothing) || (rval = Base.invokelatest(pp, val, rval))
            catch e
                err = e
            end
            # Relock the condition to update the value.
            lock(mux)
            if err === nothing
                try
                    value = Some{T}(rval)
                    setfield!(a, :value, value)
                catch e
                    # A value that cannot be converted to T is an actor error,
                    # just as an exception from the scheduled function is.
                    err = e
                end
            end
            if err !== nothing
                # There was an exception, so set the error state and stop
                # processing. The error handler runs after the actor's lock is
                # released and its task is cleared (see below).
                setfield!(a, :value, ActorException{T}(err, val, f, (queue...,)))
                empty!(queue)
                err_args = (eh, val, err)
                break
            end
            # At this point we have updated the actor's value and we still
            # hold the lock. We can go ahead and enter the next loop.
        end
        # If we reach this point, then the queue is empty, so we can exit.
    finally
        setfield!(a, :task, nothing)
        unlock(mux)
    end
    if err_args !== nothing
        eh, val, err = err_args
        # The handler runs with no lock held and with the actor's task already
        # cleared, so that it may `reset` the actor and then `send` to it and
        # have the new message processed. An exception from the handler itself
        # is not caught; it ends this task. `invokelatest` for the same reason
        # as the message and the post-processor above.
        (eh === nothing) || Base.invokelatest(eh, a, val, err)
    end
    return nothing
end
"""
    actor_start!(actor)

If the given `Actor` object is not in an error state, is not running, and does
not have an empty queue, then `actor_start!(actor)` will start `actor`'s
processing task. This should only be called by `send` when the first message to
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
function actor_start!(a::Actor{T}) where {T}
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
function actor_send!(a::Actor{T}, msg::ActorMsg) where {T}
    val = getfield(a, :value)
    isa(val, ActorException{T}) && throw(val)
    DataStructures.enqueue!(getfield(a, :queue), msg)
    (getfield(a, :task) === nothing) && actor_start!(a)
    return a
end
function actor_send!(a::Actor{T}, qq::ActorMsgQueue) where {T}
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
function actor_reset!(a::Actor{T}, val) where {T}
    val0 = getfield(a, :value)
    isa(val0, ActorException{T}) || return nothing
    setfield!(a, :value, Some{T}(val))
    empty!(getfield(a, :queue))
    # Since the queue is, as of yet, empty, there's no need to start a task.
    return val0
end

# #Volatile ####################################################################
# Volatile's, like Actor's, have a single state structure.
#
# `VolatileData` is deliberately `mutable`, and that is what makes it a heap
# reference rather than a value that has to be materialised. Reading a volatile's
# value through an unparameterised ref — which is what the commit's validation,
# `volatile_setindex!`'s comparisons and a read inside a transaction all do —
# otherwise boxes it: measured, `getfield(v, :value)` with `v::Volatile` costs
# 20.5 ns and one allocation, against 2.8 ns and none when the type is known.
# Nothing needs to be materialised for a mutable object, so every such read is
# allocation-free; the cost is one allocation wherever one is constructed, and
# reads outnumber constructions in any transaction. A write's allocations went
# 12 -> 5 and its time 849 -> 632 ns.
#
# What keeps it sound is that a write *constructs* a new `VolatileData` rather
# than changing one in place, so the `===` comparisons sprinkled through the
# staged-write logic still mean "the same change" under reference equality. Do
# not turn it back into a `struct` without re-measuring, and do not start mutating
# one in place: the concurrency tests in `test/tx_concurrency.jl` exist to catch
# that class of mistake.
mutable struct VolatileData{T}
    value::T
    # These two are `Function`-typed, which makes a filter or finalizer call
    # dynamic. Parameterizing the type would make those calls static, but it
    # would ripple into `Volatile`, `Transaction`, and every dictionary of
    # volatiles, and it would only help when a filter or finalizer is actually
    # installed — which is not the default, and when neither is set these fields
    # are never called at all. Measured, a filtered write costs about 130 ns
    # more than a plain one, most of which is the filter's own work rather than
    # the dispatch. Left as is deliberately.
    filter::Union{Nothing,Function}
    finalize::Union{Nothing,Function}
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
    # Store `v` verbatim: unlike the constructors above, the filter is *not*
    # applied. This is needed wherever the documentation requires that a value
    # not be refiltered, namely `setfilter!`, `setfinalize!`, and
    # `voldata_finalize` (the finalize/filter functions are not rerun when the
    # filter is changed, nor when the finalizer returns a value).
    function VolatileData{T}(
        v::S, flt::Union{Nothing,Function}, fin::Union{Nothing,Function}, ::Val{:raw}
    ) where {T,S}
        return new{T}(v, flt, fin)
    end
end
"""
    voldata_finalize(voldata)

Calls the finalize function of the given volatile state object (of type
`VolatileData`) if necessary. Yields the new `VolatileData` object, which may be
`voldata` if there is no finalizer or nothing has changed.

This function is a private/internal implementation detail of the `Air` library.
"""
function voldata_finalize(v::VolatileData{T}) where {T}
    (v.finalize === nothing) && return v
    u = v.finalize(v.value)
    return (u === v.value ? v : VolatileData{T}(u, v.filter, v.finalize, Val{:raw}()))
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
function Volatile(
    t::T, filterfn::Union{Nothing,Function}, finfn::Union{Nothing,Function}=nothing
) where {T}
    return Volatile{T}(t, filterfn, finfn)
end
function Base.show(io::IO, ::MIME"text/plain", v::Volatile)
    show(io, typeof(v))
    print(io, "(@")
    print(io, string(objectid(v); base=62))
    print(io, ": ")
    show(io, v[])
    return print(io, ")")
end
Base.show(io::IO, v::Volatile) = begin
    show(io, typeof(v))
    print(io, "(@")
    print(io, string(objectid(v); base=62))
    print(io, ")")
end
Base.propertynames(::Volatile) = (:value,)
Base.getproperty(v::Volatile, p::Symbol) = begin
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
    reset_value::Union{Nothing,Some{T}}
    msgs::ActorMsgQueue
    function ActorTxData{T}(w::S) where {T,S}
        return new{T}(w, w, nothing, ActorMsgQueue())
    end
end
"""
    PromiseTxData{T}

Holds a settlement staged for a `Promise` by the current transaction: a
delivery in `val` or a failure in `err`, exactly one of which is non-`nothing`.
Both are stored as `Some` wrappers, already converted, so the commit can assign
the staged slot directly.

This struct is part of the internal/private code of Air and generally should not
be used outside of the Air library.
"""
struct PromiseTxData{T}
    val::Union{Nothing,Some{T}}
    err::Union{Nothing,Some{Any}}
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

Transactions have the following properties:
* `state` is one of `:running`, `:finalizing`, `:locking`, or `:committing`;
* `rvolatiles` is the set of all `Volatile` objects that have been read during
  the transaction;
* `wvolatiles` is the set of all `Volatile` objects that have been changed; and
* `actors` is the set of all `Actor` objects to which messages have been sent
  during the transaction.
"""
mutable struct Transaction
    state::Symbol
    # `Dict`s rather than `IdDict`s, even though the keys are compared by
    # identity. A `Dict` gives identity for these keys anyway, because neither
    # `Volatile` nor `Actor` defines `hash` or `isequal`, so both fall back to
    # Base's identity versions — verified by checking that two equal-valued but
    # distinct volatiles hash differently and occupy two entries in a `Dict`, as
    # they do in an `IdDict`. What matters is that `Base.IdDict` keeps its
    # *values* in a `Vector{Any}`, which is where a transaction's 2 kB of
    # `Memory{Any}` came from; a `Dict`'s keys and values are both typed.
    reads::Dict{Volatile,VolatileData}
    writes::Dict{Volatile,NTuple{2,VolatileData}}
    actors::Dict{Actor,ActorTxData}
    # An `IdDict`, unlike the other three. `Promise` now hashes and compares by
    # identity, so a `Dict` would work here too; `IdDict` says what is meant —
    # promises are handles and are compared by identity — without depending on
    # how `Promise` defines its hash.
    promises::IdDict{Promise,PromiseTxData}
    function Transaction()
        return new(
            :running,
            Dict{Volatile,VolatileData}(),
            Dict{Volatile,NTuple{2,VolatileData}}(),
            Dict{Actor,ActorTxData}(),
            IdDict{Promise,PromiseTxData}(),
        )
    end
end
export Transaction
Base.propertynames(::Transaction) = (:state, :rvolatiles, :wvolatiles, :actors)
Base.getproperty(t::Transaction, p::Symbol) = begin
    if p == :state
        return getfield(t, :state)
    elseif p == :rvolatiles
        return Base.IdSet{Volatile}(keys(getfield(t, :reads)))
    elseif p == :wvolatiles
        return Base.IdSet{Volatile}(keys(getfield(t, :writes)))
    elseif p == :actors
        return Base.IdSet{Actor}(keys(getfield(t, :actors)))
    else
        error("type $(typeof(t)) has no field $p")
    end
end
tx_clear!(t::Transaction) = begin
    setfield!(t, :state, :running)
    empty!(getfield(t, :reads))
    empty!(getfield(t, :writes))
    empty!(getfield(t, :actors))
    empty!(getfield(t, :promises))
    return nothing
end

# The current running transaction; see currtx() below. This lives in a task-local
# slot of its own rather than in a `Var`, because a `Var`'s value is stored in
# the task's bindings dictionary, and the machinery that sets one — `withvars` —
# copies that whole dictionary twice and looks it up three times to bind a single
# value. On this path that came to seven allocations and about 130 ns, roughly
# half of what an empty transaction costs. The state is task-scoped either way, so
# a plain slot has the same meaning and none of the copying.
const _CURRTX_KEY = :AirCurrentTx

"""
    _currtx_slot()

Yields the calling task's slot holding its running `Transaction`; see `currtx`.
This is part of `Air`'s internal/private implementation details.
"""
function _currtx_slot()
    d = task_local_storage()
    s = get(d, _CURRTX_KEY, nothing)
    (s isa Ref{Union{Nothing,Transaction}}) && return s
    s = Ref{Union{Nothing,Transaction}}(nothing)
    task_local_storage(_CURRTX_KEY, s)
    return s
end

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
currtx() = _currtx_slot()[]
export currtx

# The transaction in which a `send` should be staged: the running transaction
# if there is one, and otherwise the transaction whose finalize phase this task
# is in. `currtx()` is deliberately `nothing` during finalize — callbacks must
# see committed volatile values — but a `send` made there must not escape the
# transaction, or a transaction that fails validation would still have sent its
# messages, once per attempt.
function _staging_tx()
    t = currtx()
    (t === nothing) || return t
    return _finalizing_tx()
end

"""
    TX_MAX_ATTEMPTS

The maximum number of retry attempts that a transaction will make before
aborting the transaction.
"""
# The limit is generous, and the retry policy around it has been measured rather
# than assumed. With eight tasks contending for a single volatile, a transaction
# takes 2.85 attempts on average, and a retry is *cheaper* than a first attempt:
# a failed attempt exits at validation, before the commit work, so the observed
# per-operation cost is 2827 ns against 1850 ns for one task, while attempts rise
# from 1.0 to 2.85. That is a 1.5x degradation with throughput still scaling 5.2x
# on eight threads (about 65% of ideal) — the normal price of optimistic
# concurrency on one hot ref.
#
# Two things this settles, both by measurement rather than argument. Lock
# contention is not the limit: a bare `ReentrantLock` under the same contention
# costs 33 ns per operation, flat. And a backoff would buy nothing: it would trade
# lock re-acquisition for sleeping, while the retries it removed are already the
# cheap kind.
#
# See `test/tx_concurrency.jl` for the tests this rests on.
const TX_MAX_ATTEMPTS = 2^14

# Task-local scratch buffers for the commit's lock ordering; see the note in
# `tx`. One per task, grown as needed. A transaction's commit runs once and never
# re-enters, so a given buffer is never in use twice at once — including across
# the retry loop, whose attempts are sequential.
const _VOLS_KEY = :AirLockVols
const _ACTS_KEY = :AirLockActs
const _PRMS_KEY = :AirLockPrms
const _TXPOOL_KEY = :AirTxPool
# Holds the transaction whose finalize phase is running in this task, or
# `nothing`. The task's transaction slot has already been cleared by then, so
# this serves two purposes: a nested `tx` started from a finalizer or a delivery
# callback must not take the pooled `Transaction` (it would be the very object
# the outer transaction is still using, and `tx_clear!` would wipe it
# mid-commit), and a `send` made from one of those callbacks is staged in the
# outer transaction rather than escaping it (see `_staging_tx`).
const _TXFINAL_KEY = :AirTxFinalizing
# The transaction in this task's finalize phase, if any.
_finalizing_tx() = begin
    t = get(task_local_storage(), _TXFINAL_KEY, nothing)
    return t isa Transaction ? t : nothing
end

# The transaction to reuse, held per task, so that `tx` does not construct a
# `Transaction` and its four dictionaries on every call.
#
# Reuse is safe because `tx_clear!` resets every field of a `Transaction` — the
# state and all four dictionaries — and `tx` calls it on the way out of every
# attempt: in the `catch` when the body raises, and in the `finally` when the
# attempt reaches finalize and commit. So the object is clean when an attempt
# begins, the retry path restages into empty dictionaries, and the pooled
# object holds nothing — no volatiles, actors or promises, and through them no
# staged values — between transactions. It is in fact better than safe:
# `empty!` keeps a dictionary's capacity, so a reused transaction does not
# re-grow its tables each time.
#
# A second `tx` on the same task cannot see the pooled object mid-flight: while a
# transaction is running, the task's slot holds it, and `tx` returns early
# through `fn()` rather than starting another one.
function _transaction()
    d = task_local_storage()
    t = get(d, _TXPOOL_KEY, nothing)
    (t isa Transaction) && return t
    t = Transaction()
    task_local_storage(_TXPOOL_KEY, t)
    return t
end
function _lockbuf(key::Symbol, ::Type{T}) where {T}
    d = task_local_storage()
    buf = get(d, key, nothing)
    (buf isa Vector{T}) && return buf
    buf = Vector{T}()
    task_local_storage(key, buf)
    return buf
end

"""
    tx(fn)

Runs the given function in a transaction and yields the result. The function
fn is called as `fn()` without arguments.
"""
function tx(fn::F) where {F<:Function}
    success = false
    # NOTE: there is deliberately no `res = nothing` here; the transaction's
    # result is bound inside the retry loop, for the reason given there.
    slot = _currtx_slot()
    # If there's already a transaction running, we needn't make a new one---
    # This transaction function will just get rolled up into the current one.
    (slot[] === nothing) || return fn()
    # A transaction started from a finalizer or a delivery callback runs while
    # another transaction on this task is finalizing; it must not take the
    # pooled object, which the outer transaction is still using.
    the_tx = _finalizing_tx() === nothing ? _transaction() : Transaction()
    # We need to make a new transaction for this task. We might have to do
    # this a few times if the transaction fails.
    for attempt in 1:TX_MAX_ATTEMPTS
        # A failed attempt used to retry at once, which is counterproductive under
        # contention: the retrying transaction immediately re-locks the volatiles
        # that just invalidated it, so it is likely to lose again. Yielding first
        # hands the thread to whoever holds them. Measured with several tasks
        # writing one volatile, attempts per commit fell from 3.19 to 1.96 with
        # eight tasks and from 3.05 to 2.36 with eight tasks over two volatiles,
        # and latency came down 10-20%; yielding on every retry beat yielding only
        # from the second.
        #
        # A `yield` rather than a sleep, because a transaction is short: this is a
        # scheduler hand-off, and `sleep`'s resolution is a millisecond, three
        # orders of magnitude coarser than a transaction. And it costs nothing in
        # the shape this library is built for — many volatiles, short transactions
        # touching few of them — where conflicts are rare and the retry path is
        # barely reached: at a thousand volatiles with three per transaction,
        # attempts per commit measure 1.037 both with and without it.
        (attempt > 1) && yield()
        # The transaction is already clean: every attempt clears it on the way
        # out, by the two paths below. Clearing here instead would cost a second
        # pass over the four dictionaries on every transaction.
        #
        # Bind the result with the `try` expression rather than assigning to a
        # `res` that was initialised to `nothing` before the loop: the latter
        # makes the inferred return type `Union{Nothing,T}`, and a transaction's
        # value is read on every `tx` call. The `catch` arm never falls through,
        # so the value of this expression is exactly the value of `fn`.
        res = try
            # Bind the transaction for the body by writing the task's slot
            # directly. The `withvars` binding it replaces did the same job, but
            # through the task's bindings dictionary, which it copies twice and
            # looks up three times to bind one value.
            slot[] = the_tx
            try
                fn()
            finally
                slot[] = nothing
            end
        catch e
            # Clear on the way out of the attempt: the transaction is pooled per
            # task, so anything left staged would stay reachable until this task
            # ran another transaction. A retry restages its own writes into the
            # emptied dictionaries; the finalize and commit region below clears
            # on its own exit.
            tx_clear!(the_tx)
            isa(e, TxRetryException) && continue
            rethrow()
        end
        # At this point the transaction is over and we need to ensure that the
        # current volatile values are consistent with the transaction we just
        # executed. Before we do that, we should finalize the volatiles.
        #
        # Everything below runs with `_TXFINAL_KEY` bound to this transaction:
        # the task's transaction slot has already been cleared, so a nested `tx`
        # started by a finalizer or a delivery callback must get a fresh
        # `Transaction` rather than the pooled one this transaction is still
        # using, and a `send` made by one of those callbacks is staged in this
        # transaction (see `_staging_tx`) rather than escaping it. The binding
        # is saved and restored rather than simply cleared, because a nested
        # transaction's own finalize would otherwise clear it.
        prev_final = get(task_local_storage(), _TXFINAL_KEY, nothing)
        task_local_storage(_TXFINAL_KEY, the_tx)
        try
            writes = getfield(the_tx, :writes)
            setfield!(the_tx, :state, :finalizing)
            for (v, (x0, x1)) in writes
                x = voldata_finalize(x1)
                if x !== x1
                    setfield!(the_tx, :state, :running)
                    volatile_setindex!(v, x, the_tx)
                    setfield!(the_tx, :state, :finalizing)
                end
            end
            # Settlement callbacks run after the volatile finalizers, so an
            # attempt that a finalizer aborts does not notify anyone, and
            # before any commit lock is taken, so a callback may freely take
            # other locks. They run on every attempt: a callback may therefore
            # run for a settlement that never commits, and the order among
            # several staged promises is unspecified. They are called with
            # `invokelatest`, as when a promise is settled outside a
            # transaction: a long-lived task may stage a promise whose callback
            # was written after that task was created.
            promises = getfield(the_tx, :promises)
            for (p, dat) in promises
                pd = dat::PromiseTxData
                if pd.val === nothing
                    cb = getfield(p, :onfailure)
                    (cb === nothing) || Base.invokelatest(cb, (pd.err::Some{Any}).value)
                else
                    cb = getfield(p, :ondelivery)
                    (cb === nothing) || Base.invokelatest(cb, (pd.val::Some).value)
                end
            end
            # We're now ready to attempt a commit. To do this we must lock the
            # vols, actors and promises in a sorted order. We do volatiles
            # first, since their values must all be tested, then actors, then
            # promises. Within each group the locks are taken in `objectid`
            # order; every transaction uses the same group order, so the two
            # together are a global lock order.
            reads = getfield(the_tx, :reads)
            actors = getfield(the_tx, :actors)
            setfield!(the_tx, :state, :locking)
            nw = length(writes)
            nr = length(reads)
            n = nw + nr
            m = length(actors)
            np = length(promises)
            # The lock-order buffers are task-local scratch space rather than fresh
            # vectors. The obvious version — three arrays filled from the transaction's
            # key sets — is several allocations on every commit. The commit runs once per
            # transaction and never re-enters, so a per-task buffer is safe; it does not
            # even need clearing between uses, being `resize!`d and filled from the
            # start. A delivery callback may start a nested transaction, so these are
            # populated only after the callbacks have run.
            vols = _lockbuf(_VOLS_KEY, Volatile)
            resize!(vols, n)
            i = 0
            for (v, _) in reads
                vols[i += 1] = v
            end
            for (v, _) in writes
                vols[i += 1] = v
            end
            sort!(vols; by=objectid)
            acts = _lockbuf(_ACTS_KEY, Actor)
            resize!(acts, m)
            i = 0
            for (a, _) in actors
                acts[i += 1] = a
            end
            sort!(acts; by=objectid)
            prms = _lockbuf(_PRMS_KEY, Promise)
            resize!(prms, np)
            i = 0
            for (p, _) in promises
                prms[i += 1] = p
            end
            sort!(prms; by=objectid)
            locked_vol = 0
            locked_act = 0
            locked_prm = 0
            success = false
            try
                while !success
                    checked = 0
                    # First, go through and lock all the volatiles in order.
                    for vol in vols
                        lock(getfield(vol, :mutex))
                        locked_vol += 1
                    end
                    # Now check their values. First do the read volatiles.
                    for (vol, x0) in reads
                        (x0 === getfield(vol, :value)) || break
                        checked += 1
                    end
                    (checked == nr) || break
                    # Next, check the written volatiles.
                    for (vol, (x0, x)) in writes
                        (x0 === getfield(vol, :value)) || break
                        checked += 1
                    end
                    (checked == n) || break
                    # Next we need to lock and check the actors that were read from.
                    for act in acts
                        lock(getfield(act, :mutex))
                        locked_act += 1
                    end
                    for (act, dat) in actors
                        val = getfield(act, :value)
                        (dat.start_value === val) || break
                        checked += 1
                    end
                    (checked == n + m) || break
                    # Next lock the promises. A promise is one-shot, so a
                    # promise that is no longer pending can never be delivered
                    # by this transaction. That is not a retryable conflict:
                    # throwing aborts the transaction without a retry, and this
                    # happens before any volatile or actor is committed, so the
                    # losing transaction leaves no partial writes.
                    for prm in prms
                        lock(getfield(prm, :mutex))
                        locked_prm += 1
                    end
                    for (prm, dat) in promises
                        pd = dat::PromiseTxData
                        _promise_pending(prm) ||
                            throw(ArgumentError(getfield(pd, :val) === nothing ?
                                                "given promise is already settled" :
                                                "given promise is already fulfilled"))
                    end
                    # If we reach this point, then everything checks out and we can
                    # commit all of the changes and complete the transaction.
                    setfield!(the_tx, :state, :committing)
                    for (vol, (x0, x)) in writes
                        setfield!(vol, :value, x)
                    end
                    # We need to commit the actors as well.
                    for (act, dat) in actors
                        q = dat.msgs
                        r = dat.reset_value
                        (r === nothing) || actor_reset!(act, r.value)
                        isempty(q) || actor_send!(act, q)
                    end
                    # And the promises, whose locks we still hold.
                    for prm in prms
                        dat = get(promises, prm, nothing)
                        if dat !== nothing
                            pd = dat::PromiseTxData
                            if pd.val === nothing
                                setfield!(prm, :_err, (pd.err::Some{Any}))
                            else
                                setfield!(prm, :_val, (pd.val::Some))
                            end
                            notify(getfield(prm, :cond); all=true)
                            _wake_waiters!(prm)
                        end
                    end
                    # That's all we need to do, aside from unlock.
                    success = true
                end
            finally
                for j in locked_prm:-1:1
                    unlock(getfield(prms[j], :mutex))
                end
                for j in locked_act:-1:1
                    unlock(getfield(acts[j], :mutex))
                end
                for j in locked_vol:-1:1
                    unlock(getfield(vols[j], :mutex))
                end
            end
        finally
            task_local_storage(_TXFINAL_KEY, prev_final)
            # Clear on the way out, for the reason given in the `catch` above:
            # this runs for a commit, an aborted attempt and a throwing
            # finalizer or delivery callback alike. It is also what readies the
            # pooled object for another attempt after a failed validation.
            tx_clear!(the_tx)
        end
        success && return res
    end
    # We only reach this point by exhausting the retries: a successful attempt
    # returns from inside the loop.
    error("transaction aborted after failing $TX_MAX_ATTEMPTS times")
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
function _volatile_getindex(v::Volatile{T}, t::Transaction) where {T}
    # The transaction's dictionaries are keyed by `Volatile` and hold untyped
    # `VolatileData`, so without these assertions the value read out of them is
    # inferred as `Any` and every read inside a transaction is boxed. The data
    # stored for a `Volatile{T}` is always a `VolatileData{T}`.
    w = get(getfield(t, :writes), v, nothing)
    (w === nothing) || return w[2]::VolatileData{T}
    w = get(getfield(t, :reads), v, nothing)
    (w === nothing) || return w::VolatileData{T}
    w = getfield(v, :value)
    getfield(t, :reads)[v] = w
    return w
end
_volatile_getindex(v::Volatile{T}, ::Nothing) where {T} = getfield(v, :value)
Base.getindex(v::Volatile{T}) where {T} = _volatile_getindex(v, currtx()).value
function volatile_setindex!(v::Volatile{T}, x::VolatileData{T}, t::Transaction) where {T}
    (getfield(t, :state) === :running) ||
        error("volatiles can only be set before finalizing the transaction")
    # See if it this ref has already been written to in this transaction.
    w = get(getfield(t, :writes), v, nothing)
    if w !== nothing
        (x0, x1) = w
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
function volatile_setindex!(v::Volatile{T}, x::S, ::Nothing) where {T,S}
    return error("cannot set volatile outside of transaction")
end
function Base.setindex!(v::Volatile{T}, x::S) where {T,S}
    value = getfield(v, :value)
    # The `VolatileData` constructor applies the filter, so it must not be
    # applied here as well.
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
# The `VolatileData` associated with `v` from the transaction's point of view:
# the staged write if `v` has already been written in this transaction, and the
# committed value otherwise. `setfilter!`/`setfinalize!` build on this so that
# setting both functions in a single transaction does not discard one of them.
_volatile_txdata(v::Volatile{T}, t::Transaction) where {T} = begin
    w = get(getfield(t, :writes), v, nothing)
    return w === nothing ? getfield(v, :value) : w[2]
end
"""
    setfilter!(vol, fn)

Sets the filter-function associated with the Volatile object vol. Any time that
the vol is set (`vol[] = x`) the filter-function is called and the value saved
in `vol` is instead `fn(x)`. This must be called within a transaction.
"""
function setfilter!(vol::Volatile{T}, f::Function) where {T}
    t = currtx()
    value = t === nothing ? getfield(vol, :value) : _volatile_txdata(vol, t)
    newdat = VolatileData{T}(value.value, f, value.finalize, Val{:raw}())
    return volatile_setindex!(vol, newdat, t)
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
function setfinalize!(vol::Volatile{T}, f::Function) where {T}
    t = currtx()
    value = t === nothing ? getfield(vol, :value) : _volatile_txdata(vol, t)
    newdat = VolatileData{T}(value.value, value.filter, f, Val{:raw}())
    return volatile_setindex!(vol, newdat, t)
end
export setfinalize!

# Same for actors: we can put the actor code here.
# First, we definne some private functions for handling things on the back
# end.
"""
    tx_actordata(a)

Yields the in-transaction data for the `Actor` object `a`, or `nothing` if no
transaction is running. An actor in an error state yields data whose `tx_value`
is the `ActorException`; it is `getindex` that raises on it, not this function.

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
    # The assertion is what keeps this function's result inferable: the actors
    # dictionary is keyed by `Actor` and valued by unparameterised `ActorTxData`
    # (a transaction can touch differently-typed actors), so without it the result
    # is `ActorTxData` with `T` free, and everything downstream — including a
    # plain `a[]` — infers as `Any`. Every value stored for this actor is an
    # `ActorTxData{T}`, because that is the only thing ever inserted.
    return w::ActorTxData{T}
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
    # Narrowing *to* `Some{T}` rather than testing for the exception: reading an
    # `Actor` in an error state raises, so the exception branch never returns, and
    # a positive `isa` is what leaves `x`'s type as `Some{T}` — from which `.value`
    # is `T`. Testing `isa(x, ActorException{T})` and falling through instead
    # leaves the field read on the whole two-member `Union`, which infers as `Any`.
    isa(x, Some{T}) || throw(x)
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
function actor_reset(a::Actor{T}, newval::S, t::Transaction) where {T,S}
    (getfield(t, :state) === :running) ||
        error("actors can only be reset before finalizing the transaction")
    w = tx_actordata(a, t)
    val = w.tx_value
    isa(val, ActorException{T}) || return nothing
    newval = Some{T}(newval)
    w.reset_value = newval
    w.tx_value = newval
    return val
end
function actor_reset(a::Actor{T}, newval::S, ::Nothing) where {T,S}
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
Base.reset(a::Actor{T}, s::S) where {T,S} = actor_reset(a, s, currtx())
export reset

"""
    getpostprocess(actor)

Yields the post-processing function associated with the given `Actor`, or
`nothing` if it has none. See [`Actor`](@ref) for what a post-processing
function does and when it runs.

See also: [`setpostprocess!`](@ref), [`geterrorhandler`](@ref)
"""
function getpostprocess(a::Actor)
    mux = getfield(a, :mutex)
    lock(mux)
    try
        return getfield(a, :postprocess)
    finally
        unlock(mux)
    end
end
export getpostprocess

"""
    setpostprocess!(actor, fn)

Sets the post-processing function of the given `Actor` to `fn`, or clears it
when `fn` is `nothing`, and yields the actor. The change takes effect for the
next function the actor dequeues; a function that is already running is
unaffected. This is actor configuration rather than actor state: the change is
immediate and is not rolled back by a transaction.

See also: [`getpostprocess`](@ref), [`seterrorhandler!`](@ref)
"""
function setpostprocess!(a::Actor, f::Union{Nothing,Function})
    mux = getfield(a, :mutex)
    lock(mux)
    try
        setfield!(a, :postprocess, f)
    finally
        unlock(mux)
    end
    return a
end
export setpostprocess!

"""
    geterrorhandler(actor)

Yields the error-handler function associated with the given `Actor`, or
`nothing` if it has none. See [`Actor`](@ref) for what an error handler does
and when it runs.

See also: [`seterrorhandler!`](@ref), [`getpostprocess`](@ref)
"""
function geterrorhandler(a::Actor)
    mux = getfield(a, :mutex)
    lock(mux)
    try
        return getfield(a, :errorhandler)
    finally
        unlock(mux)
    end
end
export geterrorhandler

"""
    seterrorhandler!(actor, fn)

Sets the error-handler function of the given `Actor` to `fn`, or clears it when
`fn` is `nothing`, and yields the actor. The change takes effect for the next
function the actor dequeues; a function that is already running is unaffected.
This is actor configuration rather than actor state: the change is immediate and
is not rolled back by a transaction.

See also: [`geterrorhandler`](@ref), [`setpostprocess!`](@ref)
"""
function seterrorhandler!(a::Actor, f::Union{Nothing,Function})
    mux = getfield(a, :mutex)
    lock(mux)
    try
        setfield!(a, :errorhandler, f)
    finally
        unlock(mux)
    end
    return a
end
export seterrorhandler!

"""
    actor_send(a, msg)

Adds the given `ActorMsg` message object `msg` to the queue of the given `Actor`
object `a` such that it will be processed in a separate thread. Yields
`nothing`.

This function is considered part of the internal/private interface of `Air` and
shouldn't generally be called outside of it.
"""
actor_send(a::Actor{T}, msg::ActorMsg) where {T} = actor_send(a, msg, _staging_tx())
function actor_send(a::Actor{T}, msg::ActorMsg, t::Transaction) where {T}
    # `:finalizing` is allowed: a send from a volatile finalizer or a promise
    # delivery callback is staged in the transaction that is finalizing, so it
    # is committed if and only if that transaction commits. Later phases are
    # not: the commit's lock set is already built, so a send there would be
    # discarded silently.
    st = getfield(t, :state)
    ((st === :running) || (st === :finalizing)) ||
        error("actors can only be sent functions before the transaction locks")
    w = tx_actordata(a, t)
    isa(w.tx_value, ActorException{T}) && throw(w.tx_value)
    DataStructures.enqueue!(w.msgs, msg)
    return nothing
end
function actor_send(a::Actor{T}, msg::ActorMsg, ::Nothing) where {T}
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
thread. In the actual function call that is evaluated in this separate thread,
the function is called as `fn(actor_value)`, and its return value becomes the
new stored value of the actor---unless the actor has a post-processing
function, in which case the post-processor's return value is stored instead
(see [`Actor`](@ref)).

The function is called in the latest world, so a function that was defined
after the actor's processing thread was started---a closure entered at the REPL,
for instance---can be sent to it; without that, the actor would raise a
`MethodError` saying the applicable method is too new.

Note that the value of an actor when `send()` is called does not make any
guarantee about its value when the function that is sent gets evaluated---other
functions may be processed in the interim.

If `send()` is called inside of a transaction, then the function is not
immediately queued but rather is held until the transaction successfully
completes. Once this happens, all sends are dispatched simultaneously such that
two sends to the same actor during the same transaction will always run
back-to-back in the actor's processing thread.

This includes a `send()` made from a volatile finalizer or a promise delivery
callback, which run in the transaction's finalize phase with no transaction
active (`currtx() === nothing`): the message is staged in the transaction that
is finalizing and queued only if it commits. A transaction that loses
validation, or that aborts, sends nothing.

**Important**. Due to the design pattern described in the previous paragraph,
one should never wait on the result of a function sent to an actor during a
transaction. This mistake might be made, for example, by sending an actor, during
a transaction, a function that delivers some result to a [`Promise`](@ref)
object. If the sending thread then waits on that promise while still in the same
transaction then it will deadlock because the sent function is not put in the
actor's queue until the enclosing transaction is finished. Thus the promise
will never receive the result from the actor.

TL;DR---When you touch an actor in a transaction, it's value freezes inside the
transaction. So if the transaction somehow waits for the actor to update or to
do something, that transaction will deadlock.
"""
send(f::Function, a::Actor{T}) where {T} = actor_send(a, ActorMsg(f), _staging_tx())
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
Base.write(io::TxIO, x::UInt8) = send(io.actor) do io
    write(io, x)
    return io
end
Base.write(io::TxIO, x::Char) = send(io.actor) do io
    write(io, x)
    return io
end
Base.write(io::TxIO, x::Union{String,SubString{String}}) = send(io.actor) do io
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
Base.write(io::AirOut, x::UInt8) =
    let so = stdout
        return send(io.actor) do io
            write(so, x)
            return io
        end
    end
Base.write(io::AirOut, x::Char) =
    let so = stdout
        return send(io.actor) do io
            write(so, x)
            return io
        end
    end
Base.write(io::AirOut, x::Union{String,SubString{String}}) =
    let so = stdout
        return send(io.actor) do io
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
