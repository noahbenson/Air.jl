################################################################################
# TX.jl
# Functional multi-threading tools built around atomic transactions.


import DataStructures

# #ReentrantRef ################################################################
"""
    ReentrantRef{T}

A ReentrantRef is a type of ref object that is thread-safe. There are a few
strategies for this, each of which are encoded in a different object type that
inheits from ReentrantRef. These types are as follows.

Var: Var objects act like Ref objects except that changes to them are
exclusively task-local.

Actor: Actor objects obey the actor pattern; you can send(fn, actor, args...)
where fn is a function that is, in another thread at some point, called as
fn(actor[], args...). The new value of the actor after running fn is the return
value of the call.

Volatile: Volatile objects are refs that can be changed by any thread but that
must be changed only within a synchronized transaction that ensures that all
reads and writes of volatiles, as well as reads fom and sends to actors, are
atomic: either they all happen successfully, or none of them do.
"""
abstract type ReentrantRef{T} <: Ref{T} end


# #TransactionalRef ############################################################
"""
    TransactionalRef{T}

A TransactionalRef is a reentrant ref that additional participates in
transactions. Transactional refs include Actors and Volatiles.
"""
abstract type TransactionalRef{T} <: ReentrantRef{T} end


# #Var #########################################################################
"""
    Var{T}

A Var object represents a task-local piece of data with a default value. You
can access a var with `var[]` and you can set it with `var[] = newval`. However,
the new assignment will always be task-local only. Because of this, Vars are
safe to access and update in a multi-threaded program.

All fields of a Var should be considered strictly private.
"""
mutable struct Var{T} <: ReentrantRef{T}
    mutex::ReentrantLock
    values::PTree{T}
    initval::T
    function Var{T}(initval::S) where {T, S}
        return new{T}(ReentrantLock(), PTree{T}(), initval)
    end
end
Var(initval::T) where {T} = Var{T}(initval)
Base.getindex(v::Var{T}) where {T} = get(getfield(v, :values),
                                         objectid(current_task()),
                                         getfield(v, :initval))
Base.propertynames(::Var) = (:initial_value,)
Base.getproperty(u::Var{T}, sym) where {T} = begin
    if sym == :initial_value
        return getfield(u, :initval)
    else
        throw(ArgumentError("Unknown Var property $sym"))
    end
end
# Helper function for cleaning up tasks from vars.
_var_task_cleanup(var::Var{T}, t::Task) where {T} = begin
    h = objectid(task)
    mux = getfield(var, :mutex)
    lock(mux)
    try
        vals0 = getfield(var, :values)
        vals = delete(vals0, h)
        (vals === vals0) || setfield!(var, :values, vals)
    finally
        unlock(mux)
    end
    return t
end
Base.setindex!(v::Var{T}, x::S) where {T,S} = begin
    task = current_task()
    h = objectid(task)
    mux = getfield(v, :mutex)
    lock(mux)
    try
        vals = getfield(v, :values)
        x0 = get(vals, h, vals)
        if x0 === vals
            (x === getfield(v, :initval)) && return x
            finalizer(t -> _var_task_cleanup(v, t), task)
            vals = setindex(vals, x, h)
        elseif x0 == x
            return x
        elseif x0 === getfield(v, :initval)
            vals = delete(vals, h)
        else
            vals = setindex(vals, x, h)
        end
        (getfield(v, :values) === vals) || setfield!(v, :values, vals)
    finally
        unlock(mux)
    end
    return x
end
Base.show(io::IO, ::MIME"text/plain", d::Var{T}) where {T} = begin
    print(io, "$(typeof(d))(@$(objectid(d)), init=$(getfield(d, :initval)))")
end
"""
    @var

Convenient syntax for creating a task-local Var object:
`@var name = initval` will construct a Var object with the given initial
value. `@var name = initval::T` will create a `Var{T}` object.
"""
macro var(expr::Expr)
    (expr.head === :(=)) || throw(
        ArgumentException("@var macro requires an assigmnet expression"))
    (name, initval) = expr.args
    isa(name, QuoteNode) && (name = name.value)
    isa(name, Symbol) || throw(
        ArgumentException("@var macro requires a symbol for the name"))
    s = gensym()
    if isa(initval, Expr) && initval.head === :(::)
        T = initval.args[2]
        initval = initval.args[1]
        q = :(const $name = Air.Var{$T}($initval))
    else
        q = :(const $name = Air.Var($initval))
    end
    return esc(q)
end


# #Actor #######################################################################
"""
    ActorException{T}(error, value, argno, args)

An ActorException object is thrown whenever one attempts to obtain the value
of or send a function to an actor that is in an error state. An error state
occurs when an unhandled exception is raised while an actor is processing a
sent function. In such a case the `actorerror()` and `restart()` functions may
be used. The `actorerror()` function yields an ActorException object in which
the exception that was raised is stored as `error`, the value of the actor when
the send was run is stored in `value`, and the `args` and `argno` give the
full contex of the send (i.e., the fucntion followed by the arguments in `args`
and the argument number that corresponds to the actor in `argno`.
"""
struct ActorException{T} <: Exception
    error::Any
    value::T
    argno::Int
    args::Tuple
end
Base.show(io::IO, ::MIME"text/plain", a::ActorException{T}) where {T} = begin
    print(io, "ActorException{$T}($(typeof(a.error)), $(value), ...)")
end

"""
    _ActorMsg{T}

The (private) type of a message that is queued for an actor to eventually run.
"""
struct _ActorMsg{T}
    argno::Int
    args::Tuple
end

"""
    Actor{T}

An actor is an object that represents an worker-thread to which tasks can be
scheduled. Any scheduled function is guaranteed to be evaluated at some point in
the future in some other thread, and the return value of that function will
become the new value held by the actor. Each function, when it is run, is passed
the actor's value as one of its argument.

Like with Refs, you can access an actor's current value using `actor[]`. Within
a transaction, this will be guaranteed to remain fixed for the duration of the
transaction; outside of a transaction, this may change at arbitrary times.

Actors may be initialized with post-processing functions. These functions are
called immediately after executing any scheduled function and are given as their
first argument the original agent value and as their second argument the return
value of the scheduled function. Instead of saving this return value, the
return value of the post-processing function is saved in the agent. This is
useful for agents whose job is, for example, to serialize output to a log file
being written to across many threads. If we run something like 
`send(log_agent, println, "Some log message")` the `println` function will
return nothing, which we would like to convert back into the log file so that
subsequent `send` calls can continue to use the `println` fuction.

Actors may additionally be given an error handler. Whenever an exception occurs
during a scheduled function, the error handler will be called with the arguments
of (1) the agent object, (2) the current agent value, and (3) the exception that
was caught. When an error has occurred, any attempt to schedule a function to
the agent or to read from the agent will cause an exception to be raised. The
error may be examined using `actorerror(actor)` and restarted using
`restart(actor)`.

All fields in an agent object should be considered strictly private. These
fields are likely to change between releases, and changing the values will break
your code.
"""
mutable struct Actor{T} <: TransactionalRef{T}
    _cond::Threads.Condition
    _value::Union{T, ActorException{T}}
    _task::Union{Task, Nothing}
    _msgs::DataStructures.Queue{_ActorMsg{T}}
end
Base.show(io::IO, ::MIME"text/plain", a::Actor{T}) where {T} = begin
    try
        val = a[]
        print(io, "Actor{$T}(@$(string(objectid(a), base=62)): $(val))")
    catch e
        print(io, "Actor{$T}(@$(string(objectid(a), base=62)): --error--)")
    end
end
Base.setindex!(a::Actor{T}, args...) where {T} = error(
    "Actor objects cannot be assigned--they can only be sent functions")
# The actor main function is what actually manages everything behind the scenes
# in the actor tasks. Note that this function and the task/thread it is running
# in is the only place that the _value of the actor a can be legally changed.
# Changing the _value anywhere else can result in undefined behavior. The same
# is true for the _task field.
_actor_main(a::Actor{T}) where {T} = begin
    q = a._msgs
    msg = nothing
    args = nothing
    rval = nothing
    try
        # For starters, lock the condition!
        lock(a._cond)
        # This function is a loop! We process everything in the queue until it
        # is empty, then we wait for a signal on the condition.
        while true
            try
                # See if there's anything waiting in the queue; if not we
                # should wait on a signal.
                isempty(q) && wait(a._cond)
                # (If the queue is empty now, then somethign has gone wrong.)
                msg = DataStructures.dequeue!(q)
            finally
                unlock(a._cond)
            end
            # Okay, we've grabbed the message we want to process and we've
            # unlocked so that we can process things. Do so:
            if msg.argno == 1
                f = a._value
                args = msg.args[2:end]
            else
                args = msg.args
                n = msg.argno
                f = args[1]
                args = (args[2:n-1]..., a._value, args[n+1:end]...)
            end
            try
                rval = f(args...)
            catch e
                e = ActorException{T}(e, a._value, msg.argno, msg.args)
                a._value = e
                break # break from the loop, which causes the task to exit
            end
            # Relock the condition to update the value.
            lock(a._cond)
            try
                a._value = rval
            catch e
                e = ActorException{T}(e, a._value, msg.argno, msg.args)
                a._value = e
                unlock(a._cond)
                break
            end
            # At this point we have updated the actor's value and we still
            # hold the lock. We cann go ahead and enter the next loop.
        end
    catch e
        # We've caught an exception. If it's not an actor exception we make it
        # into one...
        if !isa(e, ActorException{T})
            e = ActorException{T}(e, a._value, 0, ())
        end
        # We should set the value to this exception now.
        a._value = e
    end
    # This function/task is going to exit, so we can safely clear the task.
    a._task = nothing
    return nothing
end
# This function resets an actor.
# YOU MUST BE HOLDING THE ACTOR'S COND WHEN CALLING THIS.
_actor_reset(a::Actor{T}, u::S) where {T,S} = begin
    (a._task === nothing) || error("actor is already running")
    a._value = u
    task = Task(() -> _actor_main(a))
    a._task = task
    schedule(task)
    return u
end
Actor{T}(u::S) where {T,S} = begin
    cond = Threads.Condition()
    msgs = DataStructures.Queue{_ActorMsg{T}}()
    actr = Actor{T}(cond, u, nothing, msgs)
    _actor_reset(actr, actr._value)
    return actr
end


# #Volatile ####################################################################
struct VolatileValue{T}
    value::T
    filter::Union{Nothing, Function}
    finalize::Union{Nothing, Function}
    function VolatileValue{T}(v::S, flt::Function, fin::Function) where {T,S}
        return new{T}(flt(v), flt, fin)
    end
    function VolatileValue{T}(v::S, flt::Function, ::Nothing) where {T,S}
        return new{T}(flt(v), flt, nothing)
    end
    function VolatileValue{T}(v::S, ::Nothing, fin::Function) where {T,S}
        return new{T}(v, nothing, fin)
    end
    function VolatileValue{T}(v::S, ::Nothing, ::Nothing) where {T,S}
        return new{T}(v, nothing, nothing)
    end
    function VolatileValue{T}(v::S, flt::Function) where {T,S}
        return new{T}(flt(v), flt, nothing)
    end
    function VolatileValue{T}(v::S, ::Nothing) where {T,S}
        return new{T}(v, nothing, nothing)
    end
    function VolatileValue{T}(v::S) where {T,S}
        return new{T}(v, nothing, nothing)
    end
end
_volatile_finalize(v::VolatileValue{T}) where {T} = begin
    (v.finalize === nothing) && return v
    u = v.finalize(v.value)
    return (u === v.value ? v : VolatileValue{T}(u, v.filter, v.finalize))
end
"""
    Volatile{T}(value)
    Volatile{T}(value, filter_fn)
    volatile{T}(value, filter_fn, finalize_fn)

Volatile objects are Ref objects that must be used in conjunction with transaction
blocks (see `Air.@tx`). The value of a Volatile can be accessed at any time, and
there is no particular guarantee that a Volatile's value won't be changed by
another thread outside of a transaction. Volatile objects can only be set inside
of a transaction, however, and within a transaction, a Volatile's value is
guaranteed to be constant.

The arguments filter_fn and finnalize_fn are functions for making sure that the
value of a volatile conforms to some standard. The two functions are similar,
but are intended for slightly different use cases:
 * filter_fn is called every time the value of a volatile is set. The value
   saved in the volatile is filter_fn(value) instead of value itself.
 * finnalize_fn is called with the volatile's value immediately after the body
   of the transaction has completed, but before a commit is attempted. The
   finalize_fnn function acts like the filter_fn in that check_fn(value)
   replaces the value in the volatile.
Both filter_fn and finalize_fn can throw exceptions to abort the transaction.

For a volatile `v`, one may set `v`'s stored value via `v[] = value`. This must
be done inside a transaction (via `tx` or `@tx`). To change the filter or the
finalize functions, use the `setfilter!` and `setfinalize!` functions, both
of which also must be run in transactions as well.
"""
mutable struct Volatile{T} <: TransactionalRef{T}
    _mutex::ReentrantLock
    _value::VolatileValue{T}
    function Volatile{T}(val::S) where {T,S}
        return new{T}(ReentrantLock(), VolatileValue{T}(val, nothing, nothing))
    end
    function Volatile{T}(val::S, f::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileValue{T}(val, f, nothing))
    end
    function Volatile{T}(val::S, f::Function, ::Nothing) where {T,S}
        return new{T}(ReentrantLock(), VolatileValue{T}(val, f, nothing))
    end
    function Volatile{T}(val::S, f::Function, ff::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileValue{T}(val, f, ff))
    end
    function Volatile{T}(val::S, ::Nothing, ff::Function) where {T,S}
        return new{T}(ReentrantLock(), VolatileValue{T}(val, nothing, ff))
    end
end
Volatile(t::T) where {T} = Volatile{T}(t)
Volatile(t::T, filterfn::Function) where {T} = Volatile{T}(t, f)
Volatile(t::T, filterfn::Function, finalizefn::Function) where {T} = Volatile{T}(t, f, ff)


# #Source ######################################################################
"""
    AbtractSourceKernel{T}

An AbstractSourceKernel is an object that can be sampled to produce the next
sample for a source. Objects of types descending from AbstractSourceKernel are
not responsible for maintaining the multi-threaaded state exchange implied by
sources; rather the represent a function that will be called discretely when
samples are needed.

Source kernels use the pop! function to get values. Once popped, values cannot
be reobtained.
"""
abstract type AbstractSourceKernel{T} end
Base.pop!(k::AbstractSourceKernel{T}) where {T} = error(
    "type $(typeof(k)) <: AbstractSourceKernel{$T} must overloaded pop!")

"""
    FunctionSourceKernel{T}

A FunctionSourceKernel{T} object holds a function that is called every time a
new source sample is needed.
"""
struct FunctionSourceKernel{T} <: AbstractSourceKernel{T}
    _fn::Function
end
Base.pop!(f::FunctionSourceKernel{T}) where {T} = f()::T

"""
    Source{T,K}

Source objects work with transactions such that any samples taken from a source
during a transaction either entirely occur or do not occur. I.e., if a
transaction is rolled back or if an exception occurs during a transaction after
a sampling from a source, then the sampling is rolled back.

Sources object their samples from Source kernel objects. The type parameter K
must be the kernel type while T is the type parameter of type K (i.e., K is
something like FunctionalSourceKernel{T}).

Sources can be sampled via the get() function. Althrough get() can be called
either inside or outside of a transaction, when called inside of a transaction,
the get is only performed when the transaction succeeds. Sources are thus
safe to read from in transactions, unlike files and input streams.
"""
mutable struct Source{T,K <: AbstractSourceKernel{T}}
    _cond::Threads.Condition
    _buffer::Vector{T}
    _position::UInt128
    _state::Symbol
    _waiting::Union{Nothing, Task}
    _error::Any
    _kernel::K
end
Source{T,K}(k::K) where {T,K<:AbstractSourceKernel{T}} = begin
    cond  = Threads.Condition()
    mux   = Threads.ReentrantLock()
    buf   = Vector{T}()
    return Source{T,K}(cond, mux, buf, :ok, nothing, nothing, k)
end
Source{T}(f::Function) where {T} = begin
    k = FunctionSourceKernel{T}(f)
    return Source{T,FunctionSourceKernel{T}}(k)
end
Source(f::Function) = Source{Any}(f)
Source{K}(k::K) where {T, K <: AbstractSourceKernel{T}} = Source{T,K}(k)
Source(k::K) where {T, K <: AbstractSourceKernel{T}} = Source{T,K}(k)
Base.show(io::IO, ::MIME"text/plain", a::Source{T,K}) where {T,K} = begin
    print(io, "Source{$T}(@$(objectid(a)))")
end
# Sometimes we need to abort a source read...
"""
    SourceReadError{T,K}

Source read errors are generated when, during a transaction, is being read by
multiple threads at once and becomes invalid for the current trasnaction. This
is not really an error, and users should not encounter these exceptios; rather
they should result in transaction retries.
"""
struct SourceReadError{T,K} <: Exception
    source::Source{T,K}
    start_position::UInt128
end
# Sometimes something goes wrong and an exception has to be raised in the kernel
"""
    SourceAbortError{T,K}

A SourceAbortError is raised when a task is reading from a source (i.e., waiting
for a source-kernel to pop a value) and a reset is issued to the source. In
general this shouldn't happen because errors should only arise from the task
waiting on the kernel itself.
"""
struct SourceAbortError{T,K} <: Exception
    source::Source{T,K}
end
"""
    SourceKernelError{T,K}

A SourceKernelError is raised when an operation is attempted on a source whose
kernel has raised an exception. The error of a source s may also be extracted
via geterror(s).

For a SourceKernelError object err, err.source is the original source,
err.task is the task that was attempting to read the actor when the error was
caught, and err.error is the thrown object itself.
"""
struct SourceKernelError{T,K} <: Exception
    source::Source{T,K}
    task::Task
    error::Any
end
_source_claim_nolock(s::Source{T,K}, spos::UInt128, n::UInt128) where {T,K} = begin
    (s._state == :error) && throw(s._error)
    (s._position == start_pos) || throw(SourceReadError{T,K}(s, start_pos))
    s._position += n
    res = s._buffer[1:n]
    s._buffer = s._buffer[n+1:end]
    return res
end
"""
    _source_claim(s, start_pos, n)

Claims the first n elements after the given source position start_post from the
source s and updates the source object buffer to be free of those elements and
the source position to be position + n. If the source position has been updated
from start_pos, a SourceReadError exception is thrown.

If n values have not already been ensured, then it is an error to call this
function (i.e., it does not ensure the values, it just checks that the given
start_pos is still equal to the source position).
"""
_source_claim(s::Source{T,K}, spos::UInt128, n::UInt128) where {T,K} = lock(s._cond) do
    return _source_claim_nolock(s, spos, n)
end

# #Transaction #################################################################
const _ActorMsgQ{T} = DataStructures.Queue{_ActorMsg{T}} where {T}
mutable struct _ActorTxData{T}
    start_value::T
    tx_value::T
    msgs::_ActorMsgQ{T}
    function _ActorTxData{T}(w::S) where {T,S}
        return new{T}(w, w, _ActorMsgQ{T}())
    end
end
mutable struct _SourceTxData
    start_position::UInt128
    count::UInt128
end
"""
    Transaction

A transaction object keeps track of what is going on during a particular
transaction. These are generally low-level objects that shouldn't be
touched directly.
"""
mutable struct Transaction
    _state::Symbol
    _reads::IdDict{Volatile, VolatileValue}
    _writes::IdDict{Volatile, NTuple{2,VolatileValue}}
    _actors::IdDict{Actor, _ActorTxData}
    _sources::IdDict{Source, _SourceTxData}
    Transaction() = new(:running,
                        IdDict{Volatile, VolatileValue}(),
                        IdDict{Volatile, NTuple{2,VolatileValue}}(),
                        IdDict{Actor, _ActorTxData}(),
                        IdDict{Source, _SourceTxData}())
end

# The current running transaction.
@var _current_tx = nothing::Union{Nothing, Transaction}

"""
    current_tx()

Yields the current transaction, or nothing if there is no transaction running.
"""
current_tx() = _current_tx[]

const TX_MAX_ATTEMPTS = 2^14

_tx_lock(u::Volatile{T}) where {T} = lock(getfield(u, :_mutex))
_tx_unlock(u::Volatile{T}) where {T} = unlock(getfield(u, :_mutex))
_tx_lock(u::Actor{T}) where {T} = lock(u._cond)
_tx_unlock(u::Actor{T}) where {T} = unlock(u._cond)
_tx_lock(u::Source{T,K}) where {T,K} = lock(u._cond)
_tx_unlock(u::Source{T,K}) where {T,K} = unlock(u._cond)

"""
    tx(fn, args...)

Runs the given function in a transaction and yields the result. The function
fn is called as `fn(args...)`.
"""
function tx(fn::F, args...) where {F <: Function}
    success = false
    res = nothing
    currtx = _current_tx[]
    (currtx === nothing) || return fn(args...)
    # We need to make a new transaction for this task. We might have to do
    # this a few times if the transaction fails.
    for attempt in 1:TX_MAX_ATTEMPTS
        currtx = Transaction()
        _current_tx[] = currtx
        try
            res = fn(args...)
        catch e
            if isa(e, SourceReadError)
                continue
            else
                rethrow(e)
            end
        finally
            _current_tx[] = nothing
        end
        # At this point the transaction is over and we need to ensure that the
        # current volatile values are consistent with the transaction we just
        # executed. Before we do that, we should finalize the volatiles.
        currtx._state = :finalizing
        for (v,(x0,x1)) in currtx._writes
            x = _volatile_finalize(x1)
            if x !== x1
                currtx._state = :running
                _volatile_setindex!(v, x, currtx)
                currtx._state = :finalizing
            end
        end
        # We're now ready to attempt a commit. To do this we must lock the vols
        # and actors in a sorted order. We do volatiles first, since their
        # values must all be tested, then actors.
        currtx._state = :locking
        nw = length(currtx._writes)
        nr = length(currtx._reads)
        n  = nw + nr
        m  = length(currtx._actors)
        p  = length(currtx._sources)
        vols = Vector{Volatile}(undef, n)
        vols[1:nr] .= keys(currtx._reads)
        vols[nr+1:end] .= keys(currtx._writes)
        sort!(vols, by=objectid)
        acts = Vector{Actor}(undef, m)
        acts[1:m] .= keys(currtx._actors)
        sort!(acts, by=objectid)
        srcs = Vector{Source}(undef, p)
        srcs[1:p] .= keys(currtx._sources)
        sort!(srcs, by=objectid)
        # go through and lock; check while doinng so.
        locked = 0
        checked = 0
        success = false
        try
            while !success
                # First, go through and lock all the volatiles in order.
                for vol in vols
                    _tx_lock(vol)
                    locked += 1
                end
                # Now check their values. First do the read volatiles.
                for (vol,x0) in currtx._reads
                    (x0 === vol._value) || break
                    checked += 1
                end
                (checked == length(currtx._reads)) || break
                # Next, check the written volatiles.
                for (v,(x0,x)) in currtx._writes
                    (x0 === v._value) || break
                    checked += 1
                end
                (checked == length(vols)) || break
                # Next we need to lock and check the actors that were read from.
                for act in acts
                    _tx_lock(act)
                    locked += 1
                end
                for (act,dat) in currtx._actors
                    (dat.start_value === act._value) || break
                    checked += 1
                end
                (checked == n + m) || break
                # Next we need to lock and check the sources that were received.
                for src in srcs
                    _tx_lock(src)
                    locked += 1
                end
                for (src,dat) in currtx._sources
                    (src._state == :error) && throw(src._error)
                    (dat.start_position == src._position) || break
                    checked += 1
                end
                (checked == n + m + p) || break
                # If we reach this point, then everything checks out and we can commit
                # all of the changes and complete the transaction.
                for (vol,(x0,x)) in currtx._writes
                    vol._value = x
                end
                # We need to commit the actors as well.
                for (act,dat) in currtx._actors
                    q = dat.msgs
                    v = act._value
                    while !isempty(q)
                        msg = DataStructures.dequeue!(q)
                        DataStructures.enqueue!(act._msgs, msg)
                    end
                    # We notify the condition since there's something new to process.
                    notify(act._cond)
                end
                # And we need to claim the received items from the sources
                for (src,dat) in currtx._sources
                    _source_claim_nolock(src, dat.start_position, dat.count)
                end
                # That's all we need to do, aside from unlock.
                success = true
            end
        finally
            for mux in 1:locked
                (mux > n) && break
                _tx_unlock(vols[mux])
            end
            locked -= n
            if locked > 0
                for mux in 1:locked
                    (mux > m) && break
                    _tx_unlock(acts[mux])
                end
                locked -= m
                if locked > 0
                    for mux in 1:locked
                        #(mux > p) && break
                        _tx_unlock(src[mux])
                    end
                end
            end
        end
        success && break
    end
    # It's possible we got here because we were suceessful, but it might be that
    # we failed too many times.
    success || error("transaction aborted after failing $TX_MAX_ATTEMPTS times")
    return res
end

"""
    @tx expr

The macro @tx should be followed by an expression; that expression is run in an
atomic transaction.
"""
macro tx(expr)
    return :(tx(() -> $(esc(expr))))
end

# Now that we have definend transactions, we can define the volatile access
# methods.
_volatile_getindex(v::Volatile{T}, t::Transaction) where {T} = begin
    w = get(t._writes, v, nothing)
    (w === nothing) || return w[2]
    w = get(t._reads, v, nothing)
    (w === nothing) || return w
    w = v._value
    t._reads[v] = w
    return w
end
_volatile_getindex(v::Volatile{T}, ::Nothing) where {T} = v._value
Base.getindex(v::Volatile{T}) where {T} = _volatile_getindex(v, current_tx()).value
_volatile_setindex!(v::Volatile{T}, x::VolatileValue{T}, t::Transaction) where {T} = begin
    (t._state === :running) || error(
        "volatiles cann only be set before finalizing the transaction")
    # See if it this ref has already been written to in this transaction.
    w = get(t._writes, v, nothing)
    if w !== nothing
        (x0,x1) = w
        # If this doesn't change annything, don't do aything.
        (x === x1) && return x
        # If this changes the volatile back to its inintial value, we convert
        # this to a read operation. Otherwise we record the new value.
        if x === x0
            delete!(t._writes, v)
            t._reads[v] = x0
        else
            t._writes[v] = (x0, x)
        end
        return x
    end
    # See if this volatile has already been read.
    w = get(t._reads, v, nothing)
    if w !== nothing
        # There is a previously-read value already. If it hasn't changed, it's
        # still just a read value.
        (w === x) && return x 
        # Otherwise, we are going to promote it to a written value.
        delete!(t._reads, v)
        t._writes[v] = (w, x)
    else
        # There is no previously-read value; grab the current value.
        w = v._value
        if w === x
            t._reads[v] = w
        else
            t._writes[v] = (w, x)
        end
    end
    return x
end
_volatile_setindex!(v::Volatile{T}, x::S, ::Nothing) where {T,S} = error(
    "cannot set volatile outside of transaction")
Base.setindex!(v::Volatile{T}, x::S) where {T,S} = _volatile_setindex!(
    v, VolatileValue{T}(x, v._value.filter, v._value.finalize), current_tx())
Base.show(io::IO, ::MIME"text/plain", a::Volatile{T}) where {T} = begin
    print(io, "Volatile{$T}(@$(string(objectid(a), base=62)): $(a[]))")
end
"""
    getfilter(v)

Yields the filter-function for the volatile v.
"""
getfilter(v::Volatile{T}) where {T} = _volatile_getindex(v, current_tx()).filter
"""
    getfinalize(v)

Yields the finalize-function for the volatile v.
"""
getfinalize(v::Volatile{T}) where {T} = _volatile_getindex(v, current_tx()).finalize
"""
    setfilter!(vol, fn)

Sets the filter-function associated with the Volatile object vol. Any time that
the vol is set (`vol[] = x`) the filter-function is called and the value saved
in `vol` is instead `fn(x)`. This must be called within a transaction.
"""
setfilter!(vol::Volatile{T}, f::Function) where {T} = _volatile_setindex!(
    v, VolatileValue{T}(v._value.value, f, v._value.finalize), current_tx())

"""
    setfinalize!(vol, fn)

Sets the finalize-function associated with the Volatile object vol. Any time
that a transaction contains a change to `vol`, immediately prior to making an
attempt at committing the transaction, the finalize function is called and
the value committed to `vol` is instead `fn(x)` where `x` is the value set to
`vol` in the transaction. This funvtion must also be called within a
transaction.
"""
setfinalize!(vol::Volatile{T}, f::Function) where {T} = _volatile_setindex!(
    v, VolatileValue{T}(v._value.value, v._value.filter, f), current_tx())

# Same for actors: we can put the actor code here.
# First, we definne some private functions for handling things on the back
# end; we name these recv and send
_txdata(a::Actor{T}, t::Transaction) where {T} = begin
    # See if we've already grabbed it in the transaction.
    w = get(t._actors, a, nothing)
    if w === nothing
        w = _ActorTxData{T}(a._value)
        t._actors[a] = w
    end
    return w
end
_txdata(a::Actor{T}, ::Nothing) where {T} = nothing
_recv(a::Actor{T}, t::Transaction) where {T} = _txdata(a, t).tx_value
_recv(a::Actor{T}, ::Nothing) where {T} = a._value
_send(a::Actor{T}, argno::Int, args::Tuple, t::Transaction) where {T} = begin
    w = _txdata(a, t)
    if isa(w.tx_value, ActorException{T})
        if argno == 0
            # This is a reset!
            w.tx_value = args[1]
        else
            throw(w.tx_value)
        end
    elseif argno == 0
        error("attempt to reset actor that is not in an error state")
    end
    DataStructures.enqueue!(w.msgs, _ActorMsg{T}(argno, args))
    return nothing
end
_send(a::Actor{T}, argno::Int, args::Tuple, ::Nothing) where {T} = begin
    # We're not in a transaction, so we just lock and queue.
    lock(a._cond)
    try
        # There could be an error on the actor; if so we want to throw an
        # exception ourself.
        v = a._value
        if isa(v, ActorException{T})
            if argno == 0
                # This is a reset!
                _actor_reset(a, args[1])
            else
                throw(w.tx_value)
            end
        elseif argno == 0
            error("attempt to reset actor that is not in an error state")
        end
        # The actor isn't in an error state; we can schedule the function.
        DataStructures.enqueue!(a._msgs, _ActorMsg{T}(argno, args))
        # We notify the condition since there's something new to process.
        notify(a._cond)
        # That's all.
    finally
        unlock(a._cond)
    end
    return nothing
end
_recv(a::Actor{T}) where {T} = _recv(a, curr_tx())
_send(a::Actor{T}, argno::Int, args::Tuple) where {T} = _send(a, argno, args, curr_tx())

Base.getindex(a::Actor{T}) where {T} = begin
    x = _recv(a, current_tx())
    isa(x, ActorException{T}) && throw(x)
    return x
end
"""
   send(f, args...)

The send() method can be used to send a function to an actor. The actor may be
the function f or it may be any of the args. When send() is called, it searches
for the first occurance of an actor object in the argument list and dispatches
on that actor. In the actual function call that is evaluated, the actor object
is replaced with its stored value. The return value of the function becomes the
new stored value of the actor.

When an actor is sent a function, the function is queued to be run in some other
thread. The actor processes these threads synchronously.

If send() is called inside of a transaction, then the function is not
immediately queued but rather is held until the transaction successfully
completes. Once this happens, all sends are dispatched simultaneously such that
two sends to the same actor during the same transaction will always run
sequentially in the actor's thread.

Note that due to the design pattern described in the previous paragraph, one
should never wait on the result of a function sent to an actor during a
transaction. This might be done via a promise, for example, where the function
sent to the actor delivers to a promise also sent to the actor. Because this
function is not queued until the transaction is finished, and because the
promise will never recieve an object until the function is queued, requesting
the value of the promise will deadlock.

Note that the value of an actor when send() is does not make any guarantee
about its value when the function that is sent gets evaluated--other messages
may be processed in the interim.
"""
send(args...) = begin
    for (ii,arg) in enumerate(args)
        isa(arg, Actor) && return _send(arg, ii, args, current_tx())
    end
    error("send() must be passed at least one Actor object")
end
send(a::Actor{T}, args...) where {T} =
    _send(a, 1, (a, args...), current_tx())
send(f::F, a::Actor{T}, args...) where {F, T} =
    _send(a, 2, (f, a, args...), current_tx())
send(f::F, x::X, a::Actor{T}, args...) where {F, X, T} =
    _send(a, 3, (f, x, a, args...), current_tx())
send(f::F, x::X, y::Y, a::Actor{T}, args...) where {F, X, Y, T} =
    _send(a, 4, (f, x, y, a, args...), current_tx())
send(f::F, x::X, y::Y, z::Z, a::Actor{T}, args...) where {F, X, Y, Z, T} =
    _send(a, 4, (f, x, y, z, a, args...), current_tx())

# The actorerror and reset functions.
_actorerror(a::Actor{T}, t::Transaction) where {T} = begin
    # See if we've already grabbed it in the transaction.
    w = get(t._recvs, a, nothing)
    if w === nothing
        w = v._value
        t._recvs[a] = w
    end
    return isa(w, ActorException{T}) ? w : nothing
end
_actorerror(a::Actor{T}, ::Nothing) where {T} = begin
    v = a._value
    return (isa(v, ActorException{T}) ? v : nothing)
end
"""
    getrerror(actor)

If the given actor object is currently in an error state, then yields the
ActorException object that describes the error. Otherwise, yields nothing.
"""
geterror(a::Actor{T}) where {T} = begin
    x = _recv(a, current_tx())
    return isa(x, ActorException{T}) ? x : nothing
end
"""
    reset(actor, x)

If the given actor is in an error state, this resets it, meaning it will start
handling sent messages again, and gives it the new initial value x. Yields the
actor.
"""
reset(a::Actor{T}, s::S) where {T,S} = begin
    _send(a, 0, (s,), current_tx())
    return a
end


# Source functions.
"""
    geterror(s)

Yields the error raised by the source-kernel for the source object s
if s is in an error-state; otherwise yields nothing.
"""
geterror(s::Source{T,K}) where {T,K} = lock(s._cond) do
    return s._state == :error ? s._error : nothing
end
"""
    reset(s, k)

Resets a source object s with the new source kernel K (which must be of the same
kernel type as s, and which may be the same kernel that s already has). This
clears the error state from the source and allows it to be read again. Yields s.
"""
reset(s::Source{T,K}, newkern::K) where {T,K} = lock(s._cond) do
    (s._state   ==  :error)  || error("source object is not in an error state")
    (s._waiting === nothing) || schedule(s._waiting,
                                         SourceAbortError{T,K}(s),
                                         error=true)
    s._waiting = nothing
    s._kernel = newkern
    s._error = nothing
    s._state = :ok
    s._buffer = Vector{T}()
    s._position = 0x0
    return s
end
# Used by _source_ensure below.
_source_ensure_nolock(
    s::Source{T,K},
    n::UInt128,
    start_pos::Union{Nothing,UInt128}) where {T,K} = begin
    # Some checks to make sure that our state is fine.
    (length(s._buffer) >= n) && return (s._position, s._buffer[1:n])
    (start_pos !== nothing) && (s._position != start_pos) && throw(
        SourceReadError{T,K}(s, pos0))
    # Okay, we're goig to start reading... we loop until there are enough;
    # if something else claims the start position from underneath us during
    # this time, we can abort.
    start_pos = s._position
    buffer = s._buffer
    curr_task = current_task()
    while s._position == start_pos && length(s._buffer) < n
        if s._state === :ok
            # It's our job to read from the source next. We need to unlock the
            # cond to do this then relock onnce we've succeeded.
            s._state = :reading
            s._waiting = curr_task
            unlock(s._cond)
            try
                x = pop!(s._kernel)
            catch e
                # Uh-oh, we are now in an error-state... relock the condition.
                lock(s._cond)
                # In the unique case that what was throws was areset interrupt
                # for our task, we should not attempt to update the source; we
                # should just rethrow the error.
                if !isa(e, SourceAbortError{T,K})
                    # This is the only place allowed to set error state, so we
                    # can be confident that an error wasn't constructed by
                    # another task.
                    s._state = :error
                    s._error = SourceKernelError{T,K}(s, curr_task, e)
                    s._waiting = nothing
                    # We've raised an error on the source, so notify all the
                    # waiting tasks.
                    notify(s._cond, nothing, all=true)
                end
                # Finally, rethrow the error.
                rethrow(s._error)
            end
            # we always have to relock! We are expected to leave this function
            # locked, and it could cause deadlocks if we don't.
            lock(s._cond)
            s._state = :ok
            s._waiting = nothing
            # Regardless of where the start positoin is now, or what has
            # happend to the buffer, we need to append this new value to
            # the buffer.
            push!(s._buffer, x)
            # We've updated the source-buffer, so notify the waiters.
            notify(s._cond, nothing, all=true)
        elseif s._state === :reading
            # Secondly, someone else might be reading the kernel right now. If
            # so we want to wait on the condition instead of reading ourselves.
            wait(s._cond)
        elseif s._state === :error
            # First of all, if the buffer is in an error state, we should throw
            # the exception that was raised.
            throw(s._error)
        else
            error("invalid source state: $(s._state)")
        end
    end
    # If the start position has moved and we # need to start all over.
    (start_pos == s._position) || throw(SourceReadError{T,K}(s, start_pos))
    # Otherwise, we can return the current data!
    return (start_pos, s._buffer[1:n])
end
"""
    _source_ensure(s, n)

Ensures that the source s contains at least n items in its queue then yields a
tuple (position, first_n) where position is the current stream position and
first_n is the vector of the first n stream elements.

The condition for source s must not be locked when this is called.

If the start position moves out from underneath the source as it is waiting on
the kernel, then a SourceReadError is thrown.
"""
_source_ensure(
    s::Source{T,K},
    n::UInt128,
    start_pos::Union{Nothing,UInt128}) where {T,K} = lock(s._cond) do
    return _source_ensure_nolock(s, n, start_pos)
end
_receive(s::Source{T,K}, n::UInt128, tx::Transaction) where {T,K} = begin
    dat = get(tx.sources, s, nothing)
    if dat === nothing
        if n == 0x0
            # We don't need to lock the cond here because it doesn't matter
            # which side of any state-change this particular read of the
            # _position member lands.
            tx._sources[s] = _SourceTxData(s._position, 0x0)
            return T[]
        else
            (pos, x) = _source_ensure(s, n, nothing)
            tx._sources[s] = _SourceTxData(pos, n)
            return x
        end
    else
        if n == 0x0
            return T[]
        else
            (pos, x) = _source_ensure(s, dat.count + n, dat.start_position)
            x = x[dat.count+1:end]
            dat.count += n
            return x
        end
    end
end
_receive(s::Source{T,K}, n::UInt128, ::Nothing) where {T,K} = while true
    try
        return lock(s._cond) do
            # Basic error checks.
            (s._state == :error) && throw(s._error)
            # If the buffer is full enough, we can do a claim.
            (length(s._buffer) >= n) && return _source_claim_nolock(s, s._position, n)
            # Otherwise, we'll have to ensure it.
            (pos,res) = _source_ensure_nolock(s, n, nothing)
            return _source_claim_nolock(s, pos, n)
        end
    catch e
        isa(e, SourceReadError{T,K}) || rethrow(e)
    end
end
"""
    receive(s)
    receive(s, n)

For a Source object s, yields the next value in its queue. Reading from a
source in this way is thread-safe as long as it is done within a transaction.
Inside of a transaction, all source reads are guaranteed to be sequential
(i.e., not interrupted by other threads).

The second argument, n, may be passed, in which case, no matter what it is (even
1 or 0), a vector of the next n elements is returned; this receives each of them
so they cannot be received later. In particular, you can receive(s, 0) to ensure
that, during a transaction, no other thread processes anything from the source
without actually reading from the source.
"""
receive(s::Source{T,K}) where {T,K} = _receive(s, 0x1, current_tx())[1]
receive(s::Source{T,K}, n::UInt128) where {T,K} = _receive(s, n, current_tx())
