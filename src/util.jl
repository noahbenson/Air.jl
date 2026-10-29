################################################################################
# util.jl
#
# Utilities used in the Air library that don't depend on other components of
# Air.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

# #Delay #######################################################################
struct _DelayPending
    _mux::ReentrantLock
    _fn::Function
end
"""
    Delay{T}

`Delay` objects can be used to lazily calculate a single value the first time
it is requested. They act like `Ref`s in that you access a delay `d` via `d[]`.
`Delay` objects are thread-safe and are functionally immutable.

See also: [`@delay`](@ref), [`LazyDict`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> # Create a Delay with a long run time.
       d = Delay{Int64}(() -> (println("Running."); sleep(2); 10))
Delay{Int64}(<...>)

julia> # Start a few threads, each of which attempt to read it. The function
       # will only run once.
       for th in [(Threads.@spawn d[]) for _ in 1:5]; wait(th) end
Running.

julia> # Ensure that it produced the correct value and doesn't run again.
       d[]
10
```
"""
mutable struct Delay{T} <: Base.Ref{T}
    _val::Union{_DelayPending,Some{T}}
    function Delay{T}(t) where {T}
        return new{T}(Some{T}(t))
    end
    # We want it to be possible to share a delayed value.
    function Delay{T}(d::Delay) where {T}
        dval = d._val
        if isa(dval, _DelayPending)
            return new{T}(_DelayPending(ReentrantLock(), () -> d[]))
        else
            return new{T}(Some{T}(dval.value))
        end
    end
    function Delay{T}(d::Delay{T}) where {T}
        return d
    end
    function Delay{T}(f::Function) where {T}
        return new{T}(_DelayPending(ReentrantLock(), f))
    end
end
Delay(t) = Delay{typeof(t)}(t)
Delay(d::Delay) = d
Delay(f::Function) = Delay{Any}(f)
"""
    @delay expression

Yields a `Delay` object that matches the given expression. The expression may be
one of the following:
 1. A function of no arguments, such as `() -> 10`; in this case the delay is made
    from this function (i.e., the RHS is the `->` expression that is delayed).
 2. A function with a set of symbol arguments evaluates the RHS but uses a let
    statement to wrap all the symbols in the LHS into a closure.
 3. An expression, which is treated as equivalent to `() -> expression`.
Optionally, the expression or LHS may be tagged with a type T. In this case, a
`Delay{T}` object is yielded instead of a `Delay{Any}`.

See also: [`Delay`](@ref), [`LazyDict`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> # Create a Delay with a long run time.
       d = (@delay (println("Running."); sleep(2); 10)::Int64)
Delay{Int64}(<...>)

julia> # Start a few threads, each of which attempt to read it. The function will
       # only run once.
       for th in [(Threads.@spawn d[]) for _ in 1:5]; wait(th) end
Running.

julia> # Ensure that it produced the correct value and doesn't run again.
       d[]
10

julia> # The display now shows the realized value also.
       d
Delay{Int64}(10)

julia> # Create a Delay with a locally-bound symbol.
       d2 = (@delay (d) -> (d[] / 20.0)::Float64)
Delay{Float64}(<...>)

julia> # We can rebind d without affecting d2.
       d = 10
10

julia> d2[]
0.5
```
"""
macro delay(e)
    # First, parse the expression-- is it a function or an expression?
    if e isa Expr && e.head == :->
        syms = e.args[1]
        expr = e.args[2]
    else
        syms = ()
        expr = e
    end
    # Look for an optional trailing type annotation. The expression may be a
    # block, whose elements include `LineNumberNode`s, or it may not be an
    # `Expr` at all (e.g. a literal).
    tmp = expr
    while tmp isa Expr && tmp.head === :block
        tmp = tmp.args[end]
    end
    T = (tmp isa Expr && tmp.head === :(::)) ? tmp.args[2] : :Any
    # If there are symbols, convert them over to expressions
    isa(syms, Symbol) && (syms = :(($syms,)))
    if isa(syms, Expr)
        syms = [
            (
                if isa(sym, Symbol)
                    :($sym = $sym)
                elseif sym.head == :(::)
                    :($(sym.args[1]) = $sym)
                else
                    throw(ArgumentError("closures must be symbols or tagged-symbols"))
                end
            ) for sym in syms.args
        ]
    end
    # Generate the code
    if length(syms) == 0
        return esc(:(Air.Delay{$T}(() -> $expr)))
    else
        return esc(:(
            let $(syms...)
                Air.Delay{$T}(() -> $expr)
            end
        ))
    end
end
# Some base functions.
Base.getindex(d::Delay{T}) where {T} = begin
    v = d._val
    if isa(v, _DelayPending)
        lock(v._mux)
        try
            if d._val === v
                u = (v._fn())::T
                d._val = Some{T}(u)
            else
                u = ((d._val)::Some{T}).value
            end
            return u
        finally
            unlock(v._mux)
        end
    else
        return v.value
    end
end
Base.isready(d::Delay{T}) where {T} = isa(d._val, Some{T})
function Base.setindex!(d::Delay{T}, x...) where {T}
    return throw(ArgumentError("setindex!: Delays are immutable"))
end
# A `Delay` is a `Ref`, so Base also offers the 0-dimensional `d[i] = x` form,
# which takes a `CartesianIndex{0}`; without this method the two signatures are
# ambiguous. Dropping the index delegates to the method above, so it errors for
# the same reason.
Base.setindex!(d::Delay{T}, x, ::CartesianIndex{0}) where {T} = setindex!(d, x)
Base.isequal(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequal(a[], b[])
Base.hash(d::Delay{T}) where {T} = 0x26c850a2957fa577 + hash(d[])
function Base.show(io::IO, ::MIME"text/plain", d::Delay{T}) where {T}
    v = d._val
    if isa(v, Some)
        print(io, "$(typeof(d))($(v.value))")
    else
        print(io, "$(typeof(d))(<...>)")
    end
end
export Delay, @delay

# #memoize #####################################################################
# First, this is a helper functionthat makes sure than a function-arg
# declaration has a name.
function _memoize_fixarg(arg::Expr)
    return (
        if arg.head == :(::) && length(arg.args) == 1
            Expr(:(::), gensym(), arg.args[1])
        else
            arg
        end
    )
end
# Arguments without a type annotation, e.g. `x` in `f(x)`, arrive as plain
# symbols rather than expressions.
_memoize_fixarg(arg) = arg
# Now the memoize macro itself.
"""
    @memoize name(args...) = expr
    @memoize name(args...) where {...} = expr

`@memoize` is a macro for declaring that the function declaration that follows
should be memoized in a private dictionary and any pre-calculated value should
be returned from that dictionary instead of being recalculated. All memoization
is thread-safe: `expr` is only ever evaluated by one thread at a time, and is
only ever evaluated once per unique set of arguments.

Note that arguments are memoized according to equality, so the use of mutable
arguments can result in undefined behavior of those arguments are later changed.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> @memoize fib(n::Int) = begin
           println("Calculating fib(\$n)...")
           if n < 1
               return 0
           elseif n == 1
               return 1
           else
               return fib(n-1) + fib(n - 2)
           end
       end::Int
fib (generic function with 1 method)

julia> fib(5)
Calculating fib(5)...
Calculating fib(4)...
Calculating fib(3)...
Calculating fib(2)...
Calculating fib(1)...
Calculating fib(0)...
5

julia> fib(6)
Calculating fib(6)...
8
```
"""
macro memoize(assgn::Expr)
    (assgn.head == :(=)) ||
        throw(ArgumentError("memconst must be given an assignment expression"))
    # Parse the assignment statement.
    lhs = assgn.args[1]
    expr = assgn.args[2]
    if lhs.head == :call
        fsym = lhs.args[1]
        args = [_memoize_fixarg(a) for a in lhs.args[2:end]]
        lhs = Expr(:call, fsym, args...)
    elseif lhs.head == :where
        fsig = lhs.args[1]
        fsym = fsig.args[1]
        args = [_memoize_fixarg(a) for a in fsig.args[2:end]]
        fsig = Expr(:call, fsym, args...)
        lhs = Expr(:where, fsig, lhs.args[2:end]...)
    else
        throw(ArgumentError("memconst assignment LHS must be a call or where expression"))
    end
    # See if the expr is tagged; if so, we have a particular type we can use in
    # the memoization dict.
    MT = expr.head === :(::) ? expr.args[2] : :Any
    # Make an expression for the tuple of arguments.
    argtup = Expr(:tuple, args...)
    # Symbols we will need in the generated code.
    s_cache = gensym("cache")
    s_delay = gensym("delay")
    s_lock = gensym("lock")
    s_tmp = gensym("tmp")
    s_val = gensym("val")
    esc(
        quote
            let $s_lock = ReentrantLock(),
                $s_cache = Dict{Tuple,$MT}(),
                $s_delay = Dict{Tuple,Delay{$MT}}(),
                $s_tmp,
                $s_val

                global $lhs = begin
                    lock($s_lock)
                    try
                        $s_tmp = get($s_cache, $argtup, $s_cache)
                        if $s_tmp === $s_cache
                            $s_tmp = get!(() -> Delay{$MT}(() -> $expr), $s_delay, $argtup)
                        else
                            return $s_tmp
                        end
                    finally
                        unlock($s_lock)
                    end
                    # If we get here, we've created or grabbed a delay for the
                    # arguments; go ahead and wait on it (outside of the lock so
                    # that we don't prevent other argument tuples from computing at
                    # the same time).
                    $s_val = $s_tmp[]
                    # Now, re-grab the lock and update the dictionaries.
                    lock($s_lock)
                    try
                        # Possibly another thread updated things before we got to it.
                        $s_tmp = get($s_cache, $argtup, $s_cache)
                        if $s_tmp === $s_cache
                            # We're the first task to finish the calculation and/or 
                            # the first to grab the lock. Fix the cache.
                            $s_cache[$argtup] = $s_val
                            delete!($s_delay, $argtup)
                            return $s_val
                        else
                            return $s_tmp
                        end
                    finally
                        unlock($s_lock)
                    end
                end
                $fsym
            end
        end,
    )
end
export @memoize

# #Promise #####################################################################
"""
    Promise{T}
    Promise{T}(; ondelivery = nothing, onfailure = nothing)
    Promise(; ondelivery = nothing, onfailure = nothing)

Promise objects represent placeholders for values that are produced at some
later point, typically by another task. A promise is settled exactly once: it is
either *delivered* with a value by `put!`, or *failed* with an error payload by
[`fail!`](@ref). Until then it is pending. `take` (or `p[]` for promise `p`)
yields the delivered value, throws [`PromiseFailed`](@ref) for a failed promise,
and suspends the calling task until the promise is settled.

The `ondelivery` keyword installs a callback that is called with the delivered
value, exactly once; `onfailure` installs one that is called with the failure
payload. The callbacks are fixed when the promise is constructed and cannot be
changed later. Their return values are ignored, so they are notifications
rather than transformations, and an exception either callback throws is not
caught. A callback is called in the latest world, so a long-lived task may
settle a promise whose callback was defined after that task was created.
Outside a transaction a callback runs in the task that settled the promise,
after the settlement has been made; inside a transaction it runs in the
finalize phase, before the transaction attempts to commit.

[`take`](@ref) also accepts a timeout and a cancellation token. A `Promise`
compares by identity — two promises are equal only when they are the same
object — and its hash never waits for a settlement, so a pending promise may be
used as a dictionary key.

A `Promise` participates in transactions: a `put!` or `fail!` inside a
transaction is staged and performed only if the whole transaction commits, and
`take`, `isready`, [`isfailed`](@ref) and [`geterror`](@ref) inside that
transaction see the staged settlement. See the "Promises and transactions"
section of the manual for the details.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> p = Promise{Symbol}()
Promise{Symbol}(<...>)

julia> isready(p)
false

julia> put!(p, :value)
:value

julia> isready(p)
true

julia> take(p)
:value

julia> p[]
:value
```
"""
mutable struct Promise{T}
    # The mutex is persistent, unlike the condition in the original one-field
    # representation: a transaction's commit must be able to lock the promise
    # in order to validate and settle it, and it must be able to do so whether
    # or not the promise is still pending. `cond` is bound to `mutex`, so
    # waiting and locking are the same lock and `notify` may be called with it
    # held.
    mutex::ReentrantLock
    cond::Threads.Condition
    # At most one of `_val` and `_err` is non-`nothing`; both transitions happen
    # under `mutex`. `Some` wrappers rather than plain values, so a promise
    # delivered `nothing` or failed `nothing` is still distinguishable from a
    # pending one.
    _val::Union{Nothing,Some{T}}
    _err::Union{Nothing,Some{Any}}
    ondelivery::Union{Nothing,Function}
    onfailure::Union{Nothing,Function}
    # The waiters parked in `take`; created on the first `take` and dropped when
    # it is the last, so a promise that is never awaited allocates nothing for
    # them. Guarded by `mutex`.
    waiters::Union{Nothing,Base.IdSet{Waiter}}
    function Promise{T}(ondelivery::Union{Nothing,Function},
                        onfailure::Union{Nothing,Function}) where {T}
        mux = ReentrantLock()
        return new{T}(mux, Threads.Condition(mux), nothing, nothing, ondelivery,
                      onfailure, nothing)
    end
end
Promise{T}(; ondelivery::Union{Nothing,Function} = nothing,
           onfailure::Union{Nothing,Function} = nothing) where {T} =
    Promise{T}(ondelivery, onfailure)
Promise(; ondelivery::Union{Nothing,Function} = nothing,
        onfailure::Union{Nothing,Function} = nothing) =
    Promise{Any}(ondelivery, onfailure)

"""
    PromiseFailed(err)

An exception thrown by [`take`](@ref) (and by `p[]`) when the promise was
settled as failed with [`fail!`](@ref). The failure payload is stored as `err`.

See also: [`Promise`](@ref), [`fail!`](@ref), [`isfailed`](@ref),
[`geterror`](@ref)
"""
struct PromiseFailed <: Exception
    err::Any
end
export PromiseFailed
function Base.showerror(io::IO, e::PromiseFailed)
    print(io, "promise failed: ")
    show(io, getfield(e, :err))
    return nothing
end

"""
    ReaderTimeout(timeout)

An exception thrown by a blocking read when the given timeout expired: by
[`take`](@ref) when the promise was still pending, and by [`take`](@ref) or
`put!` when a [`Stream`](@ref) was still empty or full. The timeout, in seconds,
is stored as `timeout`.

See also: [`Promise`](@ref), [`Stream`](@ref), [`take`](@ref)
"""
struct ReaderTimeout <: Exception
    timeout::Real
end
export ReaderTimeout
function Base.showerror(io::IO, e::ReaderTimeout)
    return print(io, "reader timed out after ", getfield(e, :timeout), " seconds")
end
# Some base functions.
"""
    take(promise; timeout = nothing, token = canceltoken())

Yields the value delivered to the given `Promise` object after suspending the
current thread to wait for the value if necessary, and throws
[`PromiseFailed`](@ref) if the promise was failed. If the current transaction
has staged a settlement for the promise, that settlement is used without
waiting; see [`Promise`](@ref).

`timeout` bounds the wait in seconds: if the promise is still pending when it
expires, `take` throws [`ReaderTimeout`](@ref). A settled promise is returned or
thrown immediately, even with `timeout = 0`.

`token` is the [`CancelToken`](@ref) to observe; the default is the calling
task's own token, see [`canceltoken`](@ref). If the token has been cancelled,
`take` throws [`Cancelled`](@ref) even when the promise is already settled, so
pass `token = nothing` to wait regardless — in the cleanup of a cancelled task,
for example. Inside a transaction, a `take` that has to wait is a deadlock: the
delivery it is waiting for cannot happen until the body returns. (Waiting on a
[`Stream`](@ref) inside a transaction is not a deadlock in the same way; see its
docstring.)


```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> p = Promise()
Promise{Any}(<...>)

julia> put!(p, :done)
:done

julia> p
Promise{Any}(:done)

julia> take(p)
:done
```
"""
function take(d::Promise{T};
              timeout::Union{Nothing,Real} = nothing,
              token::Union{Nothing,CancelToken} = canceltoken()) where {T}
    # Cancellation wins over settlement, including a settlement that is already
    # present: a task whose work has been cancelled should not go on to use a
    # value it was waiting for. `token = nothing` opts out.
    (token === nothing) || checkcancel(token)
    # Read-your-writes: a settlement staged by the current transaction is
    # visible to that transaction, even though it has not committed. Without
    # this, a `take` after a `put!` in the same transaction would wait for a
    # settlement that cannot happen until the transaction's body returns.
    t = currtx()
    if t !== nothing
        dat = get(getfield(t, :promises), d, nothing)
        (dat === nothing) || return _staged_value(dat)
    end
    # The common case, and the only one that allocates nothing: an already
    # settled promise.
    r = _try_settle(d)
    (r === nothing) || return r.val
    (timeout === nothing) || (timeout >= 0) ||
        throw(ArgumentError("timeout must be non-negative"))
    w = Waiter()
    timer = nothing
    try
        _promise_register!(d, w)
        (token === nothing) || _register_waiter!(token, w)
        (timeout === nothing) || (timer = Timer(_ -> _timeout!(w), timeout))
        while true
            # Re-check now that a wake-up cannot be missed: the waiter is
            # registered before these checks, and settlement, timeout and
            # cancellation all wake it.
            (token === nothing) || checkcancel(token)
            r = _try_settle(d)
            (r === nothing) || return r.val
            if _park!(w)
                # The timer may have raced a settlement; one last check lets
                # the settlement win.
                r = _try_settle(d)
                (r === nothing) || return r.val
                throw(ReaderTimeout(timeout))
            end
        end
    finally
        _promise_unregister!(d, w)
        (token === nothing) || _unregister_waiter!(token, w)
        (timer === nothing) || close(timer)
    end
end

# A promise's settlement as `take` sees it: `nothing` means still pending, and a
# failure throws rather than being wrapped. The wrapper, rather than a bare
# value, because a promise may legitimately deliver `nothing`.
struct _Settled{T}
    val::T
end

# The settlement of `d`, or `nothing` if it is still pending. The promise's
# mutex is the only lock taken.
function _try_settle(d::Promise{T}) where {T}
    mux = getfield(d, :mutex)
    lock(mux)
    try
        v = getfield(d, :_val)
        (v !== nothing) && return _Settled{T}((v::Some{T}).value)
        e = getfield(d, :_err)
        (e === nothing) || throw(PromiseFailed((e::Some{Any}).value))
        return nothing
    finally
        unlock(mux)
    end
end

# The settlement staged for `d` by the current transaction; see `take`.
function _staged_value(dat)
    v = getfield(dat, :val)
    (v !== nothing) && return v.value
    return throw(PromiseFailed((getfield(dat, :err)::Some{Any}).value))
end

# Park `w` on `d` so that a settlement wakes it. The promise's mutex is the only
# lock taken; see the lock-order note in `cancel.jl`.
function _promise_register!(d::Promise, w::Waiter)
    mux = getfield(d, :mutex)
    lock(mux)
    try
        ws = getfield(d, :waiters)
        if ws === nothing
            ws = Base.IdSet{Waiter}()
            setfield!(d, :waiters, ws)
        end
        push!(ws, w)
    finally
        unlock(mux)
    end
    return nothing
end
function _promise_unregister!(d::Promise, w::Waiter)
    mux = getfield(d, :mutex)
    lock(mux)
    try
        ws = getfield(d, :waiters)
        if ws !== nothing
            delete!(ws, w)
            isempty(ws) && setfield!(d, :waiters, nothing)
        end
    finally
        unlock(mux)
    end
    return nothing
end
# Wake every registered waiter and drop the registry. The caller holds `d`'s
# mutex.
function _wake_waiters!(d::Promise)
    ws = getfield(d, :waiters)
    (ws === nothing) && return nothing
    setfield!(d, :waiters, nothing)
    for w in ws
        _wake!(w)
    end
    return nothing
end

Base.getindex(d::Promise{T}) where {T} = take(d)
function Base.put!(d::Promise{T}, x) where {T}
    s = Some{T}(x)
    t = currtx()
    (t === nothing) && return _promise_deliver!(d, s)
    return _promise_stage!(d, s, t)
end

"""
    fail!(promise, err)

Settles `promise` as failed with the given error payload and yields `err`; a
later [`take`](@ref) throws [`PromiseFailed`](@ref) carrying it. A promise
settles exactly once, so failing one that has already been delivered or failed
raises `ArgumentError`. Inside a transaction, the failure is staged like a
delivery and happens only if the transaction commits.

See also: [`Promise`](@ref), [`take`](@ref), [`isfailed`](@ref),
[`geterror`](@ref)
"""
function fail!(d::Promise{T}, err) where {T}
    e = Some{Any}(err)
    t = currtx()
    (t === nothing) && return _promise_fail!(d, e)
    return _promise_stage_failure!(d, e, t)
end
export fail!
# Stage a settlement in the given transaction; the commit performs it. `t` is
# deliberately unannotated: `Transaction` is defined in `TX.jl`, which is
# included after this file, and the annotation would be resolved at definition
# time.
function _promise_stage!(d::Promise{T}, s::Some{T}, t) where {T}
    _promise_check_stage!(d, t, "given promise is already fulfilled")
    getfield(t, :promises)[d] = PromiseTxData{T}(s, nothing)
    return s.value
end
function _promise_stage_failure!(d::Promise{T}, e::Some{Any}, t) where {T}
    _promise_check_stage!(d, t, "given promise is already settled")
    getfield(t, :promises)[d] = PromiseTxData{T}(nothing, e)
    return e.value
end
# The checks shared by both staging paths: the transaction must still be
# running, and the promise must be unsettled both in the transaction and in
# fact. A promise is one-shot either way, so a second settlement is a permanent
# error rather than a retryable conflict. `msg` names the settlement that was
# attempted, so `put!` keeps its historical message and `fail!` reports a
# settlement.
function _promise_check_stage!(d::Promise{T}, t, msg::String) where {T}
    (getfield(t, :state) === :running) ||
        error("promises can only be settled before finalizing the transaction")
    haskey(getfield(t, :promises), d) && throw(ArgumentError(msg))
    mux = getfield(d, :mutex)
    lock(mux)
    try
        _promise_pending(d) || throw(ArgumentError(msg))
    finally
        unlock(mux)
    end
    return nothing
end
# Whether `d` is unsettled. The caller holds `d`'s mutex.
_promise_pending(d::Promise) =
    (getfield(d, :_val) === nothing) && (getfield(d, :_err) === nothing)

# Settle immediately, outside any transaction. The settlement callback runs
# after the lock is released, so it may take or settle other promises; an
# exception from it propagates, but the promise remains settled. It is called
# with `invokelatest` for the reason the actor's functions are: the settling
# task may be long-lived, and the callback may have been written after that
# task was created.
function _promise_deliver!(d::Promise{T}, s::Some{T}) where {T}
    cb = getfield(d, :ondelivery)
    mux = getfield(d, :mutex)
    lock(mux)
    try
        _promise_pending(d) || throw(ArgumentError("given promise is already fulfilled"))
        setfield!(d, :_val, s)
        notify(getfield(d, :cond); all=true)
        _wake_waiters!(d)
    finally
        unlock(mux)
    end
    (cb === nothing) || Base.invokelatest(cb, s.value)
    return s.value
end
function _promise_fail!(d::Promise{T}, e::Some{Any}) where {T}
    cb = getfield(d, :onfailure)
    mux = getfield(d, :mutex)
    lock(mux)
    try
        _promise_pending(d) || throw(ArgumentError("given promise is already settled"))
        setfield!(d, :_err, e)
        notify(getfield(d, :cond); all=true)
        _wake_waiters!(d)
    finally
        unlock(mux)
    end
    (cb === nothing) || Base.invokelatest(cb, e.value)
    return e.value
end

function Base.isready(d::Promise{T}) where {T}
    t = currtx()
    if t !== nothing
        dat = get(getfield(t, :promises), d, nothing)
        (dat === nothing) || return true
    end
    mux = getfield(d, :mutex)
    lock(mux)
    try
        return !_promise_pending(d)
    finally
        unlock(mux)
    end
end

"""
    isfailed(promise)

Yields `true` if the given `Promise` has been failed with [`fail!`](@ref), and
`false` if it is pending or delivered. Inside the transaction that staged a
failure, it reflects the staged settlement, like `isready`.

See also: [`Promise`](@ref), [`geterror`](@ref), [`PromiseFailed`](@ref)
"""
function isfailed(d::Promise{T}) where {T}
    t = currtx()
    if t !== nothing
        dat = get(getfield(t, :promises), d, nothing)
        if dat !== nothing
            return getfield(dat, :err) !== nothing
        end
    end
    mux = getfield(d, :mutex)
    lock(mux)
    try
        return getfield(d, :_err) !== nothing
    finally
        unlock(mux)
    end
end
export isfailed

"""
    geterror(promise)

Yields the failure payload of a `Promise` that was settled as failed with
[`fail!`](@ref), and `nothing` if it is pending or delivered. This is the
non-throwing counterpart of [`take`](@ref), which throws [`PromiseFailed`](@ref)
for a failed promise. Inside the transaction that staged a failure, it reflects
the staged settlement. A promise failed with a payload of `nothing` also yields
`nothing`, so use [`isfailed`](@ref) to tell that case from a pending or
delivered promise.

See also: [`Promise`](@ref), [`isfailed`](@ref), [`fail!`](@ref)
"""
function geterror(d::Promise{T}) where {T}
    t = currtx()
    if t !== nothing
        dat = get(getfield(t, :promises), d, nothing)
        if dat !== nothing
            e = getfield(dat, :err)
            return e === nothing ? nothing : e.value
        end
    end
    mux = getfield(d, :mutex)
    lock(mux)
    try
        e = getfield(d, :_err)
        return e === nothing ? nothing : (e::Some{Any}).value
    finally
        unlock(mux)
    end
end

# Identity semantics, as for any mutable object identified by its handle: a
# promise is equal only to itself, and its hash never waits for a settlement.
# (The value-comparing forms this replaces called `take`, so an undelivered
# promise as a dictionary key hung.)
Base.isequal(a::Promise, b::Promise) = a === b
Base.hash(d::Promise, h::UInt) = hash(objectid(d), hash(0x767127451c1e402a, h))

function Base.show(io::IO, ::MIME"text/plain", d::Promise{T}) where {T}
    mux = getfield(d, :mutex)
    lock(mux)
    local v, e
    try
        v = getfield(d, :_val)
        e = getfield(d, :_err)
    finally
        unlock(mux)
    end
    if v !== nothing
        print(io, "$(typeof(d))(")
        show(io, v.value)
        print(io, ")")
    elseif e !== nothing
        print(io, "$(typeof(d))(<failed: ")
        show(io, e.value)
        print(io, ">)")
    else
        print(io, "$(typeof(d))(<...>)")
    end
    return nothing
end
export Promise, take

# #lockall #####################################################################
function _lockall(locks::Vector{T}) where {T}
    locked = 0
    try
        for l in locks
            lock(l)
            locked += 1
        end
    finally
        (locked == length(locks)) || for l in locks
            unlock(l)
            locked -= 1
            (locked == 0) && break
        end
    end
    return locked
end
_lockall(f::Function, locks::Vector{T}) where {T} = begin
    # This only gets called once we have our own copy of the vector
    sort!(locks; by=objectid)
    # This will either lock them all or raise an exception.
    _lockall(locks)
    # Now we can run the function
    try
        return f()
    finally
        for l in locks
            unlock(l)
        end
    end
end
"""
    lockall(f, l1, l2, ...)

Locks all of the lockable objects `l1`, `l2`, etc. then runs `f`, unlocks the
objects, and returns the return value of `f()`.

See also: `Threads.ReentrantLock`

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> (r1, r2, r3) = [ReentrantLock() for _ in 1:3]
       lockall(r1, r2, r3) do; :success end
:success

julia> lockall([r1, r2, r3]) do; :success end
:success

julia> lockall((r1, r2, r3)) do; :success end
:success
```
"""
function lockall end
lockall(f::Function, locks::Vector{T}) where {T} = _lockall(f, copy(locks))
lockall(f::Function, locks::NTuple{N,T}) where {T,N} = _lockall(f, [locks...])
lockall(f::Function, locks::Vararg{T,N}) where {T,N} = _lockall(f, [locks...])
export lockall

# #_to_pairs ###################################################################
"""
    _pairtypes(::Type)

Yields the key and value types `(K, V)` of a two-parameter `Pair` or `Tuple`
type, or `nothing` if the type is not one. Reading a pair type's parameters by
dispatch rather than through the `parameters` field keeps this off the runtime
reflection path.
"""
@inline _pairtypes(::Type{<:Pair{K,V}}) where {K,V} = (K, V)
@inline _pairtypes(::Type{<:Tuple{K,V}}) where {K,V} = (K, V)
@inline _pairtypes(::Type) = nothing
# A utility function for turning a list of pairs/tuples into a list of pairs.
#
# The key and value types come from the container's declared `eltype` whenever
# that is a two-parameter pair or tuple type. That is both cheaper than looking
# at the elements — it avoids a `Vector{Any}` of keys and values and a
# `typejoin` over them — and equivalent, because `eltype` of a heterogeneous
# tuple is already the join of its element types. Only a container that declares
# no element type, or too vague a one, falls through to inspecting the elements.
_to_pairs(kvs) = begin
    if Base.IteratorEltype(kvs) isa Base.HasEltype
        kv = _pairtypes(Base.eltype(kvs))
        if kv !== nothing
            K, V = kv
            return Pair{K,V}[Pair{K,V}(t[1], t[2]) for t in kvs]
        end
    end
    ks = Any[]
    vs = Any[]
    for kv in kvs
        if kv isa Pair || (kv isa Tuple && length(kv) == 2)
            push!(ks, kv[1])
            push!(vs, kv[2])
        else
            throw(ArgumentError("EquivDict: arg must be iterator of tuples or pairs"))
        end
    end
    isempty(ks) && return Pair{Any,Any}[]
    K = typejoin(map(typeof, ks)...)
    V = typejoin(map(typeof, vs)...)
    return Pair{K,V}[Pair{K,V}(k, v) for (k, v) in zip(ks, vs)]
end
