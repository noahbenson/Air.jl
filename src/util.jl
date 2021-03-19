################################################################################
# util.jl
#
# Utilities used in the Air library that don't depend on other components of
# Air.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

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
    _val::Union{_DelayPending, Some{T}}
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
macro delay(e::Expr)
    # First, parse the expression-- is it a function or an expression?
    if e.head == :->
        syms = e.args[1]
        expr = e.args[2]
    else
        syms = ()
        expr = e
    end
    tmp = expr
    while tmp.head == :block; tmp = tmp.args[2] end
    if tmp.head == :(::)
        T = tmp.args[2]
    else
        T = :Any
    end
    # If there are symbols, convert them over to expressions
    isa(syms, Symbol) && (syms = :(($syms,)))
    if isa(syms, Expr)
        syms = [
            ( isa(sym, Symbol) ? :($sym = $sym)
              : sym.head == :(::) ? :($(sym.args[1]) = $sym)
              : throw(ArgumentError("closures must be symbols or tagged-symbols")))
            for sym in syms.args]
    end
    # Generate the code
    if length(syms) == 0
        return esc(:(Air.Delay{$T}(() -> $expr)))
    else
        return esc(:(let $(syms...); Air.Delay{$T}(() -> $expr) end))
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
Base.setindex!(d::Delay{T}, x...) where {T} = throw(ArgumentError("setindex!: Delays are immutable"))
Base.isequal(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequal(a[], b[])
Base.hash(d::Delay{T}) where {T} = 0x26c850a2957fa577 + hash(d[])
Base.show(io::IO, ::MIME"text/plain", d::Delay{T}) where {T} = begin
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
_memoize_fixarg(arg::Expr) = (arg.head == :(::) && length(arg.args) == 1
                              ? Expr(:(::), gensym(), arg.args[1])
                              : arg)
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
    (assgn.head == :(=)) || throw(
        ArgumentError("memconst must be given an assignment expression"))
    # Parse the assignment statement.
    lhs  = assgn.args[1]
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
        ArgumentError(
            "memconst assignment LHS must be a call or where expression"
        ) |> throw
    end
    # See if the expr is tagged; if so, we have a particular type we can use in
    # the memoization dict.
    MT = expr.head === :(::) ? expr.head : :Any
    # Make an expression for the tuple of arguments.
    argtup = Expr(:tuple, args...)
    # Symbols we will need in the generated code.
    s_cache = gensym("cache")
    s_delay = gensym("delay")
    s_lock = gensym("lock")
    s_tmp  = gensym("tmp")
    s_val  = gensym("val")
    quote
        let $s_lock  = ReentrantLock(),
            $s_cache = Dict{Tuple, $MT}(),
            $s_delay = Dict{Tuple, Delay{$MT}}(),
            $s_tmp, $s_val
            global $lhs = begin
                lock($s_lock)
                try
                    $s_tmp = get($s_cache, $argtup, $s_cache)
                    if $s_tmp === $s_cache
                        $s_tmp = get!(() -> Delay{$MT}(() -> $expr),
                                      $s_delay, $argtup)
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
    end |> esc
end
export @memoize

# #Promise #####################################################################
"""
    Promise{T}

Promise objects represent placeholders for values that may or may not have been
delivered yet. This is effectively a `Channel` object that can only be `put!` to
a single time and all `take` calls on the promise will return that value.
Promise values can be accessed via the `take()` function (not the `take!()`
function) or via `p[]` for promise `p`. In both cases, the running thread is
suspended until a value is delivered.

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
    _val::Union{Threads.Condition, Some{T}}
end
Promise{T}() where {T} = Promise{T}(Threads.Condition())
Promise() = Promise{Any}()
# Some base functions.
"""
    take(promise)

Yields the value delivered to the given `Promise` object after suspending the
current thread to wait for the value if necessary.

See also [`Promise`](@ref), [`put!`](@ref).


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
take(d::Promise{T}) where {T} = begin
    v = d._val
    if isa(v, Threads.Condition)
        lock(v)
        try
            (d._val === v) && wait(v)
            return ((d._val)::Some{T}).value
        finally
            unlock(v)
        end
    else
        return v.value
    end
end
Base.getindex(d::Promise{T}) where {T} = take(d)
Base.put!(d::Promise{T}, x) where {T} = begin
    v = d._val
    if isa(v, Threads.Condition)
        lock(v)
        try
            (d._val === v) || throw(
                ArgumentError("given promise is already fulfilled"))
            d._val = Some{T}(x)
            notify(v, all=true)
            return d._val.value
        finally
            unlock(v)
        end
    else
        throw(ArgumentError("given promise is already fulfilled"))
    end
end
Base.isready(d::Promise{T}) where {T} = isa(d._val, Some{T})
Base.isequal(a::Promise{T}, b::Promise{S}) where {T,S} = (a === b) || isequal(a[], b[])
Base.hash(d::Promise{T}) where {T} = 0x767127451c1e402a + hash(d[])
Base.show(io::IO, ::MIME"text/plain", d::Promise{T}) where {T} = begin
    v = d._val
    if isa(v, Some)
        print(io, "$(typeof(d))(")
        show(io, v.value)
        print(io, ")")
    else
        print(io, "$(typeof(d))(<...>)")
    end
end
export Promise, take

# #lockall #####################################################################
_lockall(locks::Vector{T}) where {T} = begin
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
    sort!(locks, by=objectid)
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

See also: [`ReentrantLock`](@ref)

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
lockall(f::Function, locks::Vector{T}) where {T,N} = _lockall(f, copy(locks))
lockall(f::Function, locks::NTuple{N,T}) where {T,N} = _lockall(f, [locks...])
lockall(f::Function, locks::Vararg{T,N}) where {T,N} = _lockall(f, [locks...])
export lockall

# #_to_pairs ###################################################################
# A utility function for turning a list of pairs/tuples into a list of pairs.
_to_pairs(kvs) = begin
    if length(kvs) == 0
        K = Any
        V = Any
        if Base.IteratorEltype(kvs) isa Base.HasEltype
            ET = Base.eltype(itr)
            if isa(ET, DataType)
                if ET <: Pair
                    K = ET.parameters[1]
                    V = ET.parameters[2]
                elseif ET <: Tuple && length(ET.parameters) == 2
                    K = ET.parameters[1]
                    V = ET.parameters[2]
                end
            end
        end
        return Pair{K,V}[]
    else
        ks = []
        vs = []
        for kv in kvs
            if kv isa Pair || (kv isa Tuple && length(kv) == 2)
                push!(ks, kv[1])
                push!(vs, kv[2])
            else
                msg = "EquivDict: arg must be iterator of tuples or pairs"
                throw(ArgumentError(msg))
            end
        end
        K = typejoin(map(typeof, ks)...)
        V = typejoin(map(typeof, vs)...)
        return Pair{K,V}[Pair{K,V}(k,v) for (k,v) in zip(ks,vs)]
    end
end
