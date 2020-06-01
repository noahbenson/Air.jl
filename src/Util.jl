################################################################################
# Util.jl
# Utilities for Air that don't depend on other components of Air.
# by Noah C. Benson

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

@memoize is a macro for declaring that the function declaration that follows                                                                                                                               should be memoized in a private dictionary and any pre-calculated value should                                                                                                                             be returned from that dictionary instead of being recalculated.
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
        throw(ArgumentError("memconst assignment LHS must be a call or where expression"))
    end
    # Make an expression for the tuple of arguments.
    argtup = Expr(:tuple, args...)
    # Symbols we will need in the generated code.
    sDict = gensym()
    sLock = gensym()
    sTmp  = gensym()
    quote
        $sLock = ReentrantLock()
        $sDict = Dict{Tuple, Any}()
        $lhs = begin
            lock($sLock)
            $sTmp = get!($sDict, $argtup) do; $expr end
            unlock($sLock)
            return $sTmp
        end
        function forget(::typeof($fsym))
            global $sDict
            $sTmp = $sDict
            $sDict = Dict{Tuple, Any}()
            return $sTmp
        end
        Base.keys(::typeof($fsym)) = keys($sDict)
        $fsym
    end |> esc
end

# #Delay #######################################################################
struct _DelayPending
    _mux::ReentrantLock
    _fn::Function
end
struct _ValueRealized{T}
    _val::T
end
"""
    Delay{T}

Delay objects can be used to lazily calculate a single value the first time
it is requested. They act like Refs in that you access a delay `d` via `d[]`.
"""
mutable struct Delay{T} <: Base.Ref{T}
    _val::Union{_DelayPending, _ValueRealized{T}}
end
Delay{T}(f::Function) where {T} = Delay{T}(_DelayPending(ReentrantLock(), f))
Delay(f::Function) = Delay{Any}(f)
"""
    @delay expression

Yields a Delay object that matches the given expression. The expression may be
one of the following:
 1. A function of no arguments, such as `() -> 10`; in this case the delay is made
    from this function (i.e., the RHS is the `->` expression that is delayed).
 2. A function with a set of symbol arguments evaluates the RHS but uses a let
    statement to wrap all the symbols in the LHS into a closure.
 3. An expression, which is treated as equivalent to `() -> expression`.
Optionally, the expression or LHS may be tagged with a type T. In this case, a
`Delay{T}` object is yielded instead of a `Delay{Any}`.
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
# Delays are immutable as long as the functions are pure; they are intended for use
# this way, so we declare them immutable.
mutability(::Type{Delay}) = Immutable
mutability(::Type{Delay{T}}) where {T} = Immutable
# Some base functions.
Base.getindex(d::Delay{T}) where {T} = begin
    v = d._val
    if isa(v, _DelayPending)
        lock(v._mux)
        try
            if d._val === v
                u = v._fn()
                d._val = _ValueRealized{T}(u)
            else
                u = d._val._val
            end
        finally
            unlock(v._mux)
        end
        return u
    else
        return v._val
    end
end
Base.isready(d::Delay{T}) where {T} = isa(d._val, _ValueRealized{T})
Base.setindex!(d::Delay{T}, x...) where {T} = throw(ArgumentError("setindex!: Delays are immutable"))
Base.isequal(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequal(a[], b[])
isequiv(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequiv(a[], b[])
Base.hash(d::Delay{T}) where {T} = 0x26c850a2957fa577 + hash(d[])
equivhash(d::Delay{T}) where {T} = 0x26c850a2957fa577 + equivhash(d[])
Base.show(io::IO, ::MIME"text/plain", d::Delay{T}) where {T} = begin
    v = d._val
    if isa(v, _ValueRealized)
        print(io, "$(typeof(d))($(v._val))")
    else
        print(io, "$(typeof(d))(<...>)")
    end
end

# #Promise #####################################################################
"""
    Promise{T}

Promise objects represent placeholders for values that may or may not have been
delivered yet. This is effectively a Channel object that can only be put! to
a single time and all take calls on the promise will return that value. Promise
values can be accessed via the `take()` function (not the `take!()` function) or
via `p[]` for promise `p`. In both cases, the running thread is suspended until
a value is delivered.
"""
mutable struct Promise{T}
    _val::Union{Threads.Condition, _ValueRealized{T}}
end
Promise{T}() where {T} = Promise{T}(Threads.Condition())
Promise() = Promise{Any}()
# Promises are immutable as long as the functions are pure; they are intended
# for use this way, so we declare them immutable.
mutability(::Type{Promise}) = Immutable
mutability(::Type{Promise{T}}) where {T} = Immutable
# Some base functions.
take(d::Promise{T}) where {T} = begin
    v = d._val
    if isa(v, Threads.Condition)
        lock(v)
        try
            (d._val === v) && wait(v)
            return d._val._val
        finally
            unlock(v)
        end
    else
        return v._val
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
            d._val = _ValueRealized{T}(x)
            notify(v, all=true)
            return x
        finally
            unlock(v)
        end
    else
        throw(ArgumentError("given promise is already fulfilled"))
    end
end
isready(d::Promise{T}) where {T} = isa(d._val, _ValueRealized{T})
Base.isequal(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequal(a[], b[])
isequiv(a::Delay{T}, b::Delay{S}) where {T,S} = (a === b) || isequiv(a[], b[])
Base.hash(d::Delay{T}) where {T} = 0x767127451c1e402a + hash(d[])
equivhash(d::Delay{T}) where {T} = 0x767127451c1e402a + equivhash(d[])
Base.show(io::IO, ::MIME"text/plain", d::Promise{T}) where {T} = begin
    v = d._val
    if isa(v, _ValueRealized)
        print(io, "$(typeof(d))($(v._val))")
    else
        print(io, "$(typeof(d))(<...>)")
    end
end

