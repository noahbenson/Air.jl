################################################################################
# vars.jl
#
# Implementation of a task-local scoped Var type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson


# #Var #########################################################################
# A tricky thing about Vars is that we want each Var object to be unique (i.e.,
# to have a unique objectid). However, we also want them to be immutable. To get
# around this, we use an empty dummy struct that is heap-allocated as an
# internal field. This way nothing is mutable, but there is a unique objectid
# For each Var.
mutable struct VarID end
"""
    Var{T}

A `Var` object represents a task-local piece of data with a default value. You
can access a `Var` with the `getindex` function (`var[]`) and you can set it
with the `setindex!` function (`var[] = newval`). However, the new assignment
will always be task-local only. Because of this, `Var`s are safe to access and
update in a multi-threaded program.

All fields of a `Var` should be considered strictly private.

See also: [`@var`](@ref), [`ReentrantRef`](@ref), [`vars`](@ref),
[`withvars`](@ref), [`wrapsetvars`](@ref), [`wrapwithvars`](@ref),
 `Threads.current_task`, `Threads.task_local_storage`.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{Symbol}\\(@[0-9a-zA-Z]+: :start_sym; init=:start_sym\\)"
julia> v = Var{Symbol}(:start_sym)
Var{Symbol}(@JIm7aUS2sOl: :start_sym; init=:start_sym)

julia> v[]
:start_sym

julia> withvars(v => :new_sym) do; v[] end
:new_sym

julia> withvars(v => :new_sym) do; fetch(Threads.@spawn v[]) end
:start_sym

julia> v[]
:start_sym
```
"""
struct Var{T} <: TransactionalRef{T}
    initial_value::T
    id::VarID
    function Var{T}(initval) where {T}
        u = convert(T, initval)
        return new{T}(u, VarID())
    end
    function Var{T}(initval::S) where {T, S <: T}
        return new{T}(initval, VarID())
    end
end
Var(initval::T) where {T} = Var{T}(initval)
"""
    VarsDict

The type of a dictionary of `Var` bindings, as returned by the function `vars()`
and as is required by the function `setvars()`.
"""
const VarsDict = IdDict{Var,Any}
const VarPair = Union{Pair{<:Var,<:Any},Tuple{<:Var,<:Any}}
# These functions require a private key.
"""
    vars()

Yields an `IdDict{Var,Any}` object that contains a mapping of all `Var` objects
whose current in-task value is not its default value to that `Var`'s currently
assigned value. The dictionary returned by `vars` may be later restored using
the function `setvars`.

See also: [`setvars`](@ref), [`withvars`](@ref), [`wrapwithvars`](@ref),
[`Var`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"IdDict{Var, ?Any}\\(\\)"
julia> vars()
IdDict{Var,Any}()
```

```jldoctest; filter=r"(Var{Symbol}\\(@[0-9a-zA-Z]+: :test; init=:test\\))|(IdDict{Var, ?Any} with 1 entry:)|(  Var{Symbol}\\(@[0-9a-zA-Z]+\\) => :temp)"
julia> @var v = :test::Symbol
Var{Symbol}(@h4G6oRR9s: :test; init=:test)

julia> withvars(v => :temp) do; vars() end
IdDict{Var,Any} with 1 entry:
  Var{Symbol}(@AaE16J5Pic8) => :temp
```
"""
function vars end
"""
    withvars(f, var1 => val1, var2 => val2...)
    withvars(f, vardict)

Runs the function `f` in a context in which the given `Var`s have been bound to
the given values. In the vase of a dictionary passed as the second argument, the
keys must be `Var` objects. Once the function `f` has finished running, the
`Var` objects are reverted to their calling-frame values.

See also: [`Var`](@ref), [`@var`](@ref), [`vars`](@ref), [`setvars`](@ref),
[`wrapwithvars`](@ref), [`wrapsetvars`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{Symbol}\\(@[0-9a-zA-Z]+: :initval; init=:initval\\)"
julia> v = Var{Symbol}(:initval)
Var{Symbol}(@cFgU9Kqe8: :initval; init=:initval)

julia> v[]
:initval

julia> withvars(v => :newval) do; v[] end
:newval
```
"""
function withvars end
"""
    setvars(f, vardict)

Runs the function `f` in a context in which all `Var` objectss have been bound
to the values given in the dictionary `vardict`. The keys of this dictionary
must be `Var` objects. Once the function `f` has finished running, the `Var`
objects are reverted to their calling-frame values. The current `vardict` object
can be obtained by calling `vars()`.

See also: [`Var`](@ref), [`@var`](@ref), [`vars`](@ref), [`withvars`](@ref),
[`wrapwithvars`](@ref), [`wrapsetvars`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{Symbol}\\(@[0-9a-zA-Z]+: :initval; init=:initval\\)"
julia> v = Var{Symbol}(:initval)
Var{Symbol}(@hxC65AWl: :initval; init=:initval)

julia> v[]
:initval

julia> setvars(IdDict(v => :newval)) do; v[] end
:newval
```
"""
function setvars end
"""
    wrapwithvars(f, var1 => val1, var2 => val2...)
    wrapwithvars(f, vardict)

Equivalent to `withvars(f, args...)` except that instead of running `f` 
immediately in the context of the modified `Var`s, yields a wrapper around `f`
that, when called, passes all arguments to `f`, which is run using the current
variable bindings plus any given bindings. Notable, `wrapwithvars(f)` will
create a wrapped version of `f` that uses the `Var` bindings in the current
task.

See also: [`Var`](@ref), [`@var`](@ref), [`vars`](@ref), [`withvars`](@ref),
[`setvars`](@ref), [`wrapsetvars`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{Int64}\\(@[0-9a-zA-Z]+: 0; init=0\\)"
julia> v = Var{Int}(0)
Var{Int64}(@JIm7aUS2sOl: 0; init=0)

julia> v[]
0

julia> f = wrapwithvars(v => 2) do x; x + v[] end; f(10)
12
```
"""
function wrapwithvars end
"""
    wrapsetvars(f, vardict)

Equivalent to `setvars(f, args...)` except that instead of running `f` 
immediately in the context of the modified `Var`s, yields a wrapper around `f`
that, when called, passes all arguments to `f`, which is run using the given
variable bindings in the dictionary `vardict`.

See also: [`Var`](@ref), [`@var`](@ref), [`vars`](@ref), [`withvars`](@ref),
[`setvars`](@ref), [`wrapwithvars`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{Int64}\\(@[0-9a-zA-Z]+: 0; init=0\\)"
julia> v = Var{Int}(0)
Var{Int64}(@F6ku5d8: 0; init=0)

julia> v[]
0

julia> vs = withvars(v => 2) do; vars() end; length(vs)
1

julia> f = wrapsetvars(vs) do x; x + v[] end; f(10)
12
```
"""
function wrapsetvars end
let tls_key = gensym("Air_var_bindings")
    # Private functions that depend on the tls_key.
    getbindings() = begin
        try
            return task_local_storage(tls_key)::VarsDict
        catch e
            isa(e, KeyError) || rethrow()
            d = VarsDict()
            task_local_storage(tls_key, d)
            return d
        end
    end
    global vars() = copy(getbindings())
    global withvars(f::Function, pairs::Vararg{VarPair,N}) where {N} = begin
        b0 = getbindings()
        b1 = copy(b0)
        for pair in pairs
            push!(b1, pair)
        end
        return setvars(f, b1)
    end
    global withvars(f::Function, vs::AbstractDict{<:Var,<:Any}) = begin
        b0 = getbindings()
        b1 = copy(b0)
        for pair in vs
            push!(b1, pair)
        end
        return setvars(f, b1)
    end
    global setvars(f::Function, vs::AbstractDict{<:Var,<:Any}) =
        setvars(f, VarsDict(vs))
    global setvars(f::Function, vs::VarsDict) = begin
        b0 = getbindings()
        task_local_storage(tls_key, vs)
        try
            return f()
        finally
            task_local_storage(tls_key, b0)
        end
    end
    global wrapwithvars(f::Function, pairs::Vararg{VarPair,N}) where {N} = begin
        b0 = getbindings()
        b1 = copy(b0)
        for pair in pairs
            push!(b1, pair)
        end
        f2 = (args...; kw...) -> let b0 = getbindings()
            task_local_storage(tls_key, b1)
            try
                return f(args...; kw...)
            finally
                task_local_storage(tls_key, b0)
            end
        end
        return f2
    end
    global wrapwithvars(f::Function, vs::AbstractDict{<:Var,<:Any}) = begin
        b0 = getbindings()
        b1 = copy(b0)
        for pair in vs
            push!(b1, pair)
        end
        f2 = (args...; kw...) -> let b0 = getbindings()
            task_local_storage(tls_key, b1)
            try
                return f(args...; kw...)
            finally
                task_local_storage(tls_key, b0)
            end
        end
        return f2
    end
    global wrapsetvars(f::Function, vs::VarsDict) = begin
        # We make a copy in case it changes between now and when the function
        # gets called
        vs = copy(vs) 
        f2 = (args...; kw...) -> let b0 = getbindings()
            task_local_storage(tls_key, vs)
            try
                return f(args...; kw...)
            finally
                task_local_storage(tls_key, b0)
            end
        end
        return f2
    end
    global Base.getindex(v::Var{T}) where {T} =
        get(getbindings(), v, v.initial_value)::T
    global Base.setindex!(v::Var, u) = begin
        msg = "$(typeof(v)) objects can only be set using withvars"
        throw(DomainError(v, msg))
    end
end
Base.show(io::IO, ::MIME"text/plain", d::Var{T}) where {T} = begin
    print(io, "$(typeof(d))(@")
    print(io, string(objectid(d), base=62))
    print(io, ": ")
    show(io, d[])
    print(io, "; init=")
    show(io, getfield(d, :initial_value))
    print(io, ")")
end
Base.show(io::IO, d::Var{T}) where {T} = begin
    print(io, "$(typeof(d))(@")
    print(io, string(objectid(d), base=62))
    print(io, ")")
end
"""
    @var

Convenient syntax for creating a task-local Var object:
`@var name = initval` will construct a `Var` object with the given initial
value. `@var name = initval::T` will create a `Var{T}` object.

See also: [`Var`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"Var{(Symbol|Any)}\\(@[0-9a-zA-Z]+: :start_sym; init=:start_sym\\)"
julia> @var v = :start_sym
Var{Symbol}(@JIm7aUS2sOl: :start_sym; init=:start_sym)

julia> @var u = :start_sym::Any
Var{Any}(@JIm7aUS2sOl: :start_sym; init=:start_sym)
```
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
export Var, @var, vars, withvars, setvars, wrapwithvars, wrapsetvars
