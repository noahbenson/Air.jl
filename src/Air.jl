################################################################################
# Air.jl
#
# The Air library.
# Functional collections and tools for Julia.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

module Air

include("util.jl")
include("api.jl")
include("ptree.jl")
include("pset.jl")
include("pdict.jl")
include("lazydict.jl")
include("parray.jl")
include("pheap.jl")
include("pwdict.jl")
include("pwset.jl")
include("vars.jl")

include("TX.jl")

"""
    @p{k1 => v1, k2 => v2, ...}
    @p[x1, x2, ...]
    @p[x1 x2...]
    @p(x1, x2, ...)

Yields a persistent data structure, depending on how the macro is called. This
is a shorthand for calling the various constructors directly. The following
expressions are equivalent:
* `@p{k1 => v1, k2 => v2, ...}` and `PDict(k1 => v1, k2 => v2, ...)`
* `@p[x1, x2, ...]` and `PVector([x1, x2, ...])`
* `@p[x1 x2 ...]` and `PArray([x1 x2 ...])`
* `@p(x1, x2, ...)` and `PSet([x1, x2, ...]`

See also: [`PDict`](@ref), [`PVector`](@ref), [`PArray`](@ref), [`PSet`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PDict{Symbol, ?Int64} with 3 entries:"
julia> @p{:a => 1, :b => 2, :c => 3}
PDict{Symbol,Int64} with 3 entries:
  :c => 3
  :a => 1
  :b => 2
```

```jldoctest; filter=r"4-element (PArray{Int64, ?1}|PVector{Int64}):"
julia> @p[1, 2, 3, 4]
4-element PArray{Int64,1}:
 1
 2
 3
 4
```

```jldoctest; filter=r"2×2 (PArray{Symbol, ?2}|PMatrix{Symbol}):"
julia> @p[:q2 :q1; :q3 :q4]
2×2 PArray{Symbol,2}:
 :q2  :q1
 :q3  :q4
```

```jldoctest
julia> @p(:a, :b, :a, :c)
PSet{Symbol} with 3 elements:
  :c
  :a
  :b
```
"""
macro p(exprs...)
    return esc(:(Air.PSet([$(exprs...)])))
end
macro p(expr::Expr)
    if expr.head == :braces
        return esc(:(Air.PDict($(expr.args...))))
    elseif expr.head == :vect
        return esc(:(Air.PVector([$(expr.args...)])))
    elseif expr.head == :vcat
        return esc(:(Air.PArray($expr)))
    elseif expr.head == :hcat
        return esc(:(Air.PArray($expr)))
    elseif expr.head == :tuple
        return esc(:(Air.PSet($expr)))
    else
        return esc(:(Air.PSet([$expr])))
    end
end

export Var, @var, Volatile, Actor, Source, tx, @tx, ReentrantRef,
    TransactionalRef, AbstractSourceKernel, getfilter, getfinalize,
    setfilter!, setfinalize!, send, geterror, receive, reset, @p

end 
