################################################################################
# PWSet.jl
#
# The Persistent weighted collection PWSet which is built on top of the
# persistent heap type in PHeap.jl transient heap type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020 Noah C. Benson


# ==============================================================================
# AbstractPWSet
"""
    AbstractPWSet{T,W}

AbstractPWSet is an abstract type implemented by all persistent weighted set
types.

See also: [`PWSet`](@ref), [`AbstractPSet`](@ref).
"""
abstract type AbstractPWSet{T,W<:Number} <: AbstractPSet{T} end


# ==============================================================================
# PWSet
macro _pwset_code(name::Symbol, dicttype::Symbol)
    return quote
        struct $name{T,W} <: AbstractPWSet{T,W}
            heap::PHeap{T,W,typeof(>),$dicttype{T,Int}}
        end
        function $name{T,W}(itr::AbstractArray) where {T,W<:Number}
            return $name{T,W}(itr...)
        end
        function $name{T,W}(itr::AbstractSet) where {T,W<:Number}
            return $name{T,W}(itr...)
        end
        function $name{T,W}(itr::AbstractDict) where {T,W<:Number}
            return $name{T,W}(itr...)
        end
        function $name{T,W}() where {T,W<:Number}
            return $name{T,W}(PHeap{T,W,typeof(>),$dicttype{T,Int}}(>))
        end
        function $name{T,W}(d::$name{T,W}) where {T,W}
            return d
        end
        function $name{T,W}(ps::Union{Tuple,Pair}...) where {T,W<:Number}
            return reduce(push, ps, init=$name{T,W}())
        end
        function $name{T}(itr) where {T}
            return $name{T}(itr...)
        end
        function $name{T}(::Tuple{}) where {T,W<:Number}
            return $name{T,Float64}(PHeap{T,Float64,typeof(>),$dicttype{T,Int}}(>))
        end
        function $name{T}() where {T,W<:Number}
            return $name{T,Float64}(PHeap{T,Float64,typeof(>),$dicttype{T,Int}}(>))
        end
        function $name{T}(d::$name{T,W}) where {T,W}
            return d
        end
        function $name{T}(ps::Union{Tuple,Pair}...) where {T}
            return $name{T,Float64}(ps...)
        end
        function $name(itr)
            return $name(itr...)
        end
        function $name(d::$name{T,W}) where {T,W}
            return d
        end
        function $name(ps::Union{Tuple,Pair}...)
            return $name{Any,Float64}(ps...)
        end
        function $name(::Tuple{})
            return $name{Any,Float64}()
        end
        function $name()
            return $name{Any,Float64}(PHeap{Any,Float64,typeof(>),$dicttype{T,Int}}(>))
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = equalfn($(dicttype))
        hashfn(::Type{T}) where {T <: $name} = hashfn($(dicttype))
        equalfn(::$name) = equalfn($(dicttype))
        hashfn(::$name) = hashfn($(dicttype))
        # Generic construction.
        # Base methods.
        Base.length(s::$name) = length(s.heap)
        Base.iterate(u::$name{T,W}) where {T,W} = iterate(u.heap)
        Base.iterate(u::$name{T,W}, x) where {T,W} = iterate(u.heap, x)
        Base.in(s::S, u::$name{T,W}, eqfn::Function) where {S,T,W} = in(s, u.heap, eqfn)
        Base.in(s::S, u::$name{T,W}) where {S,T,W} = in(s, u.heap)
        # And the Air methods.
        Air.push(u::$name{T,W}, x::Tuple{S,X}) where {T,W,S,X<:Number} = begin
            heap = push(u.heap, x)
            (heap === u.heap) && return u
            return $name{T,W}(heap)
        end
        Air.push(u::$name{T,W}, x::Pair{S,X}) where {T,W,S,X<:Number} = push(u, (x[1],x[2]))
        Air.pop(u::$name{T,W}) where {T,W} = begin
            (length(u) == 0) && throw(ArgumentError("$name must be non-empty"))
            return $name{T,W}(pop(u.heap))
        end
        Base.first(u::$name{T,W}) where {T,W} = begin
            (length(u) == 0) && throw(ArgumentError("$name must be non-empty"))
            return first(u.heap)
        end
        Air.delete(u::$name{T,W}, s::S) where {T,W,S} = begin
            heap = delete(u.heap, s)
            (heap === u.heap) && return u
            return $name{T,W}(heap)
        end
        Air.getweight(u::$name{T,W}, x::J) where {T,W,J} = getweight(u.heap, x)
        Air.setweight(u::$name{T,W}, x::J, w::X) where {T,W,J,X<:Number} = begin
            heap = setweight(u.heap, x, w)
            (heap === u.heap) && return u
            return $name{T,W}(heap)
        end
        Random.rand(u::$name{T,W}) where {T,W} = begin
            (length(u) > 0) || throw(ArgumentError("PWSet must be non-empty"))
            return rand(u.heap)
        end
    end |> esc
end
# Declare the types.    
@_pwset_code PWSet PDict
@_pwset_code PWIdSet PIdDict
# Document the types.
@doc """
    PWSet{K,V}

A persistent set with weighted elements. As such, a `PWSet` supports the
typical operations of a `PSet` as well as the following:
* The `first` function yields the element with the highest weight.
* The `pop` function yields a copy of the set without the element that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PSet`](@ref), [`getweight`](@ref), [`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> PWSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWSet{Symbol,Float64} with 3 elements:
  :c
  :b
  :a
```
""" PWSet
@doc """
    PWIdSet{T}

A persistent set with weighted elements. As such, a `PWIdSet` supports the
typical operations of a `PIdSet` as well as the following:
* The `first` function yields the element with the highest weight.
* The `pop` function yields a copy of the set without the element that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PIdSet`](@ref), [`PWSet`](@ref), [`getweight`](@ref),
[`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> PWIdSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWIdSet{Symbol,Float64} with 3 elements:
  :c
  :b
  :a
```
""" PWIdSet
# A few functions needed here instead of inside the macro.
Base.empty(s::PWSet{T,W}, ::Type{S}=T, ::Type{X}=W) where {T,W,S,X} = PWSet{S,X}()
Base.empty(s::PWIdSet{T,W}, ::Type{S}=T, ::Type{X}=W) where {T,W,S,X} = PWIdSet{S,X}()
Base.isequal(s::ST, t) where {ST<:AbstractPWSet} = false
Base.isequal(t, s::ST) where {ST<:AbstractPWSet} = false
Base.isequal(s::SS, t::ST) where {
    S,T,
    SS <: AbstractPWSet{S},
    ST <: AbstractPWSet{T}} = begin
    (length(t) == length(s)) || return false
    for x in s
        (x in t) || return false
        (getweight(s, x) == getweight(t, x)) || return false
    end
    (equalfn(SS) === equalfn(ST)) && return true
    for x in t
        (x in s) || return false
    end
    return true
end
Base.hash(u::PWS) where {T,W,PWS<:AbstractPWSet{T,W}} = begin
    h = length(u)
    for s in u.heap
        w = getweight(u.heap, s)
        h += hash(s) ⊻ hash(w)
    end
    return h
end

# Export the relevant symbols.
export PWSet, PWIdSet


