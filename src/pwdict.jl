################################################################################
# pwdict.jl
#
# The Persistent weighted collection PWDict which is built on top of the
# persistent heap type in PHeap.jl# transient heap type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020 Noah C. Benson


# ==============================================================================
# AbstractPWDict
"""
    AbstractPWDict

AbstractPWDict is a subtype of AbstractPDict that is overloaded by classes that
represent weighted persistent dictionaries.

See also: [`AbstractPDict`](@ref), [`PWDict`](@ref), [`AbstractPWSet`](@ref).
"""
abstract type AbstractPWDict{K,V,W<:Number} <: AbstractPDict{K,V} end


# ==============================================================================
# PWDict
_to_pwtup(p::Tuple{K,V,W}) where {K,V,W<:Number} = p
_to_pwtup(p::Tuple{K,Tuple{V,W}}) where {K,V,W<:Number} = begin
    (k,(v,w)) = p
    return (k,v,w)
end
_to_pwtup(p::Pair{K,Tuple{V,W}}) where {K,V,W<:Number} = begin
    (k,(v,w)) = p
    return (k,v,w)
end
_to_pwtup(p::Pair{K,Pair{V,W}}) where {K,V,W<:Number} = begin
    (k,(v,w)) = p
    return (k,v,w)
end
_to_pwtup(p::Tuple{Tuple{K,V},W}) where {K,V,W<:Number} = begin
    ((k,v),w) = p
    return (k,v,w)
end
_to_pwtup(p::Tuple{Pair{K,V},W}) where {K,V,W<:Number} = begin
    ((k,v),w) = p
    return (k,v,w)
end
_to_pwtup(p::Pair{Pair{K,V},W}) where {K,V,W<:Number} = begin
    ((k,v),w) = p
    return (k,v,w)
end
const PWDictHeap{K,W,D} = PHeap{K,W,typeof(>),D} where {K,W<:Number,D<:AbstractPDict{K,Int}}
macro _pwdict_code(name::Symbol, dicttype::Symbol)
    eq = gensym()
    h = gensym()
    return quote
        struct $name{K,V,W<:Number} <: AbstractPWDict{K,V,W}
            heap::PWDictHeap{K,W,$dicttype{K,Int}}
            dict::$dicttype{K,V}
        end
        function $name{K,V,W}() where {K,V,W<:Number}
            return $name{K,V,W}(PWDictHeap{K,W,$dicttype{K,Int}}(>),
                                $dicttype{K,V}())
        end
        function $name{K,V,W}(d::$name{K,V,W}) where {K,V,W}
            return d
        end
        function $name{K,V,W}(ps::Union{Tuple,Pair}...) where {K,V,W<:Number}
            return reduce(push, ps, init=$name{K,V,W}())
        end
        function $name{K,V}(::Tuple{}) where {K,V}
            return $name{K,V,Float64}(PWDictHeap{K,Float64,$dicttype{K,Int}}(>),
                                      $dicttype{K,V}())
        end
        function $name{K,V}() where {K,V}
            return $name{K,V,Float64}(PWDictHeap{K,Float64,$dicttype{K,Int}}(>),
                                      $dicttype{K,V}())
        end
        function $name{K,V}(d::$name{K,V,W}) where {K,V,W<:Number}
            return d
        end
        function $name{K,V}(ps::Union{Tuple,Pair}...) where {K,V}
            return $name{K,V,Float64}(ps...)
        end
        function $name(d::$name{K,V,W}) where {K,V,W<:Number}
            return d
        end
        function $name(ps::Union{Tuple,Pair}...)
            return $name{Any,Any,Float64}(ps...)
        end
        function $name(::Tuple{})
            return $name{Any,Any,Float64}()
        end
        function $name()
            return $name{Any,Any,Float64}(PWDictHeap{Any,Float64,$dicttype{Any,Int}}(),
                                          $dicttype{Any,Any}())
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = equalfn($dicttype)
        hashfn(::Type{T}) where {T <: $name} = hashfn($dicttype)
        equalfn(::$name) = equalfn($dicttype)
        hashfn(::$name) = hashfn($dicttype)
        # Generic construction.
        # Base methods.
        #Base.empty(s::$name{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = $name{J,U,X}()
        Base.length(s::$name) = length(s.heap)
        Base.iterate(u::$name{K,V,W}) where {K,V,W} = begin
            (length(u) == 0) && return nothing
            return iterate(u, u.heap)
        end
        Base.iterate(u::$name{K,V,W}, x) where {K,V,W} = begin
            nxt = iterate(u.heap, x)
            (nxt === nothing) && return nothing
            (k,b) = nxt
            return (k => u.dict[k], b)
        end
        Base.get(u::$name{K,V,W}, k, df) where {K,V,W} = get(u.dict, k, df)
        Base.in(kv::Pair, u::$name{K,V,W}, eqfn::Function) where {K,V,W} = in(kv, u.dict, eqfn)
        Base.in(kv::Pair, u::$name{K,V,W}) where {K,V,W} = in(kv, u.dict)
        # And the Air methods.
        Air.push(u::$name{K,V,W}, x::Pair{J,U}) where {K,V,W,J,U} = begin
            (k,v,w) = _to_pwtup(x)
            heap = push(u.heap, (k,w))
            dict = push(u.dict, k => v)
            (heap === u.heap && dict === u.dict) && return u
            return $name{K,V,W}(heap, dict)
        end
        Air.push(u::$name{K,V,W}, x::Tuple) where {K,V,W,J,U} = begin
            (k,v,w) = _to_pwtup(x)
            heap = push(u.heap, (k,w))
            dict = push(u.dict, k => v)
            (heap === u.heap && dict === u.dict) && return u
            return $name{K,V,W}(heap, dict)
        end
        Air.pop(u::$name{K,V,W}) where {K,V,W} = begin
            (length(u) == 0) && throw(ArgumentError("PWDict must be non-empty"))
            x = first(u.heap)
            dict = delete(u.dict, x)
            return $name{K,V,W}(pop(u.heap), dict)
        end
        Air.first(u::$name{K,V,W}) where {K,V,W} = begin
            (length(u) == 0) && throw(ArgumentError("PWDict must be non-empty"))
            k = first(u.heap)
            return k => u.dict[k]
        end
        Air.setindex(u::$name{K,V,W}, v, k::J) where {K,V,W,J} = begin
            (get(u.heap._index, k, 0) == 0) && throw(
                ArgumentError("when inserting into PWDict, value and weight are required"))
            d = setindex(u.dict, v, k)
            (d === u.dict) && return u
            return $name{K,V,W}(u.heap, d)
        end
        Air.setindex(u::$name{K,V,W}, v::Tuple{VV,WW}, k::J) where {K,V,W,VV,WW<:Number,J} = begin
            return push(u, k => v)
        end
        Air.setindex(u::$name{K,V,W}, v::Pair{VV,WW}, k::J) where {K,V,W,VV,WW<:Number,J} = begin
            return push(u, k => v)
        end
        Air.delete(u::$name{K,V,W}, k::J) where {K,V,W,J} = begin
            heap = delete(u.heap, k)
            (heap === u.heap) && return u
            dict = delete(u.dict, k)
            return $name{K,V,W}(heap, dict)
        end
        Air.getweight(u::$name{K,V,W}, x::J) where {K,V,W,J} = getweight(u.heap, x)
        Air.setweight(u::$name{K,V,W}, x::J, w::X) where {K,V,W,J,X<:Number} = begin
            heap = setweight(u.heap, x, w)
            (heap === u.heap) && return u
            return $name{K,V,W}(heap, u.dict)
        end
        Random.rand(u::$name{K,V,W}) where {K,V,W} = begin
            (length(u) > 0) || throw(ArgumentError("PWDict must be non-empty"))
            k = rand(u.heap)
            return k => u.dict[k]
        end
    end |> esc
end
# Declare the types.    
@_pwdict_code PWDict PDict
@_pwdict_code PWIdDict PIdDict
# And document them.
@doc """
    PWDict{K,V}

A persistent dictionary with weighted pairs. As such, a `PWDict` supports the
typical operations of a `PDict` as well as the following:
* The `first` function yields the key-value pair with the highest weight.
* The `pop` function yields a copy of the dictionary without the pair that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PDict`](@ref), [`getweight`](@ref), [`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PWDict{Symbol, ?Int64, ?Float64} with 3 entries:"
julia> PWDict{Symbol,Int}(:a => (1,0.1), :b => (2, 0.2), :c => (3, 0.3))
PWDict{Symbol,Int64,Float64} with 3 entries:
  :c => 3
  :b => 2
  :a => 1
```
""" PWDict
@doc """
    PWIdDict{K,V}

A persistent dictionary with weighted pairs. As such, a `PWIdDict` supports the
typical operations of a `PIdDict` as well as the following:
* The `first` function yields the key-value pair with the highest weight.
* The `pop` function yields a copy of the dictionary without the pair that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PIdDict`](@ref), [`PWDict`](@ref), [`getweight`](@ref),
[`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PWIdDict{Symbol, ?Int64, ?Float64} with 3 entries:"
julia> PWIdDict{Symbol,Int}(:a => (1,0.1), :b => (2, 0.2), :c => (3, 0.3))
PWIdDict{Symbol,Int64,Float64} with 3 entries:
  :c => 3
  :b => 2
  :a => 1
```
""" PWIdDict
# For some reason, these need to be out of the macro above.
Base.empty(s::PWDict{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = PWDict{J,U,X}()
Base.empty(s::PWIdDict{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = PWIdDict{J,U,X}()
# Equality functions also.
Base.isequal(s::DS, t) where {DS <: AbstractPWDict} = false
Base.isequal(t, s::DS) where {DS <: AbstractPWDict} = false
Base.isequal(s::DS, t::DT) where {
    KT,VT,KS,VS,
    DS <: AbstractPWDict{KS,VS},
    DT <: AbstractPWDict{KT,VT}} = begin
    (length(t) == length(s)) || return false
    for (k,v) in s
        tt = get(t, k, t)
        (tt === t) && return false
        isequal(tt, v) || return false
        (getweight(s, k) == getweight(t, k)) || return false
    end
    (equalfn(DS) === equalfn(DT)) && return true
    for (k,v) in t
        ss = get(s, k, s)
        (ss === s) && return false
        isequal(ss, v) || return false
    end
    return true
end
Base.hash(u::PWD) where {K,V,W,PWD<:AbstractPWDict{K,V,W}} = begin
    h = length(u)
    for (k,v) in u.dict
        w = getweight(u.heap, k)
        h += hash(k) ⊻ hash(v) ⊻ hash(w)
    end
    return h
end

# Export the relevant symbols.
export PWDict, PWIdDict
