################################################################################
# PWDict.jl
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
const PWDictHeap{K,W,D} = PHeap{K,W,typeof(>),D} where {K,W<:Number}
macro _pwdict_code(name::Symbol, dicttype::Symbol)
    eq = gensym()
    h = gensym()
    return quote
        begin
            local $eq = equalfn($dicttype)
            local $h = hashfn($dicttype)
            struct $name{K,V,W<:Number} <: AbstractPWDict{K,V,W}
                heap::PWDictHeap{K,W,$dicttype{K,Int}}
                dict::$dicttype{K,V}
            end
            function $name{K,V,W}() where {K,V,W<:Number}
                return $name{K,V,W}(PWDictHeap{K,W,dicttype}(),
                                    dicttype{K,V}())
            end
            function $name{K,V,W}(d::$name{K,V,W}) where {K,V,W}
                return d
            end
            function $name{K,V,W}(ps::Union{Tuple,Pair}...) where {K,V,W<:Number}
                return reduce(push, ps, init=$name{K,V,W}())
            end
            function $name{K,V}(::Tuple{}) where {K,V,W<:Number}
                return $name{K,V,Float64}(PWDictHeap{K,Float64,dicttype}(),
                                          dicttype{K,V}())
            end
            function $name{K,V}() where {K,V,W<:Number}
                return $name{K,V,Float64}(PWDictHeap{K,Float64,dicttype}(),
                                          dicttype{K,V}())
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
                return $name{Any,Any,Float64}(PWDictHeap{Any,Float64,dicttype}(),
                                              dicttype{Any,Any}())
            end
            # Document the equal/hash types.
            equalfn(::Type{T}) where {T <: $name} = $eq
            hashfn(::Type{T}) where {T <: $name} = $h
            equalfn(::$name) = $eq
            hashfn(::$name) = $h
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
            Air.push(u::$name{K,V,W}, x::Union{Pair,Tuple}) where {K,V,W} = begin
                (k,v,w) = _to_pwtup(x)
                heap = push(u.heap, (k,w))
                dict = push(u.dict, k => v)
                (heap === u.heap && dict === u.dict) && return u
                return $name{K,V,W}(heap, dict)
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
                (lengt(u) > 0) || throw(ArgumentError("PWDict must be non-empty"))
                k = rand(u.heap)
                return k => u.dict[k]
            end
        end
    end |> esc
end
    
@_pwdict_code PWDict PDict
@_pwdict_code PWIdDict PIdDict
@_pwdict_code PWEqualDict PEqualDict

Base.empty(s::PWDict{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = PWDict{J,U,X}()
Base.empty(s::PWIdDict{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = PWIdDict{J,U,X}()
Base.empty(s::PWEqualDict{K,V,W}, ::Type{J}=K, ::Type{U}=V, ::Type{X}=W) where {K,V,W,J,U,X} = PWEqualDict{J,U,X}()

mutability(::Type{PWDict}) = Immutable
mutability(::Type{PWDict{K,V}}) where {K,V} = Immutable
mutability(::Type{PWIdDict}) = Immutable
mutability(::Type{PWIdDict{K,V}}) where {K,V} = Immutable
mutability(::Type{PWEqualDict}) = Immutable
mutability(::Type{PWEqualDict{K,V}}) where {K,V} = Immutable

_isequiv(::Immutable, ::Immutable, s::DS, t) where {DS <: AbstractPWDict} = false
_isequiv(::Immutable, ::Immutable, t, s::DS) where {DS <: AbstractPWDict} = false
_isequiv(::Immutable, ::Immutable, s::DS, t::DT) where {
    KT,VT,KS,VS,
    DS <: AbstractPWDict{KS,VS},
    DT <: AbstractPWDict{KT,VT}} = begin
    (length(t) == length(s)) || return false
    for (k,v) in s
        tt = get(t, k, t)
        (tt === t) && return false
        isequiv(tt, v) || return false
        (getweight(s, k) == getweight(t, k)) || return false
    end
    (equalfn(DS) === equalfn(DT)) && return true
    for (k,v) in t
        ss = get(s, k, s)
        (ss === s) && return false
        isequiv(ss, v) || return false
    end
    return true
end
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
_equivhash(::Immutable, u::PWD) where {K,V,W,PWD<:AbstractPWDict{K,V,W}} = begin
    h = length(u)
    for (k,v) in u.dict
        w = getweight(u.heap, k)
        h += equivhash(k) ⊻ equivhash(v) ⊻ equivhash(w)
    end
    return h
end
Base.hash(u::PWD) where {K,V,W,PWD<:AbstractPWDict{K,V,W}} = begin
    h = length(u)
    for (k,v) in u.dict
        w = getweight(u.heap, k)
        h += hash(k) ⊻ hash(v) ⊻ hash(w)
    end
    return h
end
