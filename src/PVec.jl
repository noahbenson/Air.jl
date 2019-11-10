################################################################################
# PVec
# The PVec type is a simple small immutable vector type that is formed out of
# julia structs. In this sense they are 100% immutable, unlike alternative
# implementations of persistent vector types backed by a julia array.
# The PVec type is intended as 

import Base.getindex, Base.setindex!, Base.length, Base.size, Base.iterate,
       Base.show, Base.axes, Base.similar, Base.similar, Base.convert
import Printf, Printf.@sprintf

# PVec types are based on this abstract type:
abstract type PVec{T} <: AbstractArray{T,1} end
# How we print pvecs:
_print_pvec(io::IO, u::PVec{T}, head::String) where {T} = begin
    print(io, "$(head)[")
    n = length(u)
    if n < 50
        for ii in 1:n
            if ii > 1 print(io, ", $(repr(u[ii]))") else print(io, "$(repr(u[ii]))") end
        end
    else
        for ii in 1:20
            if ii > 1 print(io, ", $(repr(u[ii]))") else print(io, "$(repr(u[ii]))") end
        end
        print(io, " ... ")
        for ii in n-20:n
            if ii < n print(io, "$(repr(u[ii])), ") else print(io, "$(repr(u[ii]))") end
        end
    end
    print(io, "]")
end
Base.show(io::IO, ::MIME"text/plain", pv::PVec{T}) where {T} =
    _print_pvec(io, pv, "PVec")
Base.setindex!(u::PVec{T}, v::U, args...) where {T,U<:T} = error(
    "setindex! immutable type PVec cannot be changed")
Base.getindex(u::PVec, k1::Integer, k2::Integer, args::Vararg{<:Integer}) =
    getindex(getindex(u, k1), k2, args...)
_iseq(t::PVec{T}, s::AbstractArray{S,1}, eqfn::Function) where {T, S} = begin
    (length(t) == length(s)) || return false
    for (it,is) in zip(t, s)
        eqfn(it, is) || return false
    end
    return true
end
isequiv(t::PVec{T}, s::AbstractArray{S,1}) where {T, S} = _iseq(t, s, iequiv)
isequiv(s::AbstractArray{S,1}, t::PVec{T}) where {T, S} = _iseq(t, s, iequiv)
Base.isequal(t::PVec{T}, s::AbstractArray{S,1}) where {T, S} = _iseq(t, s, isequal)
Base.isequal(s::AbstractArray{S,1}, t::PVec{T}) where {T, S} = _iseq(t, s, isequal)

#Base.getindex(u::PVec, k::Integer) = throw(BoundsError(u, [k]))
function Base.iterate(u::PVec{T}, k::Integer)::Union{Nothing,Tuple{T,Integer}} where {T}
    if 1 <= k <= length(u)
        return (u[k], k+1)
    else
        return nothing
    end
end
function Base.iterate(u::PVec{T})::Union{Nothing,Tuple{T,Integer}} where {T}
    return iterate(u, 1)
end
Base.size(u::PVec) = (length(u),)
#TODO: axes(::PVec) (and how to handle arrays?)
# The various persistent operations:
dissoc(u::PVec, k) = u
remove(u::PVec, k) = throw(ArgumentError("invalid index: $k of type $(typeof(k))"))
struct DissocError <: Exception
    object::Any
    key::Any
end
Base.showerror(io::IO, e::DissocError) = let k=repr(e.key), t=repr(typeof(e.object))
    print(io, "cannot dissociate key $k from association of type $t")
end
# dissoc and remove can be defined generically for PVec:
dissoc(u::PVec{T}, k::Integer) where {T} = let n = length(u)
    if k == n
        return pop(u)
    elseif k > n || k < 1
        return u
    else
        throw(DissocError(u, k))
    end
end
remove(u::PVec{T}, k::Integer) where {T} = let n = length(u)
    if k == n
        return pop(u)
    elseif k > n || k < 1
        throw(BoundsError(u, [k]))
    elseif Missing <: T
        return assoc(u, k, missing)
    else
        return assoc(convert(PVec{Union{Missing,T}}, u), k, missing)
    end
end
# convert can also be defined...
Base.convert(::Type{PVec{T}}, u::AbstractArray{S,1}) where {T, S} = PVec{T}(
    [convert(T, k) for k in u])
# we can build slice-indexing out of what we have also:
struct PSubVec{T,I} <: PVec{T} where {I <: Integer}
    _vector::PVec{T}
    _index::AbstractArray{I,1}
    function PSubVec{T,I}(u::PVec{T}, ii::AbstractArray{I,1})::PSubVec{T,I} where {T, I<:Integer}
        let n = length(u)
            # we need to check every element to make sure it's ok...
            for k in ii
                if k < 1 || k > n
                    throw(BoundsError(u, [k]))
                end
            end
            # okay, they're all fine...
            return new(u, PVec{I}(ii))
        end
    end
    # For certain instances we can do less work...
    function PSubVec{T,I}(u::PVec{T}, ii::StepRange{I,J})::PSubVec{T,I} where {T, I<:Integer, J<:Integer}
        let n = length(u), f = ii.start, l = ii.stop, s = ii.step
            if (f > l && s > 0) || (f < l && s < 0)
                return new(PVec0{T}(), PVec0{I}())
            elseif f < 1 || f > n
                throw(BoundsError(u, [f]))
            elseif l < 1 || l > n
                throw(BoundsError(u, [l]))
            else
                return new(u, ii)
            end
        end
    end
    function PSubVec{T,I}(u::PVec{T}, ii::UnitRange{I})::PSubVec{T,I} where {T, I<:Integer}
        let n = length(u), f = ii.start, l = ii.stop
            if f > l
                return new(PVec0{T}(), PVec0{I}())
            elseif f < 1 || f > n
                throw(BoundsError(u, [f]))
            elseif l < 1 || l > n
                throw(BoundsError(u, [l]))
            else
                return new(u, ii)
            end
        end
    end
end
Base.length(u::PSubVec) = length(u._index)
Base.getindex(u::PSubVec{T,I}, k::Integer) where {T,I<:Integer} = u._vector[u._index[k]]
# getindex for any pvec and slice/subarray:
Base.getindex(u::PVec{T}, ii::AbstractArray{I,1}) where {T, I<:Integer} = PSubVec{T,I}(u, ii)
push(u::PSubVec{T,I}, x::S) where {T, S<:T, I<:Integer} = push(PVec{T}(u), x)
pop(u::PSubVec{T,I}) where {T, I<:Integer} = pop(PVec{T}(u))
assoc(u::PSubVec{T,I}, x::S) where {T, S<:T, I<:Integer} = assoc(PVec{T}(u), x)
dissoc(u::PSubVec{T,I}, k::J) where {T, I<:Integer, J<:Integer} = dissoc(PVec{T}(u), x)

################################################################################
# #PSmallVec
# The small vectors are for finite-sized vectors of length 0-32. After size 32,
# we switch over to PBigVec objects, which use PSmallVec objects to store larger
# vectors in a nested fashion.

# The simplest concrete type is the empty vector
struct PVec0{T} <: PVec{T} end
const _PVEC0_SIZE = (0,)
PVec0(u::Vector{T}) where {T} = PVec0{T}()
PVec0{T}(u::Vector{T}) where {T} = PVec0{T}()
# basic functions for the empty PVec:
Base.length(::PVec0) = 0
Base.size(::PVec0) = _PVEC0_SIZE
Base.in(::PVec0, x) = false
push(::PVec0{T}, x::S) where {T, S<:T} = PVec1{T}(x)
pop(e::PVec0) = e
assoc(::PVec0{T}, ::Type{Val{1}}, v::S) where {T, S<:T} = PVec1{T}(v)

# we create other pvec types with this macro:
abstract type PSmallVec{T} <: PVec{T} end
macro psmallvectype(name::Symbol, n, lname::Symbol, uname::Symbol)
    return esc(
        quote
            # Declare the struct:
            struct $name{T} <: PSmallVec{T}
                _elements::NTuple{$n, T}
            end
            function $name{T}(u::AbstractArray{S,1}) where {T,S<:T}
                (length(u) >= $n) || error("wrong size given to small pvec type")
                return $name{T}(NTuple{$n, T}(u))
            end
            function push(u::$name{T}, v::S) where {T, S<:T}
                $(uname === :PBigVec
                  ? :($uname{T}(33, PVec2{PSmallVec{T}}((u, PVec1{T}((v,))))))
                  : :($uname{T}(NTuple{$(n+1),T}([u._elements..., v]))))
            end
            function pop(u::$name{T}) where {T}
                return $lname{T}(u._elements[1:end-1])
            end
            function assoc(u::$name{T}, k::Integer, v::S) where {T, S<:T}
                if k == $(n+1)
                    return push(u, v)
                elseif k > $(n+1) || k < 1
                    throw(BoundsError(u, [k]))
                else
                    a = T[u._elements...]
                    a[k] = v
                    return $name{T}(a)
                end
            end
            # The getindex operator...
            function Base.getindex(u::$name{T}, k::Integer) where {T}
                if k < 1 || k > $n
                    throw(BoundsError(u, [k]))
                else
                    return u._elements[k]
                end
            end
            function Base.length(u::$name{T}) where {T}
                return $n
            end
            function Base.in(x::S, u::$name{T}, eqfn::Function) where {T, S<:T}
                return in(x, u._elements, eqfn)
            end
            function Base.in(x::S, u::$name{T}) where {T, S<:T}
                return in(x, u._elements)
            end
            # The in operator...
            # equality and equivalence
            #$(esc(:(
            function _iseq(t::$name{T}, s::$name{S}, eqfn::Function) where {S,T}
                for (a,b) in zip(t._elements, s._elements)
                    eqfn(a, b) || return false
                end
                return true
            end
            function _iseq(t::$name{T}, s::AbstractArray{S,1}, eqfn::Function) where {S,T}
                (length(s) == $n) || return false
                for (a,b) in zip(t._elements, s)
                    eqfn(a, b) || return false
                end
                return true
            end
            function equivhash(t::$name{T}) where {T}
                let h = objectid(T)
                    for (k,x) in enumerate(t)
                        h += equivhash(x) * k
                    end
                    return h
                end
            end
            # Arithmetic operatrs
            function Base.map(f::Function, u::$name{T}) where {T}
                return $name(map(f, u._elements))
            end
            function Base.:+(a::$name{T}) where {T}
                return a
            end
            function Base.:-(a::$name{T}) where {T}
                return $name{T}([-x for x in a._elements])
            end
        end)
end

# Okay, make pvec types for 1-32:
@psmallvectype PVec1  1  PVec0  PVec2 
@psmallvectype PVec2  2  PVec1  PVec3 
@psmallvectype PVec3  3  PVec2  PVec4 
@psmallvectype PVec4  4  PVec3  PVec5 
@psmallvectype PVec5  5  PVec4  PVec6 
@psmallvectype PVec6  6  PVec5  PVec7 
@psmallvectype PVec7  7  PVec6  PVec8 
@psmallvectype PVec8  8  PVec7  PVec9 
@psmallvectype PVec9  9  PVec8  PVec10
@psmallvectype PVec10 10 PVec9  PVec11
@psmallvectype PVec11 11 PVec10 PVec12
@psmallvectype PVec12 12 PVec11 PVec13
@psmallvectype PVec13 13 PVec12 PVec14
@psmallvectype PVec14 14 PVec13 PVec15
@psmallvectype PVec15 15 PVec14 PVec16
@psmallvectype PVec16 16 PVec15 PVec17
@psmallvectype PVec17 17 PVec16 PVec18
@psmallvectype PVec18 18 PVec17 PVec19
@psmallvectype PVec19 19 PVec18 PVec20
@psmallvectype PVec20 20 PVec19 PVec21
@psmallvectype PVec21 21 PVec20 PVec22
@psmallvectype PVec22 22 PVec21 PVec23
@psmallvectype PVec23 23 PVec22 PVec24
@psmallvectype PVec24 24 PVec23 PVec25
@psmallvectype PVec25 25 PVec24 PVec26
@psmallvectype PVec26 26 PVec25 PVec27
@psmallvectype PVec27 27 PVec26 PVec28
@psmallvectype PVec28 28 PVec27 PVec29
@psmallvectype PVec29 29 PVec28 PVec30
@psmallvectype PVec30 30 PVec29 PVec31
@psmallvectype PVec31 31 PVec30 PVec32
@psmallvectype PVec32 32 PVec31 PBigVec

const _psmallvec_bysize = PVec32{Type}(
    [PVec1,  PVec2,  PVec3,  PVec4,  PVec5,  PVec6,
     PVec7,  PVec8,  PVec9,  PVec10, PVec11, PVec12,
     PVec13, PVec14, PVec15, PVec16, PVec17, PVec18,
     PVec19, PVec20, PVec21, PVec22, PVec23, PVec24,
     PVec25, PVec26, PVec27, PVec28, PVec29, PVec30,
     PVec31, PVec32])
PSmallVec{T}(n::Integer, u::AbstractArray{S,1}) where {T,S<:T} = let nu = length(u)
    if     n <= 0  return PVec0{T}()
    elseif n > 32  error("PSmallVec given vector of size > 32")
    elseif n > nu  error("PSmallVec given incorrect size")
    elseif n == nu return _psmallvec_bysize[n]{T}(u)
    else           return _psmallvec_bysize[n]{T}(u[1:n])
    end
end
PSmallVec{T}(u::AbstractArray{S,1}) where {T,S<:T} = PSmallVec{T}(length(u), u)

# next we need vectors that can be composed of smaller vectors...
struct PBigVec{T} <: PVec{T}
    _n::Integer
    _elements::PVec{PSmallVec{T}}
end
_pbigvec_iisplit(ii::Integer) = let kk=ii-1
    ((kk >> 5) + 1, (kk & 0b11111) + 1)
end
Base.length(u::PBigVec) = u._n
Base.getindex(u::PBigVec{T}, k::Integer) where {T, P<:PVec} = begin
    (k < 1 || k > u._n) && throw(BoundsError(u, [k]))
    let (k1,k2) = _pbigvec_iisplit(k)
        return u._elements[k1][k2]
    end
end
Base.in(u::PBigVec{T}, x::S) where {T, S<:T} = begin
    for y in u._elements
        if x in u return true end
    end
    return false
end
push(u::PBigVec{T}, v::S) where {T, S<:T} = begin
    let n0 = length(u), n = n0+1, (k1,k2) = _pbigvec_iisplit(n)
        return PBigVec{T}(n,
                          k1 > length(u._elements)
                          ? push(u._elements, PVec1{T}(v))
                          : assoc(u._elements, k1, push(u._elements[k1], v)))
    end
end
pop(u::PBigVec{T}) where {T} = begin
    u._n < 33 && return pop(u._elements[1])
    let (k1,k2) = _pbigvec_iisplit(u._n), ee = u._elements[k1]
        return PBigVec{T}(u._n - 1,
                          length(ee) == 1
                          ? pop(u._elements)
                          : assoc(u._elements, k1, pop(ee)))
    end
end
assoc(u::PBigVec{T}, k::Integer, v::S) where {T, S<:T} = begin
    k == u._n + 1 && return push(u, v)
    (k > u._n || k < 1) && throw(BoundsError(u, [k]))
    let (k1,k2) = _pbigvec_iisplit(k), el = u._elements[k1]
        return PBigVec{T}(u._n, assoc(u._elements, k1, assoc(el, k2, v)))
    end
end

PVec{T}(u::AbstractArray{S,1}) where {T,S<:T} = let n = length(u)
    if typeof(u) <: PVec  return u
    elseif n <= 0  return PVec0{T}()
    elseif n <= 32 return _psmallvec_bysize[n]{T}(NTuple{n,T}(u))
    else
        let parts = [PSmallVec{T}(u[ii:min(ii+32-1, n)]) for ii in 1:32:n]
            return PBigVec{T}(n, PVec{PSmallVec{T}}(parts))
        end
    end
end
PVec(u::AbstractArray{T, 1}) where {T} = PVec{T}(u)

# The important thing in the operations below is to operate in chunks whenever
# possible.
function Base.map(f::Function, u::PBigVec{T}) where {T}
    let v   = [map(f, x) for x in u._elements],
        tts = [typeof(x).parameters[1] for x in v],
        U   = typejoin(tts...),
        v   = [(tt.parameters[1] === U ? x : convert(PVec{U}, x))
               for (x,tt) in zip(v,tts)]
        return PBigVec{U}(u._n, PVec{PSmallVec{U}}(v))
    end
end
function Base.map(f::Function, uu::Vararg{PBigVec{T}}) where {T}
    let v   = [map(f, args...) for args in zip(uu...)],
        tts = [typeof(x).parameters[1] for x in v],
        U   = typejoin(tts...),
        v   = [(tt.parameters[1] === U ? x : convert(PVec{U}, x))
               for (x,tt) in zip(v,tts)]
        return PBigVec{U}(u._n, PVec{PSmallVec{U}}(v))
    end
end
function Base.broadcast(f::Function, u::PBigVec{T}) where {T}
    let v   = [map(f, x) for x in u._elements],
        tts = [typeof(x).parameters[1] for x in v],
        U   = typejoin(tts...),
        v   = [(tt.parameters[1] === U ? x : convert(PVec{U}, x))
               for (x,tt) in zip(v,tts)]
        return PBigVec{U}(u._n, PVec{PSmallVec{U}}(v))
    end
end
function Base.broadcast(f::Function, uu::Vararg{PBigVec{T}}) where {T}
    let v   = [map(f, args...) for args in zip([u._elements for u in uu]...)],
        tts = [typeof(x).parameters[1] for x in v],
        U   = typejoin(tts...),
        v   = [(tt.parameters[1] === U ? x : convert(PVec{U}, x))
               for (x,tt) in zip(v,tts)]
        return PBigVec{U}(u._n, PVec{PSmallVec{U}}(v))
    end
end
Base.:+(a::PBigVec{T}) where {T} = a
Base.:-(a::PBigVec{T}) where {T} = begin
    return PBigVec{T}(a._n, PVec{PSmallVec{T}}([-x for x in a._elements]))
end
function Base.:+(a::PBigVec{T}, b::PBigVec{S}) where {T,S}
    let ss  = [aa + bb for (aa,bb) in zip(a._elements, b._elements)],
        tts = [typeof(s).parameters[1] for s in ss],
        U   = typejoin(tts...),
        ss  = [(tt === U ? s : convert(PVec{U}, s))
               for (s,tt) in zip(ss, tts)]
        return PBigVec{U}(a._n, PVec{PSmallVec{U}}(ss))
    end
end
function Base.:-(a::PBigVec{T}, b::PBigVec{S}) where {T,S}
    let ss  = [aa - bb for (aa,bb) in zip(a._elements, b._elements)],
        tts = [typeof(s).parameters[1] for s in ss],
        U   = typejoin(tts...),
        ss  = [(tt === U ? s : convert(PVec{U}, s))
               for (s,tt) in zip(ss, tts)]
        return PBigVec{U}(a._n, PVec{PSmallVec{U}}(ss))
    end
end

