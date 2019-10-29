################################################################################
# PVec
# The PVec type is a simple small immutable vector type that is formed out of
# julia structs. In this sense they are 100% immutable, unlike alternative
# implementations of persistent vector types backed by a julia array.
# The PVec type is intended as 

import Base.getindex, Base.setindex!, Base.length, Base.size, Base.iterate,
       Base.show, Base.axes, Base.similar, Base.similar, Base.convert
import Printf.@sprintf

# PVec types are based on this abstract type:
abstract type PVec{T} <: AbstractArray{T,1} end
# How we print pvecs:
_print_pvec(io::IO, u::PVec{T}, head::String) where {T} = begin
    print(io, "$(head)[")
    n = length(u)
    if n < 50
        for ii in 1:n
            if ii > 1 print(io, ", $(u[ii])") else print(io, "$(u[ii])") end
        end
    else
        for ii in 1:20
            if ii > 1 print(io, ", $(u[ii])") else print(io, "$(u[ii])") end
        end
        print(io, " ... ")
        for ii in n-20:n
            if ii < n print(io, "$(u[ii]), ") else print(io, "$(u[ii])") end
        end
    end
    print(io, "]")
end
Base.show(io::IO, ::MIME"text/plain", pv::PVec{T}) where {T} =
    _print_pvec(io, pv, "PVector")
Base.setindex!(u::PVec{T}, v::U, args...) where {T,U<:T} = error(
    "setindex! immutable type PVec cannot be changed")
Base.getindex(u::PVec, k1::Integer, k2::Integer, args::Vararg{<:Integer}) =
    getindex(getindex(u, k1), k2, args...)
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
#TODO: similar(::PVec)
#TODO: axes(::PVec) (and how to handle arrays?)
# The various persistent operations:
assoc(u::PVec, k, v) = throw(ArgumentError("invalid index: $k of type $(typeof(k))"))
assoc(u::PVec, k::Integer, v) = throw(BoundsError(u, [k]))
dissoc(u::PVec, k) = u
remove(u::PVec, k) = throw(ArgumentError("invalid index: $k of type $(typeof(k))"))
remove(u::PVec, k::Integer) = throw(BoundsError(u, [k]))
# dissoc and remove can be defined generically for PVec:
dissoc(u::PVec{T}, k::Integer) where {T} = let n = length(u)
    if k == n
        return pop(u)
    elseif k > n || k < 1
        return u
    elseif Missing <: T
        return assoc(u, k, missing)
    else
        return assoc(convert(PVec{Union{Missing,T}}, u), k, missing)
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
    function PSubVec{T,I}(u::PVec{T}, ii::StepRange{I})::PSubVec{T,I} where {T, I<:Integer}
        let n = length(u), f = ii[1], l = ii[end]
            if f < 1 || f > n
                throw(BoundsError(u, [f]))
            elseif l < 1 || l > n
                throw(BoundsError(u, [l]))
            else
                return new(u, ii)
            end
        end
    end
    function PSubVec{T,I}(u::PVec{T}, ii::UnitRange{I})::PSubVec{T,I} where {T, I<:Integer}
        let n = length(u), f = ii[1], l = ii[end]
            if f < 1 || f > n
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
push(::PVec0{T}, x::S) where {T, S<:T} = PVec1{T}(x)
pop(e::PVec0) = e
assoc(::PVec0{T}, ::Val{1}, v::S) where {T, S<:T} = PVec1{T}(v)

# we create other pvec types with this macro:
abstract type PSmallVec{T} <: PVec{T} end
macro psmallvectype(name::Symbol, n, lname::Symbol, uname::Symbol)
    let vals = [Symbol("_element_" * @sprintf("%03d", k)) for k in 1:n]
        return quote
            # Declare the struct:
            struct $(esc(name)){$(esc(:T))} <: PSmallVec{$(esc(:T))}
                $([:($(val)::$(esc(:T))) for val in vals]...)
            end
            $(esc(:push))(u::$(esc(name)){T}, v::S) where {T, S<:T} =
                $(uname === :PBigVec
                  ? :($(esc(uname)){T}(33, PVec2{PSmallVec{T}}(u, PVec1{T}(v))))
                  : :($(esc(uname)){T}($([:(u.$(esc(vals[k]))) for k in 1:n]...), v)))
            $(esc(:pop))(u::$(esc(name)){T}) where {T} = $(esc(lname))(
                $([:(u.$(esc(vals[k]))) for k in 1:n-1]...))
            $([quote
               $(esc(:assoc))($(esc(:u))::$(esc(name)){T}, ::$(esc(:(Base.Type))){$(esc(:(Base.Val))){$k}}, $(esc(:v))::S) where {T, S<:T} = $(esc(name)){T}(
                   $((let q = [esc(:(u.$(vals[kk]))) for kk in 1:n]
                          q[k] = esc(:v)
                          q
                      end)...))
               end
               for k in 1:n]...)
            $(esc(:assoc))(u::$(esc(name)){T}, ::$(esc(:(Base.Type))){$(esc(:(Base.Val))){$(n+1)}}, v::S) where {T, S<:T} = push(u, v)
            $(esc(:assoc))(u::$(esc(name)){T}, k::Integer, v::S) where {T, S<:T} = assoc(u, $(esc(Base.Val)){k}, v)
                                   
            # The getindex operator...
            $([:(_pvec_getindex($(esc(:x))::$(esc(name)), ::$(esc(:(Base.Type))){$(esc(:(Base.Val))){$k}}) = $(esc(:(x.$(vals[k]))))) for k in 1:n]...)
            $(esc(:(Base.getindex)))(x::$(esc(name)), ii::Integer) = _pvec_getindex(x, $(esc(:(Base.Val))){ii})
            $(esc(:(Base.length)))(x::$(esc(name))) = $n
        end
    end
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
    PVec1,  PVec2,  PVec3,  PVec4,  PVec5,  PVec6,
    PVec7,  PVec8,  PVec9,  PVec10, PVec11, PVec12,
    PVec13, PVec14, PVec15, PVec16, PVec17, PVec18,
    PVec19, PVec20, PVec21, PVec22, PVec23, PVec24,
    PVec25, PVec26, PVec27, PVec28, PVec29, PVec30,
    PVec31, PVec32)
PSmallVec{T}(n::Integer, u::AbstractArray{S,1}) where {T,S<:T} = let nu = length(u)
    if     n <= 0  return PVec0{T}()
    elseif n > 32  error("PSmallVec given vector of size > 32")
    elseif n > nu  error("PSmallVec given incorrect size")
    elseif n == nu return _psmallvec_bysize[n]{T}(u...)
    else           return _psmallvec_bysize[n]{T}(u[1:n]...)
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
push(u::PBigVec{T}, v::S) where {T, S<:T} = begin
    let n0 = length(u._elements), (k1,k2) = _pbigvec_iisplit(n+1)
        return PBigVec{T}(u._n+1,
                          k1 > n0
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
    elseif n <= 32 return _psmallvec_bysize[n]{T}(u[1:n]...)
    else
        let parts = [PSmallVec{T}(u[ii:min(ii+32-1, n)]) for ii in 1:32:n]
            return PBigVec{T}(n, PVec{PSmallVec{T}}(parts))
        end
    end
end
PVec(u::AbstractArray{T, 1}) where {T} = PVec{T}(u)