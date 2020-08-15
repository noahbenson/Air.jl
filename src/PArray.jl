################################################################################
# PArrray.jl
# The Persistent Array type and related types such as TArray (the transient
# array type).
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2019 Noah C. Benson


# ==============================================================================
# #TODO List

# psparse() method for making PArrays similar to SparseArrays.sparse().
# pzeros() similar to zeros() and spzeros()
# prand() similar to rand() and sprand()
# prandn() similar to randn() and sprandn()
# pdiagm() similar to spdiagm()
# blockdiag() method instance.
# permute() method
# broadcast() method
# Many others: overload arithmetic operators? reshape, others.


# ==============================================================================
# PArray definition.

# We expand on some sparse array methods with PArrays (which are implemented as
# efficitn sparse arrays anyway.
import SparseArrays

"""
    AbstractPArray{T,N}

The Abstract persistent array type represents any Array-like type that is
persistent. Reified types include PArray{T,N} and LazyArray{T,N}.
"""
abstract type AbstractPArray{T,N} <: AbstractArray{T,N} end

"""
    PArray{T,N}

The PArray type is a persistent/immutable corrolary to the Array{T,N} type. Like
Array, PArray can store n-dimensional non-ragged arrays. However, unlike Arrays,
PArrays can create duplicates of themselves with finite edits in log-time.

PArrays have a similar interface as Arrays, but instead of the functions
`push!`, `pop!`, and `setindex!`, PArrays use `push`, `pop`, and `setindex`.
PArrays also have efficient implementations of `pushfirst` and `popfirst`.
"""
struct PArray{T,N} <: AbstractPArray{T,N}
    # The initial element of the list; the tree is basically an enormous
    # circular buffer, so it's okay for this to roll around at UInt max.
    _i0::_PTREE_KEY_T
    # The dimensions of the array and the indices into the tree.
    _index::LinearIndices{N,NTuple{N,Base.OneTo{Int}}}
    # The data in a long array.
    _tree::PTree{T}
    # The default value, if any (for sparse arrays).
    _default::Union{Nothing, Tuple{T}}
end
mutability(::Type{PArray}) = Immutable()
mutability(::Type{PArray{T,N}}) where {T,N} = Immutable()


# ==============================================================================
# PArray Constructors

PArray{T,N}(default::S, size::NTuple{N,<:Integer}) where {T,N,S} = begin
    default = isa(default, UndefInitializer) ? nothing : Tuple{T}(default)
    return PArray{T,N}(0x0, LinearIndices(size), PTree{T}(), default)
end
PArray{T,N}(default::S, size::Vararg{<:Integer,N}) where {T,N,S} = begin
    return PArray{T,N}(dflt, NTuple{N,Int}(size))
end
PArray{T,N}(dflt::S, size::Vector{<:Integer}) where {T,N,S} = begin
    return PArray{T,N}(dflt, NTuple{N,Int}(size))
end
PArray(default::T, size::NTuple{N,<:Integer}) where {T,N} = PArray{T,N}(default, size)
PArray(default::T, size::Vararg{<:Integer,N}) where {T,N} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray(default::T, size::Vector{<:Integer}) where {T,N} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray{T,N}(a::AbstractArray{S,N}) where {T,N,S} = begin
    # #TODO: rewrite using TArray or TTree
    tree = PTree{T}()
    k = _PTREE_KEY_T(0x0)
    for x in a
        k += 0x1
        tree = setindex(tree, x, k)
    end
    return PArray{T,N}(0x0, _lindex(size(a)...), tree, nothing)
end
PArray{T,N}(p::PArray{T,N}) where {T,N} = p
PArray{T,N}() where {T,N} = PArray{T,N}(undef, (0, [1 for _ in 2:N]...))
PArray(a::AbstractArray{T,N}) where {T,N} = PArray{T,N}(a)
PArray(p::PArray{T,N}) where {T,N} = p
PArray() = PArray{Any,1}()
# #TODO: psparse() method for making PArrays similar to SparseArrays.sparse().


# ==============================================================================
# PArray aliases

# Declare some aliases of PArray; this is most efficient using a macro to
# generate the basic code block.
macro _declare_parray_alias(name::Symbol, d::Int)
    sname = String(name)
    docstr = ["    $(sname){T}",
              "",
              "$(sname){T} is an alias for PArray{T,$(d)}"]
    docstr = join(docstr, "\n")
    docstr = join(["\n", docstr, "\n"])
    return quote
        $docstr
        const $(name){T} = PArray{T,$(d)} where {T}
        function $(name)(a::AbstractArray{T,$(d)}) where {T}
            return PArray{T,$(d)}(a)
        end
        function $(name)(a::PArray{T,$(d)}) where {T}
            return a
        end
        function $(name)()
            return PArray{Any,$(d)}()
        end
    end |> esc
end

@_declare_parray_alias P1Tensor 1
@_declare_parray_alias P2Tensor 2
@_declare_parray_alias P3Tensor 3
@_declare_parray_alias P4Tensor 4
@_declare_parray_alias P5Tensor 5
@_declare_parray_alias P6Tensor 6
#@_declare_parray_alias PVector 1
#@_declare_parray_alias PMatrix 2
const PVector{T} = PArray{T,1} where {T}
const PMatrix{T} = PArray{T,2} where {T}


# ==============================================================================
# SparseArrays methods.

"""
    nnz(p::PArray)

Yields the number of explicitly set values in the persistent array p, regardless
of the number that are zero. This is different from the sparse-array library
only in that persistent arrays support arbitrary default values instead of
supporting only the value zero. Thus this counts explicit values instead of
non-zero values.
"""
SparseArrays.nnz(u::PArray{T,N}) where {T,N} = length(u._tree)
_dropzeros(p::PTree{T}, ::Nothing) where {T} = p
_dropzeros(p::PTree{T}, df::Tuple{T}) where {T} = begin
    df = df[1]
    for (k,v) in p
        if v == df
            p = delete(p, k)
        end
    end
    return p
end
"""
    dropzeros(p::PArray)

Drops explicit values of the given array p that are equal to the array's default
value. This differs from the SparseArrays implementation of dropzeros() only in
that PArrays allow arbitrary default values, while SparseArrays allows only the
default value of zero.
"""
SparseArrays.dropzeros(p::PArray{T,N}) where {T,N} = begin
    t = _dropzeros(p._tree, p._default)
    (t === p._tree) && return p
    return PArray{T,N}(p._i0, p._index, t, p._default)
end
SparseArrays.dropzeros!(p::PArray{T,N}) where {T,N} = error(
    "dropzeros!: object of type $(typeof(p)) is immutable")
"""
    findnz(p::PArray)

Yields the explicitly set elements of the given persistent array p. This method
is identical to the typical SparseArrays implementation of findnz() except that
it respects the arbitrary default-value that persistent arrays are allowed to
have rather than assuming that this value is a zero, as is done in the
SparseArrays library.
"""
SparseArrays.findnz(p::PArray{T,N}) where {T,N} = begin
    sz = size(p._index)
    ndims = length(sz)
    n = nnz(p)
    iilists = [Vector{Int}(undef, n) for _ in 1:N]
    vals = Vector{T}(undef, n)
    cis = CartesianIndices(Tuple([p._size[k] for k in p._dimorder]...))
    for (elno,(k,v)) in enumerate(p._tree)
        ci = cis[Int(k) + 1]               
        for (iilist,oo) in zip(iilists,ci)
            iilist[elno] = oo
        end
        vals[elno] = v
    end
    return Tuple(iilists) + (vals,)
end
"""
    nonzeros(p::PArray)

Yields the explicitly set values of the given persistent array p. This method
is identical to the typical SparseArrays implementation of nonzeros() except
that it respects the arbitrary default-value that persistent arrays are allowed
to have rather than assuming that this value is a zero, as is done in the
SparseArrays library.
"""
SparseArrays.nonzeros(p::PArray{T,N}) where {T,N} = T[v for (k,v) in p._tree]
"""
    issparse(p::PArray)

Yields true, as PArrays are implicitly sparse.
"""
SparseArrays.issparse(::PArray) = true


# ==============================================================================
# Iterator methods

Base.IteratorSize(::Type{PArray{T,N}}) where {T,N} = HasSize{N}()
Base.IteratorEltype(::Type{PArray{T,N}}) where {T,N} = HasEltype()

Base.length(u::PArray{T,N}) where {T,N} = length(u._index)
Base.size(u::PArray{T,N}) where {T,N} = size(u._index)
Base.eltype(u::PArray{T,N}) where {T,N} = T

Base.iterate(u::PArray{T,N}, k::Int) where {T,N} = begin
    return k > length(u) ? nothing : (u[k], k+1)
end
Base.iterate(u::PArray{T,N}) where {T,N} = iterate(u, 1)


# ==============================================================================
# Indexing methods

IndexStyle(::Type{PArray{T,N}}) where {T,N} = IndexLinear()
_parray_get(u::PTree{T}, ii::_PTREE_KEY_T, ::Nothing) where {T} = begin
    x = get(u, ii, nothing)
    (x === nothing) || return x
    isa(nothing, T) && in(ii => nothing, u) && return x
    #error("PArray has unset values and no default")
    return Array{T}(undef, 1)[1]
end
_parray_get(u::PTree{T}, ii::_PTREE_KEY_T, d::Tuple{T}) where {T} = begin
    return get(u, ii, d[1])
end
Base.getindex(u::PArray{T,N}, k::Vararg{Int,N}) where {T,N} = begin
    if N > 1
        k = u._index[k...]
    else
        k = k[1]
    end
    (k < 1) && throw(BoundsError(u, k))
    (k > length(u)) && throw(BoundsError(u, k))
    return _parray_get(u._tree, _PTREE_KEY_T(k) + u._i0, u._default)
end
Base.setindex!(u::PArray{T,N}, v, k::Int) where {T,N} = error(
    "setindex!: object of type $(typeof(u)) is immutable; see setindex()")
Base.firstindex(u::PArray{T,N}) where {T,N} = 1
Base.lastindex(u::PArray{T,N}) where {T,N} = length(u)


# ==============================================================================
# AbstractArray methods

Base.push!(u::PArray{T,N}, v) where {T,N} = error(
    "push!: object of type $(typeof(u)) is immutable; see push()")
Base.pop!(u::PArray{T,N}, v) where {T,N} = error(
    "pop!: object of type $(typeof(u)) is immutable; see last() and pop()")
#Base.pushfist!(u::PArray{T,N}, v) where {T,N} = error(
#    "pushfirst!: object of type $(typeof(u)) is immutable; see first() pushfirst()")
#Base.popfirst!(u::PArray{T,N}, v) where {T,N} = error(
#    "popfist!: object of type $(typeof(u)) is immutable; see first() and popfist()")
_lindex(u::Vararg{Int,N}) where {N} = LinearIndices{N,NTuple{N,Base.OneTo{Int}}}(
    ((Base.OneTo{Int}.(u))...,))
Base.permutedims(u::PArray{T,N}, dims::NTuple{N,Int}) where {T,N} = begin
    # This is actually pretty easy, code-wise:
    # #TODO: use TArray here
    v = PVector{T}(0x0, LinearIndices(()), PTree{T}(), u._default)
    for k in permutedims(u._index, dims)
        v = push(v, u[k])
    end
    sz = size(u)
    idx = _lindex([sz[k] for k in dims]...)
    return PArray{T,N}(v._i0, idx, v._tree, u._default)
end
#Base.broadcast(fn::F, u::PArray{T,N}, args...) where {F<:Function,T,N} = begin
#
#end

# ==============================================================================
# Persistent array methods

_defaultvalue(::Nothing) = undef
_defaultvalue(u::Tuple{T}) where {T} = u[1]
_eqdefault(::Nothing, x) = false
_eqdefault(dflt::Tuple{T}, x::S) where {T,S} = (dflt[1] == x)
defaultvalue(u::PArray{T,N}) where {T,N} = _defaultvalue(u._default)
setindex(u::PArray{T,N}, v::S, ci::CartesianIndex{N}) where {T,N,S} = begin
    return setindex(u, v, u._index[ci])
end
setindex(u::PArray{T,N}, v::S, k::Vararg{Idx,N}) where {T,N,S,Idx<:Integer} = begin
    if N == 1
        k = k[1]
    else
        k = u._index[k...]
    end
    (k < 1) && throw(BoundsError(u,k))
    n = length(u)
    (k > n + 1) && throw(BoundsError(u,k))
    (N == 1) && (k > n) && return push(u, v)
    ii = u._i0 + _PTREE_KEY_T(k)
    t = (_eqdefault(u._default, v) ? delete(u._tree, ii)
                                   : setindex(u._tree, v, ii))
    return t === u._tree ? u : PArray{T,N}(u._i0, u._index, t, u._default)
end
setindex(u::PArray{T,N}, v::S, ii...) where {T,N,S} = begin
    idcs = getindex(u._index, ii...)
    pp = broadcast(Pair{Int,T}, idcs, v)
    if isa(pp, Pair)
        return setindex(u, pp[2], pp[1])
    else
        # #TODO: Use TArray for this
        for (k,v) in pp
            u = setindex(u, v, k)
        end
        return u
    end
end
# Push and pop methods are only defined for vectors
push(u::PVector{T}, x::S) where {T,S} = begin
    n = length(u)
    if _eqdefault(u._default, x)
        tree = u._tree
    else
        ii   = u._i0 + _PTREE_KEY_T(n + 1)
        tree = setindex(u._tree, x, ii)
    end
    index = _lindex(n+1)
    return PVector{T}(u._i0, index, tree, u._default)
end
pushfirst(u::PVector{T}, x::S) where {T,S} = begin
    n = length(u)
    if _eqdefault(u._default, x)
        tree = u._tree
    else
        tree = setindex(u._tree, x, u._i0)
    end
    index = _lindex(n+1)
    return PVector{T}(u._i0 - 0x1, index, tree, u._default)
end
pop(u::PVector{T}) where {T} = begin
    n = length(u)
    (n == 0) && throw(ArgumentError("PArray must be non-empty"))
    ii = u._i0 + _PTREE_KEY_T(n)
    tree = delete(u._tree, ii)
    return PVector{T}(u._i0, _lindex(n-1), tree, u._default)
end
popfirst(u::PVector{T}) where {T} = begin
    n = length(u)
    (n == 0) && throw(ArgumentError("PArray must be non-empty"))
    tree = delete(u._tree, u._i0 + 0x1)
    return PVector{T}(u._i0 + 0x1, _lindex(n-1), tree, u._default)
end
isequiv(u::AS, v::AT) where {
    S,T,N,
    AS <: AbstractPArray{S,N},
    AT <: AbstractPArray{T,N}
} = begin
    (u === v) && return true
    (size(u) == size(v)) || return false
    for (a,b) in zip(u, v)
        isequiv(a, b) || return false
    end
    return true
end
equivhash(u::AA) where {T,N,AA<:AbstractPArray{T,N}} = begin
    h = equivhash(size(u))
    for x in u
        h *= 0x1f
        h += equivhash(x)
    end
    return h
end
