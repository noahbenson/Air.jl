################################################################################
# parray.jl
# The Persistent Array type, composed using the PTree type.
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
    PArray{T,N}

The PArray type is a persistent/immutable corrolary to the Array{T,N} type. Like
Array, PArray can store n-dimensional non-ragged arrays. However, unlike Arrays,
PArrays can create duplicates of themselves with finite edits in log-time.

PArrays have a similar interface as Arrays, but instead of the functions
`push!`, `pop!`, and `setindex!`, PArrays use `push`, `pop`, and `setindex`.
PArrays also have efficient implementations of `pushfirst` and `popfirst`.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"0-element (PArray{Any, ?1}|PVector{Any})"
julia> PArray()
0-element PArray{Any,1}
```

```jldoctest; filter=r"[23]-element (PArray{Int64, ?1}|PVector{Int64})"
julia> u = PArray{Int,1}([1,2])
2-element PArray{Int64,1}:
 1
 2

julia> push(u, 3)
3-element PArray{Int64,1}:
 1
 2
 3
```

```jldoctest; filter=r"2×3 (PArray{Symbol, 2}|PMatrix{Symbol})"
julia> PArray{Symbol,2}(:abc, (2,3))
2×3 PArray{Symbol,2}:
 :abc  :abc  :abc
 :abc  :abc  :abc
```
"""
struct PArray{T,N} <: AbstractPArray{T,N}
    # The initial element of the list; the tree is basically an enormous
    # circular buffer, so it's okay for this to roll around at UInt max.
    _i0::HASH_T
    # The dimensions of the array and the indices into the tree.
    _index::LinearIndices{N,NTuple{N,Base.OneTo{Int}}}
    # The data in a long array.
    _tree::PTree{T}
    # The default value, if any (for sparse arrays).
    _default::Union{Nothing, Tuple{T}}
end


# ==============================================================================
# PArray Constructors

PArray{T,N}(default, size::NTuple{N,<:Integer}) where {T,N} = begin
    default = Tuple{T}((default,))
    return PArray{T,N}(0x0, LinearIndices(size), PTree{T}(), default)
end
PArray{T,N}(::UndefInitializer, size::NTuple{N,<:Integer}) where {T,N} =
    return PArray{T,N}(0x0, LinearIndices(size), PTree{T}(), nothing)
PArray{T,N}(default::S, size::Vararg{<:Integer,N}) where {T,N,S} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray{T,N}(default::S, size::Vector{<:Integer}) where {T,N,S} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray(default::T, size::NTuple{N,<:Integer}) where {T,N} = PArray{T,N}(default, size)
PArray(default::T, size::Vararg{<:Integer,N}) where {T,N} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray(default::T, size::Vector{<:Integer}) where {T,N} = begin
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray{T,N}(a::AbstractArray{S,N}) where {T,N,S} = begin
    tree = PTree{T}(a)
    return PArray{T,N}(HASH_T(0x0), _lindex(size(a)...), tree, nothing)
end
PArray{T,N}(p::PArray{T,N}) where {T,N} = p
PArray{T,N}() where {T,N} = PArray{T,N}(undef, (0, [1 for _ in 2:N]...))
PArray(a::AbstractArray{T,N}) where {T,N} = PArray{T,N}(a)
PArray(p::PArray{T,N}) where {T,N} = p
PArray() = PArray{Any,1}()
# Convert function also.
Base.convert(::Type{PArray{T,N}}, x) where {T,N} = PArray{T,N}(x)
Base.convert(::Type{PArray{T,N}}, x::PArray{T,N}) where {T,N} = x


# ==============================================================================
# PArray aliases

"""
    PVector{T}

An alias for `PArray{T,1}`, representing a persistent vector.
"""
const PVector{T} = PArray{T,1} where {T}
PVector(default::T, len::II) where {T,II<:Integer} = PArray{T,1}(default, (len,))
PVector(default::T, len::Tuple{II}) where {T,II<:Integer} =
    PArray{T,1}(default, len)
PVector(a::AbstractArray{T,1}) where {T} = PArray{T,1}(a)
PVector(p::PArray{T,1}) where {T,N} = p
PVector() = PArray{Any,1}()
export PVector

"""
    PMatrix{T}

An alias for `PArray{T,2}`, representing a persistent matrix.
"""
const PMatrix{T} = PArray{T,2} where {T}
PMatrix(args...) = PArray{Any,2}(args...)
PMatrix(default::T, rs::II, cs::JJ) where {T,II<:Integer,JJ<:Integer} =
    PArray{T,2}(default, (rs,cs))
PMatrix(default::T, sz::Tuple{<:Integer,<:Integer}) where {T} =
    PArray{T,2}(default, len)
PMatrix(a::AbstractArray{T,2}) where {T} = PArray{T,2}(a)
PMatrix(p::PArray{T,2}) where {T,N} = p
PMatrix() = PArray{Any,2}()
export PMatrix

# ==============================================================================
# SparseArrays methods.

"""
    nnz(p::PArray)

Yields the number of explicitly set values in the persistent array p, regardless
of the number that are zero. This is different from the sparse-array library
only in that persistent arrays support arbitrary default values instead of
supporting only the value zero. Thus this counts explicit values instead of
non-zero values.

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"4-element (PArray{Int64, 1}|PVector{Int64}):"
julia> u = PVector{Int}(0, (4,))
4-element PArray{Int64,1}:
 0
 0
 0
 0

julia> v = setindex(u, 2, 3)
4-element PArray{Int64,1}:
 0
 0
 2
 0

julia> nnz(v)
1

julia> nnz(u)
0
```
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

Drops explicit values of the given array `p` that are equal to the array's
default value. This differs from the SparseArrays implementation of dropzeros()
only in that `PArray`s allow arbitrary default values, while `SparseArray`s
allow only the default value of zero.

Note that under most circumstances, a `PArrray` will not encode explicit zeros,
so this function typically returns the object `p` untouched.

See also: `SparseArrays.nnz`, `SparseArrays.findnz`, [`PArray`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"4-element (PArray{Int64, ?1}|PVector{Int64}):"
julia> u = PVector{Int}([0,1,2,3])
4-element PArray{Int64,1}:
 0
 1
 2
 3

julia> dropzeros(u) === u
true

julia> v = setindex(u, 0, 3)
4-element PArray{Int64,1}:
 0
 1
 0
 3

julia> dropzeros(v) === v
true
```
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

Yields the explicitly set elements of the given persistent array `p`. This
method is identical to the typical SparseArrays implementation of findnz()
except that it respects the arbitrary default-value that persistent arrays are
allowed to have rather than assuming that this value is a zero, as is done in
the `SparseArray`s library.

Note that under most circumstances, a `PArray` will not encode explicit zeros,
so this function typically returns indices and values for all values that aren't
equal to the default value of the array `p` (which is zero by default).

See also: `SparseArrays.nnz`, [`PArray`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"4-element (PArray{Int64, ?1}|PVector{Int64}):"
julia> u = PVector([0,10,20,30])
4-element PArray{Int64,1}:
  0
 10
 20
 30

julia> findnz(u)
([1, 2, 3, 4], [0, 10, 20, 30])
```

```jldoctest; filter=r"4-element (PArray{Float64, ?1}|PVector{Float64}):"
julia> u = setindex(PVector(0.0, 4), 20, 2)
4-element PArray{Float64,1}:
  0.0
 20.0
  0.0
  0.0

julia> findnz(u)
([2], [20.0])
```
"""
SparseArrays.findnz(p::PArray{T,N}) where {T,N} = begin
    sz = size(p._index)
    ndims = length(sz)
    n = SparseArrays.nnz(p)
    iilists = [Vector{Int}(undef, n) for _ in 1:N]
    vals = Vector{T}(undef, n)
    cis = CartesianIndices(size(p))
    for (elno,(k,v)) in enumerate(p._tree)
        ci = cis[Int(k - p._i0) + 1]
        for (iilist,oo) in zip(iilists,Tuple(ci))
            iilist[elno] = oo
        end
        vals[elno] = v
    end
    return (iilists..., vals)
end
"""
    nonzeros(p::PArray)

Yields the explicitly set values of the given persistent array `p`. This method
is identical to the typical `SparseArrays` implementation of `nonzeros()` for 
its sparse array classes except that it returns a persistent array of values and
that it respects the arbitrary default-value that persistent arrays are allowed
to have rather than assuming that this value is a zero, as is done in the
`SparseArrays` library.

Note that because `PArray`s don't typically store values equal to their default
value explicitly, this will typically yield a vector of every non-default value
in the array.

See also: `SparseArrays.findnz`, `SparseArrays.nnz`, [`PArray`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"4-element (PArray{Int64, ?1}|PVector{Int64}):"
julia> u = PVector([0,10,20,30])
4-element PArray{Int64,1}:
  0
 10
 20
 30

julia> nonzeros(u)
4-element PArray{Int64,1}:
  0
 10
 20
 30
```

```jldoctest; filter=r"[14]-element (PArray{Float64, ?1}|PVector{Float64}):"
julia> u = setindex(PVector(0.0, 4), 20, 2)
4-element PArray{Float64,1}:
  0.0
 20.0
  0.0
  0.0

julia> nonzeros(u)
1-element PArray{Float64,1}:
 20.0
```
"""
SparseArrays.nonzeros(p::PArray{T,N}) where {T,N} =
    PVector{T}(T[v for (k,v) in p._tree])
# PArrays are always considered sparse.
SparseArrays.issparse(::PArray) = true


# ==============================================================================
# Iterator methods

Base.IteratorSize(::Type{PArray{T,N}}) where {T,N} = Base.HasShape{N}()
Base.IteratorEltype(::Type{PArray{T,N}}) where {T,N} = Base.HasEltype()

Base.length(u::PArray{T,N}) where {T,N} = length(u._index)
Base.size(u::PArray{T,N}) where {T,N} = size(u._index)
Base.eltype(::Type{PArray{T,N}}) where {T,N} = T
Base.eltype(u::PArray{T,N}) where {T,N} = T

Base.iterate(u::PArray{T,N}, k::Int) where {T,N} = begin
    return k > length(u) ? nothing : (u[k], k+1)
end
Base.iterate(u::PArray{T,N}) where {T,N} = iterate(u, 1)


# ==============================================================================
# Indexing methods

IndexStyle(::Type{PArray{T,N}}) where {T,N} = IndexLinear()
_parray_get(u::PTree{T}, ii::HASH_T, ::Nothing) where {T} = begin
    x = get(u, ii, nothing)
    (x === nothing) || return x
    isa(nothing, T) && in(ii => nothing, u) && return x
    #error("PArray has unset values and no default")
    return Array{T}(undef, 1)[1]
end
_parray_get(u::PTree{T}, ii::HASH_T, d::Tuple{T}) where {T} = begin
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
    return _parray_get(u._tree, HASH_T(k - 1) + u._i0, u._default)
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
    ii = u._i0 + HASH_T(k - 1)
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
        ii   = u._i0 + HASH_T(n + 1 - 1)
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
        tree = setindex(u._tree, x, u._i0 - 0x1)
    end
    index = _lindex(n+1)
    return PVector{T}(u._i0 - 0x1, index, tree, u._default)
end
pop(u::PVector{T}) where {T} = begin
    n = length(u)
    (n == 0) && throw(ArgumentError("PArray must be non-empty"))
    ii = u._i0 + HASH_T(n - 1)
    tree = delete(u._tree, ii)
    return PVector{T}(u._i0, _lindex(n-1), tree, u._default)
end
popfirst(u::PVector{T}) where {T} = begin
    n = length(u)
    (n == 0) && throw(ArgumentError("PArray must be non-empty"))
    tree = delete(u._tree, u._i0)
    return PVector{T}(u._i0 + 0x1, _lindex(n-1), tree, u._default)
end

# pzeros, pones, and other useful utility functions.
"""
    pzeros(dims...)
    pzeros(T, dims...)

Yields a persistent array (`PArray`) of zeros exactly as done by the `zeros()`
function.

See also [`pones`](@ref), [`pfill`](@ref), `zeros`

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"3-element (PArray{Float64, ?1}|PVector{Float64}):"
julia> pzeros(Integer, 3)
3-element PArray{Integer,1}:
 0
 0
 0
```

```jldoctest; filter=r"1×2 (PArray{Bool, ?2}:|PMatrix{Bool}):"
julia> pzeros(Bool, (1,2))
1×2 PArray{Bool,2}:
 0  0
```

```jldoctest; filter=r"2×3 (PArray{Float64, ?2}:|PMatrix{Float64}):"
julia> pzeros(2, 3)
2×3 PArray{Float64,2}:
 0.0  0.0  0.0
 0.0  0.0  0.0
```

```jldoctest; filter=r"1×1×1×1 PArray{Float64, ?4}:"
julia> pzeros((1, 1, 1, 1))
1×1×1×1 PArray{Float64,4}:
[:, :, 1, 1] =
 0.0
```
"""
pzeros(::Type{T}, dims::Vararg{Integer}) where {T} = PArray{T,length(dims)}(0, dims...)
pzeros(::Type{T}, dims::Tuple) where {T} = PArray{T,length(dims)}(0, dims...)
pzeros(dims::Tuple) = PArray{Float64,length(dims)}(0.0, dims)
pzeros(dims::Vararg{Integer}) = PArray{Float64,length(dims)}(0.0, dims...)

"""
    pones(dims...)
    pones(T, dims...)

Yields a persistent array (`PArray`) of ones exactly as done by the `ones()`
function.

See also [`pzeros`](@ref), [`pfill`](@ref), `ones`

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"3-element (PArray{Integer, ?1}|PVector{Integer}):"
julia> pones(Integer, 3)
3-element PArray{Integer,1}:
 1
 1
 1
```

```jldoctest; filter=r"1×2 (PArray{Bool, ?2}|PMatrix{Bool}):"
julia> pones(Bool, (1,2))
1×2 PArray{Bool,2}:
 1  1
```

```jldoctest; filter=r"2×3 (PArray{Float64, ?2}|PMatrix{Float64}):"
julia> pones(2, 3)
2×3 PArray{Float64,2}:
 1.0  1.0  1.0
 1.0  1.0  1.0
```

```jldoctest; filter=r"1×1×1×1 PArray{Float64, ?4}:"
julia> pones((1, 1, 1, 1))
1×1×1×1 PArray{Float64,4}:
[:, :, 1, 1] =
 1.0
```
"""
pones(::Type{T}, dims...) where {T} = PArray{T,length(dims)}(1, dims...)
pones(::Type{T}, dims::Tuple) where {T} = PArray{T,length(dims)}(1, dims...)
pones(dims::Tuple) = PArray{Float64,length(dims)}(1.0, dims)
pones(dims::Vararg{Integer}) = PArray{Float64,length(dims)}(1.0, dims...)

"""
    pfill(val, dims...)

Yields a persistent array (`PArray`) of values exactly as done by the `fill()`
function.

See also [`pones`](@ref), [`pzeros`](@ref), `fill`

# Examples

```@meta
DocTestSetup = quote
    using Air, SparseArrays
end
```

```jldoctest; filter=r"3-element (PArray{Float64, ?1}|PVector{Float64}):"
julia> pfill(NaN, 3)
3-element PArray{Float64,1}:
 NaN
 NaN
 NaN
```

```jldoctest; filter=r"2×3 (PArray{Symbol, ?2}|PMatrix{Float64}):"
julia> pfill(:abc, 2, 3)
2×3 PArray{Symbol,2}:
 :abc  :abc  :abc
 :abc  :abc  :abc
```
"""
pfill(val::T, dims::Vararg{Integer}) where {T} = PArray{T,length(dims)}(val, dims...)
pfill(val::T, dims::Tuple) where {T} = PArray{T,length(dims)}(val, dims)

export PArray, PVector, pzeros, pones, pfill
