################################################################################
# parray.jl
# The Persistent Array type, composed using the PTree type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2019-2026 Noah C. Benson

# ==============================================================================
# #TODO List

# This list is a snapshot; several of the entries it used to carry are done.
# Done: `pzeros`/`pones`/`pfill`; `broadcast` (see broadcast.jl); `reshape` and
# the elementwise arithmetic operators (see arrayops.jl); and every array
# operation that used to fall back to a mutable `Array` — `copy`, `map`,
# `filter`, `reverse`, indexed selection and `vcat`/`hcat` (also arrayops.jl).
# Done: `psparse`, `prand`, `prandn`, `pdiagm`, `blockdiag` and `permutedims` —
# see the generators below. (A `permute` method is not needed: `permutedims` is
# the name Julia uses for the non-mutating operation and `permute!` for the
# mutating one, which an immutable array cannot offer.)
# Matrix multiplication and the rest of linear algebra are deliberately left to
# Base, which answers with a mutable `Array` — the right answer for them.

# ==============================================================================
# PArray definition.

# We expand on some sparse array methods with PArrays (which are implemented as
# efficitn sparse arrays anyway.
using SparseArrays: SparseArrays
# `Random` is brought into scope by `pheap.jl`, which is included later, and the
# generators below name `Random.AbstractRNG` in a signature — evaluated when this
# file is defined, not when it is called.
using Random: Random

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

```jldoctest; filter=r"0-element (PArray{Any, ?1}|PVector{Any})|Any\\[\\]"
julia> PArray()
0-element PArray{Any,1}
```

```jldoctest; filter=r"[23]-element (PArray{Int64, ?1}|PVector{Int64}):"
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

```jldoctest; filter=r"2×3 (PArray{Symbol, ?2}|PMatrix{Symbol}):"
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
    _default::Union{Nothing,Tuple{T}}
end

# ==============================================================================
# PArray Constructors

function PArray{T,N}(default, size::NTuple{N,<:Integer}) where {T,N}
    default = Tuple{T}((default,))
    return PArray{T,N}(0x0, LinearIndices(size), PTree{T}(), default)
end
function PArray{T,N}(::UndefInitializer, size::NTuple{N,<:Integer}) where {T,N}
    return PArray{T,N}(0x0, LinearIndices(size), PTree{T}(), nothing)
end
function PArray{T,N}(default::S, size::Vararg{<:Integer,N}) where {T,N,S}
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
function PArray{T,N}(default::S, size::Vector{<:Integer}) where {T,N,S}
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
PArray(default::T, size::NTuple{N,<:Integer}) where {T,N} = PArray{T,N}(default, size)
function PArray(default::T, size::Vararg{<:Integer,N}) where {T,N}
    return PArray{T,N}(default, NTuple{N,Int}(size))
end
function PArray(default::T, size::Vector{<:Integer}) where {T}
    return PArray(default, NTuple{length(size),Int}(size))
end
function PArray{T,N}(a::AbstractArray{S,N}) where {T,N,S}
    tree = PTree{T}(a)
    return PArray{T,N}(HASH_T(0x0), _lindex(size(a)...), tree, nothing)
end
PArray{T,N}(p::PArray{T,N}) where {T,N} = p
PArray{T,N}() where {T,N} = PArray{T,N}(undef, (0, [1 for _ in 2:N]...))
PArray(a::AbstractArray{T,N}) where {T,N} = PArray{T,N}(a)
PArray(p::PArray{T,N}) where {T,N} = p
PArray() = PArray{Any,1}()
# Convert function also.
# The argument is `x::AbstractArray` rather than `x`: with `Any` this overlaps
# `LinearAlgebra`'s `convert(::Type{T<:AbstractArray}, ::AbstractQ)` and its
# `Factorization` variant, and neither `AbstractQ` nor `Factorization` is a
# subtype of `AbstractArray`, so declaring the array restriction makes both
# intersections empty. No constructor accepted a non-array single argument
# anyway, so nothing is lost by saying so.
Base.convert(::Type{PArray{T,N}}, x::AbstractArray) where {T,N} = PArray{T,N}(x)
Base.convert(::Type{PArray{T,N}}, x::PArray{T,N}) where {T,N} = x

# ==============================================================================
# PArray aliases

"""
    PVector{T}

An alias for `PArray{T,1}`, representing a persistent vector.
"""
const PVector{T} = PArray{T,1} where {T}
PVector(default::T, len::II) where {T,II<:Integer} = PArray{T,1}(default, (len,))
PVector(default::T, len::Tuple{II}) where {T,II<:Integer} = PArray{T,1}(default, len)
PVector(a::AbstractArray{T,1}) where {T} = PArray{T,1}(a)
PVector(p::PArray{T,1}) where {T} = p
PVector() = PArray{Any,1}()
export PVector

"""
    PMatrix{T}

An alias for `PArray{T,2}`, representing a persistent matrix.
"""
const PMatrix{T} = PArray{T,2} where {T}
PMatrix(args...) = PArray{Any,2}(args...)
function PMatrix(default::T, rs::II, cs::JJ) where {T,II<:Integer,JJ<:Integer}
    return PArray{T,2}(default, (rs, cs))
end
PMatrix(default::T, sz::Tuple{<:Integer,<:Integer}) where {T} = PArray{T,2}(default, sz)
PMatrix(a::AbstractArray{T,2}) where {T} = PArray{T,2}(a)
PMatrix(p::PArray{T,2}) where {T} = p
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

```jldoctest; filter=r"4-element (PArray{Int64, ?1}|PVector{Int64}):"
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
    for (k, v) in p
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
function SparseArrays.dropzeros(p::PArray{T,N}) where {T,N}
    t = _dropzeros(p._tree, p._default)
    (t === p._tree) && return p
    return PArray{T,N}(p._i0, p._index, t, p._default)
end
function SparseArrays.dropzeros!(p::PArray{T,N}) where {T,N}
    return error("dropzeros!: object of type $(typeof(p)) is immutable")
end
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
function SparseArrays.findnz(p::PArray{T,N}) where {T,N}
    sz = size(p._index)
    ndims = length(sz)
    n = SparseArrays.nnz(p)
    iilists = [Vector{Int}(undef, n) for _ in 1:N]
    vals = Vector{T}(undef, n)
    cis = CartesianIndices(size(p))
    for (elno, (k, v)) in enumerate(p._tree)
        ci = cis[Int(k - p._i0) + 1]
        for (iilist, oo) in zip(iilists, Tuple(ci))
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
SparseArrays.nonzeros(p::PArray{T,N}) where {T,N} = PVector{T}(T[v for (k, v) in p._tree])
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

function Base.iterate(u::PArray{T,N}, k::Int) where {T,N}
    return k > length(u) ? nothing : (u[k], k+1)
end
Base.iterate(u::PArray{T,N}) where {T,N} = iterate(u, 1)

# ==============================================================================
# Indexing methods

# A `PArray` only defines N-dimensional `getindex`, so it is a Cartesian array.
# Claiming `IndexLinear` here meant that linear indexing was left to Base's
# fallback; since Julia 1.13 that fallback raises a `CanonicalIndexError`
# instead, so `p[1]` (and anything built on it, such as array equality) failed.
Base.IndexStyle(::Type{PArray{T,N}}) where {T,N} = IndexCartesian()
function _parray_get(u::PTree{T}, ii::HASH_T, ::Nothing) where {T}
    x = get(u, ii, nothing)
    (x === nothing) || return x
    isa(nothing, T) && in(ii => nothing, u) && return x
    # An array with no default should hold a value at every position, so reaching
    # here means one was never set. This used to `return Array{T}(undef, 1)[1]`,
    # which yields *fresh* uninitialized memory on every read: nondeterministic
    # (two reads of the same position disagree), one allocation per read, and for
    # a `T` holding references a pointer that must never be dereferenced. Raising
    # is what the author's commented-out line intended, and it matches what Julia
    # does for a reference-holding `Array` built with `undef`.
    error(
        "PArray has unset values and no default: an array built with `undef` " *
        "must have every position set before it is read",
    )
end
_parray_get(u::PTree{T}, ii::HASH_T, d::Tuple{T}) where {T} = begin
    return get(u, ii, d[1])
end
function Base.getindex(u::PArray{T,N}, k::Vararg{Int,N}) where {T,N}
    if N > 1
        k = u._index[k...]
    else
        k = k[1]
    end
    (k < 1) && throw(BoundsError(u, k))
    (k > length(u)) && throw(BoundsError(u, k))
    return _parray_get(u._tree, HASH_T(k - 1) + u._i0, u._default)
end
function Base.setindex!(u::PArray{T,N}, v, k::Int) where {T,N}
    return error("setindex!: object of type $(typeof(u)) is immutable; see setindex()")
end
Base.firstindex(u::PArray{T,N}) where {T,N} = 1
Base.lastindex(u::PArray{T,N}) where {T,N} = length(u)

# ==============================================================================
# AbstractArray methods

function Base.push!(u::PArray{T,N}, v) where {T,N}
    return error("push!: object of type $(typeof(u)) is immutable; see push()")
end
function Base.pop!(u::PArray{T,N}, v) where {T,N}
    return error("pop!: object of type $(typeof(u)) is immutable; see last() and pop()")
end
#Base.pushfist!(u::PArray{T,N}, v) where {T,N} = error(
#    "pushfirst!: object of type $(typeof(u)) is immutable; see first() pushfirst()")
#Base.popfirst!(u::PArray{T,N}, v) where {T,N} = error(
#    "popfist!: object of type $(typeof(u)) is immutable; see first() and popfist()")
function _lindex(u::Vararg{Int,N}) where {N}
    return LinearIndices{N,NTuple{N,Base.OneTo{Int}}}(((Base.OneTo{Int}.(u))...,))
end
# There used to be a `permutedims(u::PArray, dims::NTuple{N,Int})` here. It could
# not have worked: it built an empty `PVector` with a `LinearIndices{0}` for an
# index, which its own constructor rejects — and being more specific than the
# general method below, it shadowed it. `permutedims` for a `PArray` is defined
# further down, where the generators are.
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
function setindex(u::PArray{T,N}, v::S, ci::CartesianIndex{N}) where {T,N,S}
    return setindex(u, v, u._index[ci])
end
function setindex(u::PArray{T,N}, v::S, k::Vararg{Idx,N}) where {T,N,S,Idx<:Integer}
    if N == 1
        k = k[1]
    else
        k = u._index[k...]
    end
    (k < 1) && throw(BoundsError(u, k))
    n = length(u)
    (k > n + 1) && throw(BoundsError(u, k))
    (N == 1) && (k > n) && return push(u, v)
    ii = u._i0 + HASH_T(k - 1)
    t = (_eqdefault(u._default, v) ? delete(u._tree, ii) : setindex(u._tree, v, ii))
    return t === u._tree ? u : PArray{T,N}(u._i0, u._index, t, u._default)
end
function setindex(u::PArray{T,N}, v::S, ii...) where {T,N,S}
    idcs = getindex(u._index, ii...)
    pp = broadcast(Pair{Int,T}, idcs, v)
    if isa(pp, Pair)
        return setindex(u, pp[2], pp[1])
    else
        for (k, v) in pp
            u = setindex(u, v, k)
        end
        return u
    end
end
# Push and pop methods are only defined for vectors
function push(u::PVector{T}, x::S) where {T,S}
    n = length(u)
    if _eqdefault(u._default, x)
        tree = u._tree
    else
        ii = u._i0 + HASH_T(n + 1 - 1)
        tree = setindex(u._tree, x, ii)
    end
    index = _lindex(n+1)
    return PVector{T}(u._i0, index, tree, u._default)
end
function pushfirst(u::PVector{T}, x::S) where {T,S}
    n = length(u)
    if _eqdefault(u._default, x)
        tree = u._tree
    else
        tree = setindex(u._tree, x, u._i0 - 0x1)
    end
    index = _lindex(n+1)
    return PVector{T}(u._i0 - 0x1, index, tree, u._default)
end
function pop(u::PVector{T}) where {T}
    n = length(u)
    (n == 0) && throw(ArgumentError("PArray must be non-empty"))
    ii = u._i0 + HASH_T(n - 1)
    tree = delete(u._tree, ii)
    return PVector{T}(u._i0, _lindex(n-1), tree, u._default)
end
function popfirst(u::PVector{T}) where {T}
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

```jldoctest; filter=r"3-element (PArray{Integer, ?1}|PVector{Integer}):"
julia> pzeros(Integer, 3)
3-element PArray{Integer,1}:
 0
 0
 0
```

```jldoctest; filter=r"1×2 (PArray{Bool, ?2}|PMatrix{Bool}):"
julia> pzeros(Bool, (1,2))
1×2 PArray{Bool,2}:
 0  0
```

```jldoctest; filter=r"2×3 (PArray{Float64, ?2}|PMatrix{Float64}):"
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

```jldoctest; filter=r"2×3 (PArray{Symbol, ?2}|PMatrix{Symbol}):"
julia> pfill(:abc, 2, 3)
2×3 PArray{Symbol,2}:
 :abc  :abc  :abc
 :abc  :abc  :abc
```
"""
pfill(val::T, dims::Vararg{Integer}) where {T} = PArray{T,length(dims)}(val, dims...)
pfill(val::T, dims::Tuple) where {T} = PArray{T,length(dims)}(val, dims)

# ==============================================================================
# Generators

# The raw default of an operand, or `nothing` — as `_parray_from` wants it, rather
# than the `Tuple{T}` the representation uses.
_air_rawdefault(u) = begin
    d = _air_default(u)
    return d === nothing ? nothing : _defaultvalue(d)
end

# Everything here that builds an array from a list of positions goes through this:
# it drops entries equal to the default, so the result is as sparse as its inputs
# were and never stores a value it would only have to read back as the default.
function _parray_from(
    ::Type{T}, ::Val{N}, dims::NTuple{N,Int}, default, entries
) where {T,N}
    tree = PTree{T}()
    dflt = default === nothing ? nothing : Tuple{T}((T(default),))
    li = LinearIndices(dims)
    for (idx, val) in entries
        v = T(val)
        _eqdefault(dflt, v) && continue
        tree = setindex(tree, v, HASH_T(li[idx...] - 1))
    end
    return PArray{T,N}(HASH_T(0x0), _lindex(dims...), tree, dflt)
end

"""
    psparse(I, J, V, m, n)
    psparse(A; default=zero(eltype(A)))

Yields the `PArray` holding `V[k]` at position `(I[k], J[k])`, of size `m` by `n`,
as `SparseArrays.sparse` does — or, from an array, one holding the entries of `A`
that are not zero, again as `sparse(A)` does.

Unlike a `SparseArray`, the array it yields may have a default value of its own;
`default` names it, and an entry equal to it is not stored at all. So
`psparse(A; default=NaN)` yields an array that reads as `NaN` wherever `A` is
`NaN`, rather than storing those positions — which is the only way to express
what a `SparseArray` cannot.

Repeated positions are an error rather than being summed, since a persistent
array cannot distinguish "written twice" from "written once".

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"2×3 (PArray{Float64, ?2}|PMatrix{Float64}):"
julia> psparse([1, 2], [2, 3], [10.0, 20.0], 2, 3)
2×3 PMatrix{Float64}:
 0.0  10.0   0.0
 0.0   0.0  20.0
```
"""
function psparse(
    I::AbstractVector{<:Integer}, J::AbstractVector{<:Integer},
    V::AbstractVector, m::Integer, n::Integer
)
    length(I) == length(J) == length(V) ||
        throw(DimensionMismatch("psparse: I, J and V must be the same length"))
    T = eltype(V)
    return _parray_from(
        T, Val(2), (Int(m), Int(n)), zero(T),
        (((Int(I[k]), Int(J[k])), V[k]) for k in eachindex(V)),
    )
end
psparse(I::AbstractVector, J::AbstractVector, V::AbstractVector) =
    psparse(I, J, V, maximum(I; init = 0), maximum(J; init = 0))
# The default is zero, so that this does what `sparse(A)` does and drops the
# zeros; `default` overrides it for the arrays a `SparseArray` cannot express.
psparse(A::AbstractArray{T,N}; default = zero(T)) where {T,N} =
    _parray_from(T, Val(N), size(A), default, _pnz(A))
# The entries of an array that are not its default; from a `PArray` that is
# exactly its tree, and from anything else every position.
function _pnz(p::PArray{T,N}) where {T,N}
    i0 = getfield(p, :_i0)
    ci = CartesianIndices(size(p))
    return (
        (Tuple(ci[Int(ii - i0) + 1]), v) for (ii, v) in getfield(p, :_tree)
    )
end
_pnz(a::AbstractArray) = ((Tuple(I), a[I]) for I in CartesianIndices(a))

"""
    prand(dims...)
    prand(m, n, p)

Yields a `PArray` of random values: of the given dimensions, as `rand` would, or
of size `m` by `n` with each position set with probability `p`, as `sprand`
would. The sparse form has a default of zero, so only the positions it fills are
stored.
"""
prand(dims::Vararg{Integer}) = PArray(rand(dims...))
prand(dims::Tuple) = PArray(rand(dims))
prand(rng::Random.AbstractRNG, dims::Vararg{Integer}) = PArray(rand(rng, dims...))
function prand(rng::Random.AbstractRNG, m::Integer, n::Integer, p::Real)
    return _parray_from(
        Float64, Val(2), (Int(m), Int(n)), 0.0,
        (((i, j), rand(rng)) for (i, j) in _prand_positions(rng, m, n, p))
    )
end
prand(m::Integer, n::Integer, p::Real) = prand(Random.default_rng(), m, n, p)

"""
    prandn(dims...)
    prandn(m, n, p)

As [`prand`](@ref), with normally distributed values.
"""
prandn(dims::Vararg{Integer}) = PArray(randn(dims...))
prandn(dims::Tuple) = PArray(randn(dims))
prandn(rng::Random.AbstractRNG, dims::Vararg{Integer}) = PArray(randn(rng, dims...))
function prandn(rng::Random.AbstractRNG, m::Integer, n::Integer, p::Real)
    return _parray_from(
        Float64, Val(2), (Int(m), Int(n)), 0.0,
        (((i, j), randn(rng)) for (i, j) in _prand_positions(rng, m, n, p))
    )
end
prandn(m::Integer, n::Integer, p::Real) = prandn(Random.default_rng(), m, n, p)

# The positions a sparse random array fills: each position with probability `p`,
# which is what `sprand` means and, incidentally, cannot produce a duplicate.
# Drawing a count from a binomial would be faster for a dense-in-the-limit `p`,
# but the draw has to be the definition rather than merely fast.
function _prand_positions(rng::Random.AbstractRNG, m::Integer, n::Integer, p::Real)
    ci = CartesianIndices((Int(m), Int(n)))
    return (Tuple(ci[ix]) for ix in 1:(Int(m) * Int(n)) if rand(rng) < p)
end

"""
    pdiagm(v)
    pdiagm(k => v, ...)

Yields the `PArray` with the vector `v` on its diagonal — or, for the second form,
on the `k`-th diagonal, as `SparseArrays.spdiagm` does. Entries equal to zero are
not stored, so the result is diagonal in the sparse sense as well as the
structural one.
"""
pdiagm(v::AbstractVector) = pdiagm(0 => v)
function pdiagm(offsets::Pair{<:Integer,<:AbstractVector}...)
    n = maximum(length(v) + abs(Int(k)) for (k, v) in offsets; init = 0)
    return _parray_from(
        promote_type(map(v -> eltype(last(v)), offsets)...), Val(2), (n, n), 0,
        _pdiagm_entries(offsets),
    )
end
function _pdiagm_entries(offsets)
    return (
        ((i, j), diag[i]) for (k, diag) in offsets for (i, j) in
        ((j + max(-Int(k), 0), j + max(Int(k), 0)) for j in 1:length(diag))
    )
end

"""
    permutedims(u::PArray, perm)

Yields the `PArray` whose dimensions are those of `u` rearranged by `perm`. As
with `reshape` this only moves the stored entries — a permutation of the axes
changes which linear position each one has, so the tree is rebuilt, but the
default is carried along and a sparse array stays sparse.

This is the persistent counterpart of `Base.permutedims`.
"""
function Base.permutedims(u::PArray{T,N}, perm) where {T,N}
    p = _perm_check(perm, Val(N))
    dims = ntuple(d -> size(u, p[d]), N)
    li = LinearIndices(dims)
    return _parray_from(
        T, Val(N), dims, _air_rawdefault(u),
        ((ntuple(d -> ci[p[d]], N), u[Tuple(ci)...]) for ci in CartesianIndices(size(u))),
    )
end
Base.permutedims(u::PArray{T,N}) where {T,N} =
    permutedims(u, ntuple(d -> N - d + 1, N))

function _perm_check(perm, ::Val{N}) where {N}
    p = ntuple(d -> Int(perm[d]), N)
    (sort(collect(p)) == collect(1:N)) ||
        throw(ArgumentError("permutedims: $perm is not a permutation of 1:$N"))
    return p
end

"""
    blockdiag(a::PArray, others...)

Yields the `PArray` holding the given arrays on its diagonal and its default
everywhere else. This is the persistent counterpart of
`SparseArrays.blockdiag`.
"""
# (`PArray` only: `TArray` is defined in `transient.jl`, which is included after
# this file, so it cannot appear in a signature here.)
function SparseArrays.blockdiag(a::PArray, others::PArray...)
    args = (a, others...)
    all(u -> ndims(u) == 2, args) ||
        throw(DimensionMismatch("blockdiag: every argument must be a matrix"))
    T = promote_type(map(eltype, args)...)
    rs = cumsum([size(u, 1) for u in args])
    cs = cumsum([size(u, 2) for u in args])
    dims = (rs[end], cs[end])
    return _parray_from(
        T, Val(2), dims, _air_rawdefault(a),
        _blockdiag_entries(args, rs, cs, T),
    )
end
function _blockdiag_entries(args, rs, cs, ::Type{T}) where {T}
    return (
        (((ro + ci[1], co + ci[2])), u[ci]) for
        (u, ro, co) in zip(args, ([0; rs[1:(end - 1)]]), ([0; cs[1:(end - 1)]])) for
        ci in CartesianIndices(size(u))
    )
end
export psparse, prand, prandn, pdiagm

export PArray, PVector, pzeros, pones, pfill
