################################################################################
# arrayops.jl
#
# Array operations that answer a persistent array with a persistent array.
#
# These are the operations the `AbstractArray` interface provides that have no
# reason to leave the collection's own hierarchy: copying, mapping, selecting,
# reordering, and indexing with a vector of positions. Left to the generic
# fallbacks they allocate a mutable `Array` through `similar`, so an operation
# over a persistent collection quietly stops being persistent — `copy(u)`, for
# one, returned a `Vector`.
#
# Each of them preserves the operand's *default* as well as its entries, which is
# what keeps a sparse array sparse: selecting or reordering a `PArray` moves its
# stored entries and carries its default along, exactly as a broadcast does. An
# element equal to the default is not stored in the result, so `filter`ing a
# sparse array does not materialise the positions it discards.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

# A tree rebuilt from an operand's, with each stored entry moved from one linear
# position to another. `f` maps an old position to a new one, or to `nothing` to
# drop it.
function _air_retree(::Type{T}, u, f) where {T}
    tree = PTree{T}()
    tr = _air_tree(u)
    (tr === nothing) && return tree
    i0 = _air_i0(u)
    for (ii, v) in tr
        d = ii - i0
        j = f(Int(d) + 1)
        (j === nothing) && continue
        tree = setindex(tree, v, HASH_T(j - 1))
    end
    return tree
end

# ==============================================================================
# Copying

# A `PArray` is immutable, so a copy may be the original — as in `StaticArrays`,
# where `copy(::SArray)` returns the same array. What matters is that the result
# is the same *kind*.
Base.copy(u::PArray{T,N}) where {T,N} = u

# A transient, by contrast, is mutable: a copy has to be an independent one. That
# is `copy(t::TArray)`, in `transient.jl`.

# ==============================================================================
# Mapping

# `map` over arrays is elementwise, so it is the broadcast — which already
# answers a `PArray` with a `PArray` and a transient with a transient, and
# already keeps a sparse result sparse.
Base.map(f::F, u::Union{PArray,TArray}, rest...) where {F} = broadcast(f, u, rest...)

# ==============================================================================
# Selecting and reordering

"""
    filter(f, u::PVector)

Yields the `PVector` of the elements of `u` for which `f` is true, in order. The
result has the same default as `u`, so an element equal to the default is not
stored in it — filtering a sparse vector does not materialise the positions it
keeps.

This is the persistent counterpart of `Base.filter`, which would otherwise
return a mutable `Array`.
"""
function Base.filter(f, u::PVector{T}) where {T}
    tree = PTree{T}()
    dflt = getfield(u, :_default)
    j = 0
    for k in 1:length(u)
        x = _air_value(u, k)
        f(x) || continue
        j += 1
        _eqdefault(dflt, x) && continue
        tree = setindex(tree, x, HASH_T(j - 1))
    end
    return PVector{T}(HASH_T(0x0), _lindex(j), tree, dflt)
end

"""
    reverse(u::PVector)

Yields the `PVector` of the elements of `u` in the opposite order, with the same
default. Only the stored entries move, so the cost is the number of entries
rather than the length.
"""
function Base.reverse(u::PVector{T}) where {T}
    n = length(u)
    tree = _air_retree(T, u, k -> n + 1 - k)
    return PVector{T}(HASH_T(0x0), _lindex(n), tree, getfield(u, :_default))
end

# ==============================================================================
# Indexing with a vector of positions

# `u[[1, 3]]` and `u[1:2]`: the result is a `PVector` with `u`'s default, so a
# selected element equal to the default is not stored. The generic fallback
# builds a `Vector` through `similar`.
function Base.getindex(u::PVector{T}, I::AbstractVector{<:Integer}) where {T}
    tree = PTree{T}()
    dflt = getfield(u, :_default)
    for (j, k) in enumerate(I)
        x = u[k]                                    # bounds-checks `k`
        _eqdefault(dflt, x) && continue
        tree = setindex(tree, x, HASH_T(j - 1))
    end
    return PVector{T}(HASH_T(0x0), _lindex(length(I)), tree, dflt)
end

# The same for any number of dimensions, where a scalar index drops its dimension
# exactly as it does for an `Array`: `m[1:2, [1, 3]]` is a 2×2 `PArray` while
# `m[1:3, 3]` is a `PVector`. As in the vector case above, the result carries the
# operand's default, so a selected element equal to the default is not stored —
# selecting a sparsely-set matrix yields a sparsely-set result.
#
# A call whose indices are all integers is handled by the `Vararg{Int,N}` method
# above, which is the narrower of the two, so this one sees at least one vector
# index and the result is always an array rather than a scalar.
function Base.getindex(
    u::PArray{T,N}, I::Vararg{Union{Integer,AbstractVector{<:Integer}},N}
) where {T,N}
    dims = Tuple(length(i) for i in I if !(i isa Integer))
    tree = PTree{T}()
    dflt = getfield(u, :_default)
    li = LinearIndices(dims)
    for ci in CartesianIndices(dims)
        x = u[_air_select(I, ci)...]                # bounds-checks each index
        _eqdefault(dflt, x) && continue
        tree = setindex(tree, x, HASH_T(li[ci] - 1))
    end
    return PArray{T,length(dims)}(HASH_T(0x0), _lindex(dims...), tree, dflt)
end

# The index tuple of `u` that a result position selects: a vector index takes its
# position from `ci`, and a scalar index is used as it stands.
function _air_select(I::Tuple, ci::CartesianIndex)
    sel = Vector{Int}(undef, length(I))
    d = 0
    for (n, i) in enumerate(I)
        if i isa Integer
            sel[n] = i
        else
            d += 1
            sel[n] = i[ci[d]]
        end
    end
    return sel
end

# ==============================================================================
# Concatenation

# An operand's default, converted to the result's element type. The first
# operand's default becomes the result's, so that a concatenation of sparse
# arrays stays sparse.
function _air_default_as(u, ::Type{T}) where {T}
    d = _air_default(u)
    return d === nothing ? nothing : Tuple{T}((T(_defaultvalue(d)),))
end

# A transient operand makes a transient result, as it does for a broadcast.
_air_transient(args) = Val(any(a -> a isa TArray, args))

# The concatenations differ only in how an operand's linear positions map into
# the result, so they share this. An operand whose default is not the result's —
# or which has no default at all — is materialised position by position, because
# its default value is an explicit value of the result. An operand that shares
# the result's default only has to move the entries it stores.
function _air_concat(
    ::Type{T}, ::Val{N}, sz, dflt, args, posmap, ::Val{Transient}
) where {T,N,Transient}
    tree = PTree{T}()
    for (m, a) in enumerate(args)
        ad = _air_default(a)
        materialise =
            (dflt === nothing) || (ad === nothing) ||
            (_defaultvalue(ad) != _defaultvalue(dflt))
        if materialise
            for k in 1:length(a)
                v = T(_air_value(a, k))
                _eqdefault(dflt, v) && continue
                tree = _air_set(tree, v, HASH_T(posmap(m, k) - 1), Val(Transient))
            end
        else
            tr = _air_tree(a)
            (tr === nothing) && continue
            i0 = _air_i0(a)
            for (ii, v) in tr
                k = Int(ii - i0) + 1
                tree = _air_set(
                    tree, T(v), HASH_T(posmap(m, k) - 1), Val(Transient)
                )
            end
        end
    end
    return _air_pack(T, Val(N), tree, sz, dflt, Val(Transient))
end

"""
    vcat(u::PVector, rest::AbstractVector...)

Yields the `PVector` of the given vectors, one after another. The result's default
is the first operand's, so concatenating sparse vectors stays sparse; an operand
whose default differs from it is stored explicitly, as its default value is an
explicit value of the result.

This replaces the `SparseArrays` method that would otherwise be reached and
return a mutable `SparseVector`, which is why the element types are restricted to
`Number`: `SparseArrays` defines its own `vcat` over `AbstractVecOrMat{<:Number}`,
and a method over the same arguments but without that restriction is not more
specific than it, so the call would be *ambiguous* rather than intercepted. With
it, this method is strictly narrower in both arguments. A `PVector` of a
non-`Number` element type is left to `Base`, and its `vcat` is an `Array`.
"""
function Base.vcat(
    u::Union{PVector{T},TVector{T}}, us::AbstractVector{S}...
) where {T<:Number,S<:Number}
    args = (u, us...)
    Tt = promote_type(T, map(eltype, us)...)
    n = sum(length, args)
    offs = Vector{Int}(undef, length(args))
    acc = 0
    for (m, a) in enumerate(args)
        offs[m] = acc
        acc += length(a)
    end
    return _air_concat(
        Tt, Val(1), (n,), _air_default_as(u, Tt), args,
        (m, k) -> offs[m] + k, _air_transient(args),
    )
end

"""
    hcat(u::PMatrix, rest::AbstractMatrix...)

Yields the `PArray{Float64,2}`-style matrix of the given matrices, side by side.
The operands must have the same number of rows. As with `vcat`, the result's
default is the first operand's.

This replaces the `SparseArrays` method that would otherwise be reached and
return a mutable `SparseMatrixCSC`. As with `vcat`, the element types are
restricted to `Number` so that this method is strictly narrower than
`SparseArrays`' rather than ambiguous with it.
"""
function Base.hcat(
    u::Union{PMatrix{T},TArray{T,2}}, us::AbstractMatrix{S}...
) where {T<:Number,S<:Number}
    args = (u, us...)
    Tt = promote_type(T, map(eltype, us)...)
    rows = size(u, 1)
    all(a -> size(a, 1) == rows, us) || throw(
        DimensionMismatch("hcat: the matrices must have the same number of rows")
    )
    cols = sum(a -> size(a, 2), args)
    # An operand occupies a block of columns, so its linear positions shift by a
    # constant: the number of columns before it, times the number of rows.
    offs = Vector{Int}(undef, length(args))
    acc = 0
    for (m, a) in enumerate(args)
        offs[m] = acc * rows
        acc += size(a, 2)
    end
    return _air_concat(
        Tt, Val(2), (rows, cols), _air_default_as(u, Tt), args,
        (m, k) -> offs[m] + k, _air_transient(args),
    )
end

"""
    vcat(u::PMatrix, rest::AbstractMatrix...)

Yields the matrix of the given matrices, one above another. The operands must
have the same number of columns. As with the vectors, the result's default is the
first operand's, and the element types are restricted to `Number` for the same
reason — to be strictly narrower than `SparseArrays`' method rather than
ambiguous with it.
"""
function Base.vcat(
    u::Union{PMatrix{T},TArray{T,2}}, us::AbstractMatrix{S}...
) where {T<:Number,S<:Number}
    args = (u, us...)
    Tt = promote_type(T, map(eltype, us)...)
    cols = size(u, 2)
    all(a -> size(a, 2) == cols, us) || throw(
        DimensionMismatch("vcat: the matrices must have the same number of columns")
    )
    rows = sum(a -> size(a, 1), args)
    r = [size(a, 1) for a in args]
    off = Vector{Int}(undef, length(args))
    acc = 0
    for m in eachindex(args)
        off[m] = acc
        acc += r[m]
    end
    # Stacking rows is not a constant shift: column `j` of an operand lands in
    # the wider column of the result, so the mapping has to go through the
    # operand's own row count.
    posmap = function (m, k)
        i = mod1(k, r[m])
        j = div(k - 1, r[m]) + 1
        return (off[m] + i - 1) + (j - 1) * rows + 1
    end
    return _air_concat(
        Tt, Val(2), (rows, cols), _air_default_as(u, Tt), args, posmap,
        _air_transient(args),
    )
end

# ==============================================================================
# The elementwise arithmetic operators

# Each of these is the broadcast, which already answers a `PArray` with a
# `PArray` (a `TArray` for a transient operand), already maps the default as well
# as the stored entries, already drops an entry that lands on the new default, and
# already handles a scalar or a plain-array operand. Defining them is what stops
# `u + v` from falling back to `Base`, which allocates a mutable `Array`.
#
# The set of signatures is deliberate. `Base` defines `+` for `(AbstractArray,
# AbstractArray)`, `(AbstractArray, Number)` and `(Number, AbstractArray)`; a
# method whose argument is a `Union` of those alternatives is *not* narrower than
# them, so the call would be ambiguous rather than intercepted. Hence one method
# per combination of Air type and Base kind, plus the Air-Air case, which is the
# narrowest of them all and resolves the overlaps between the others.
#
# Matrix multiplication is deliberately not here. `A * B` is not elementwise, and
# `Base`'s method — which answers with a mutable `Matrix` — is the right answer
# for it. Only multiplication and division *by a number* are included.
const _AirArray = Union{PArray,TArray}

Base.:+(u::_AirArray, v::_AirArray) = broadcast(+, u, v)
Base.:+(u::_AirArray, v::Number) = broadcast(+, u, v)
Base.:+(u::_AirArray, v::AbstractArray) = broadcast(+, u, v)
Base.:+(u::Number, v::_AirArray) = broadcast(+, u, v)
Base.:+(u::AbstractArray, v::_AirArray) = broadcast(+, u, v)

Base.:-(u::_AirArray, v::_AirArray) = broadcast(-, u, v)
Base.:-(u::_AirArray, v::Number) = broadcast(-, u, v)
Base.:-(u::_AirArray, v::AbstractArray) = broadcast(-, u, v)
Base.:-(u::Number, v::_AirArray) = broadcast(-, u, v)
Base.:-(u::AbstractArray, v::_AirArray) = broadcast(-, u, v)
Base.:-(u::_AirArray) = broadcast(-, u)

# Multiplication and division by a number are elementwise, so they are included;
# array-by-array multiplication is `Base`'s (see the note above).
Base.:*(a::Number, u::_AirArray) = broadcast(*, a, u)
Base.:*(u::_AirArray, a::Number) = broadcast(*, u, a)
Base.:/(u::_AirArray, a::Number) = broadcast(/, u, a)

# ==============================================================================
# Reshape

# `reshape` preserves the linear order, and a `PArray`'s tree is keyed by that
# order, so the tree does not change at all: only the shape does. That makes this
# O(1) rather than a copy, and it keeps a sparse array sparse for free.
#
# A `TArray` is left to `Base`, which answers with a `ReshapedArray` — a mutable
# view over the transient, which is the right answer for a mutable array, just as
# a new `PArray` is for an immutable one.

# The dimensions a `:` stands for, following `Base`: at most one colon, and it
# takes up whatever the others leave. A zero among the others pins it to zero,
# and only if there is nothing to hold.
function _air_reshape_dims(n::Int, dims::Tuple)
    ncolon = count(d -> d === Colon(), dims)
    (ncolon > 1) &&
        throw(DimensionMismatch("reshape: at most one dimension may be `:`"))
    (ncolon == 0) && return dims
    known = prod(d -> d === Colon() ? 1 : Int(d), dims; init=1)
    if known == 0
        (n == 0) ||
            throw(DimensionMismatch("reshape: the dimensions do not multiply to $n"))
        return map(d -> d === Colon() ? 0 : Int(d), dims)
    end
    (n % known == 0) ||
        throw(DimensionMismatch("reshape: $n is not divisible by $known"))
    return map(d -> d === Colon() ? n ÷ known : Int(d), dims)
end

"""
    reshape(u::PArray, dims)

Yields the `PArray` of the same elements in the given shape. Since a `PArray`'s
tree is keyed by linear position and reshaping preserves the linear order, this
changes only the shape: it is O(1), and a sparse array stays sparse.

A `:` may stand for one of the dimensions, as it does for an `Array`, and takes
up whatever the other dimensions leave.

See also: `Base.reshape`, [`PArray`](@ref).
"""
function _air_reshape(u::PArray{T}, dims::Tuple) where {T}
    d = _air_reshape_dims(length(u), dims)
    (prod(d) == length(u)) ||
        throw(DimensionMismatch("reshape: the dimensions do not multiply to $(length(u))"))
    return PArray{T,length(d)}(
        getfield(u, :_i0), _lindex(d...), getfield(u, :_tree), getfield(u, :_default)
    )
end
Base.reshape(u::PArray{T}, dims::Tuple{Vararg{Union{Int,Colon}}}) where {T} =
    _air_reshape(u, dims)
Base.reshape(u::PArray, dims::Vararg{Union{Int,Colon}}) = _air_reshape(u, dims)
# The all-integer case, and the two `Colon` forms that Base defines for a vector,
# each have to be named explicitly: a signature over `Union{Int,Colon}` is not
# narrower than one over all `Int` (or over a bare `Colon`), so those calls would
# be ambiguous rather than reaching the methods above.
Base.reshape(u::PArray{T}, dims::Dims{M}) where {T,M} = _air_reshape(u, dims)
Base.reshape(u::PVector{T}, dims::Tuple{Colon}) where {T} = _air_reshape(u, dims)
Base.reshape(u::PVector, ::Colon) = _air_reshape(u, (Colon(),))

