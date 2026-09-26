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
