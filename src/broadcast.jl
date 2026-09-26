################################################################################
# broadcast.jl
#
# Broadcasting for the persistent arrays.
#
# A `PArray` is a default value plus a set of explicit entries, so broadcasting
# maps both halves of it: the default becomes the function applied to the
# operands' defaults, and an explicit entry becomes the function applied to the
# operands' entries at the same position. An entry that lands on the new default
# is dropped, which preserves the representation's invariant that nothing stored
# is equal to the default. So `pvec .+ 1` increments the sparse value as well as
# the explicit ones, and a batch of explicit values that happen to equal the
# sparse value simply stores nothing.
#
# The result is a `PArray` rather than a mutable `Array`: the point of a
# persistent collection is that operations defined on it stay persistent.
# `similar` is deliberately left alone and so still returns an `Array` — it is
# documented to return scratch space for the caller to write into, which an
# immutable array cannot be. (The same division holds for `StaticArrays`, where
# `SArray .+ 1` is an `SArray` while `similar(::SArray)` is mutable.)
#
# Two situations force a dense result, in which every position is computed and
# stored. A positional operand with no default of its own — a plain `Array`, or a
# `PArray` built with `undef` — leaves no default to map. And operands whose
# shapes differ can place several result positions at one stored entry, so the
# sparse shortcut is not available.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

# ==============================================================================
# Broadcast style

"""
    AirArrayStyle{N}

The broadcast style of the persistent arrays. Declaring one is what keeps a
broadcast over a `PArray` a `PArray`: without it the operation falls back to
`DefaultArrayStyle` and produces a mutable `Array`.

This type is part of `Air`'s internal/private implementation details.
"""
struct AirArrayStyle{N} <: Base.Broadcast.AbstractArrayStyle{N} end
Base.BroadcastStyle(::Type{<:PArray{T,N}}) where {T,N} = AirArrayStyle{N}()
Base.BroadcastStyle(::AirArrayStyle{M}, ::AirArrayStyle{N}) where {M,N} =
    AirArrayStyle{max(M, N)}()
# A broadcast involving a `PArray` is a `PArray` broadcast, so this style wins
# over the default one in either order. A scalar contributes
# `DefaultArrayStyle{0}`, for which `max` gives `M`, so it is covered here too.
Base.BroadcastStyle(
    ::AirArrayStyle{M}, ::Base.Broadcast.DefaultArrayStyle{N}
) where {M,N} = AirArrayStyle{max(M, N)}()
Base.BroadcastStyle(
    ::Base.Broadcast.DefaultArrayStyle{N}, ::AirArrayStyle{M}
) where {M,N} = AirArrayStyle{max(M, N)}()

# ==============================================================================
# Reading operands at a position and at the default

# Whether an operand occupies positions in the result at all. A scalar is
# broadcast to every position, so it affects the mapped default but is never
# stored.
_air_ispositional(::AbstractArray) = true
_air_ispositional(::Any) = false
# An operand's own default, in the representation's `Tuple{T}` form. `nothing`
# means the operand is dense — every position is explicit.
_air_default(u::PArray) = getfield(u, :_default)
_air_default(::AbstractArray) = nothing
_air_hasdefault(a) = _air_ispositional(a) && (_air_default(a) !== nothing)
# An operand's value at linear position `k`. The sparse path only runs when every
# positional operand has the result's shape, so linear indexing agrees with the
# broadcast's own index mapping.
function _air_value(u::PArray, k::Int)
    return _parray_get(
        getfield(u, :_tree),
        getfield(u, :_i0) + HASH_T(k - 1),
        getfield(u, :_default),
    )
end
_air_value(u::AbstractArray, k::Int) = u[k]
_air_value(x, ::Int) = x

# ==============================================================================
# The result's default and its explicit entries

# The result's default, or `nothing` if some positional operand is dense. The
# function is applied to the operands' defaults because the default is the value
# the array has at every position it does not store — so it has to be mapped just
# as an explicit entry is.
function _air_result_default(f, args, ::Type{Tr}) where {Tr}
    vals = Vector{Any}(undef, length(args))
    for (i, a) in enumerate(args)
        if _air_ispositional(a)
            d = _air_default(a)
            (d === nothing) && return nothing
            vals[i] = _defaultvalue(d)
        else
            vals[i] = a
        end
    end
    return Tuple{Tr}((Tr(f(vals...)),))
end

# The linear positions where at least one operand stores an entry explicitly.
# Only those can differ from the result's default, so only those need visiting:
# this is what keeps a broadcast over a sparse array proportional to the number
# of entries rather than to its length.
function _air_positions(args, n::Int)
    ks = Set{Int}()
    for a in args
        (a isa PArray) || continue
        i0 = getfield(a, :_i0)
        for (ii, _) in getfield(a, :_tree)
            d = ii - i0
            (d < HASH_T(n)) || continue
            push!(ks, Int(d) + 1)
        end
    end
    return ks
end

# ==============================================================================
# The broadcast itself

function _air_broadcast(
    bc::Base.Broadcast.Broadcasted{AirArrayStyle{N}}, ::Type{Tr}
) where {N,Tr}
    f = getfield(bc, :f)
    args = getfield(bc, :args)
    ax = axes(bc)
    sz = ntuple(d -> length(ax[d]), N)
    n = prod(sz)
    dflt = _air_result_default(f, args, Tr)
    # Every position can be derived from the stored entries only when each
    # positional operand has the result's shape and a default of its own.
    sparse = all(
        a -> !_air_ispositional(a) || (size(a) == sz && _air_hasdefault(a)), args
    )
    tree = PTree{Tr}()
    if sparse
        for k in _air_positions(args, n)
            v = Tr(f((_air_value(a, k) for a in args)...))
            (dflt !== nothing) && (v == _defaultvalue(dflt)) && continue
            tree = setindex(tree, v, HASH_T(k - 1))
        end
    else
        # A dense result: compute and store every position. Evaluating through
        # the broadcasted object, rather than argument by argument, is what makes
        # this correct when the operands' shapes differ.
        for I in CartesianIndices(ax)
            v = Tr(Base.Broadcast._broadcast_getindex(bc, I))
            k = LinearIndices(sz)[I]
            tree = setindex(tree, v, HASH_T(k - 1))
        end
    end
    return PArray{Tr,N}(HASH_T(0x0), _lindex(sz...), tree, dflt)
end

function Base.copy(bc::Base.Broadcast.Broadcasted{AirArrayStyle{N}}) where {N}
    bc = Base.Broadcast.flatten(bc)
    Tr = Base.Broadcast.combine_eltypes(getfield(bc, :f), getfield(bc, :args))
    return _air_broadcast(bc, Tr)
end
