################################################################################
# broadcast.jl
#
# Broadcasting for the persistent arrays and their transients.
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
# persistent collection is that operations defined on it stay persistent. A
# broadcast involving a `TArray` is a transient broadcast and yields a `TArray`,
# whose tree the result owns. `similar` is deliberately left alone for a
# `PArray`, and so still returns an `Array` — it is documented to return scratch
# space for the caller to write into, which an immutable array cannot be. (The
# same division holds for `StaticArrays`, where `SArray .+ 1` is an `SArray`
# while `similar(::SArray)` is mutable.) `similar` of a `TArray` does return a
# `TArray`, since a transient is mutable and so has no such problem.
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

The broadcast style of the persistent arrays and their transients. Declaring one
is what keeps a broadcast over a `PArray` a `PArray`, and a broadcast over a
`TArray` a `TArray`: without it the operation falls back to `DefaultArrayStyle`
and produces a mutable `Array`.

This type is part of `Air`'s internal/private implementation details.
"""
struct AirArrayStyle{N} <: Base.Broadcast.AbstractArrayStyle{N} end
Base.BroadcastStyle(::Type{<:PArray{T,N}}) where {T,N} = AirArrayStyle{N}()
Base.BroadcastStyle(::Type{<:TArray{T,N}}) where {T,N} = AirArrayStyle{N}()
Base.BroadcastStyle(::AirArrayStyle{M}, ::AirArrayStyle{N}) where {M,N} =
    AirArrayStyle{max(M, N)}()
# A broadcast involving either kind is one of ours, so this style wins over the
# default one in either order. A scalar contributes `DefaultArrayStyle{0}`, for
# which `max` gives `M`, so scalars are covered here too.
Base.BroadcastStyle(
    ::AirArrayStyle{M}, ::Base.Broadcast.DefaultArrayStyle{N}
) where {M,N} = AirArrayStyle{max(M, N)}()
Base.BroadcastStyle(
    ::Base.Broadcast.DefaultArrayStyle{N}, ::AirArrayStyle{M}
) where {M,N} = AirArrayStyle{max(M, N)}()
# A `PArray` and a `TArray` broadcast together to a `TArray`, so that the result
# is the mutable kind whenever either operand was; this is decided by the
# operands, in `_air_broadcast`.

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
_air_default(u::TArray) = getfield(u, :default)
_air_default(::AbstractArray) = nothing
_air_hasdefault(a) = _air_ispositional(a) && (_air_default(a) !== nothing)
# The tree of an operand that stores entries, or `nothing` if it stores none.
# The persistent and transient types name their fields differently.
_air_tree(u::PArray) = getfield(u, :_tree)
_air_tree(u::TArray) = getfield(u, :tree)
# Anything else — a plain array, or a scalar — stores no entries. This has to be
# `Any` rather than `AbstractArray`: a scalar operand is not an array.
_air_tree(::Any) = nothing
_air_i0(u::PArray) = getfield(u, :_i0)
_air_i0(u::TArray) = getfield(u, :i0)
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
function _air_value(u::TArray, k::Int)
    return _parray_get(
        getfield(u, :tree), getfield(u, :i0) + HASH_T(k - 1), getfield(u, :default)
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
        tree = _air_tree(a)
        (tree === nothing) && continue
        i0 = _air_i0(a)
        for (ii, _) in tree
            d = ii - i0
            (d < HASH_T(n)) || continue
            push!(ks, Int(d) + 1)
        end
    end
    return ks
end

# ==============================================================================
# Building the result

# `Val`s rather than a run-time branch: which of these applies is fixed by the
# operands' types, so the loop below compiles to one of the two.
_air_set(tree::PTree, v, ii::HASH_T, ::Val{false}) = setindex(tree, v, ii)
# The transient form claims the nodes it builds, which is what makes the result a
# transient that owns its own tree rather than one sharing another's.
_air_set(tree::PTree, v, ii::HASH_T, ::Val{true}) = _ptree_tsetindex(tree, v, ii)
function _air_pack(
    ::Type{Tr}, ::Val{N}, tree::PTree{Tr}, sz, dflt, ::Val{false}
) where {Tr,N}
    return PArray{Tr,N}(HASH_T(0x0), _lindex(sz...), tree, dflt)
end
function _air_pack(
    ::Type{Tr}, ::Val{N}, tree::PTree{Tr}, sz, dflt, ::Val{true}
) where {Tr,N}
    return TArray{Tr,N}(HASH_T(0x0), _lindex(sz...), tree, dflt, prod(sz))
end

function _air_broadcast(
    bc::Base.Broadcast.Broadcasted{AirArrayStyle{N}}, ::Type{Tr}
) where {N,Tr}
    f = getfield(bc, :f)
    args = getfield(bc, :args)
    ax = axes(bc)
    sz = ntuple(d -> length(ax[d]), N)
    n = prod(sz)
    dflt = _air_result_default(f, args, Tr)
    # A transient operand makes a transient result, so that the result of a batch
    # update is still a batch update.
    transient = Val(any(a -> a isa TArray, args))
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
            tree = _air_set(tree, v, HASH_T(k - 1), transient)
        end
    else
        # A dense result: compute and store every position. Evaluating through
        # the broadcasted object, rather than argument by argument, is what makes
        # this correct when the operands' shapes differ.
        for I in CartesianIndices(ax)
            v = Tr(Base.Broadcast._broadcast_getindex(bc, I))
            k = LinearIndices(sz)[I]
            tree = _air_set(tree, v, HASH_T(k - 1), transient)
        end
    end
    return _air_pack(Tr, Val(N), tree, sz, dflt, transient)
end

function Base.copy(bc::Base.Broadcast.Broadcasted{AirArrayStyle{N}}) where {N}
    bc = Base.Broadcast.flatten(bc)
    Tr = Base.Broadcast.combine_eltypes(getfield(bc, :f), getfield(bc, :args))
    return _air_broadcast(bc, Tr)
end
