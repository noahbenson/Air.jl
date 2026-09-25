################################################################################
# transient.jl
#
# Transient (mutable) versions of the persistent collections, for batch updates.
#
# A transient is a way to make many updates and pay for the persistence once,
# rather than on every step. It is *not* a copy: `transient` is O(1) and shares
# every node of the collection it is made from, and `persistent!` is O(1) and
# hands that structure back. What changes is what an update costs. A persistent
# update cannot touch anything it was given, so it must copy every node on the
# path from the root to the leaf it changes. A transient first *claims* those
# nodes — replaces each with a copy that it alone owns (see PTREE_OWNED_FLAG in
# ptree.jl) — and thereafter changes the ones it owns in place, so that after the
# first update of a batch the ancestors stop being copied and their cell vectors
# stop being reallocated.
#
# The ownership flag is what makes that safe rather than merely fast: an owned
# node is reachable only from the transient that made it, so mutating it cannot
# be observed by any collection the caller still holds. A node with the flag
# clear is never mutated, only copied.
#
# Contract, as in Clojure: after `persistent!`, the transient must not be used
# again. The collection it hands back shares nodes with the transient, and the
# transient's remaining nodes have the flag set, so a later update through it
# would change a value the caller believes is immutable.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

# # Tree operations ============================================================

"""
    _ptree_claim(u)

Yields a copy of the tree node `u` that the calling transient owns, or `u`
itself if it is already owned. The copy's cells vector is fresh, so it may be
changed in place; `u` is left untouched.
"""
@inline function _ptree_claim(u::PTree{T}) where {T}
    # The canonical empty node has no cells to change; there is nothing to own.
    (getfield(u, :numel) == 0) && return u
    id = getfield(u, :id)
    ptree_owned(id) && return u
    cells = getfield(u, :cells)
    return PTree{T}(
        ptree_owned(id, true),
        getfield(u, :bits),
        getfield(u, :numel),
        cells === nothing ? nothing : copy(cells),
    )
end
"""
    _ptree_claimcells(u, cells)

Yields a cells vector that the calling transient may change in place: `cells`
itself if it already owns the node `u`, and a copy of it otherwise.
"""
@inline _ptree_claimcells(u::PTree{T}, cells::Vector{S}) where {T,S} =
    (ptree_owned(getfield(u, :id)) ? cells : copy(cells))
# An owned twig holding a single leaf, which is how a transient introduces new
# subtrees.
@inline function _ptree_ttwig(::Type{T}, k::HASH_T, v) where {T}
    id = ptree_owned(ptree_id(k & ~PTREE_DEPTH_MASK, PTREE_TWIG_DEPTH), true)
    bits = PTREE_BITS_T(0x1) << (k & PTREE_DEPTH_MASK)
    return PTree{T}(id, bits, 1, T[v])
end
"""
    _ptree_tsetindex(u, v, k)

`PTree.setindex`, except that the nodes it touches end up owned by the calling
transient and their cells vectors are reused rather than copied. Yields the new
tree, sharing with `u` everything the operation did not touch.
"""
function _ptree_tsetindex(u::PTree{T}, v::V, k::HASH_T) where {T,V}
    if getfield(u, :numel) == 0
        return _ptree_ttwig(T, k, v)
    end
    id = getfield(u, :id)
    if !ptree_isbeneath(id, k)
        # `k` is outside this node's range: build an owning ancestor holding both,
        # exactly as `setindex` does, but with the flag set.
        kv = _ptree_ttwig(T, k, v)
        topmiss = ptree_highbitdiff(id, k)
        if topmiss <= HASH_BITCOUNT - PTREE_ROOT_SHIFT
            nodelev = div(topmiss - PTREE_TWIG_SHIFT, PTREE_NODE_SHIFT)
            bit0 = PTREE_TWIG_SHIFT + nodelev * PTREE_NODE_SHIFT
            shift = PTREE_NODE_SHIFT
            newdepth = PTREE_LEVELS - 2 - nodelev
        else
            newdepth = 0
            bit0 = HASH_BITCOUNT - PTREE_ROOT_SHIFT
            shift = PTREE_ROOT_SHIFT
        end
        newminleaf = id & ~((HASH_ONE << bit0) - HASH_ONE)
        newid = ptree_owned(ptree_id(newminleaf, newdepth), true)
        newmask = (HASH_ONE << shift) - HASH_ONE
        ii_k = newmask & (id >> bit0)
        ii_u = newmask & (k >> bit0)
        newbits = (BITS_ONE << ii_k) | (BITS_ONE << ii_u)
        newchs = ii_u < ii_k ? PTree{T}[kv, u] : PTree{T}[u, kv]
        return PTree{T}(newid, newbits, 1 + getfield(u, :numel), newchs)
    end
    bits = getfield(u, :bits)
    numel = getfield(u, :numel)
    (bit0, shift) = ptree_bitshift(id)
    idx = (k >> bit0) & lowmask(shift)
    flag = BITS_ONE << idx
    idx = 1 + count_ones(bits & lowmask(idx))
    oid = ptree_owned(id, true)
    if ptree_depth(id) == PTREE_TWIG_DEPTH
        cells = getfield(u, :cells)::Vector{T}
        cells = _ptree_claimcells(u, cells)
        if bits & flag == flag
            @inbounds cells[idx] = v
            return PTree{T}(oid, bits, numel, cells)
        else
            insert!(cells, idx, v)
            return PTree{T}(oid, bits | flag, numel + 1, cells)
        end
    else
        cells = getfield(u, :cells)::Vector{PTree{T}}
        if bits & flag == flag
            oldc = @inbounds cells[idx]
            oldn = getfield(oldc, :numel)
            newc = _ptree_tsetindex(oldc, v, k)
            (newc === oldc) && return u
            newn = getfield(newc, :numel)
            cells = _ptree_claimcells(u, cells)
            @inbounds cells[idx] = newc
            return PTree{T}(oid, bits, numel - oldn + newn, cells)
        else
            cells = _ptree_claimcells(u, cells)
            insert!(cells, idx, _ptree_ttwig(T, k, v))
            return PTree{T}(oid, bits | flag, numel + 1, cells)
        end
    end
end

"""
    _ptree_tdelete(u, k)

`PTree.delete`, except that the nodes it touches end up owned by the calling
transient and their cells vectors are shortened in place rather than copied. The
minimal-tree collapse `delete` does — a branch left with one child becomes that
child — is preserved; the node handed back for a collapsed branch is a child,
which the parent then claims in the ordinary way.
"""
function _ptree_tdelete(u::PTree{T}, k::HASH_T) where {T}
    (getfield(u, :numel) == 0) && return u
    id = getfield(u, :id)
    bits = getfield(u, :bits)
    (inq, bitidx, idx) = ptree_cellindex(id, bits, k)
    inq || return u
    numel = getfield(u, :numel)
    oid = ptree_owned(id, true)
    if ptree_depth(id) == PTREE_TWIG_DEPTH
        # The last leaf leaves the canonical empty node, which holds no cells and
        # so has nothing for a transient to own.
        (numel == 1) && return PTree{T}()
        cells = getfield(u, :cells)::Vector{T}
        cells = _ptree_claimcells(u, cells)
        deleteat!(cells, idx)
        return PTree{T}(oid, bits & ~(BITS_ONE << bitidx), numel - 1, cells)
    else
        cells = getfield(u, :cells)::Vector{PTree{T}}
        oldc = @inbounds cells[idx]
        oldn = getfield(oldc, :numel)
        newc = _ptree_tdelete(oldc, k)
        (newc === oldc) && return u
        newn = getfield(newc, :numel)
        numel += newn - oldn
        (numel == 0) && return newc
        if newn == 0
            cells = _ptree_claimcells(u, cells)
            deleteat!(cells, idx)
            (length(cells) == 1) && return @inbounds cells[1]
            return PTree{T}(oid, bits & ~(BITS_ONE << bitidx), numel, cells)
        else
            (length(cells) == 1) && return newc
            cells = _ptree_claimcells(u, cells)
            @inbounds cells[idx] = newc
            return PTree{T}(oid, bits, numel, cells)
        end
    end
end

"""
    _ptree_clean(u)

Yields a tree equal to `u` with every node's ownership flag cleared, visiting
only the nodes a transient owns.

This is what makes `persistent!` sound. An owned node is one a transient may
change in place, so a tree handed to the caller must not contain any: otherwise
a later `transient` of that tree would see the flag, take the node for its own,
and write through it — changing a collection the caller believes is immutable.

The walk is cheap for the same reason the transient is. Every node on the path a
transient updates is rebuilt as owned, so the owned nodes form a connected spine
from the root down to each claimed leaf: an owned node's ancestors are all owned
too. So the walk descends into a child only when the child is owned, and stops
everywhere else — it visits exactly the nodes the batch created, no others, and
never reaches an untouched subtree. For an append-heavy batch that is a handful
of ancestors plus one twig per 64 appends, since consecutive appends claim the
same rightmost twig. It rebuilds a node struct per visited node (the id is
immutable) but copies no cells vector, since it owns the ones it visits.
"""
function _ptree_clean(u::PTree{T}) where {T}
    id = getfield(u, :id)
    ptree_owned(id) || return u
    clean = id & ~PTREE_OWNED_FLAG
    cells = getfield(u, :cells)
    if cells === nothing || ptree_depth(id) == PTREE_TWIG_DEPTH
        # A twig's cells are leaves, and an empty node has none: nothing below to
        # visit either way.
        return PTree{T}(clean, getfield(u, :bits), getfield(u, :numel), cells)
    end
    cs = cells::Vector{PTree{T}}
    changed = false
    for i in eachindex(cs)
        c = @inbounds cs[i]
        ptree_owned(getfield(c, :id)) || continue
        nc = _ptree_clean(c)
        # The visited vectors are all owned, so they can be changed in place.
        @inbounds cs[i] = nc
        changed = true
    end
    return PTree{T}(clean, getfield(u, :bits), getfield(u, :numel), cs)
end

# #TArray =====================================================================
"""
    TArray{T,N}

A mutable counterpart of `PArray{T,N}`, for batch updates. See the file comment
for the contract; use [`transient`](@ref) to make one and
[`persistent!`](@ref) to get a `PArray` back.

So far only the vector operations are provided — `push!` and `pop!`, both on
`TArray{T,1}` — so a `TArray` is made from and yields a `PVector`.
N-dimensional updates, which change entries without changing the shape, will
follow.
"""
mutable struct TArray{T,N}
    i0::HASH_T
    index::LinearIndices{N,NTuple{N,Base.OneTo{Int}}}
    tree::PTree{T}
    default::Union{Nothing,Tuple{T}}
    n::Int
end
"""
    TVector{T}

An alias for `TArray{T,1}`, the transient counterpart of a [`PVector`](@ref).
"""
const TVector{T} = TArray{T,1} where {T}
export TArray, TVector

"""
    transient(coll)

Yields a transient (mutable) counterpart of `coll`, which shares its structure
and is O(1). Make updates through the transient, then take the persistent result
with `persistent!`. `coll` is not modified, and must not be used afterwards if
you keep updating the transient.
"""
transient(u::PArray{T,N}) where {T,N} =
    TArray{T,N}(getfield(u, :_i0), getfield(u, :_index), getfield(u, :_tree),
                getfield(u, :_default), length(u))

"""
    persistent!(t)

Yields a persistent collection holding what the transient `t` currently holds,
in O(1). `t` must not be used afterwards.
"""
function persistent!(t::TArray{T,N}) where {T,N}
    # A vector's shape follows its length, which the transient has been tracking;
    # any other shape is unchanged, since nothing outside the vector operations
    # can change a transient's length.
    index = N == 1 ? _lindex(t.n) : t.index
    return PArray{T,N}(t.i0, index, _ptree_clean(t.tree), t.default)
end

Base.length(t::TArray) = t.n
function Base.push!(t::TArray{T,1}, x::S) where {T,S}
    # As in `push` for a PVector: a value equal to the array's default needs no
    # entry in the tree at all.
    if !_eqdefault(t.default, x)
        t.tree = _ptree_tsetindex(t.tree, x, t.i0 + HASH_T(t.n))
    end
    t.n += 1
    return t
end
function Base.pop!(t::TArray{T,1}) where {T}
    (t.n == 0) && throw(ArgumentError("PArray must be non-empty"))
    t.tree = _ptree_tdelete(t.tree, t.i0 + HASH_T(t.n - 1))
    t.n -= 1
    return t
end

# #TDict ======================================================================
# The persistent dictionaries are generated from a macro so that the hash and
# identity flavours stay in step; the transients mirror them the same way. The
# tree operations are the same `_ptree_tsetindex` the array transient uses — it
# is generic in the leaf type, which here is a collision bucket.
macro _tdict_code(tname::Symbol, pname::Symbol, hashfn::Symbol, dicttype::Symbol)
    return esc(
        quote
            mutable struct $tname{K,V} <: AbstractPDict{K,V}
                count::Int
                root::PTree{$dicttype{K,V}}
            end
            Base.length(t::$tname) = getfield(t, :count)
            transient(u::$pname{K,V}) where {K,V} =
                $tname{K,V}(getfield(u, :count), getfield(u, :root))
            persistent!(t::$tname{K,V}) where {K,V} =
                $pname{K,V}(getfield(t, :count), _ptree_clean(getfield(t, :root)))
            Base.get(t::$tname{K,V}, k, df) where {K,V} = begin
                ld = get(getfield(t, :root), $hashfn(k), nothing)
                return ld === nothing ? df : get(ld, k, df)
            end
            Base.getindex(t::$tname{K,V}, k) where {K,V} = begin
                ld = get(getfield(t, :root), $hashfn(k), nothing)
                (ld === nothing) && throw(KeyError(k))
                return ld[k]
            end
            Base.haskey(t::$tname{K,V}, k) where {K,V} = begin
                ld = get(getfield(t, :root), $hashfn(k), nothing)
                return ld === nothing ? false : haskey(ld, k)
            end
            Base.push!(t::$tname{K,V}, kv::Pair) where {K,V} = setindex!(t, kv.second, kv.first)
            function Base.setindex!(t::$tname{K,V}, v, k) where {K,V}
                hh = $hashfn(k)
                root = getfield(t, :root)
                ld0 = get(root, hh, nothing)
                if ld0 === nothing
                    t.count += 1
                    t.root = _ptree_tsetindex(root, $dicttype{K,V}(Pair{K,V}(k, v)), hh)
                else
                    ld1 = setindex(ld0, v, k)
                    (ld1 === ld0) && return t
                    t.count += length(ld1) - length(ld0)
                    t.root = _ptree_tsetindex(root, ld1, hh)
                end
                return t
            end
            function Base.delete!(t::$tname{K,V}, k) where {K,V}
                hh = $hashfn(k)
                root = getfield(t, :root)
                ld0 = get(root, hh, nothing)
                (ld0 === nothing) && return t
                ld1 = delete(ld0, k)
                (ld1 === ld0) && return t
                t.count -= 1
                # An emptied bucket goes away; otherwise the bucket replaces the
                # old one. (The bucket itself is copied rather than claimed: a
                # bucket has no ownership flag, and one is usually touched once
                # per batch anyway.)
                #
                # TODO (after the 1.0 release): a `TLinearDict`/`TLinearSet` —
                # a transient form of the collision bucket itself — would let a
                # bucket be claimed and changed in place, as the tree nodes are
                # here. Test whether the gain is worth it first: buckets are
                # small, and a batch usually touches a given bucket once, so the
                # ceiling on this is the per-bucket copy of two short vectors.
                t.root =
                    length(ld1) == 0 ? _ptree_tdelete(root, hh) :
                    _ptree_tsetindex(root, ld1, hh)
                return t
            end
        end,
    )
end
@_tdict_code TDict PDict hash PLinearDict
@_tdict_code TIdDict PIdDict objectid PIdLinearDict
export TDict, TIdDict
@doc """
    TDict{K,V}

A mutable counterpart of `PDict{K,V}`, for batch updates: make many updates
through it, then take the persistent result with `persistent!`. See
[`transient`](@ref) and the `transient.jl` file comment for the contract;
[`TIdDict`](@ref) is the identity-keyed counterpart of `PIdDict`.

# When a transient is worth it

A batch update costs `N + 1` allocations fewer through a `TDict` than through a
`PDict` for every node on the path from the root that the transient has already
claimed, but it also pays for making the transient and for the walk
`persistent!` does to clear the ownership flags. So the gain is proportional to
the *batch* and constant in the dictionary, and small batches lose:

| entries in the dictionary | entries updated | `PDict` | `TDict` | |
|---|---|---|---|---|
| 1000 | 1 | 19 allocs | 24 allocs | 0.79x |
| 1000 | 10 | 181 | 184 | 0.98x |
| 100 | 100 | 1612 | 1276 | 1.26x |
| 1000 | 1000 | 18451 | 13673 | 1.35x |
| 10000 | 1000 | 19948 | 15665 | 1.27x |

The crossover is around **ten updates**: below it a transient costs more than it
saves (about 20% more allocations for a single update), and above it the gain
grows slowly to roughly a third fewer allocations by a few hundred updates. The
time gain is smaller than the allocation gain — about 1.1x at the sizes above —
because the transient does a little more work per update.

Note that this is the *overwrite* case, which is what a batch update usually
means. Two things bound the gain: a collision bucket has no ownership flag, so
bucket updates are copied as before (a bucket is usually touched once per batch
anyway), and the ceiling is the tree path, which is most but not all of an
update's allocations. For appending to a `TArray` the gain is much larger
(2.9x fewer allocations), because consecutive appends claim the *same* twig.
""" TDict
@doc """
    TIdDict{K,V}

A mutable counterpart of `PIdDict{K,V}`, otherwise as [`TDict`](@ref). Keys are
compared by identity, as in `IdDict`.
""" TIdDict

export transient, persistent!
