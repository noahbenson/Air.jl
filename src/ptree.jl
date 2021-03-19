################################################################################
# ptree.jl
#
# Code implementing a persistent hash array mapped trie tyoe on which the
# persistent collection library of Air is built. This type is not intended for
# use outside of the library.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020 Noah C. Benson

# ==============================================================================
# Private
#
# This section contains private constants and function definitions intendef for
# use only in this library or occasionally within the library, but not outside
# of the current module.

# # Constants ==================================================================
#   This is the type of the unsigned integer used for the leaves.
const HASH_T              = typejoin(typeof(hash(nothing)),
                                     typeof(objectid(nothing)))
const HASH_ZERO           = HASH_T(0x00)
const HASH_ONE            = HASH_T(0x01)
const HASH_MAX            = ~HASH_ZERO
# Number of bits in the total addressable hash-space of a PTree.
const HASH_BITCOUNT       = sizeof(HASH_T) * 8
const HASH_TOPBIT         = HASH_T(0x01) << (HASH_BITCOUNT - 1)
#   We use a constant shift of 5 throughout except at the root node (which is
#   Due to the implemeentation details, the twig shifts *MUST* be >= 3 and all
#   shifts must be <= 7.
const PTREE_NODE_SHIFT    = 6
const PTREE_TWIG_SHIFT    = 6
const PTREE_ROOT_SHIFT    = (let tmp = rem(HASH_BITCOUNT - PTREE_TWIG_SHIFT,
                                           PTREE_NODE_SHIFT)
                                 (tmp == 0 ? PTREE_NODE_SHIFT : tmp)
                             end)
const PTREE_ROOT_FIRSTBIT = HASH_BITCOUNT - PTREE_ROOT_SHIFT
const PTREE_ROOT_NCHILD   = (1 << PTREE_ROOT_SHIFT)
const PTREE_NODE_NCHILD   = (1 << PTREE_NODE_SHIFT)
const PTREE_TWIG_NCHILD   = (1 << PTREE_TWIG_SHIFT)
const PTREE_BITS_BITCOUNT = max(PTREE_ROOT_NCHILD,
                                PTREE_NODE_NCHILD,
                                PTREE_TWIG_NCHILD)
const PTREE_BITS_T        = (if     (PTREE_BITS_BITCOUNT <= 8)  UInt8
                             elseif (PTREE_BITS_BITCOUNT <= 16) UInt16
                             elseif (PTREE_BITS_BITCOUNT <= 32) UInt32
                             elseif (PTREE_BITS_BITCOUNT <= 64) UInt64
                             else                               UInt128
                             end)
const BITS_ZERO           = PTREE_BITS_T(0x0)
const BITS_ONE            = PTREE_BITS_T(0x1)
const PTREE_NODE_BITS     = HASH_BITCOUNT - PTREE_ROOT_SHIFT - PTREE_TWIG_SHIFT
const PTREE_NODE_LEVELS   = div(PTREE_NODE_BITS, PTREE_NODE_SHIFT)
const PTREE_LEVELS        = PTREE_NODE_LEVELS + 2
const PTREE_ROOT_DEPTH    = 0
const PTREE_TWIG_DEPTH    = PTREE_ROOT_DEPTH + PTREE_NODE_LEVELS + 1
const PTREE_LEAF_DEPTH    = PTREE_TWIG_DEPTH + 1
# Functions for making masks
"""
    lowmask(bitno, T)
    lowmask(bitno)

Yields a mask of (unsigned integer) type T with all bits above the given bit
number set to false and all bits below that number set to true. The bit itself
is set to false. Bits are indexed starting at 0.

The default value of T is the PTree hash type (HASH_T).

lowmask(bitno) is equal to ~highmask(bitno).
"""
lowmask(bitno::K, ::Type{T}) where {T<:Unsigned,K<:Integer} =
    ((T(0x01) << bitno) - 0x01)
lowmask(bitno::K) where {K <: Integer} = lowmask(bitno, HASH_T)
"""
    highmask(bitno, T)
    highmask(bitno)

Yields a mask of (unsigned integer) type T with all bits above the given bit
number set to true and all bits below that number set to false. The bit itself
is set to true. Bits are indexed starting at 0.

The default value of T is the PTree hash type (HASH_T).

highmask(bitno) is equal to ~lowmask(bitno).
"""
highmask(bitno::K, ::Type{T}) where {T<:Unsigned,K<:Integer} =
    ~lowmask(bitno, T)
highmask(bitno::K) where {K <: Integer} = highmask(bitno, HASH_T)
# The mask of the bits that represent the depth of the node-id in a PTree.
const PTREE_DEPTH_MASK = lowmask(PTREE_TWIG_SHIFT)
# The pair type for the dict interface.
const HASHPAIR_T{T} = Pair{HASH_T, T} where {T}

# # Functions
"""
    ptree_depth(nodeaddr)

Yields the depth of the node with the given node address. This depth is in the
theoretical complete tree, not in the reified tree represented with memory.
"""
ptree_depth(nodeid::HASH_T) = nodeid & PTREE_DEPTH_MASK
"""
    depth_to_bitshift(depth)

Yields the tuple (B0,S) of the first bit index and the shift for the given
depth.
"""
depth_to_bitshift(depth::Integer) = depth_to_bitshift(Int(depth))
depth_to_bitshift(depth::Int) = begin
    if depth == PTREE_TWIG_DEPTH
        return (0, PTREE_TWIG_SHIFT)
    elseif depth == 0
        return (PTREE_ROOT_FIRSTBIT, PTREE_ROOT_SHIFT)
    else
        return (PTREE_ROOT_FIRSTBIT - depth*PTREE_NODE_SHIFT,
                PTREE_NODE_SHIFT)
    end
end
"""
    ptree_bitshift(nodeid)

Yields the bitshift for the given ptree node id.
"""
ptree_bitshift(nodeid::HASH_T) = depth_to_bitshift(ptree_depth(nodeid))
"""
    ptree_shift(nodeid)

Yields the shift for the given ptree node id.
"""
ptree_shift(nodeid::HASH_T) = ptree_bitshift(nodeid)[2]
"""
    ptree_firstbit(nodeid)

Yields the first-bit for the given ptree node id.
"""
ptree_firstbit(nodeid::HASH_T) = ptree_bitshift(nodeid)[1]
"""
    ptree_minleaf(nodeid)

Yields the minimum child leaf index associated with the given nodeid.
"""
ptree_minleaf(nodeid::HASH_T) = nodeid & ~PTREE_DEPTH_MASK
"""
    ptree_maxleaf(nodeid)

Yields the maximum child leaf index assiciated with the given nodeid.
"""
ptree_maxleaf(nodeid::HASH_T) = begin
    (b0,s) = ptree_bitshift(nodeid)
    mask = (HASH_ONE << (b0 + s)) - HASH_ONE
    return nodeid | mask
end
"""
    ptree_minmaxleaf(nodeid)

Yields the (min, max) child leaf index assiciated with the given nodeid.
"""
ptree_minmaxleaf(nodeid::HASH_T) = begin
    mn = ptree_minleaf(nodeid)
    (bit0,shift) = ptree_bitshift(nodeid)
    mask = (HASH_ONE << (bit0+shift)) - HASH_ONE
    return (nodeid & ~mask, nodeid | mask)
end
"""
    ptree_id(minleaf, depth)

Yields the node-id for the node whose minimum leaf and depth are given.
"""
ptree_id(minleaf::Integer,depth::Integer) = 
    ptree_id(HASH_T(minleaf), HASH_T(depth))
ptree_id(minleaf::HASH_T, depth::HASH_T) = minleaf | depth
"""
    ptree_parentid(nodeid0)

Yields the node-id of the parent of the given node. Note that node 0 (the tree's
theoretical root) has no parent. If given a node id of 0, this function will
return an arbitrary large number.
"""
ptree_parentid(nodeid0::HASH_T) = begin
    d = ptree_depth(nodeid0) - 0x1
    node = nodeid0 & ~((-HASH_T(1)) >> (d*PTREE_NODE_SHIFT))
    return node | d
end
"""
    ptree_isbeneath(nodeid, leafid)

Yields true if the given leafid can be found beneath the given node-id.
"""
ptree_isbeneath(nodeid::HASH_T, leafid::HASH_T) = begin
    (mn,mx) = ptree_minmaxleaf(nodeid)
    return leafid <= mx && leafid >= mn
end
"""
    ptree_highbitdiff(id1, id2)

Yields the highest bit that is different between id1 and id2.
"""
ptree_highbitdiff(id1::HASH_T, id2::HASH_T) =
    HASH_BITCOUNT - leading_zeros(xor(id1, id2)) - 1


# ==============================================================================
# Structures / Types
#
# Definitions of structures and types that should come before the various
# methods and constructors that operate on them.
# #PTree =======================================================================
"""
    PTree{T}

PTree is a persistent tree/array hybrid that maps arbitrary unsigned integers
(type HASH_T) to values. It can efficiently be used for arrays or for
hash-maps, and supports efficient lookup, association, and dissociation.

All fields of a PTree should be considered strictly private, as modification
of the fields may result in a kernel crash. The properties of a PTree are both
immutable and safe to inpect.
"""
struct PTree{T} <: AbstractDict{HASH_T, T}
    id::HASH_T
    bits::PTREE_BITS_T
    numel::Int
    cells::Union{Nothing, Vector{T}, Vector{PTree{T}}}
end
const PTREE_NODE_CELL_T{T} = Vector{PTree{T}} where {T}
const PTREE_TWIG_CELL_T{T} = Vector{T} where {T}
const PTREE_EMPTY_CELL_T = Nothing
const PTREE_CELLS_T{T} = Union{PTREE_EMPTY_CELL_T,
                               PTREE_NODE_CELL_T{T},
                               PTREE_TWIG_CELL_T{T}} where {T}
# We will occasionally want to do something across all the children of a node;
# we can exploit the bits data to speed this up in sparse cases. The code for
# doing that is in these two macros.
macro _ptree_forcells(tree, var::Symbol, expr::Expr)
    b = gensym("bits")
    return quote
        let $b = $(tree).bits, $var
            while true
                $var = trailing_zeros($b)
                if $var == PTREE_BITS_BITCOUNT
                    break
                else
                    $expr
                end
                $b &= ~($(HASH_T(0x01)) << $var)
            end
        end
    end |> esc
end
macro _ptree_forcells_from(tree, var::Symbol, from, expr::Expr)
    b = gensym("bits")
    f = gensym("from")
    return quote
        let $b = $(tree).bits & highmask($from), $var
            while true
                $var = trailing_zeros($b)
                if $var == PTREE_BITS_BITCOUNT
                    break
                else
                    $expr
                end
                b &= ~($(HASH_T(0x01)) << $var)
            end
        end
    end |> esc
end
"""
    ptree_cellindex(ptree, leafid)
    ptree_cellindex(nodeid, bits, leafid)
Yields a tuple (present, bitindex, cellindex) in which [1] present is a boolean
indicating whether a ptree containing the given leafid or the leaf itself is a
child of the given ptree; [2] bitindex is the index into ptree's bits integer 
for the given leaf, and cellindex is the inndex innto the ptree's cell array 
where that child is or would be found. If the leafid is outside of the given
ptree (i.e., it cannot exist beneath this ptree) then the cellindex returned
returned is 0, but the bitindex will still match the appropriate shift for the
ptree's depth.
"""
ptree_cellindex(id::HASH_T, bits::PTREE_BITS_T, leafid::HASH_T) = begin
    # Check that the leaf is below this leaf.
    (bit0, shift) = ptree_bitshift(id)
    # Grab the index out of the leaf id.
    idx = (leafid >> bit0) & lowmask(shift)
    ptree_isbeneath(id, leafid) || return (false, idx, 0)
    # See if the bit is set?
    flag = HASH_ONE << idx
    chno = 1 + count_ones(bits & lowmask(idx))
    return (bits & flag == flag, idx, chno)
end
ptree_cellindex(u::PTree{T}, leafid::HASH_T) where {T} =
    ptree_cellindex(getfield(u, :id), getfield(u, :bits), leafid)
"""
    ptree_cellindex!(ptree, leafid)
    ptree_cellindex!(nodeid, bits, leafid)
Yields the index into the ptree's cell vector of the child containing the leaf
with the given id. If the ptree does not contain the given leafid because the
appropriate bit is not set then 0 is returned. However, unlike the function
ptree_cellindex(), this function does not check whether or not the given leafid
is in the set of possible children of the tree.
"""
ptree_cellindex!(id::HASH_T, bits::PTREE_BITS_T, leafid::HASH_T) = begin
    # Check that the leaf is below this leaf.
    (bit0, shift) = ptree_bitshift(id)
    # Grab the index out of the leaf id.
    idx = (leafid >> bit0) & lowmask(shift)
    # See if the bit is set?
    flag = HASH_ONE << idx
    chno = 1 + count_ones(bits & lowmask(idx))
    return (bits & flag == flag, idx, chno)
end
ptree_cellindex!(ptree::PTree{T}, leafid::HASH_T) where {T} =
    ptree_cellindex!(getfield(ptree, :id), getfield(ptree, :bits), leafid)
"""
    ptree_cellkey(ptree, childidx)
Yields the leafid (a HASH_T value) of the key that goes with the particular
child index that is given. This only works correctly for twig nodes.
"""
ptree_cellkey(id::HASH_T, k::HASH_T) = begin
    mn = ptree_minleaf(id)
    return mn | k
end
ptree_cellkey(id::HASH_T, k::II) where {II<:Integer} =
    ptree_cellkey(id, (HASH_T(k)))
ptree_cellkey(u::PTree{T}, k::II) where {T,II<:Integer} =
    ptree_cellkey(getfield(u, :id), HASH_T(k))
ptree_cellkey(u::PTree{T}, k::HASH_T) where {T} =
    ptree_cellkey(getfield(u, :id), k)

# # Constructors ===============================================================
#   First, constructors for making empty Trees.
PTree{T}() where {T} = PTree{T}(0x0, 0x0, 0, nothing)
#   Now for copying trees of the same type.
PTree{T}(u::PTree{T}) where {T}   = u
PTree{T}(u::PTree{U}) where {T,U} = begin
    (id, bits, n, cells) = (getfield(u, :id),    getfield(u, :bits),
                            getfield(u, :numel), getfield(u, :cells))
    d = ptree_depth(id)
    if d == PTREE_TWIG_DEPTH
        cells::Vector{U}
        return PTree{T}(id, bits, n, convert(Vector{T}, cells))
    else
        cells::Vector{PTree{U}}
        return PTree{T}(id, bits, n, PTree{T}[PTree{T}(c) for c in cells])
    end
end
PTree{T}(d::AbstractDict{HASH_T,S}) where {T,S} = begin
    u = PTree{T}()
    for (k,v) in d
        u = setindex(u, v, k)
    end
    return u
end
#   We want to be able to construct the twig for a pair:
PTree{T}(kv::Pair{HASH_T,V}) where {T,V} = begin
    id = ptree_id(kv[1] & ~PTREE_DEPTH_MASK, PTREE_TWIG_DEPTH)
    bits = PTREE_BITS_T(0x1) << (kv[1] & PTREE_DEPTH_MASK)
    return PTree{T}(id, bits, 1, T[kv[2]])
end
PTree(kv::Pair{HASH_T,V}) where {V} = PTree{V}(kv)
#   We also construct from arrays; they are considered to be 0-indexed when used
#   with trees.
PTree{T}(u::AbstractArray{S,N}) where {T,S,N} = begin
    depth = PTREE_TWIG_DEPTH
    height = 1
    n = length(u)
    (n == 0) && return PTree{T}()
    (n == 1) && return PTree{T}(ptree_id(0x0, depth), 0x1, 1, T[u[1]])
    # Start by making all the twigs
    maxidx = fld(n + PTREE_TWIG_NCHILD - 1, PTREE_TWIG_NCHILD)
    bits = ~PTREE_BITS_T(0x0)
    nodes = Vector{PTree{T}}(undef, maxidx)
    for ii in 1:maxidx 
        kk0 = (ii - 1)*PTREE_TWIG_NCHILD
        nn = min(PTREE_TWIG_NCHILD, n - kk0)
        nid = ptree_id(PTREE_TWIG_NCHILD * (ii-1), depth)
        ch = Vector{T}(undef, PTREE_TWIG_NCHILD)
        if nn == PTREE_TWIG_NCHILD
            for kk in 1:PTREE_TWIG_NCHILD
                @inbounds ch[kk] = u[kk + kk0]
            end
            nodes[ii] = PTree{T}(nid, bits, Int(PTREE_TWIG_NCHILD), ch)
        else
            for kk in 1:nn
                @inbounds ch[kk] = u[kk + kk0]
            end
            nodes[ii] = PTree{T}(nid, lowmask(nn, PTREE_BITS_T), Int(nn), ch)
        end
    end
    # Then iterate over the nodes doing the same thing.
    while true
        # Even if we are given an array with the max number of elements, we must
        # eventually have a list of only one node; that's the node we return.
        (maxidx == 1) && return nodes[1]
        # Increment
        depth -= 1
        height += 1
        n = length(nodes)
        oldnodes = nodes
        maxidx = div(n + PTREE_NODE_NCHILD - 1, PTREE_NODE_NCHILD)
        nodes = Vector{PTree{T}}(undef, maxidx)
        hnp = PTREE_NODE_NCHILD^(height-1) # number per node at the prev level
        hn  = Int(PTREE_NODE_NCHILD*hnp) # number per node at this level
        bits = ~PTREE_BITS_T(0x0)
        for ii in 1:maxidx
            kk0 = (ii - 1)*PTREE_NODE_NCHILD
            nn = min(PTREE_NODE_NCHILD, n - kk0)
            nid = ptree_id(hn * (ii-1), depth)
            ch = Vector{PTree{T}}(undef, nn)
            if nn == PTREE_NODE_NCHILD
                # split this out for faster looping
                for kk in 1:PTREE_NODE_NCHILD
                    @inbounds ch[kk] = (@inbounds oldnodes[kk + kk0])
                end
                nodes[ii] = PTree{T}(nid, bits, hn - hnp + getfield(ch[end], :numel), ch)
            else
                for kk in 1:nn
                    @inbounds ch[kk] = (@inbounds oldnodes[kk + kk0])
                end
                nodes[ii] = PTree{T}(nid, (PTREE_BITS_T(0x1) << nn) - 0x1,
                                     (nn - 1)*hnp + getfield(ch[end], :numel), ch)
            end
        end
    end
end
# Finally, we want to duplicate all of these constructors using the convert
Base.convert(::Type{PTree{T}}, u::PTree{T}) where {T} = u
Base.convert(::Type{PTree{T}}, u::PTree{U}) where {T,U} = PTree{T}(u)
Base.convert(::Type{PTree{T}}, u::AbstractDict{HASH_T,U}) where {T,U} =
    PTree{T}(u)
Base.convert(::Type{PTree{T}}, u::AbstractArray{1,U}) where {T,U} =
    PTree{T}(u)

# ==============================================================================
# Methods
#
# Definitions of methods of the above type.
Base.empty(u::PTree{T}) where {T} = PTree{T}(r)
Base.empty(u::PTree{T}, S::Type) where {T} = PTree{S}()
Base.isempty(u::PTree{T}) where {T} = (getfield(u, :numel) == 0)
Base.length(u::PTree{T}) where {T} = getfield(u, :numel)
Base.isequal(t::PTree{T}, s::PTree{S}) where {T,S} = begin
    bits = getfield(t, :bits)
    (bits == getfiield(s, :bits)) || return false
    (getfield(t, :id) == getfield(s, :id)) || return false
    tcells = getfield(t, :cells)
    scells = getfield(s, :cells)
    for k in 1:count_ones(bits)
        t = @inbounds tcells[k]
        s = @inbounds scells[k]
        isequal(t, s) || return false
    end
    return true
end
Base.copy(t::PTree{T}) where {T} = t
Base.propertynames(u::PTree{T}, private::Bool=false) where {T} =
    (:leafcount,
     :children,
     :minleaf,
     :maxleaf,
     :depth,
     :address)
Base.getproperty(u::PTree{T}, symbol::Symbol) where {T} = begin
    if symbol == :leafcount
        return getfield(u, :numel)
    elseif symbol == :children
        id = getfield(u, :id)
        (getfield(u, :numel) == 0) && return ()
        if ptree_depth(id) == PTREE_TWIG_DEPTH
            return tuple(u...)
        else
            cells = getfield(u, :cells)::Vector{PTree{T}}
            return tuple(cells...)
        end
    elseif symbol == :minleaf
        return ptree_minleaf(getfield(u, :id))
    elseif symbol == :maxleaf
        return ptree_maxleaf(getfield(u, :id))
    elseif symbol == :depth
        return ptree_depth(getfield(u, :id))
    elseif symbol == :address
        return getfield(u, :id)
    else
        throw(ArgumentError("Type PTree has no property $symbol"))
    end
end
getpair(u::PTree{T}, k::HASH_T) where {T} = begin
    # Start by making sure that the address spaces are compatible:
    # we want to make sure (once) that k is beneath u.
    id = getfield(u, :id)
    ptree_isbeneath(id, k) || return missing
    # Also, if we are empty, we need to return right away.
    (getfield(u, :numel) == 0) && return missing
    # Okay, let's try to find the child.
    while true
        # Is it in this Tree?
        (inq, bitidx, idx) = ptree_cellindex!(u, k)
        inq || return missing
        # Are we a twig?
        #if isa(getfield(u, :cells), Vector{T})
        if ptree_depth(id) == PTREE_TWIG_DEPTH
            qq = getfield(u, :cells)::Vector{T}
            return k => (@inbounds qq[idx])
        else
            q = getfield(u, :cells)::Vector{PTree{T}}
            u = (@inbounds q[idx])
        end
        id = getfield(u, :id)
    end
end
Base.get(u::PTree{T}, k::HASH_T, df) where {T} = begin
    # Start by making sure that the address spaces are compatible:
    # we want to make sure (once) that k is beneath u.
    id = getfield(u, :id)
    ptree_isbeneath(id, k) || return df
    # Also, if we are empty, we need to return right away.
    (getfield(u, :numel) == 0) && return df
    # Okay, let's try to find the child.
    while true
        # Is it in this Tree?
        (inq, bitidx, idx) = ptree_cellindex(u, k)
        inq || return df
        # Are we a twig?
        #if isa(getfield(u, :cells), Vector{T})
        if ptree_depth(id) == PTREE_TWIG_DEPTH
            qq = getfield(u, :cells)::Vector{T}
            return (@inbounds qq[idx])
        else
            q = getfield(u, :cells)::Vector{PTree{T}}
            u = (@inbounds q[idx])
        end
        id = getfield(u, :id)
    end
end
Base.in(kv::Pair{HASH_T,T}, u::PTree{T}, f::F) where {T,F<:Function} = begin
    id = getfield(u, :id)
    ptree_isbeneath(id, k) || return df
    # Also, if we are empty, we need to return right away.
    (getfield(u, :numel) == 0) && return df
    # Okay, let's try to find the child.
    while true
        # Is it in this Tree?
        (inq, bitidx, idx) = ptree_cellindex!(u, k)
        inq || return df
        # Are we a twig?
        #if isa(getfield(u, :cells), Vector{T})
        if ptree_depth(id) == PTREE_TWIG_DEPTH
            qq = getfield(u, :cells)::Vector{T}
            return f(kv[2], (@inbounds qq[idx])) == true
        else
            q = getfield(u, :cells)::Vector{PTree{T}}
            u = (@inbounds q[idx])
        end
        id = getfield(u, :id)
    end
end
Base.in(kv::Pair{HASH_T,T}, u::PTree{T}) where {T} = in(kv, u, (===))

Base.iterate(u::PTree{T}) where {T} = begin
    (getfield(u, :numel) == 0) && return nothing
    id = getfield(u, :id)
    d = ptree_depth(id)
    while d < PTREE_TWIG_DEPTH
        cells = getfield(u, :cells)::Vector{PTree{T}}
        u = cells[1]
        id = getfield(u, :id)
        d = ptree_depth(id)
    end
    leaves = getfield(u, :cells)::Vector{T}
    bits = getfield(u, :bits)
    bitno = trailing_zeros(bits)
    k = ptree_cellkey(id, bitno)
    v = leaves[1] # @inbounds leaves[1]
    return (Pair{HASH_T,T}(k,v), k)
end
iterkeys(u::PTree{T}) where {T} = begin
    (getfield(u, :numel) == 0) && return nothing
    id = getfield(u, :id)
    d = ptree_depth(id)
    while d < PTREE_TWIG_DEPTH
        cells = getfield(u, :cells)::Vector{PTree{T}}
        u = cells[1]
        id = getfield(u, :id)
        d = ptree_depth(id)
    end
    leaves = getfield(u, :cells)::Vector{T}
    bits = getfield(u, :bits)
    bitno = trailing_zeros(bits)
    k = ptree_cellkey(id, bitno)
    v = leaves[1] # @inounds leaves[1]
    return (k, k)
end
itervals(u::PTree{T}) where {T} = begin
    (getfield(u, :numel) == 0) && return nothing
    id = getfield(u, :id)
    d = ptree_depth(id)
    while d < PTREE_TWIG_DEPTH
        cells = getfield(u, :cells)::Vector{PTree{T}}
        u = cells[1]
        id = getfield(u, :id)
        d = ptree_depth(id)
    end
    leaves = getfield(u, :cells)::Vector{T}
    bits = getfield(u, :bits)
    bitno = trailing_zeros(bits)
    k = ptree_cellkey(id, bitno)
    v = leaves[1] # @inounds leaves[1]
    return (v, k)
end

macro _ptree_iterate_gencode(rtype::Symbol)
    cids = [gensym("cid") for _ in 1:PTREE_LEVELS-1]
    cels = [gensym("cel") for _ in 1:PTREE_LEVELS-1]
    bits = [gensym("bit") for _ in 1:PTREE_LEVELS-1]
    if rtype == :keys
        rexpr = :(k)
        fn = :iterkeys
    elseif rtype == :vals
        rexpr = :(v)
        fn = :itervals
    else
        rexpr = :(Pair{HASH_T,T}(k,v))
        fn = :(Base.iterate)
    end
    twigexpr = quote
        # We must be a twig at this point.
        leaves = getfield(u, :cells)::Vector{T}
        id = getfield(u, :id)
        bits = getfield(u, :bits)
        (inq, bitno, idx) = ptree_cellindex(id, bits, k0)
        mask = (BITS_ONE << (bitno + 1)) - BITS_ONE
        nextbitno = trailing_zeros(bits & ~mask)
        if nextbitno < PTREE_BITS_BITCOUNT
            k = ptree_cellkey(id, nextbitno)
            idx = count_ones(bits & mask) + 1
            v = (@inbounds leaves[idx])
            return ($rexpr, k)
        end
        # If we reach here, we did not find a next leaf.
    end
    # At the lowest level, the try and twig expressions are the same.
    tryexpr = twigexpr
    for ii in PTREE_LEVELS-1:-1:1
        cel = cels[ii]
        cid = cids[ii]
        bit = bits[ii]
        tryexpr = quote
            if d == PTREE_TWIG_DEPTH
                $twigexpr
            else
                $cel = getfield(u, :cells)::Vector{PTree{T}}
                $bit = getfield(u, :bits)
                (inq, bitno, $cid) = ptree_cellindex(id, $bit, k0)
                u = $cel[$cid]
                id = getfield(u, :id)
                d = ptree_depth(id)
                $tryexpr
                # If we reach this point, we haven't found a next pair yet.
                n = count_ones($bit)
                if $cid < n
                    u = $cel[$cid + 1]
                    break
                end
            end
        end
    end
    # Put these together into the functions.
    return quote
        $fn(u::PTree{T}, k0::HASH_T) where {T} = begin
            id = getfield(u, :id)
            d = ptree_depth(id)
            while true
                $tryexpr
                # If we reach the end of the loop, there's nothing past k0.
                return nothing
            end
            # Upon finding a cell that hasn't been iterated, the code breaks
            # from the abbove loop and brings us here.
            id = getfield(u, :id)
            d = ptree_depth(id)
            while d < PTREE_TWIG_DEPTH
                cells = getfield(u, :cells)::Vector{PTree{T}}
                u = (@inbounds cells[1])
                id = getfield(u, :id)
                d = ptree_depth(id)
            end
            leaves = getfield(u, :cells)::Vector{T}
            bits = getfield(u, :bits)
            nextbitno = trailing_zeros(bits)
            k = ptree_cellkey(id, nextbitno)
            v = (@inbounds leaves[1])
            return ($rexpr, k)
        end
    end |> esc
end
# Run the macro to generate the iteration functions:
(@_ptree_iterate_gencode keys)
(@_ptree_iterate_gencode vals)
(@_ptree_iterate_gencode pairs)

setindex(u::PTree{T}, v::V, k::HASH_T) where {T,V} = begin
    # First of all, if this node is empty, we just return a new node
    (getfield(u, :numel) == 0) && return PTree{T}(k => v)
    # Next, if k is a child of this node, we return a new ancestor node that
    # contains us both and is the highest possible ancestor that contains us
    # both.
    id = getfield(u, :id)
    if !ptree_isbeneath(id, k)
        # Make the new twig that will hold k => v.
        kv = PTree{T}(k => v)
        # What's the lowest bit at which the mismatch in address starts?
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
        newid = ptree_id(newminleaf, newdepth)
        newmask = (HASH_ONE << shift) - HASH_ONE
        ii_k = newmask & (id >> bit0)
        ii_u = newmask & (k >> bit0)
        newbits = (BITS_ONE << ii_k) | (BITS_ONE << ii_u)
        newchs = ii_u < ii_k ? PTree{T}[kv, u] : PTree{T}[u, kv]
        return PTree{T}(newid, newbits, 1 + getfield(u, :numel), newchs)
    end
    # Otherwise, this key goes beneath this node, so find out where.
    bits = getfield(u, :bits)
    numel = getfield(u, :numel)
    (bit0, shift) = ptree_bitshift(id)
    # Grab the index out of the leaf id.
    idx = (k >> bit0) & lowmask(shift)
    # We'll want to see if this bit is set and then insert/edit at this index:
    flag = BITS_ONE << idx
    idx = 1 + count_ones(bits & lowmask(idx))
    # However, what we do next depends on whether we are in a twig node or not.
    if ptree_depth(id) == PTREE_TWIG_DEPTH    
        cells = getfield(u, :cells)::Vector{T}
        if bits & flag == flag
            # We are replacing an item
            return PTree{T}(id, bits, numel, setindex(cells, v, idx))
        else
            # If we get here, we are adding a new node, which means we are
            # inserting a leaf.
            return PTree{T}(id, bits | flag, numel + 1, insert(cells, idx, v))
        end
    else
        cells = getfield(u, :cells)::Vector{PTree{T}}
        if bits & flag == flag
            # We pass control along then return a copy of ourselves on the
            # way back up the control stack.
            oldc = cells[idx]
            oldn = getfield(oldc, :numel)
            newc = setindex(cells[idx], v, k)
            (newc === cells[idx]) && return u
            newn = getfield(newc, :numel)
            newcells = setindex(cells, newc, idx)
            return PTree{T}(id, bits, numel - oldn + newn, newcells)
        else
            # We add a twig beneath this node.
            newcells = insert(cells, idx, PTree{T}(k => v))
            return PTree{T}(id, bits | flag, numel + 1, newcells)
        end
    end
end
delete(u::PTree{T}, k::HASH_T) where {T} = begin
    id = getfield(u, :id)
    bits = getfield(u, :bits)
    (inq, bitidx, idx) = ptree_cellindex(id, bits, k)
    # If the node lies outside of this node or the bit isn't set, there's
    # nothing to do.
    inq || return u
    # There's something to delete; start by determining if we are or are not a
    # twig node:
    numel = getfield(u, :numel)
    # Check if we're a twig or something else.
    if ptree_depth(id) == PTREE_TWIG_DEPTH
        # If we are deleting the only element, we return a new empty node.
        (numel == 1) && return PTree{T}()
        cells = getfield(u, :cells)::Vector{T}
        newcells = delete(cells, idx)
        return PTree{T}(id, bits & ~(BITS_ONE << bitidx), numel-1, newcells)
    else
        cells = getfield(u, :cells)::Vector{PTree{T}}
        oldc = cells[idx]
        oldn = getfield(oldc, :numel)
        newc = delete(oldc, k)
        newn = getfield(newc, :numel)
        numel += newn - oldn
        (numel == 0) && return newc
        if newn == 0
            newcells = delete(cells, idx)
            return PTree{T}(id, bits & ~(BITS_ONE << bitidx), numel, newcells)
        else
            newcells = setindex(cells, newc, idx)
            return PTree{T}(id, bits, numel, newcells)
        end
    end
end
