################################################################################
# PTree.jl
# Code implementing a persistent tree type on which the persistent collection
# library is built. This type is not intended for use outside of the library.
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
#   We use a constant shift of 4 throughout; we pick this across all array
#   dimensionalities and for all dimensions because it simplifies slicing not to
#   have to worry about arrays having mismatched shifts.
#   Due to the implemeentation details, the shift *MUST* be >= 4 and <= 7.
const _PTREE_SHIFT    = 4
const _PTREE_COUNT    = (1 << _PTREE_SHIFT)
const _PTREE_BITS_T   = (_PTREE_COUNT <= 8  ? UInt8  :
                         _PTREE_COUNT <= 16 ? UInt16 :
                         _PTREE_COUNT <= 32 ? UInt32 :
                         _PTREE_COUNT <= 64 ? UInt64 :
                         UInt128)
#   This is the type of the unsigned integer used for the leaves.
const _PTREE_KEY_T    = typejoin(typeof(hash(nothing)),
                                 typeof(equivhash(nothing)),
                                 typeof(objectid(nothing)))
const _PTREE_KEY_MAX  = _PTREE_KEY_T(0x00) - _PTREE_KEY_T(0x01)
# Number of bits in the total addressable hash-space of a PTree.
const _PTREE_KEY_BITS = sizeof(_PTREE_KEY_T) * 8
# The depth of the leaves and twigs in the tree.
const _PTREE_LEAF_DEPTH = _PTREE_KEY_T(div(_PTREE_KEY_BITS + _PTREE_SHIFT - 1, _PTREE_SHIFT))
const _PTREE_TWIG_DEPTH = _PTREE_KEY_T(_PTREE_LEAF_DEPTH - 1)
# The mask of the bits that represent the depth of the node-id in a PTree.
const _PTREE_DEPTH_MASK = ((_PTREE_KEY_T(0x1) << (_PTREE_SHIFT)) - 1)
# The pair type for the dict interface.
const _PTREE_PAIR_T{T} = Pair{_PTREE_KEY_T, T} where {T}

# # Functions
"""
    _node_depth(nodeid)
Yields the depth of the node with the given nodeid. This depth is in the
theoretical complete tree, not in the reified tree represented with memory.
"""
_node_depth(nodeid::_PTREE_KEY_T) = nodeid & _PTREE_DEPTH_MASK
"""
    _node_min(nodeid)
Yields the minimum child leaf index associated with the given nodeid.
"""
_node_min(nodeid::_PTREE_KEY_T) = nodeid & ~_PTREE_DEPTH_MASK
"""
    _node_max(nodeid)
Yields the maximum child leaf index assiciated with the given nodeid.
"""
_node_max(nodeid::_PTREE_KEY_T) = begin
    h = _PTREE_LEAF_DEPTH - _node_depth(nodeid)
    mask = (_PTREE_KEY_T(1) << (_PTREE_SHIFT * h)) - 0x1
    return nodeid | mask
end
"""
    _nodeid(minleaf, depth)
Yields the node-id for the node whose minimum leaf and depth are given.
"""
_nodeid(minleaf::_PTREE_KEY_T, depth::_PTREE_KEY_T) = begin
    return (depth & _PTREE_DEPTH_MASK) | ((minleaf + 0x1) & ~_PTREE_DEPTH_MASK)
end
"""
    _node_child(nodeid0, childno)
Yields the node-id of the child with the child-number `childno` of the node 
whose node-id is `nodeid0`.
"""
_node_child(nodeid0::_PTREE_KEY_T, childno::UInt) = begin
    d = _node_depth(nodeid0) + 0x1
    return _nodeid(_node_min(nodeid0) | (childno << d), d)
end
"""
    _node_parent(nodeid0)
Yields the node-id of the parent of the given node. Note that node 0 (the tree's
theoretical root) has no parent. If given a node id of 0, this function will
return an arbitrary large number.
"""
_node_parent(nodeid0::_PTREE_KEY_T) = begin
    d = _node_depth(nodeid0) - 0x1
    node = nodeid0 & ~((-_PTREE_KEY_T(1)) >> (d*_PTREE_SHIFT))
    return node | d
end
"""
    _node_index(nodeid, leafid)
Yields the index into the children of the given nodeid at which the given leafid
can be found; if it cannot be found in the node's children, yields 0.
"""
_node_index(nodeid::_PTREE_KEY_T, leafid::_PTREE_KEY_T) = begin
    d = _node_depth(nodeid)
    mn = _node_min(nodeid)
    (leafid < mn) && return 0
    h = _PTREE_LEAF_DEPTH - d
    one = _PTREE_KEY_T(0x1)
    mask = (one << (_PTREE_SHIFT * h)) - 0x1
    mx = nodeid | mask
    (leafid > mx) && return 0
    mask = (one << _PTREE_SHIFT) - 1
    shift = _PTREE_SHIFT * (_PTREE_TWIG_DEPTH - d)
    return 1 + Int((leafid & (mask << shift)) >> shift)
end

# ==============================================================================
# Structures / Types
#
# Definitions of structures and types that should come before the various
# methods and constructors that operate on them.
#
# #_PNought ====================================================================
"""
    _PNought()
_PNought is an empty type that is used to indicate that an element of a PTree is
empty, internally. It shouldn't be used outside of the PTree library code.
"""
struct _PNought end
const _nought = _PNought()
# The types for the leaves and cells of PTree and TTree. Because we're
# declaring these before the struct, we force the stuct to be a type parameter.
# This is just hackerty in order to be able to delcare the types of these only
# once; we declare aliases for them below.
const _PTREE_CELL_TYPE{PT,T} = Union{_PNought, T, PT} where {PT,T}
const _TTREE_CELL_TYPE{PT,TT,T} = Union{_PNought,T,TT,PT} where {PT,TT,T}
# # Types of the containers.
#   The container is the type (e.g., Array or Tuple) that holds the list of
#   nodes for the tree. There are advantages to making this a somewhat flexible
#   feature of the implementation for now, so I've parameterized it here.
const _PTREE_CELLS_TYPE{PT, T} =
    Vector{_PTREE_CELL_TYPE{PT,T}} where {PT,T}
const _TTREE_CELLS_TYPE{PT, TT, T} =
    Vector{_TTREE_CELL_TYPE{PT,TT,T}} where {PT,TT,T}
# This is the NTuple that we convert into every empty PTree's _data field. Using
# just this one object is an easy way to reduce the memory footprint of empty
# trees.
const _ptree_noughts =
    NTuple{_PTREE_COUNT, _PNought}([_nought for _ in 1:_PTREE_COUNT])
# #PTree =======================================================================
"""
    PTree{T}
PTree is a persistent tree/array hybrid that maps arbitrary unsigned integers
(type _PTREE_KEY_T) to values. It can efficiently be used for arrays or for
hash-maps, and supports efficient lookup, association, and dissociation.

`PTree{T}()` yields an empty tree.
"""
struct PTree{T} <: AbstractDict{_PTREE_KEY_T, T}
    _n::_PTREE_KEY_T
    _id::_PTREE_KEY_T
    _bits::_PTREE_BITS_T
    _data::_PTREE_CELLS_TYPE{PTree{T}, T}
end
# #TTee ========================================================================
"""
    TTree{T}
TTree is a transient tree/array hybrid that maps arbitrary unsigned integers
(type _PTREE_KEY_T) to values. It is identical in most ways to PTree but that
it is mutable. However, TTree and PTree have similar underlying structures
such that a TTree can be constructed from a PTree in O(1) time, and a PTree
can be constructed from a TTree in O(n log n) time where n is the number of
changes made to the TTree since it was instantiated from a (potentially empty)
PTree.
"""
struct TTree{T} <: AbstractDict{_PTREE_KEY_T, T}
    _n::_PTREE_KEY_T
    _id::_PTREE_KEY_T
    _bits::_PTREE_BITS_T
    _data::_TTREE_CELLS_TYPE{PTree{T}, TTree{T}, T}
end
# These simply alias the _TYPE types above, for convenience now that we have
# defined both PTree and TTree.
const _PTREE_CELL_T{T}  = _PTREE_CELL_TYPE{PTree{T},T}
const _PTREE_CELLS_T{T} = _PTREE_CELLS_TYPE{PTree{T},T}
const _TTREE_CELL_T{T}  = _TTREE_CELL_TYPE{PTree{T},TTree{T},T}
const _TTREE_CELLS_T{T} = _TTREE_CELLS_TYPE{PTree{T},TTree{T},T}
# These are private functions used by PTree and TTree to manage its children.
_ptree_empty_cells(::Type{T}) where {T} = _PTREE_CELL_T{T}[]
_ttree_empty_cells(::Type{T}) where {T} = _TTREE_CELL_T{T}[]
_ptree_cell_assoc(u::_PTREE_CELLS_T{T}, k::Int, x::_PTREE_CELL_T{S}) where {T,S<:T} = begin
    if length(u) < k
        n = length(u)
        u = _PTREE_CELL_T{T}[u..., [_nought for _ in 1:(k - n - 1)]..., x]
        return u
    elseif u[k] === x
        return u
    else
        u = copy(u)
        u[k] = x
        return u
    end
end
_ptree_cell_dissoc(u::_PTREE_CELLS_T{T}, k::Int) where {T} = begin
    n = length(u)
    if n < k || u[k] === _nought
        return u
    elseif n == k
        return u[1:end-1]
    else
        u = copy(u)
        u[k] = _nought
        return u
    end
end
_ttree_cell_assoc(u::_TTREE_CELLS_T{T}, k::Int, x::_PTREE_CELL_T{S}) where {T,S<:T} = begin
    if length(u) < k
        n = length(u)
        u = _TTREE_CELL_T{T}[u..., [_nought for _ in 1:(n - k - 1)]..., x]
        return u
    elseif u[k] === x
        return u
    else
        u = copy(u)
        u[k] = x
        return u
    end
end
_ttree_cell_dissoc(u::_TTREE_CELLS_T{T}, k::Int) where {T} = begin
    u[k] = _nought
    return u
end
# We will occasionally want to do something across all the children of a node;
# we can exploit the bits data to speed this up in sparse cases. The code for
# doing that is in this macro, which should work with either PTrees or TTrees.
macro _forcells(tree, var::Symbol, expr::Expr)
    quote
        let q = _PTREE_BITS_T(0x1), bits = $(esc(tree))._bits
            for $(esc(var)) in 1:_PTREE_COUNT
                ((bits & q) > 0) && $(esc(expr))
                q <<= 1
            end
        end
    end
end
macro _forcells_from(tree, var::Symbol, from, expr::Expr)
    quote
        let ff = $(esc(from)),
            q = _PTREE_BITS_T(0x1) << (ff - 1),
            bits = $(esc(tree))._bits
            for $(esc(var)) in ff:_PTREE_COUNT
                ((bits & q) > 0) && $(esc(expr))
                q <<= 1
            end
        end
    end
end
# # Constructors ===============================================================
#   First, constructors for making empty Trees.
PTree{T}() where {T} = PTree{T}(0x0, 0x0, 0x0, _ptree_empty_cells(T))
TTree{T}() where {T} = TTree{T}(0x0, 0x0, 0x0, _ttree_empty_cells(T))
PTree() = PTree{Any}()
TTree() = TTree{Any}()
#   Now for copying trees of the same type.
PTree{T}(u::PTree{T}) where {T} = u
PTree{T}(u::PTree{S}) where {T,S} = PTree{T}(u._n, u._id, u._bits, _PTREE_CELLS_T{T}(u._data))
TTree{T}(u::TTree{S}) where {T,S} = begin
    newdata = _TTREE_CELL_T{T}[x isa TTree{S} ? TTree{T}(x) : x for x in u._data]
    return TTree{T}(u._n, u._id, u._bits, newdata)
end
#   Now for converting back and forth between PType and TType.
PTree{T}(u::TTree{S}) where {T,S} = begin
    newdata = _PTREE_CELL_T{T}[(x isa PTree{S} ? PTree{T}(x) :
                                x isa TTree{S} ? PTree{T}(x) : x)
                               for x in u._data]
    return PTree{T}(u._n, u._id, u._bits, _PTREE_CELLS_T{T}(newdata))
end
TTree{T}(u::PTree{S}) where {T,S} = begin
    return TTree{T}(u._n, u._id, u._bits,
                    _TTREE_CELL_T{T}[u._data...])
end
#   Now for constructors from abstract dictionaries.
PTree{T}(d::AbstractDict{_PTREE_KEY_T, S}) where {T,S} = let u = TTree{T}()
    for (k,v) in d
        u[k] = v
    end
    return PTree{T}(u)
end
TTree{T}(d::AbstractDict{_PTREE_KEY_T, S}) where {T,S} = let u = TTree{T}()
    for (k,v) in d
        u[k] = v
    end
    return u
end
#   We also construct from arrays; they are considered to be 0-indexed when used
#   with trees.
PTree{T}(a::AbstractArray{1,S}) where {T,S} = begin
    depth = _PTREE_TWIG_DEPTH
    height = 1
    n = length(a)
    if n == 0
        return PTree{T}()
    elseif n == 1
        ch = _ptree_cell_assoc(_ptree_empty_cells(T), 1, a[1])
        return PTree{T}(0x1, _nodeid(0x0, depth), 0x1, ch)
    else
        while n > 1
            # First, divide the previous elements into grouped tree nodes.
            u = a
            a = PTree{T}[]
            maxidx = div(n + _PTREE_COUNT - 1, _PTREE_COUNT)
            hnp = _PTREE_COUNT^(h-1) # number per node at the prev level
            hn  = _PTREE_COUNT*hnp
            bits = ~_PTREE_BITS_T(0x0)
            for ii in 1:maxidx - 1
                nid = _nodeid(hnp * (ii-1), depth)
                ch = _PTREE_CELLS_T{T}(u[1 + (ii-1)*_PTREE_COUNT:ii*_PTREE_COUNT])
                push!(a, PTree{T}(hn, nid, bits, ch))
            end
            # The last element we handle specially because we need to be
            # careful about its n.
            ii = maxidx
            nid = _nodeid(hnp * (ii-1), depth)
            ch = u[1 + (ii-1)*_PTREE_COUNT:end]
            nch = length(ch)
            ch = _PTREE_CELLS_T{T}([ch..., _ptree_noughts[1 + nch:end]...])
            bids = (_PTREE_BITS_T(0x1) << nch) - 0x1
            t = PTree{T}((nch-1)*hnp + ch[end-1]._n, nid, bits, ch)
            push!(a, t)
            # Increment
            depth -= 1
            height += 1
            n = length(a)
        end
        # we're down to one node, so we can return it.
        return a[1]
    end
end
TTree{T}(a::AbstractArray{1,S}) where {T,S} = begin
    depth = _PTREE_TWIG_DEPTH
    height = 1
    n = length(a)
    if n == 0
        return TTree{T}()
    elseif n == 1
        ch = _ttree_cell_assoc(_ttree_empty_cells(T), 1, a[1])
        return TTree{T}(0x1, _nodeid(0x0, depth), 0x1, ch)
    else
        while n > 1
            # First, divide the previous elements into grouped tree nodes.
            u = a
            a = TTree{T}[]
            maxidx = div(n + _PTREE_COUNT - 1, _PTREE_COUNT)
            hnp = _PTREE_COUNT^(h-1) # number per node at the prev level
            hn  = _PTREE_COUNT*hnp
            bits = ~_PTREE_BITS_T(0x0)
            for ii in 1:maxidx - 1
                nid = _nodeid(hnp * (ii-1), depth)
                ch = _TTREE_CELLS_T{T}(u[1 + (ii-1)*_PTREE_COUNT:ii*_PTREE_COUNT])
                push!(a, TTree{T}(hn, nid, bits, ch))
            end
            # The last element we handle specially because we need to be
            # careful about its n.
            ii = maxidx
            nid = _nodeid(hnp * (ii-1), depth)
            ch = u[1 + (ii-1)*_PTREE_COUNT:end]
            nch = length(ch)
            ch = _TTREE_CELLS_T{T}([ch..., _ptree_noughts[1 + nch:end]...])
            bids = (_PTREE_BITS_T(0x1) << nch) - 0x1
            t = TTree{T}((nch-1)*hnp + ch[end-1]._n, nid, bits, ch)
            push!(a, t)
            # Increment
            depth -= 1
            height += 1
            n = length(a)
        end
        # we're down to one node, so we can return it.
        return a[1]
    end
end
# Finally, we want to duplicate all of these constructors using the convert
# function. It's not 100% clear to me why there is/should be both constructors
# and convert functions, or whether there is a good way to avoid reproducing
# all of the code required for that, but for now I'm just duplicating it.
# #TODO: duplicate above code for convert() once it's been tested.


# ==============================================================================
# Methods
#
# Definitions of methods of the above type.
Base.empty(u::PTree{T}) where {T} = PTree{T}()
Base.empty(u::PTree{T}, S::Type) where {T} = PTree{S}()
Base.isempty(u::PTree) = (u._n == 0x0)
Base.length(u::PTree) = Int(u._n)
Base.isequal(t::PTree{T}, s::PTree{S}) where {T,S} = begin
    (t._n == s._n) || return false
    (t._bits == s._bits) || return false
    (t._nodeid == s._nodeid) || return false
    @_forcells t k ((t._data[k] == s._data[k]) || return false)
    return true
end
isequiv(t::PTree{T}, s::PTree{S}) where {T,S} = begin
    (t._n == s._n) || return false
    (t._bits == s._bits) || return false
    (t._nodeid == s._nodeid) || return false
    @_forcells t k isequiv(t._data[k], s._data[k]) || return false
    return true
end
Base.hash(t::PTree{T}) where {T} = let h = hash(PTree)
    for (k,v) in t
        h += hash(v) * (k + 0x1f)
    end
    return h
end
equivhash(t::PTree{T}) where {T} = let h = hash(PTree)
    for (k,v) in t
        h += equivhash(v) * (k + 0x1f)
    end
    return h
end
Base.convert(::Type{PTree{T}}, t::PTree{S}) where {T,S} = PTree{T}(t)
Base.copy(t::PTree{T}) where {T} = t
Base.get(u::PTree{T}, k::_PTREE_KEY_T, df) where {T} = begin
    # Is it in this Tree?
    idx = _node_index(u._id, k)
    (idx == 0) && return df
    ((u._bits & (_PTREE_BITS_T(0x1) << (idx-1))) > 0) || return df
    # See what we have at that position.
    uu = u._data[idx]
    (uu === _nought) && return df
    # If we're a twig-node, this must be the item to return.
    (_node_depth(u._id) == _PTREE_TWIG_DEPTH) && return uu
    # Otherwise, this must be a tree.
    @assert isa(uu, PTree{T}) "Invalid PTree child at depth $(d)"
    return get(uu, k, df)
end
Base.in(kv::Pair{_PTREE_KEY_T,T}, m::PTree{T}) where {T} = begin
    v = get(m, kv[1], _nought)
    return v !== _nought && v === kv[2]
end
Base.iterate(t::PTree{T}) where {T} = begin
    if _node_depth(t._id) == _PTREE_TWIG_DEPTH
        mnk = _node_min(t._id)
        @_forcells t k begin
            return (_PTREE_PAIR_T{T}(_PTREE_KEY_T(k-1) + mnk, t._data[k]), k)
        end
    else
        @_forcells t k begin
            val = iterate(t._data[k])
            (val === nothing) || return (val[1], (k, val[2]))
        end
    end
    return nothing
end
Base.iterate(t::PTree{T}, ks) where {T} = begin
    if _node_depth(t._id) == _PTREE_TWIG_DEPTH
        mnk = _node_min(t._id)
        start = ks + 1
        @_forcells_from t k start begin
            return (_PTREE_PAIR_T{T}(_PTREE_KEY_T(k-1) + mnk, t._data[k]), k)
        end
    else
        (k0, it) = ks
        val = iterate(t._data[k0], it)
        (val === nothing) || return (val[1], (k0, val[2]))
        k0 += 1
        @_forcells_from t k k0 begin
            val = iterate(t._data[k])
            (val === nothing) || return (val[1], (k, val[2]))
        end
    end
    return nothing
end
_maketwig(k::_PTREE_KEY_T, v::T) where {T} = begin
    mnk = k - rem(k, _PTREE_COUNT)
    kk = 1 + (k - mnk)
    bits = _PTREE_BITS_T(0x1) << (kk - 1)
    dat = _ptree_empty_cells(T)
    dat = _ptree_cell_assoc(dat, Int(kk), v)
    return PTree{T}(0x1, _nodeid(mnk, _PTREE_TWIG_DEPTH), bits, dat)
end
setindex(t::PTree{T}, u::S, k::_PTREE_KEY_T) where {T,S} = begin
    idx = _node_index(t._id, k)
    idxbit = _PTREE_BITS_T(1) << (idx - 1)
    if idx == 0
        # The node lies outside of this node, so we need to make a new tree.
        # We need to make the node that would hold both this node (t) and
        # the new twig node. We need to find the nearest common ancestor, which
        # we can find from the bits that are equal.
        nid = t._id
        while idx == 0
            # Go up one parent level
            nid = _node_parent(nid)
            idx = _node_index(nid, k)
        end
        # idx is no longer 0 so we can create this node and insert the new twig.
        twig = _maketwig(k, u::T)
        one = _PTREE_BITS_T(0x1)
        tidx = _node_index(nid, _node_min(t._id))
        bits = (one << (idx - 1)) | (one << (tidx - 1))
        ch = _ptree_empty_cells(T)
        ch = _ptree_cell_assoc(ch, idx, twig)
        ch = _ptree_cell_assoc(ch, tidx, t)
        # We return this parent node, which now contains us plus the new twig.
        return PTree{T}(t._n + 1, nid, bits, ch)
    elseif t._bits & idxbit > 0
        # We might be a twig node.
        if _node_depth(t._id) == _PTREE_TWIG_DEPTH
            ch = _ptree_cell_assoc(t._data, idx, u)
            return PTree{T}(t._n, t._id, t._bits, ch)
        else
            # The node lies below us in an existing node, so we can just pass along
            # the responsibility of inserting it.
            ch0 = t._data[idx]
            ch = setindex(ch0, u, k)
            (ch === ch0) && return t
            (t._n == 0) && return ch
            ch = _ptree_cell_assoc(t._data, idx, ch)
            return PTree{T}(t._n+1, t._id, t._bits, ch)
        end
    else
        # We might be a twig node.
        if _node_depth(t._id) == _PTREE_TWIG_DEPTH
            ch = _ptree_cell_assoc(t._data, idx, u)
            return PTree{T}(t._n+1, t._id, t._bits | idxbit, ch)
        else
            # The node lies beneath us but the lower node doesn't exist; we can just
            # make a twig node for this and put it there.
            twig = _maketwig(k, u::T)
            (t._n == 0) && return twig
            # Now us with that modification
            dat = _ptree_cell_assoc(t._data, idx, twig)
            return PTree{T}(t._n+1, t._id, t._bits | idxbit, dat)
        end
    end
end
delete(t::PTree{T}, k::_PTREE_KEY_T) where {T} = begin
    idx = _node_index(t._id, k)
    idxbit = _PTREE_BITS_T(1) << (idx - 1)
    if idx == 0
        # The node lies outside of this node, so nothing to do.
        return t
    elseif t._bits & idxbit > 0
        (length(t._data) < idx) && return t
        v = t._data[idx]
        (v === _nought) && return t
        # We might be a twig node.
        if _node_depth(t._id) == _PTREE_TWIG_DEPTH
            ch = _ptree_cell_assoc(t._data, idx, _nought)
            return PTree{T}(t._n - 1, t._id, t._bits & ~idxbit, ch)
        else
            # The node lies below us in an existing node, so we can just pass along
            # the responsibility of inserting it.
            ch0 = t._data[idx]
            ch = delete(ch0, k)
            (ch === ch0) && return t
            if ch._n == 0
                ch = _ptree_cell_assoc(t._data, idx, _nought)
                bits = t._bits & ~idxbit
            else
                ch = _ptree_cell_assoc(t._data, idx, ch)
                bits = t._bits
            end
            return PTree{T}(t._n - 1, t._id, bits, ch)
        end
    else
        # We're dropping something that isn't there
        return t
    end
end
