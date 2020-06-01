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
# The types for the leaves and cells of PTree and TTree. Because we're
# declaring these before the struct, we force the stuct to be a type parameter.
# This is just hackerty in order to be able to delcare the types of these only
# once; we declare aliases for them below.
const _PTREE_CELLS_TYPE{PT, T} = Union{Vector{PT}, Vector{T}} where {PT,T}
const _TTREE_CELLS_TYPE{PT, TT, T} = Union{Vector{Union{PT,TT}},
                                           Vector{T}} where {PT,TT,T}
# #PTree =======================================================================
"""
    PTree{T}

PTree is a persistent tree/array hybrid that maps arbitrary unsigned integers
(type _PTREE_KEY_T) to values. It can efficiently be used for arrays or for
hash-maps, and supports efficient lookup, association, and dissociation.

`PTree{T}()` yields an empty tree.
"""
struct PTree{T} <: AbstractDict{_PTREE_KEY_T, T}
    _id::_PTREE_KEY_T
    _bits::_PTREE_BITS_T
    _data::_PTREE_CELLS_TYPE{PTree{T}, T}
end
mutability(::Type{PTree}) = Immutable
mutability(::Type{PTree{T}}) where {T} = Immutable
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
const _PTREE_CELLS_T{T} = _PTREE_CELLS_TYPE{PTree{T},T} where {T}
const _TTREE_CELLS_T{T} = _TTREE_CELLS_TYPE{PTree{T},TTree{T},T} where {T}
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
PTree{T}() where {T} = PTree{T}(0x0, 0x0, PTree{T}[])
TTree{T}() where {T} = TTree{T}(0x0, 0x0, 0x0, PTree{T}[])
#   Now for copying trees of the same type.
PTree{T}(u::PTree{T}) where {T} = u
PTree{T}(u::PTree{S}) where {T,S} = PTree{T}(u._id, u._bits, _PTREE_CELLS_T{T}(u._data))
#   Now for converting back and forth between PType and TType.
PTree{T}(u::TTree{S}) where {T,S} = begin
    newdata = _PTREE_CELL_T{T}[(x isa PTree{S} ? PTree{T}(x) :
                                x isa TTree{S} ? PTree{T}(x) : x)
                               for x in u._data]
    return PTree{T}(u._id, u._bits, _PTREE_CELLS_T{T}(newdata))
end
#   Now for constructors from abstract dictionaries.
PTree{T}(d::AbstractDict{_PTREE_KEY_T, S}) where {T,S} = let u = TTree{T}()
    for (k,v) in d
        u[k] = v
    end
    return PTree{T}(u)
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
        return PTree{T}(_nodeid(0x0, depth), 0x1, T[a[1]])
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
                push!(a, PTree{T}(nid, bits, ch))
            end
            # The last element we handle specially because we need to be
            # careful about its n.
            ii = maxidx
            nid = _nodeid(hnp * (ii-1), depth)
            ch = u[1 + (ii-1)*_PTREE_COUNT:end]
            nch = length(ch)
            ch = _PTREE_CELLS_T{T}(T[ch...])
            bids = (_PTREE_BITS_T(0x1) << nch) - 0x1
            t = PTree{T}(nid, bits, ch)
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
Base.isempty(u::PTree) = (u._bits == 0x0)
_lencount(u::PTree{T}, v::Vector{T}) where {T} = begin
    n = 0
    @_forcells u k n += 1
    return n
end
_lencount(u::PTree{T}, v::Vector{PTree{T}}) where {T} = begin
    n = 0
    @_forcells u k n += length(v[k])
    return n
end
Base.length(u::PTree) = _lencount(u, u._data)
Base.isequal(t::PTree{T}, s::PTree{S}) where {T,S} = begin
    (t._bits == s._bits) || return false
    (t._nodeid == s._nodeid) || return false
    @_forcells t k isequal(t._data[k], s._data[k]) || return false
    return true
end
isequiv(t::PTree{T}, s::PTree{S}) where {T,S} = begin
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
    one = _PTREE_BITS_T(0x1)
    while true
        # Is it in this Tree?
        idx = _node_index(u._id, k)
        (idx == 0) && return df
        ((u._bits & (one << (idx-1))) > 0) || return df
        # Are we a twig?
        if _node_depth(u._id) == _PTREE_TWIG_DEPTH
        #if isa(u._data, Vector{T})
            q = u._data::Vector{T}
            return u._data[idx]
        else
            q = u._data::Vector{PTree{T}}
            u = q[idx]
        end
    end
end
Base.in(kv::Pair{_PTREE_KEY_T,T}, u::PTree{T}) where {T} = begin
    k = kv[1]
    while true
        # Is it in this Tree?
        idx = _node_index(u._id, k)
        (idx == 0) && return false
        ((u._bits & (_PTREE_BITS_T(0x1) << (idx-1))) > 0) || return false
        # Are we a twig?
        #if _node_depth(u._id) == _PTREE_TWIG_DEPTH
        if isa(u._data, Vector{PTree{T}})
            u = u._data[idx]::PTree{T}
        else
            return u._data[idx] === kv[2]
        end
    end
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
_maketwig(::Type{T}, k::_PTREE_KEY_T, v::T) where {T} = begin
    mnk = k - rem(k, _PTREE_COUNT)
    kk = k - mnk
    bits = _PTREE_BITS_T(0x1) << kk
    kk = Int(kk + 1)
    dat = Vector{T}(undef, kk)
    dat[kk] = v
    return PTree{T}(_nodeid(mnk, _PTREE_TWIG_DEPTH), bits, dat)
end
_cells_assoc(t::PTree, u::Vector{T}, k::Int, x) where {T} = begin
    dat = Vector{T}(undef, max(length(u), k))
    @_forcells t ii dat[ii] = u[ii]
    dat[k] = x
    return dat
end
function _ptree_setindex(data::Vector{T}, idx::Int, idxbit::_PTREE_BITS_T,
                         t::PTree{T}, u::S, k::_PTREE_KEY_T) where {T,S}
    bits = t._bits & idxbit > 0 ? t._bits : t._bits | idxbit
    ch = _cells_assoc(t, data, idx, u)
    return PTree{T}(t._id, bits, ch)
end
function _ptree_setindex(data::Vector{PTree{T}}, idx::Int, idxbit::_PTREE_BITS_T,
                         t::PTree{T}, u::S, k::_PTREE_KEY_T) where {T,S}
    if t._bits & idxbit > 0
        # The node lies below us in an existing node, so we can just pass along
        # the responsibility of inserting it.
        ch0 = data[idx]
        ch = setindex(ch0, u, k)
        (ch === ch0) && return t
        (t._bits == 0) && return ch
        ch = _cells_assoc(t, data, idx, ch)
        return PTree{T}(t._id, t._bits, ch)
    else
        # The node lies beneath us but the lower node doesn't exist; we can just
        # make a twig node for this and put it there.
        twig = _maketwig(T, k, u)
        (t._bits == 0) && return twig
        # Now update with that modification
        dat = _cells_assoc(t, t._data, idx, twig)
        return PTree{T}(t._id, t._bits | idxbit, dat)
    end
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
        twig = _maketwig(T, k, u)
        one = _PTREE_BITS_T(0x1)
        tidx = _node_index(nid, _node_min(t._id))
        bits = (one << (idx - 1)) | (one << (tidx - 1))
        ch = Vector{PTree{T}}(undef, max(idx, tidx))
        ch[idx] = twig
        ch[tidx] = t
        # We return this parent node, which now contains us plus the new twig.
        return PTree{T}(nid, bits, ch)
    else
        return _ptree_setindex(t._data, idx, idxbit, t, u, k)
    end
end
delete(t::PTree{T}, k::_PTREE_KEY_T) where {T} = begin
    idx = _node_index(t._id, k)
    idxbit = _PTREE_BITS_T(1) << (idx - 1)
    if idx == 0
        # The node lies outside of this node, so nothing to do.
        return t
    elseif t._bits & idxbit == 0
        # Nothing there
        return t
    elseif isa(t._data, Vector{T})
        # We'e a twig node, so we just clear the bit
        return PTree{T}(t._id, t._bids & ~idxbit, t._data)
    else
        # We recurse down...
        ch0 = t._data[idx]
        ch = delete(ch0, k)
        (ch === ch0) && return t
        if ch._bits == 0
            # Just clear the bit
            ch = t._data
            bits = t._bits & ~idxbit
        else
            ch = _cells_assoc(t, t._data, idx, ch)
            bits = t._bits
        end
        return PTree{T}(t._id, bits, ch)
    end
end
