################################################################################
# Ptrie.jl
#
# **NOTE**
# PTrie is not used in or really part of Air; however, I wrote it as an
# experiment into whether or not giving the compiler more information about
# the structure of the hash-mapped trie while simultaneously explicitly
# representing all levels of the tree, would lead to better performance due to
# substantially more type-safety. The result is that it does not: the number of
# allocations dominates the time constraints instead, and you still have to do
# more array-lookups, which ends up mattering more for time than the type-safety
# issue.
#
# Code implementing a persistent trie type on which the persistent collection
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
"""
    DEFAULT_TRIEROOT_SHIFT

DEFAULT_TRIEROOT_SHIFT (4) is the bit-shift size of the root node of the
PTrie and TTrie types.
"""
const DEFAULT_TRIEROOT_SHIFT = 4
"""
    DEFAULT_TRIENODE_SHIFT

DEFAULT_TRIENODE_SHIFT (5) is the bit-shift size of the non-terminal nodes of 
PTrie and TTrie types.
"""
const DEFAULT_TRIENODE_SHIFT = 5
"""
    DEFAULT_TRIELEAF_SHIFT

DEFAULT_TRIELEAF_SHIFT (5) is the default shift used by leaves if the
type-parameter is not provided.
"""
const DEFAULT_LEAF_SHIFT = 5
"""
    TRIEKEY_T

TRIEKEY_T is the type of a PTrie or TTrie key.
"""
const TRIEKEY_T = typejoin(typeof(hash(nothing)),
                           typeof(equivhash(nothing)),
                           typeof(objectid(nothing)))
"""
    TRIEKEY_MIN

TRIEKEY_MIN is the minimum value of a TRIEKEY_T.
"""
const TRIEKEY_MIN = TRIEKEY_T(0x00)
"""
    TRIEKEY_MAX

TRIEKEY_MAX is the maximum value of a TRIEKEY_T.
"""
const TRIEKEY_MAX = TRIEKEY_T(0x00) - TRIEKEY_T(0x01)
"""
    TRIEKEY_BITCOUNT

TRIEKEY_BITCOUNT is the number of bits in a TRIKEY_T, ergo, the number of bits
in the total addressable hash-space of a PTrie or TTrie.
"""
const TRIEKEY_BITCOUNT = sizeof(TRIEKEY_T) * 8
"""
    TRIEPAIR_T

The pair type for the dict interface of the PTrie and TTrie types.
"""
const TRIEPAIR_T{T} = Pair{TRIEKEY_T, T} where {T}


# # Low-level Utilities ========================================================
#   The functions below are designed to be exploited in low-level code optimized
#   for the compiler. They deal largely with the arithmetic of the Trie nodes.
"""
    trie_count(shiftt)

trie_count(shift) is the number of elements in a single node of a PTrie or TTrie
that uses the given bit-shift.
"""
trie_count(::Val{shift}) where {shift} = (1 << shift)
"""
    triebits_t(shift)

triebits_t(shift) yields the type of the entry bits object of a PTrie or TTrie
node that uses the given shift.
"""
triebits_t(s::Val{shift}) where {shift} = begin
    n = trie_count(s)
    return (n <= 8  ? UInt8  :
            n <= 16 ? UInt16 :
            n <= 32 ? UInt32 :
            n <= 64 ? UInt64 : UInt128)
end
"""
    trieleaf_depth(shift)

trileaf_depth(shift) yields the depth of a leaf (terminal node) in a trie that
uses the given bit shift.
"""
trieleaf_depth(u::Val{shift}) where {shift} =
    TRIEKEY_T(div(TRIEKEY_BITCOUNT + shift - 1, shift))
"""
    trietwig_depth(shift)

trietwig_depth(shift) yields the depth of a twig-node (parent of terminal node)
in a trie that uses the given bit-shift.
"""
trietwig_depth(s::Val{shift}) where {shift} =
    TRIEKEY_T(trieleaf_depth(s) - 1)
"""
    trie_depth_mask(shift)

trie_depth_mask(shift) yields the mask of the bits that represent the depth of
the node-id in a trie that uses the given shift.
"""
trie_depth_mask(::Val{shift}) where {shift} =
    ((TRIEKEY_T(0x1) << (shift)) - 1)
"""
    shift_to_bitstype(shift)

Yields the UInt type large enough to hold the given shift, or Nothing if there is
no such type.
"""
shift_to_bitstype(s::Val{shift}) where {shift} = begin
    count = trie_count(s)
    return (count == 1   ? Bool    :
            count <= 8   ? UInt8   :
            count <= 16  ? UInt16  :
            count <= 32  ? Uint32  :
            count <= 64  ? UInt64  :
            count <= 128 ? UInt128 : Nothing)
end
"""
    shift_to_mask(shift)
    shift_to_mask(start_bit, shift)

Yields the mask of the given shift parameters. The default start_bit is 0.
"""
shift_to_mask(::Val{start_bit}, ::Val{shift}) where {start_bit, shift} =
    ((TRIEKEY_T(0x1) << shift) - 0x1) << start_bit
shift_to_mask(s::Val{shift}) where {shift} = shift_to_mask(Val{0}(), s)
"""
    key_address(S0, S, k)

Calculates the address of the given key k for the shift of S bits starting
at bit S0. Yields a TRIEKEY_T (unsigned integer).
"""
key_address(::Val{S0}, ::Val{S}, k::TRIEKEY_T) where {S0, S} = begin
     return k & ~((TRIEKEY_T(0x1) << (S0 + S)) - 0x1)
end
key_address(s::Val{S}, k::TRIEKEY_T) where {S} = key_address(Val{0}(), s, k)
"""
    key_shiftindex(S0, S, k)
    key_shiftindex(S, k)

Calculates the shift-index of the given key k for the shift of S bits starting
at bit S0. Yields a UInt 0-based index into a the bit-flags for a trie node. 
The default value of S0 is 0.
"""
key_shiftindex(s0::Val{S0}, s::Val{S}, k::TRIEKEY_T) where {S0, S} = begin
    mask = shift_to_mask(s0, s)
    return (k & mask) >> S0
end
key_shiftindex(s::Val{S}, k::TRIEKEY_T) where {S} =
    key_shiftindex(Val{0}(), s, k)


# # The Low-level Trie Types ===================================================
#    Definitions of the structures of the Trie low-level machinery and types
#    that should come before the various methods and constructors that operate
#    on them.
"""
    AbstractTrie{T}

An AbstractTrie represents a hash of unsigned integer (e.g., hash-code or
objectid) keys mapped to values of type T. The hash is implemented using a
bit-wise trie.
"""
abstract type AbstractTrie{T} <: AbstractDict{TRIEKEY_T, T} end


# ==============================================================================
# TrieTwig
#
# The TrieTwig type is a Trie object that stores only the final bits of the
# trie-mapped hash.
"""
    TrieTwig{T, S, B}

A persistent trie leaf-node that uses the first S bits. The parameter B
must be the bitstype appropriate for the given shift (see the function
shift_to_bitstype).

Type Parameters:
 * T: The type of the values stored in the leaves.
 * S: The shift, which must be an integer such that 1 <= S < 8. The shift
      determines the number of children of each node (2^shift).
 * B: The type of the bits value. This must be an unsigned integer type
      that has at least 2^shift bits (which is why the max S is 7: UInt128).
It is up to whoever constructs the PTrieLeaf object to ensure that B has
enough bits for the shift S.

Note that TrieTwig is in fact a mutable struct. This is because all Trie nodes
are intended as either transient nodes or persistent nodes based on their
persistent flag. This means that TrieTwigs are only safe as long as their
fields are not modified outside of the file in which they are defined.
"""
mutable struct TrieTwig{T,S,B} <: AbstractTrie{T}
    _bits::B
    _address::TRIEKEY_T
    _elements::Vector{T}
end
# Note: The last several bits of address are never used! (the last B bits
# to be precise). Accordingly, we can assign them flags:
const TRIE_TRANSIENT_FLAG = TRIEKEY_T(0x1)
_trieaddr(addr::TRIEKEY_T) = addr & ~TRIE_TRANSIENT_FLAG
_triepers(addr::TRIEKEY_T) = addr & TRIE_TRANSIENT_FLAG == 0x0
_trietrns(addr::TRIEKEY_T) = addr & TRIE_TRANSIENT_FLAG == 0x1
_trietopers(addr::TRIEKEY_T) = addr & ~TRIE_TRANSIENT_FLAG
_trietotrns(addr::TRIEKEY_T) = addr | TRIE_TRANSIENT_FLAG
"""
    trieaddress(trie)

Yields the address of the given trie object. The address is the key of the
first value the given trie node could possibly store.
"""
trieaddress(addr::TRIEKEY_T) where {S} = _trieaddr(addr)
trieaddress(u::TrieTwig{T,S,B}) where {T,S,B} = _trieaddr(u._address)
"""
    ispersistent(trie)

Yields true if the given trie node or twig is persistent and false otherwise.
"""
ispersistent(u::TrieTwig{T,S,B}) where {T,S,B} =
    _triepers(u._address)
"""
    istransient(trie)

Yields true if the given trie node or twig is transient and false otherwise.
"""
istransient(u::TrieTwig{T,S,B}) where {T,S,B} =
    _trietrns(u._address)
"""
    addresscompare(trie, key)

Yields -1, 0, or 1 if the given key is before, within, or after the address
space of the given trie node. In particular, a return value of 0 indicates that
the given key is in the address space of trie.
"""
addresscompare(trie::TrieTwig{T,S,B}, key::TRIEKEY_T) where {T,S,B} = begin
    mask = ~shift_to_mask(Val{S}())
    taddr = trieaddress(trie) & mask
    kaddr = key & mask
    return (kaddr == taddr ? 0 : kaddr < taddr ? -1 : 1)
end
Base.IteratorSize(u::TrieTwig) = Base.HasLength()
Base.IteratorEltype(u::TrieTwig) = Base.HasEltype()
Base.length(u::TrieTwig{T,S,B}) where {T,S,B} = count_ones(u._bits)
Base.eltype(u::TrieTwig{T,S,B}) where {T,S,B} = T
Base.eltype(::Type{TrieTwig{T,S,B}}) where {T,S,B} = T
Base.show(io::IO, ::MIME"text/plain", u::TrieTwig{T,S,B}) where {T,S,B} =
    print(io, "TrieTwig{$T,$S,$B}($(u._address), $(u._bits))")
Base.show(io::IO, u::TrieTwig{T,S,B}) where {T,S,B} =
    print(io, "TrieTwig{$T,$S,$B}($(u._address), $(u._bits))")
_private_iterate(u::TrieTwig{T,S,B}, k0::TRIEKEY_T) where {T,S,B} = begin
    mask = shift_to_mask(Val{S}())
    idx = k0 & mask
    bitmaskx = ~((B(0x1) << idx) - 0x1)
    b = u._bits & bitmaskx
    nextidx = trailing_zeros(b)
    (nextidx >= (B(0x1) << S)) && return nothing
    k = (k0 & ~mask) | TRIEKEY_T(nextidx)
    return (k => u._elements[nextidx + 1], k + 0x1)
end
Base.iterate(u::TrieTwig{T,S,B}, k0::TRIEKEY_T) where {T,S,B} = begin
    cmp = addresscompare(u, k0)
    if cmp < 0
        k0 = trieaddress(u)
    elseif cmp > 0
        return nothing
    end
    return _private_iterate(u, k0)
end
Base.iterate(u::TrieTwig{T,S,B}) where {T,S,B} =
    _private_iterate(u, trieaddress(u))
_private_haskey(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = B(0x1) << idx
    return u._bits & flag == flag
end
Base.haskey(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    (addresscompare(u, k) == 0) || return false
    return _private_haskey(u, k)
end
_private_get(u::TrieTwig{T,S,B}, k::TRIEKEY_T, failval) where {T,S,B} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = B(0x1) << idx
    (u._bits & flag == flag) || return failval
    return u._elements[idx + 1]
end
Base.get(u::TrieTwig{T,S,B}, k::TRIEKEY_T, failval) where {T,S,B} = begin
    (addresscompare(u, k) == 0) || return failval
    return _private_get(u, k, failval)
end
_private_getindex(u::TrieTwig{T,S,B}, k::TRIEKEY_T, u0) where {T,S,B} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = B(0x1) << idx
    (u._bits & flag == flag) || throw(IndexError(u0, k))
    return u._elements[idx + 1]
end
Base.getindex(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    (addresscompare(u, k) == 0) || throw(IndexError(u, k))
    return _private_getindex(u, k, u)
end
Base.firstindex(u::TrieTwig{T,S,B}) where {T,S,B} =
    trailing_zeros(u._bits)
Base.lastindex(u::TrieTwig{T,S,B}) where {T,S,B} =
    (B(0x1) << S) - leading_zeros(u._bits)
Base.isequal(t::TrieTwig{T,S,B}, u::TrieTwig{U,S,B}) where {T,S,B,U} = begin
    # Identical objects are always equal.
    (t === u) && return true
    # Otherwise, equality depends only on the same keys being mapped to equal
    # values.
    (trieaddress(t) == tieaddress(u) && t._bits == u._bits) || return false
    bits = t._bits
    while bits != 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        isequal(t._elements[idx], u._elements[idx]) || return false
        bits &= ~flag
    end
    return true
end
isequiv(t::TrieTwig{T,S,B}, u::TrieTwig{U,S,B}) where {T,S,B,U} = begin
    # Identical objects are always equivalent.
    (t === u) && return true
    # Otherwise, equivalence depends on persistence.
    (istransient(t) || istransient(u)) && return false
    # Equavalence depends on the same keys being mapped to equivalent values.
    (t._address == u._address && t._bits == u._bits) || return false
    bits = t._bits
    while bits != 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        isequiv(t._elements[idx], u._elements[idx]) || return false
        bits &= ~flag
    end
    return true
end

# Constructors
TrieTwig{T,S,B}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {T,S,B} =
    TrieTwig{T,S,B}(B(0x0),
                    persistent ? _trietopers(addr) : _trietotrns(addr),
                    Vector{T}(undef, 2^S))
TrieTwig{T,S}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {T,S} =
    TrieTwig{T,S,shift_to_bitstype(Val{S}())}(addr, persistent)
TrieTwig{T}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {T} =
    TrieTwig{T,DEFAULT_TRIELEAF_SHIFT}(addr, persistent)
TrieTwig(addr::TRIEKEY_T=0x0, persistent::Bool=true) =
    TrieTwig{Any}(addr, persistent)
# We want to be able to construct with a key-valye pair, in order to make it
# easy to add new nodes into a trie.
TrieTwig{T,S,B}(kv::Pair{TRIEKEY_T,V}, persistent::Bool=true) where {T,S,B,V} = begin
    k = kv.first
    # We start by figuring out what address this is at this level.
    addr = key_address(Val{S}(), k)
    idx = key_shiftindex(Val{S}(), k)
    u = TrieTwig{T,S,B}(addr, persistent)
    u._elements[idx + 1] = kv.second
    u._bits |= B(0x1) << idx
    return u
end
# Other duplicators and updaters
Base.copy(u::TrieTwig{T,S,B}) where {T,S,B} = begin
    ispersistent(u) && return u
    return TrieTwig{T,S,B}(u._bits, u._address, copy(u._elements))
end
Base.convert(::Type{TrieTwig{T,S,B}}, u::TrieTwig{U,S,B}) where {T,U,S,B} =
    TrieTwig{T,S,B}(u._bits, u._address, convert(Vector{T}, u._elements))
TrieTwig{T,S,B}(u::TrieTwig{U,S,B}) where {T,U,S,B} =
    TrieTwig{T,S,B}(u._bits, u._address, convert(Vector{T}, u._elements))
TrieTwig{T,S,B}(u::TrieTwig{T,S,B}) where {T,S,B} = begin
    ispersistent(u) && return u
    return TrieTwig{T,S,B}(u._bits, u._address, copy(u._elements))
end
TrieTwig{T}(u::TrieTwig{U,S,B}) where {T,U,S,B} =
    TrieTwig{T,S,B}(u._bits, u._address, convert(Vector{T}, u._elements))
TrieTwig{T}(u::TrieTwig{T,S,B}) where {T,S,B} = TrieTwig{T,S,B}(u)
TrieTwig(u::TrieTwig{T,S,B}) where {T,S,B} = TrieTwig{T,S,B}(u)

"""
    TrieNode{T,S0,S,B}

A node in a trie hash. The type parameter S is the bit-shift while the
parameter B is the (Unsigned) bits type. The parameter T must be
exactly the type of the children. The parameter S0 is the initial bit
that the node tests for. Note that the eltype of a TrieNode is in
fact not T but the the eltype of T.
"""
mutable struct TrieNode{T,S0,S,B} <: AbstractTrie{T}
    _bits::B
    _length::Int
    _address::TRIEKEY_T
    _elements::Vector{T}
end
trieaddress(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    _trieaddr(u._address)
ispersistent(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    _triepers(u._address)
istransient(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    _trietrns(u._address)
addresscompare(trie::TrieNode{T,S0,S,B}, key::TRIEKEY_T) where {T,S0,S,B} = begin
    mask = ~shift_to_mask(Val{S + S0}())
    taddr = trieaddress(trie) & mask
    kaddr = key & mask
    return (kaddr == taddr ? 0 : kaddr < taddr ? -1 : 1)
end
Base.IteratorSize(u::TrieNode) = Base.HasLength()
Base.IteratorEltype(u::TrieNode) = Base.HasEltype()
Base.length(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = u._length
Base.eltype(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = eltype(T)
Base.eltype(::Type{TrieNode{T,S0,S,B}}) where {T,S0,S,B} = eltype(T)
Base.show(io::IO, ::MIME"text/plain", u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    print(io, "TrieNode{$T,$S0,$S,$B}($(u._address), $(u._bits))")
Base.show(io::IO, u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    print(io, "TrieNode{$T,$S0,$S,$B}($(u._address), $(u._bits))")
_private_iterate(u::TrieNode{T,S0,S,B}, k0::TRIEKEY_T) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k0)
    # Start by checking if k0 is a child...
    bitmask = B(0x1) << idx
    if u._bits & bitmask == bitmask
        el = u._elements[idx + 1]
        next = _private_iterate(el, k0)
        (next === nothing) || return next
    end
    # At this point, there was not a value for k0; look for the next
    bitmask_hi = ~((bitmask << 1) - 0x1)
    b = u._bits & bitmask_hi
    nextidx = trailing_zeros(b)
    (nextidx >= (B(0x1) << S)) && return nothing
    # Since this is the fist time coming to nextidx, we know there must
    # be something for it to return:
    kmask_hi = ~((TRIEKEY_T(0x1) << (S0 + S)) - 0x1)
    k = (k0 & kmask_hi) | (TRIEKEY_T(nextidx) << S0)
    return _private_iterate(u._elements[nextidx + 1], k)
end
Base.iterate(u::TrieNode{T,S0,S,B}, k0::TRIEKEY_T) where {T,S0,S,B} = begin
    cmp = addresscompare(u, k0)
    if cmp < 0
        k0 = trieaddress(u)
    elseif cmp > 0
        return nothing
    end
    return _private_iterate(u, k0)
end
Base.iterate(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} =
    _private_iterate(u, trieaddress(u))
_private_haskey(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    return ((u._bits & flag == flag) &&
            _private_haskey(u._elements[idx], k))
end
Base.haskey(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    adr = key_address(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    return ((adr == trieaddress(u)) &&
            (u._bits & flag == flag) &&
            _private_haskey(u._elements[idx], k))
end
_private_get(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T, failval, u0) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    if (u._bits & flag) == flag
        return _private_get(u._elements[idx + 1], k, failval, u0)
    else
        return failval
    end
end
Base.get(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T, failval) where {T,S0,S,B} = begin
    adr = key_address(Val{S0}(), Val{S}(), k)
    (adr == trieaddress(u)) || return failval
    return _private_get(u, k, failval, k)
end
_private_getindex(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T, u0) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    (u._bits & flag == flag) || throw(IndexError(u0, k))
    return _private_getindex(u._elements[idx], k, u0)
end
Base.getindex(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    adr = key_address(Val{S0}(), Val{S}(), k)
    (adr == trieaddress(u)) || throw(IndexError(u, k))
    return _private_getindex(u, k, u)
end
Base.firstindex(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = begin
    idx = trailing_zeros(u._bits) + 1
    if idx > (B(0x1) << S)
        return TRIEKEY_T(0x0) - 0x1
    else
        return firstindex(u._elements[idx])
    end
end
Base.lastindex(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = begin
    maxbit = sizeof(B) * 8
    idx = maxbit - trailing_zeros(u._bits)
    if idx > (B(0x1) << S)
        return TRIEKEY_T(0x0)
    else
        return lastindex(u._elements[idx])
    end
end
Base.isequal(t::TrieNode{T,S0,S,B}, u::TrieNode{U,S0,S,B}) where {T,S0,S,B,U} = begin
    # Identical objects are always equal.
    (t === u) && return true
    # Otherwise, equality depends only on the same keys being mapped to equal
    # values.
    (trieaddress(t) == tieaddress(u)  &&
     t._bits == u._bits               &&
     t._length == u._length         ) || return false
    bits = t._bits
    while bits != 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        isequal(t._elements[idx], u._elements[idx]) || return false
        bits &= ~flag
    end
    return true
end
isequiv(t::TrieNode{T,S0,S,B}, u::TrieNode{U,S0,S,B}) where {T,S0,S,B,U} = begin
    # Identical objects are always equivalent.
    (t === u) && return true
    # Otherwise, equivalence depends on persistence.
    (istransient(t) || istransient(u)) && return false
    # Equavalence depends on the same keys being mapped to equivalent values.
    (t._address == u._address  &&
     t._bits == u._bits        &&
     t._length == u._length  ) || return false
    bits = t._bits
    while bits != 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        isequiv(t._elements[idx], u._elements[idx]) || return false
        bits &= ~flag
    end
    return true
end

# Constructors (empty and from a single pair)
TrieNode{T,S0,S,B}(addr::TRIEKEY_T=TRIEKEY_T(0x0), persistent::Bool=true) where {T,S0,S,B} =
    TrieNode{T,S0,S,B}(B(0x0), 0,
                       persistent ? _trietopers(addr) : _trietotrns(addr),
                       Vector{T}(undef, 2^S))
TrieNode{T,S0,S}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {T,S0,S} =
    TrieNode{T,S0,S,shift_to_bitstype(Val{S}())}(addr, persistent)
TrieNode{T,S0}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {T,S0} =
    TrieNode{T,S0,DEFAULT_TRIENODE_SHIFT}(addr, persistent)
TrieNode{S0}(addr::TRIEKEY_T=0x0, persistent::Bool=true) where {S0} =
    TrieNode{Any,S0,DEFAULT_TRIENODE_SHIFT}(addr, persistent)
TrieNode{T,S0,S,B}(kv::Pair{TRIEKEY_T,V},
                   persistent::Bool=true) where {T,S0,S,B,V} = begin
    k = kv.first
    # We start by figuring out what address this is at this level.
    addr = key_address(Val{S0}(), Val{S}(), k)
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    u = TrieNode{T,S0,S,B}(addr, persistent)
    u._bits |= B(0x1) << idx
    el = T(kv, persistent)
    u._elements[idx + 1] = el
    u._length = 1
    return u
end
TrieNode{T,S0,S,B}(u::TrieNode{T,S0,S,B}, persistent::Bool=true) where {T,S0,S,B} = begin
    els = persistent && ispersistent(u) ? u._elements : copy(u._elements)
    addr = trieaddress(u)
    return TrieNode{T,S0,S,B}(u._bits, u._length,
                              persistent ? _trietopers(addr) : _trietotrns(addr),
                              els)
end
# Other duplicators and updaters
Base.copy(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = begin
    ispersistent(u) && return u
    u = TrieNode{T,S0,S,B}(u._bits, u._length, u._address, copy(u._elements))
    bits = u._bits
    while bits > 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        bits &= ~flag
        uu = u._elements[idx]
        ispersistent(uu) || (u._elements[idx] = copy(uu))
    end
    return u
end
Base.convert(::Type{TrieNode{T,S0,S,B}},
             u::TrieNode{U,S0,S,B}) where {T,U,S0,S,B} = begin
    els = Vector{T}(undef, 2^S)
    bits = u._bits
    while bits > 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        bits &= ~flag
        uu = u._elements[idx]
        els[idx] = convert(T, uu)
    end
    u = TrieNode{T,S0,S,B}(u._bits, u._length, u._address, els)
end
Base.convert(::Type{TrieNode{T,S0,S,B}}, u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = u
TrieNode{T,S0,S,B}(u::TrieNode{U,S0,S,B}) where {T,U,S0,S,B} =
    convert(TrieNode{T,S0,S,B}, u)
TrieNode{T}(u::TrieNode{U,S0,S,B}) where {T,U,S0,S,B} = TrieNode{T,S0,S,B}(u)
TrieNode(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = TrieNode{T,S0,S,B}(u)

# Aliases for TrieNode types with a particular shift.
const TrieNode3{T,S0} = TrieNode{T,S0,3,UInt8} where {T,S0}
const TrieNode4{T,S0} = TrieNode{T,S0,4,UInt16} where {T,S0}
const TrieNode5{T,S0} = TrieNode{T,S0,5,UInt32} where {T,S0}
const TrieNode6{T,S0} = TrieNode{T,S0,6,UInt64} where {T,S0}
const TrieNode7{T,S0} = TrieNode{T,S0,7,UInt128} where {T,S0}
const TrieTwig3{T} = TrieTwig{T,3,UInt8} where {T}
const TrieTwig4{T} = TrieTwig{T,4,UInt16} where {T}
const TrieTwig5{T} = TrieTwig{T,5,UInt32} where {T}
const TrieTwig6{T} = TrieTwig{T,6,UInt64} where {T}
const TrieTwig7{T} = TrieTwig{T,7,UInt128} where {T}

# PTrie objects use a particular preset nesting of node objects;
# the total shifts must sum to TRIEKEY_BITCOUNT
if TRIEKEY_BITCOUNT == 8
    const TrieL0{T} = TrieTwig5{T}             where {T}
    const Trie{T}   = TrieNode3{TrieL0{T}, 5}  where {T}
elseif TRIEKEY_BITCOUNT == 16
    const TrieL0{T} = TrieTwig5{T}             where {T}
    const TrieL1{T} = TrieNode5{TrieL0{T},  5} where {T}
    const TrieL2{T} = TrieNode3{TrieL1{T}, 10} where {T}
    const Trie{T}   = TrieNode3{TrieL2{T}, 13} where {T}
elseif TRIEKEY_BITCOUNT == 32
    const TrieL0{T} = TrieTwig5{T}             where {T}
    const TrieL1{T} = TrieNode5{TrieL0{T},  5} where {T}
    const TrieL2{T} = TrieNode5{TrieL1{T}, 10} where {T}
    const TrieL3{T} = TrieNode5{TrieL2{T}, 15} where {T}
    const TrieL4{T} = TrieNode5{TrieL3{T}, 20} where {T}
    const TrieL5{T} = TrieNode4{TrieL4{T}, 25} where {T}
    const Trie{T}   = TrieNode3{TrieL5{T}, 29} where {T}
elseif TRIEKEY_BITCOUNT == 64
    #const TrieL0{T} = TrieTwig4{T}             where {T}
    #const TrieL1{T} = TrieNode4{TrieL0{T},  4} where {T}
    #const TrieL2{T} = TrieNode4{TrieL1{T},  8} where {T}
    #const TrieL3{T} = TrieNode4{TrieL2{T}, 12} where {T}
    #const TrieL4{T} = TrieNode4{TrieL3{T}, 16} where {T}
    #const TrieL5{T} = TrieNode4{TrieL4{T}, 20} where {T}
    #const TrieL6{T} = TrieNode4{TrieL5{T}, 24} where {T}
    #const TrieL7{T} = TrieNode4{TrieL6{T}, 28} where {T}
    #const TrieL8{T} = TrieNode4{TrieL7{T}, 32} where {T}
    #const TrieL9{T} = TrieNode4{TrieL8{T}, 36} where {T}
    #const TrieLA{T} = TrieNode4{TrieL9{T}, 40} where {T}
    #const TrieLB{T} = TrieNode4{TrieLA{T}, 44} where {T}
    #const TrieLC{T} = TrieNode4{TrieLB{T}, 48} where {T}
    #const TrieLD{T} = TrieNode4{TrieLC{T}, 52} where {T}
    #const TrieLE{T} = TrieNode4{TrieLD{T}, 56} where {T}
    #const Trie{T}   = TrieNode4{TrieLE{T}, 60} where {T}
    const TrieL0{T} = TrieTwig5{T}             where {T}
    const TrieL1{T} = TrieNode5{TrieL0{T},  5} where {T}
    const TrieL2{T} = TrieNode5{TrieL1{T}, 10} where {T}
    const TrieL3{T} = TrieNode5{TrieL2{T}, 15} where {T}
    const TrieL4{T} = TrieNode5{TrieL3{T}, 20} where {T}
    const TrieL5{T} = TrieNode5{TrieL4{T}, 25} where {T}
    const TrieL6{T} = TrieNode5{TrieL5{T}, 30} where {T}
    const TrieL7{T} = TrieNode5{TrieL6{T}, 35} where {T}
    const TrieL8{T} = TrieNode5{TrieL7{T}, 40} where {T}
    const TrieL9{T} = TrieNode5{TrieL8{T}, 45} where {T}
    const TrieLA{T} = TrieNode5{TrieL9{T}, 50} where {T}
    const TrieLB{T} = TrieNode5{TrieLA{T}, 55} where {T}
    const Trie{T}   = TrieNode4{TrieLB{T}, 60} where {T}
    #const TrieL0{T} = TrieTwig6{T}             where {T}
    #const TrieL1{T} = TrieNode6{TrieL0{T},  6} where {T}
    #const TrieL2{T} = TrieNode6{TrieL1{T}, 12} where {T}
    #const TrieL3{T} = TrieNode6{TrieL2{T}, 18} where {T}
    #const TrieL4{T} = TrieNode6{TrieL3{T}, 24} where {T}
    #const TrieL5{T} = TrieNode6{TrieL4{T}, 30} where {T}
    #const TrieL6{T} = TrieNode6{TrieL5{T}, 36} where {T}
    #const TrieL7{T} = TrieNode6{TrieL6{T}, 42} where {T}
    #const TrieL8{T} = TrieNode6{TrieL7{T}, 48} where {T}
    #const TrieL9{T} = TrieNode6{TrieL8{T}, 54} where {T}
    #const Trie{T}   = TrieNode4{TrieL9{T}, 60} where {T}
    #const TrieL0{T} = TrieTwig7{T}             where {T}
    #const TrieL1{T} = TrieNode7{TrieL0{T},  7} where {T}
    #const TrieL2{T} = TrieNode7{TrieL1{T}, 14} where {T}
    #const TrieL3{T} = TrieNode7{TrieL2{T}, 21} where {T}
    #const TrieL4{T} = TrieNode7{TrieL3{T}, 28} where {T}
    #const TrieL5{T} = TrieNode7{TrieL4{T}, 35} where {T}
    #const TrieL6{T} = TrieNode7{TrieL5{T}, 42} where {T}
    #const TrieL7{T} = TrieNode7{TrieL6{T}, 49} where {T}
    #const TrieL8{T} = TrieNode4{TrieL7{T}, 56} where {T}
    #const Trie{T}   = TrieNode4{TrieL6{T}, 60} where {T}
elseif TRIEKEY_BITCOUNT == 128
    const TrieL00{T} = TrieTwig5{T}               where {T}
    const TrieL01{T} = TrieNode5{TrieL00{T},   5} where {T}
    const TrieL02{T} = TrieNode5{TrieL01{T},  10} where {T}
    const TrieL03{T} = TrieNode5{TrieL02{T},  15} where {T}
    const TrieL04{T} = TrieNode5{TrieL03{T},  20} where {T}
    const TrieL05{T} = TrieNode5{TrieL04{T},  25} where {T}
    const TrieL06{T} = TrieNode5{TrieL05{T},  30} where {T}
    const TrieL07{T} = TrieNode5{TrieL06{T},  35} where {T}
    const TrieL08{T} = TrieNode5{TrieL07{T},  40} where {T}
    const TrieL09{T} = TrieNode5{TrieL08{T},  45} where {T}
    const TrieL0A{T} = TrieNode5{TrieL09{T},  50} where {T}
    const TrieL0B{T} = TrieNode5{TrieL0A{T},  55} where {T}
    const TrieL0C{T} = TrieNode5{TrieL0B{T},  60} where {T}
    const TrieL0D{T} = TrieNode5{TrieL0C{T},  65} where {T}
    const TrieL0E{T} = TrieNode5{TrieL0D{T},  70} where {T}
    const TrieL0F{T} = TrieNode5{TrieL0E{T},  75} where {T}
    const TrieL10{T} = TrieNode5{TrieL0F{T},  80} where {T}
    const TrieL11{T} = TrieNode5{TrieL10{T},  85} where {T}
    const TrieL12{T} = TrieNode5{TrieL11{T},  90} where {T}
    const TrieL13{T} = TrieNode5{TrieL12{T},  95} where {T}
    const TrieL14{T} = TrieNode5{TrieL13{T}, 100} where {T}
    const TrieL14{T} = TrieNode5{TrieL14{T}, 105} where {T}
    const TrieL16{T} = TrieNode5{TrieL15{T}, 110} where {T}
    const TrieL17{T} = TrieNode5{TrieL16{T}, 115} where {T}
    const TrieL18{T} = TrieNode4{TrieL17{T}, 120} where {T}
    const Trie{T}    = TrieNode4{TrieL18{T}, 124} where {T}
else
    error("PTrie hash type has $(TRIEKEY_BITCOUNT) bits;" *
          " this size is not supported")
end

# We can process the Trie{T} type to determine what the bit-pieces of the
# TRIEKEY_T type are.
const trie_split = let layers = [],
                       TT = Trie{Symbol}
    while TT <: TrieNode
        (T, S0, S, B) = TT.parameters
        mask = shift_to_mask(Val(S0), Val(S))
        push!(layers, (UInt8(S0), TRIEKEY_T(mask)))
        TT = T
    end
    (T, S, B) = TT.parameters
    mask = shift_to_mask(Val(S))
    push!(layers, (UInt8(0), TRIEKEY_T(mask)))
    Tuple(layers)
end

# Mutators
_private_setindex!(u::TrieTwig{T,S,B}, val::V, k::TRIEKEY_T, ispers::Bool=false) where {T,S,B,V} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = (B(0x1) << idx)
    u._elements[idx] = val
    u._bits |= flag
    return u
end
_private_setindex!(u::TrieNode{T,S0,S,B}, val::V, k::TRIEKEY_T, ispers::Bool=false) where {T,S0,S,B,V} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    idx += 1
    if u._bits & flag == flag
        # There's already a child here!
        # If this is a persistent child, we make a new (transient) child;
        # otherwise we edit it inplace
        els = u._elements
        el = els[idx]
        if !ispers && ispersistent(el)
            el = T(el, ispers)
            els[idx] = el
        end
        l0 = length(el)
        _private_setindex!(el, val, k, ispers)
        l1 = length(el)
        (l0 == l1) || (u._length += l1 - l0)
        return u
    else
        # We make a new child.
        u._bits |= flag
        u._elements[idx] = T(k => val, ispers)
        u._length += 1
    end
    return u
end
_private_setindex(u::TrieTwig{T,S,B}, val::V, k::TRIEKEY_T) where {T,S,B,V} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = (B(0x1) << idx)
    bits = u._bits
    els = copy(u._elements)
    addr = u._address
    els[idx + 1] = val
    return TrieTwig{T,S,B}(bits | flag, addr, els)
end
_private_setindex(u::TrieNode{T,S0,S,B}, val::V, k::TRIEKEY_T) where {T,S0,S,B,V} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    idx += 1
    bits = u._bits
    els = copy(u._elements)
    addr = u._address
    len = u._length
    if bits & flag == flag
        uu = els[idx]
        l0 = length(uu)
        uu = _private_setindex(uu, val, k)
        l1 = length(uu)
        len += l1 - l0
    else
        uu = T(k => val, true)
        len += 1
    end
    els[idx] = uu
    return TrieNode{T,S0,S,B}(bits | flag, len, addr, els)
end
setindex(u::TrieTwig{T,S,B}, val::V, k::TRIEKEY_T) where {T,S,B,V} = begin
    # Start by comparing addresses.
    cmp = addresscompare(u, k)
    # If the key is in this twig's address space, we just add the val to this
    # twig.
    (cmp == 0) && return _private_setindex(u, val, k)
    addr = trieaddress(u)
    # Otherwise, we need to return a node of a higher type. We can just create
    # a valid upper node then pass control to it.
    UpS0 = S
    UpS = min(DEFAULT_NODE_SHIFT, TRIEKEY_BITCOUNT - S)
    UpT = TrieTwig{T,S,B}
    UpB = shift_to_bitstype(Val{S}())
    upels = Vector{TT}(undef, 2^UpS)
    idx = key_shiftindex(Val{UpS0}(), Val{UpS}(), addr)
    upels[idx + 1] = u
    upmask = shift_to_mask(Val{UpS0 + UpS}())
    upaddr = addr & ~upmask
    upbits = UpB(0x1) << idx
    upper = TrieNode{UpT,UpS0,UpS,UpB}(length(u), upbits, upaddr, upels)
    # If upper now has the key in scope, we use the private setindex; otherwise
    # just call on up.
    if addresscompare(upaddr, k) == 0
        return _private_setindex!(upper, val, k, true)
    else
        return setindex(upper, val, k)
    end
end
setindex(u::TrieNode{T,S0,S,B}, val::V, k::TRIEKEY_T) where {T,S0,S,B,V} = begin
    # Start by comparing addresses.
    cmp = addresscompare(u, k)
    # If the key is in this twig's address space, we just add the val to this
    # twig.
    (cmp == 0) && return _private_setindex(u, val, k)
    addr = trieaddress(u)
    # Otherwise, we need to return a node of a higher type. We can just create
    # a valid upper node then pass control to it.
    UpS0 = S + S0
    UpS = min(DEFAULT_NODE_SHIFT, TRIEKEY_BITCOUNT - UpS0)
    UpT = TrieNode{T,S0,S,B}
    UpB = shift_to_bitstype(Val{UpS}())
    upels = Vector{UpT}(undef, 2^UpS)
    idx = key_shiftindex(Val{UpS0}(), Val{UpS}(), addr)
    upels[idx] = u
    upmask = shift_to_mask(Val{UpS0 + UpS}())
    upaddr = addr & ~upmask
    upbits = UpB(0x1) << idx
    upper = TrieNode{UpT,UpS0,UpS,UpB}(length(u), B(0x1) << idx, upaddr, upels)
    if addresscompare(upaddr, k) == 0
        return _private_setindex!(upper, val, k, true)
    else
        return setindex(upper, val, k)
    end
end
setindex!(u::TrieTwig{T,S,B}, val::V, k::TRIEKEY_T) where {T,S,B,V} = begin
    # Check for transience first.
    if ispersistent(u)
        throw(ArgumentError("setindex! illegal for persistent TrieTwig"))
    end
    # Start by comparing addresses.
    cmp = addresscompare(u, k)
    # If the key is in this twig's address space, we just add the val to this
    # twig.
    if cmp == 0
        return _private_setindex!(u, val, k, false)
    else
        throw(ArgumentError("setindex! called outside of TrieTwig's range"))
    end
end
setindex!(u::TrieNode{T,S0,S,B}, val::V, k::TRIEKEY_T) where {T,S0,S,B,V} = begin
    # Check for transience first.
    if ispersistent(u)
        throw(ArgumentError("setindex! illegal for persistent TrieNode"))
    end
    # Start by comparing addresses.
    cmp = addresscompare(u, k)
    # If the key is in this node's address space, we just add the val to this
    # node; otherwise, error.
    if cmp == 0
        return _private_setindex!(u, val, k)
    else
        throw(ArgumentError("setindex! called outside of TrieNode's range"))
    end
end
_private_delete(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = (B(0x1) << idx)
    bits = u._bits
    (bits & flag == 0x0) && return u
    els = u._elements
    addr = u._address
    # no actual need to delete since we deref the flag
    return TrieTwig{T,S,B}(bits & ~flag, addr, els)
end
_private_delete(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    idx += 1
    bits = u._bits
    (bits & flag == 0x0) && return u
    els = copy(u._elements)
    addr = u._address
    el0 = els[idx]
    el1 = _private_delete(el0, k)
    (el1 === el0) && return u
    if length(el1) == 0
        bits &= ~flag
    end        
    # no actual need to delete since we deref the flag
    return TrieTwig{T,S,B}(bits, addr, els)
end
delete(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    # Check the address.
    cmp = addresscompare(u, k)
    return (cmp == 0 ? _private_delete(u, k) : u)
end
delete(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    # Check the address.
    cmp = addresscompare(u, k)
    return (cmp == 0 ? _private_delete(u, k) : u)
end
_private_delete!(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    idx = key_shiftindex(Val{S}(), k)
    flag = (B(0x1) << idx)
    # no actual need to delete since we deref the flag
    u._bits &= ~flag
    return u
end
_private_delete!(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    idx = key_shiftindex(Val{S0}(), Val{S}(), k)
    flag = (B(0x1) << idx)
    idx += 1
    bits = u._bits
    (bits & flag == 0x0) && return u
    els = u._elements
    addr = u._address
    el = els[idx]
    n0 = length(el)
    el = _private_delete!(el, k)
    n1 = length(e1)
    (n1 === n0) && return u
    (n1 == 0) && (u._bits &= ~flag)
    # no actual need to delete since we deref the flag
    return u
end
delete!(u::TrieTwig{T,S,B}, k::TRIEKEY_T) where {T,S,B} = begin
    # If we're persistent we can't delete
    ispersistent(u) && throw(
        ArgumentError("cannot delete! from persistent trie"))
    # Check the address.
    cmp = addresscompare(u, k)
    return (cmp == 0 ? _private_delete!(u, k) : u)
end
delete!(u::TrieNode{T,S0,S,B}, k::TRIEKEY_T) where {T,S0,S,B} = begin
    # If we're persistent we can't delete
    ispersistent(u) && throw(
        ArgumentError("cannot delete! from persistent trie"))
    # Check the address.
    cmp = addresscompare(u, k)
    return (cmp == 0 ? _private_delete!(u, k) : u)
end
# Persisting a Trie sets all the transient bits to false. This should always be
# permanent.
"""
    persist!(t)

Yields t after setting all the transient bits in the given Trie t to be false.
If t is already persistent, yields t without changing anything.
"""
persist!(u::TrieTwig{T,S,B}) where {T,S,B} = begin
    ispersistent(u) && return u
    # Clear the transient flag.
    u._address = _trietopers(u._address)
    return u
end
persist!(u::TrieNode{T,S0,S,B}) where {T,S0,S,B} = begin
    ispersistent(u) && return u
    # Persist all the children first.
    bits = u._bits
    while bits > 0x0
        idx = trailing_zeros(bits)
        flag = B(0x1) << idx
        idx += 1
        bits &= ~flag
        persist!(u._elements[idx])
    end
    # Clear the transient flag.
    u._address = _trietopers(u._address)
    return u
end


# Public =======================================================================
# Below is the public interface of the Trie API (not everything below here is
# intended to be public, per se, but the main public interface is defined
# below).

# # Types ======================================================================
"""
    TTrie{T}

TTrie{T} is a transient trie that maps unsigned integers (i.e., hash values) to
objects of type T. Transient Tries are like persistent Tries in they're
structured in a way that makes it fairly easy/fast to (1) create a TTrie
duplicate of a PTrie [O(1)], (2) update a TTrie in a single thread, and (3)
create a PTrie duplicate of the TTrie [O(n) technically, but very fast in
practice].
"""
mutable struct TTrie{T} <: AbstractDict{TRIEKEY_T, T}
    var"_!private!unsafe_root"::Trie{T}
end
mutability(::Type{TTrie}) = Mutable()
mutability(::Type{TTrie{T}}) where {T} = Mutable()
equalfn(::Type{TTrie}) = (===)
equalfn(::Type{TTrie{T}}) where {T} = (===)
"""
    PTrie{T}

PTrie{T} is a persistent trie that maps unsigned integers (i.e., hash values) to
objects of type T. Persistent Tries are like transient Tries in they're
structured in a way that makes it fairly easy/fast to (1) create a TTrie
duplicate of a PTrie [O(1)], (2) update a TTrie in a single thread, and (3)
create a PTrie duplicate of the TTrie [O(n) technically, but very fast in
practice].
"""
struct PTrie{T} <: AbstractDict{TRIEKEY_T, T}
    var"_!private!unsafe_root"::Trie{T}
end
mutability(::Type{PTrie}) = Immutable()
mutability(::Type{PTrie{T}}) where {T} = Immutable()
equalfn(::Type{PTrie}) = (===)
equalfn(::Type{PTrie{T}}) where {T} = (===)

# # Constructors ===============================================================
#   First, constructors for making empty Tries.
PTrie{T}() where {T} = PTrie{T}(Trie{T}(TRIEKEY_T(0x0), true))
TTrie{T}() where {T} = TTrie{T}(Trie{T}(TRIEKEY_T(0x0), false))
PTrie() = PTrie{Any}()
TTrie() = TTrie{Any}()
#   Now for copying trees of the same type.
PTrie{T}(u::PTrie{S}) where {T,S} =
    PTrie{T}(convert(Trie{T}, u.var"_!private!unsafe_root"))
PTrie{T}(u::PTrie{T}) where {T} = u
PTrie(u::PTrie{T}) where {T} = u
TTrie{T}(u::TTrie{S}) where {T,S} =
    TTrie{T}(convert(Trie{S}, u.var"_!private!unsafe_root"))
TTrie{T}(u::TTrie{T}) where {T} =
    TTrie{T}(copy(u.var"_!private!unsafe_root"))
TTrie(u::TTrie{T}) where {T} = TTrie{T}(u)
"""
    PTrie(t::TTrie)

A PTrie can be constructed efficiently from a TTrie. In this construction of
such a transformation, the TTrie object t is not modified internally, so the
construction of the PTrie is a read-only operation on t. In contrast to the
persist!(t) function, the constructor version of TTrie persistence is slightly
slower and requires more memory allocation than the persist! version, but,
unlike persist!, it can safely be done concurrently with other read-only
operations.
"""
PTrie{T}(t::TTrie{U}) where {T,U} = 
    PTrie{T}(persist!(convert(Trie{T}, t.var"_!private!unsafe_root")))
PTrie{T}(t::TTrie{T}) where {T} = 
    PTrie{T}(persist!(copy(t.var"_!private!unsafe_root")))
PTrie(t::TTrie{T}) where {T} = PTrie{T}(t)
"""
    persist!(t::TTrie)

A PTrie can be constructed efficiently from a TTrie. In this construction of
such a transformation, the TTrie object t is modified internally; though its
external representation will not change at all. This is important to know,
however, as constructing a PTrie from a TTrie using the persist! function is a
write operation that cannot be safely be distributed across threads at the same
time as cross-thread read operations or other write operations. For a read-only
version of the same operation, one can use PTrie(t).
"""
persist!(t::TTrie{T}) where {T} =
    PTrie{T}(persist!(t.var"_!private!unsafe_root"))
# It is perfectly fine for a TTrie to point at a persistent Trie, so conversion
# in the other direction is even easier.
"""
    TTrie(p::PTrie)

A TTrie can be constructed from a PTrie in O(1) time.
"""
TTrie{T}(p::PTrie{T}) where {T} =
    TTrie{T}(p.var"_!private!unsafe_root")
TTrie{T}(p::PTrie{U}) where {T,U} =
    TTrie{T}(convert(Trie{T}, p.var"_!private!unsafe_root"))
TTrie(p::PTrie{T}) where {T} = TTrie{T}(p.var"_!private!unsafe_root")
#   Now for constructors from abstract dictionaries.
TTrie{T}(p::Pair{K,V}) where {T,K<:AbstractVector{TRIEKEY_T},V<:AbstractVector} = begin
    u = TTrie{T}()
    for (k,v) in zip(p.first, p.second)
        u[k] = v
    end
    return u
end
TTrie(p::Pair{K,V}) where {K<:AbstractVector{TRIEKEY_T},V<:AbstractVector} = begin
    T = typejoin(map(typeof, p.second))
    return TTrie{T}(p)
end
PTrie{T}(p::Pair{J,U}) where {T,K,J<:AbstractVector{K},V,U<:AbstractVector{V}} = begin
    return persist!(TTrie{T}(p))
end
PTrie(p::Pair{J,U}) where {K,J<:AbstractVector{K},V,U<:AbstractVector{V}} = begin
    return persist!(TTrie(p))
end
const TrieArgType = Union{Pair{TRIEKEY_T,<:Any}, Tuple{TRIEKEY_T,<:Any}}
TTrie{T}(kvs::Vararg{TrieArgType,N}) where {T,N} = begin
    u = TTrie{T}()
    for (k,v) in kvs
        u[k] = v
    end
    return u
end
TTrie(kvs::Vararg{TrieArgType,N}) where {N} = begin
    T = typejoin([typeof(v) for (k,v) in kvs])
    return TTrie{T}(kvs...)
end
PTrie{T}(kvs::Vararg{TrieArgType,N}) where {T,N} = persist!(TTrie{T}(kvs...))
PTrie(kvs::Vararg{TrieArgType,N}) where {T,N} = persist!(TTrie(kvs...))


# ==============================================================================
# Methods
#
# Definitions of methods of the above type.
Base.empty(u::PTrie{T}) where {T} = PTrie{T}()
Base.empty(u::PTrie{T}, S::Type) where {T} = PTrie{S}()
Base.empty(u::TTrie{T}) where {T} = TTrie{T}()
Base.empty(u::TTrie{T}, S::Type) where {T} = TTrie{S}()
Base.isempty(u::PTrie) = (length(u) == 0)
Base.isempty(u::TTrie) = (length(u) == 0)
Base.length(u::PTrie) = length(u.var"_!private!unsafe_root")
Base.length(u::TTrie) = length(u.var"_!private!unsafe_root")
Base.convert(::Type{PTrie{T}}, t::PTrie{S}) where {T,S} = PTrie{T}(t)
Base.convert(::Type{TTrie{T}}, t::TTrie{S}) where {T,S} = TTrie{T}(t)
Base.copy(t::PTrie{T}) where {T} = t
Base.copy(t::TTrie{T}) where {T} = TTrie{T}(t)
Base.get(u::PTrie{T}, k::TRIEKEY_T, df) where {T} =
    get(u.var"_!private!unsafe_root", k, df)
Base.get(u::TTrie{T}, k::TRIEKEY_T, df) where {T} =
    get(u.var"_!private!unsafe_root", k, df)
Base.in(kv::Pair{TRIEKEY_T,U}, u::PTrie{T}) where {U,T} =
    in(kv, u.var"_!private!unsafe_root")
Base.in(kv::Pair{TRIEKEY_T,U}, u::TTrie{T}) where {U,T} =
    in(kv, u.var"_!private!unsafe_root")
Base.in(kv::Pair{TRIEKEY_T,U}, u::PTrie{T}, f) where {U,T} =
    in(kv, u.var"_!private!unsafe_root", f)
Base.in(kv::Pair{TRIEKEY_T,U}, u::TTrie{T}, f) where {U,T} =
    in(kv, u.var"_!private!unsafe_root", f)
Base.iterate(t::PTrie{T})   where {T} = iterate(t.var"_!private!unsafe_root")
Base.iterate(t::PTrie{T},k) where {T} = iterate(t.var"_!private!unsafe_root",k)
Base.iterate(t::TTrie{T})   where {T} = iterate(t.var"_!private!unsafe_root")
Base.iterate(t::TTrie{T},k) where {T} = iterate(t.var"_!private!unsafe_root",k)
Base.haskey(p::PTrie{T}, k::TRIEKEY_T) where {T} =
    haskey(p.var"_!private!unsafe_root", k)
Base.haskey(t::TTrie{T}, k::TRIEKEY_T) where {T} =
    haskey(t.var"_!private!unsafe_root", k)
setindex(t::PTrie{T}, u::S, k::TRIEKEY_T) where {T,S} = begin
    root0 = t.var"_!private!unsafe_root"
    root1 = setindex(root0, u, k)
    (root0 === root1) && return t
    return PTrie{T}(root1)
end
Base.setindex!(t::TTrie{T}, u::S, k::TRIEKEY_T) where {T,S} =
    setindex!(t.var"_!private!unsafe_root", u, k)
delete(t::PTrie{T}, k::TRIEKEY_T) where {T} = begin
    idx = triekey_child_index(t._id, k)
    idxbit = TRIEBITS_T(1) << (idx - 1)
    if idx == 0
        # The node lies outside of this node, so nothing to do.
        return t
    elseif t._bits & idxbit == 0
        # Nothing there
        return t
    elseif isa(t._data, Vector{T})
        # We'e a twig node, so we just clear the bit
        return PTrie{T}(t._id, t._bits & ~idxbit, t._data)
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
        return PTrie{T}(t._id, bits, ch)
    end
end
