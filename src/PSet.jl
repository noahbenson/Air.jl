################################################################################
# PSet
# A persistent set type that uses PVec and PMap types.

import Base.IdSet

# Before we declare our own sets, let's use a clever trick to make a mutable set
# type based on equivalence.

"""
    EquivRef{T}
An EquivRef wraps an object and mimics its equivhash() and equiv() functions on
the hash() and isequal() functions.
"""
struct EquivRef{T} <: Ref{T}
    object::T
end
Base.hash(r::EquivRef) = equivhash(r.object)
Base.isequal(r::EquivRef, s::EquivRef) = isequiv(r.object, s.object)
Base.:(==)(r::EquivRef, s::EquivRef) = isequiv(r.object, s.object)
Base.getindex(r::EquivRef{T}) where {T} = r.object
"""
    EquivSet{T}
The EquivSet type is like the Set and IdSet types, except that it works on the
isequiv() and equivhash() functions instead of the isequal() and hash()
functions or the === and objectid() functions.
"""
struct EquivSet{T} <: AbstractSet{T}
    set::Set{EquivRef{T}}
end
EquivSet{T}() where {T} = EquivSet{T}(Set{EquivRef{T}}())
EquivSet() = EquivSet{Any}()
EquivSet{T}(arr::AbstractArray{S,1}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    return EquivSet{T}(Set{EquivRef{T}}(ers))
end
EquivSet{T}(arr::AbstractSet{S}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    return EquivSet{T}(Set{EquivRef{T}}(ers))
end
EquivSet(arr::AbstractArray{S,1}) where {T,S<:T} = begin
    let T = typejoin([typeof(x) for x in arr]...), ers = [EquivRef{T}(x) for x in arr]
        return EquivSet{T}(Set{EquivRef{T}}(ers))
    end
end
EquivSet(arr::AbstractSet{S}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    let T = typejoin([typeof(x) for x in arr]...), ers = [EquivRef{T}(x) for x in arr]
        return EquivSet{T}(Set{EquivRef{T}}(ers))
    end
end
Base.length(s::EquivSet) = length(s.set)
Base.in(x::S, s::EquivSet{T}) where {T, S<:T} = in(EquivRef{T}(x), s.set)
Base.in(x::S, s::EquivSet{T}, f) where {T, S<:T} = in(EquivRef{T}(x), s.set, f)
Base.iterate(s::EquivSet{T}) where {T} = let u = iterate(s.set)
    (u === nothing) && return nothing
    return (u[1].object, u[2])
end
Base.iterate(s::EquivSet{T}, it) where {T} = let u = iterate(s.set, it)
    (u === nothing) && return nothing
    return (u[1].object, u[2])
end
Base.push!(s::EquivSet{T}, x::S) where {T, S<:T} = begin
    push!(s.set, EquivRef{T}(x))
    return s
end
Base.delete!(s::EquivSet{T}, x::S) where {T, S<:T} = begin
    delete!(s.set, EquivRef{T}(x))
    return s
end
Base.show(io::IO, m::MIME"text/plain", s::EquivSet{T}) where {T} = begin
    print(io, "EquivSet(")
    (T === Any) || print(io, repr(m, T))
    let ii = 1, st = iterate(s)
        while st !== nothing && ii < 20
            (ii > 1) && print(io, ", ")
            print(io, repr(m, st[1].object))
            st = iterate(s, st[2])
            ii += 1
        end
    end
    (st === nothing) || print(io, " ...")
    print(io, "])")
end
equalfn(::Type{EquivSet}) = isequiv
equalfn(::Type{EquivSet{T}}) where {T} = isequiv
hashfn(::Type{EquivSet}) = equivhash
hashfn(::Type{EquivSet{T}}) where {T} = equivhash
        

# Okay, now declare the persistent types:
abstract type PSet{T} <: AbstractSet{T} end
"""
    PIdSet{T}
A persistent version of IdSet{T}.
"""
abstract type PIdSet{T} <: PSet{T} end
"""
    PEqualSet{T}
A persistent form of Set{T}.
"""
abstract type PEqualSet{T} <: PSet{T} end
"""
    PEquivSet{T}
A persistent set based on object equivalence rather than equality or identity
(see isequiv and equivhash). This is the default kind of PSet.
"""
abstract type PEquivSet{T} <: PSet{T} end

# setup our equality functions functions
equalfn(::Type{T}) where {T <: PIdSet} = (===)
equalfn(::Type{T}) where {T <: PEqualSet} = isequal
equalfn(::Type{T}) where {T <: PEquivSet} = isequiv
hashfn(::Type{T}) where {T <: PIdSet} = objectid
hashfn(::Type{T}) where {T <: PEqualSet} = hash
hashfn(::Type{T}) where {T <: PEquivSet} = equivhash
equalfn(::PIdSet) = (===)
equalfn(::PEqualSet) = isequal
equalfn(::PEquivSet) = isequiv
hashfn(::PIdSet) = objectid
hashfn(::PEqualSet) = hash
hashfn(::PEquivSet) = equivhash

_pset_show(io::IO, u::PSet{T}, name::String) where {T} = begin
    print(io, "$(name)($(repr(T))[")
    if length(u) > 40
        for (i,x) in enumerate(u)
            (x > 40) && break
            print(io, (i == 1 ? "" : ", ") * repr(x))
        end
        print(io, " ... ")
    else
        for (i,x) in enumerate(u)
           print(io, (i == 1 ? "" : ", ") * repr(x))
        end 
    end
    print(io, "])")
end
Base.show(io::IO, ::MIME"text/plain", u::PSet{T}) where {T} = _pset_show(io, u, "PSet")
Base.show(io::IO, ::MIME"text/plain", u::PIdSet{T}) where {T} = _pset_show(io, u, "PIdSet")
Base.show(io::IO, ::MIME"text/plain", u::PEqualSet{T}) where {T} = _pset_show(io, u, "PEqualSet")
Base.show(io::IO, ::MIME"text/plain", u::PEquivSet{T}) where {T} = _pset_show(io, u, "PSet")

_iseq(t::PSet{T}, s::AbstractSet{S}, eqfn::Function) where {T,S} = begin
    (length(t) == length(s)) || return false
    # we want the general types to match...
    (equalfn(typeof(t)) === equalfn(typeof(s))) || return false
    for x in s
        in(x, t, eqfn) || return false
    end
    return true
end
Base.isequal(t::PSet{T}, s::AbstractSet{S}) where {T,S} = _iseq(t, s, isequal)
Base.isequal(s::AbstractSet{S}, t::PSet{T}) where {T,S} = _iseq(t, s, isequal)
isequiv(t::PSet{T}, s::AbstractSet{S}) where {T,S} = _iseq(t, s, isequiv)
isequiv(s::AbstractSet{S}, t::PSet{T}) where {T,S} = _iseq(t, s, isequiv)

# PVectorSet objects don't attenpt to differentiate objects; they just act in
# a naive linear fashion; this makes them good if, for example, you need to
# keep track of objects that have already had a hash collision.
struct PIdVectorSet{T} <: PIdSet{T}
    _items::PVec{T}
    # When given a PVec or IdSet, PVectorSet-types don't have to check for correctness!
    PIdVectorSet{T}(u::IdSet{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    PIdVectorSet{T}(u::PIdSet{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    # Otherwise, they do:
    PIdVectorSet{T}(u::AbstractSet{S}) where {T,S<:T} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === (===) && hfn === objectid
                return new{T}(PVec{T}(T[u...]))
            else
                return new{T}(PVec{T}([IdSet{T}(u)...]))
            end
        end
    end
    PIdVectorSet{T}(u::AbstractArray{S,1}) where {T,S<:T} = new{T}(PVec{T}(T[IdSet{T}(u)...]))
    # However, it is too expensive to check every possible iteration of a set for
    # correctness (e.g. as we build up or edit a set), so we do not check when
    # we are passed a PVec that is of the appropriate type:
    PIdVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
end
PIdVectorSet{T}(u::PIdVectorSet{T}) where {T} = u
struct PEqualVectorSet{T} <: PEqualSet{T}
    _items::PVec{T}
    PEqualVectorSet{T}(u::Set{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    PEqualVectorSet{T}(u::PEqualSet{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    PEqualVectorSet{T}(u::AbstractSet{S}) where {T,S<:T} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === isequal && hfn === hash
                return new{T}(PVec{T}(T[u...]))
            else
                return new{T}(PVec{T}([Set{T}(u)...]))
            end
        end
    end
    PEqualVectorSet{T}(u::AbstractArray{S,1}) where {T,S<:T} = new{T}(PVec{T}(T[Set{T}(u)...]))
    PEqualVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
end
PEqualVectorSet{T}(u::PEqualVectorSet{T}) where {T} = u
struct PVectorSet{T} <: PEquivSet{T}
    _items::PVec{T}
    PVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
    PVectorSet{T}(u::EquivSet{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    PVectorSet{T}(u::PSet{S}) where {T,S<:T} = new{T}(PVec{T}(T[u...]))
    PVectorSet{T}(u::AbstractSet{S}) where {T,S<:T} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === isequiv && hfn === equivhash
                return new{T}(PVec{T}(T[u...]))
            else
                return new{T}(PVec{T}([EquivSet{T}(u)...]))
            end
        end
    end
    PVectorSet{T}(u::AbstractArray{S,1}) where {T,S<:T} = new{T}(PVec{T}(T[EquivSet{T}(u)...]))
end
const PVECSETS = (Union{PIdVectorSet{T}, PEqualVectorSet{T}, PVectorSet{T}} where {T})
Base.convert(::Type{PVec{T}}, u::PVectorSet{T}) where {T} = u._items
Base.convert(::Type{PVec{T}}, u::PIdVectorSet{T}) where {T} = u._items
Base.convert(::Type{PVec{T}}, u::PEqualVectorSet{T}) where {T} = u._items
Base.length(u::PVectorSet) = length(u._items)
Base.length(u::PIdVectorSet) = length(u._items)
Base.length(u::PEqualVectorSet) = length(u._items)
Base.iterate(u::PVectorSet) = iterate(u._items)
Base.iterate(u::PIdVectorSet) = iterate(u._items)
Base.iterate(u::PEqualVectorSet) = iterate(u._items)
Base.iterate(u::PVectorSet, x) = iterate(u._items, x)
Base.iterate(u::PIdVectorSet, x) = iterate(u._items, x)
Base.iterate(u::PEqualVectorSet, x) = iterate(u._items, x)
Base.in(x::S, u::PVECSETS{T}, eqfn::Function) where {T, S<:T} = begin
    for y in convert(PVec{T}, u)
        eqfn(y, x) && return true
    end
    return false
end
Base.in(x::S, u::PVECSETS{T}) where {T, S<:T} = in(x, u, equalfn(u))
"""
    push(set, x)
Yields the set that is formed of the given set conjoined with
the object x.
"""
push(u::PVECSETS{T}, x::S) where {T,S<:T} = let U = typeof(u)
    (x in u) && return u
    return U(push(convert(PVec{T}, u), x))
end
"""
    delete(set, x)
Yields the set that is formed of the given set disjoined from
the object x.
"""
delete(u::PVECSETS{T}, x::S) where {T, S<:T} = begin
    let k, eq = equalfn(u), v = convert(PVec{T}, u), U = typeof(u)
        for (k,y) in enumerate(v)
            eq(y, x) && return U(PVec{T}(vcat(v[1:k-1], v[k+1:end])))
        end
        return u
    end
end

# PHashSet is a a set that uses persistent maps to keep track of items in and
# not in the set. The expectation is that the sets that make up the map values
# are PVectorSet type sets.
struct PHashSet{T} <: PEquivSet{T}
    _n::Integer
    _table::PMap{PSet{T}}
end
struct PIdHashSet{T} <: PIdSet{T}
    _n::Integer
    _table::PMap{PSet{T}}
end
struct PEqualHashSet{T} <: PEqualSet{T}
    _n::Integer
    _table::PMap{PSet{T}}
end
_pvecset_type(::Type{PHashSet}) = PVectorSet
_pvecset_type(::Type{PIdHashSet}) = PIdVectorSet
_pvecset_type(::Type{PEqualHashSet}) = PEqualVectorSet
_pvecset_type(::Type{PHashSet{T}}) where {T} = PVectorSet{T}
_pvecset_type(::Type{PIdHashSet{T}}) where {T} = PIdVectorSet{T}
_pvecset_type(::Type{PEqualHashSet{T}}) where {T} = PEqualVectorSet{T}
_pvecset_type(::Type{PVectorSet}) = PVectorSet
_pvecset_type(::Type{PIdVectorSet}) = PIdVectorSet
_pvecset_type(::Type{PEqualVectorSet}) = PEqualVectorSet
_pvecset_type(::Type{PVectorSet{T}}) where {T} = PVectorSet{T}
_pvecset_type(::Type{PIdVectorSet{T}}) where {T} = PIdVectorSet{T}
_pvecset_type(::Type{PEqualVectorSet{T}}) where {T} = PEqualVectorSet{T}
_pvecset_type(::Type{EquivSet}) = PVectorSet
_pvecset_type(::Type{IdSet}) = PIdVectorSet
_pvecset_type(::Type{Set}) = PEqualVectorSet
_pvecset_type(::Type{EquivSet{T}}) where {T} = PVectorSet{T}
_pvecset_type(::Type{IdSet{T}}) where {T} = PIdVectorSet{T}
_pvecset_type(::Type{Set{T}}) where {T} = PEqualVectorSet{T}
_pvecset_type(::Type{AbstractSet}) = PEqualVectorSet
_pvecset_type(::Type{AbstractSet{T}}) where {T} = PEqualVectorSet{T}
_phashset_type(::Type{PVectorSet}) = PHashSet
_phashset_type(::Type{PIdVectorSet}) = PIdHashSet
_phashset_type(::Type{PEqualVectorSet}) = PEqualHashSet
_phashset_type(::Type{PVectorSet{T}}) where {T} = PHashSet{T}
_phashset_type(::Type{PIdVectorSet{T}}) where {T} = PIdHashSet{T}
_phashset_type(::Type{PEqualVectorSet{T}}) where {T} = PEqualHashSet{T}
_phashset_type(::Type{PHashSet}) = PHashSet
_phashset_type(::Type{PIdHashSet}) = PIdHashSet
_phashset_type(::Type{PEqualHashSet}) = PEqualHashSet
_phashset_type(::Type{PHashSet{T}}) where {T} = PHashSet{T}
_phashset_type(::Type{PIdHashSet{T}}) where {T} = PIdHashSet{T}
_phashset_type(::Type{PEqualHashSet{T}}) where {T} = PEqualHashSet{T}
_phashset_type(::Type{EquivSet}) = PVectorSet
_phashset_type(::Type{IdSet}) = PIdHashSet
_phashset_type(::Type{Set}) = PEqualHashSet
_phashset_type(::Type{EquivSet{T}}) where {T} = PVectorSet{T}
_phashset_type(::Type{IdSet{T}}) where {T} = PIdHashSet{T}
_phashset_type(::Type{Set{T}}) where {T} = PEqualHashSet{T}
_phashset_type(::Type{AbstractSet}) = PEqualHashSet
_phashset_type(::Type{AbstractSet{T}}) where {T} = PEqualHashSet{T}
# For hash-type sets:
_pset_init(P::DataType, arr::AbstractArray{S,1}) where {S} = begin
    let P = (isa(P, UnionAll) ? P{S} : P), T = P.parameters[1],
        V = _pvecset_type(P), u = TMap{PSet{T}}(), eqfn = equalfn(P), hfn = hashfn(P), h, x
        for a in arr
            h = hfn(a)
            x = get(u, h, nothing)
            if (x === nothing)
                u[h] = V([a])
            else
                u[h] = push(x, a)
            end
        end
        return P(length(u), freeze(u))
    end
end
_pset_init(P::DataType, s::AbstractSet{S}) where {S} = begin
    let P = (isa(P, UnionAll) ? P{S} : P), T = P.parameters[1]
        return _pset_init(P, T[s...])
    end
end
# PHash type methods:
const PHASHSETS = (Union{PHashSet{T}, PIdHashSet{T}, PEqualHashSet{T}} where {T})
Base.convert(::Type{PMap{PSet{T}}}, u::PHashSet{T}) where {T} = u._table
Base.convert(::Type{PMap{PSet{T}}}, u::PIdHashSet{T}) where {T} = u._table
Base.convert(::Type{PMap{PSet{T}}}, u::PEqualHashSet{T}) where {T} = u._table
Base.length(u::PHashSet)       = u._n
Base.length(u::PIdHashSet)     = u._n
Base.length(u::PEqualHashSet)  = u._n
Base.iterate(u::PHASHSETS{T}) where {T} = let t = convert(PMap{PSet{T}}, u), x = iterate(t)
    if x  === nothing
        return x
    else
        let ((k,vs),addr1) = x, (v,addr2) = iterate(vs)
            return (v, (t, vs, addr1, addr2))
        end
    end
end
Base.iterate(u::PHASHSETS{T}, addr::Tuple) where {T} = let (t,vs,ad1,ad2) = addr, x
    x = iterate(vs, ad2)
    if x === nothing
        x = iterate(t, ad1)
        (x === nothing) && return nothing
        let ((k,vs),addr1) = x, (v,addr2) = iterate(vs)
            return (v, (t, vs, addr1, addr2))
        end
    else
        let (v,addr2) = x
            return (v, (t, vs, ad1, addr2))
        end
    end
end
Base.in(x::S, u::PHASHSETS{T}, eqfn::Function) where {T, S<:T} = begin
    let h = UInt(hashfn(u)(x)), t = convert(PMap{PSet{T}}, u), ys = get(t, h, t)
        (ys === t) && return false
        return in(x, ys, eqfn)
    end
end
Base.in(x::S, u::PHASHSETS{T}) where {T, S<:T} = in(x, u, equalfn(u))
push(u::PHASHSETS{T}, x::S) where {T, S<:T} = begin
    let U = typeof(u), t = convert(PMap{PSet{T}}, u), eqfn = equalfn(u),
        h = UInt(hashfn(u)(x)), v = get(t, h, nothing)
        if v === nothing
            return U(u._n + 1, assoc(t, h, _pvecset_type(U)(PVec1{T}((x,)))))
        elseif in(x, v, eqfn)
            return u
        else
            let vv = push(v, x)
                if vv === v
                    return u
                else
                    return U(u._n + 1, assoc(t, h, vv))
                end
            end
        end
    end
end
delete(u::PHASHSETS{T}, x::S) where {T, S<:T} = begin
    let U = typeof(u), t = convert(PMap{PSet{T}}, u), eqfn = equalfn(u),
        h = UInt(hashfn(u)(x)), v = get(t, h, nothing)
        (v === nothing || !in(x, v)) && return u
        return U(u._n - 1,
                 length(v) == 1
                 ? dissoc(t, h)
                 : assoc(t, h, delete(v, x)))
    end
end
# Constructors:
PHashSet{T}() where {T} = PHashSet{T}(0, PMap{PSet{T}}())
PHashSet{T}(s::PHashSet{T}) where {T} = s
PHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(PHashSet{T}, s)
PHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(PHashSet{T}, arr)
PIdHashSet{T}() where {T} = PIdHashSet{T}(0, PMap{PSet{T}}())
PIdHashSet{T}(s::PIdHashSet{T}) where {T} = s
PIdHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(PIdHashSet{T}, s)
PIdHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(PIdHashSet{T}, arr)
PEqualHashSet{T}() where {T} = PEqualHashSet(0, PMap{PSet{T}}())
PEqualHashSet{T}(s::PEqualHashSet{T}) where {T} = s
PEqualHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(PEqualHashSet{T}, s)
PEqualHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(PEqualHashSet{T}, arr)

PSet{T}() where {T} = PHashSet{T}()
PSet{T}(s::AbstractSet{S}) where {T, S<:T} = PHashSet{T}(s)
PSet{T}(s::AbstractArray{S,1}) where {T, S<:T} = PHashSet{T}(s)
PIdSet{T}() where {T} = PIdHashSet{T}()
PIdSet{T}(s::AbstractSet{S}) where {T, S<:T} = PIdHashSet{T}(s)
PIdSet{T}(s::AbstractArray{S,1}) where {T, S<:T} = PIdHashSet{T}(s)
PEqualSet{T}() where {T} = PEqualHashSet{T}()
PEqualSet{T}(s::AbstractSet{S}) where {T, S<:T} = PEqualHashSet{T}(s)
PEqualSet{T}(s::AbstractArray{S,1}) where {T, S<:T} = PEqualHashSet{T}(s)

_pset_construct(PS::Type, s::AbstractArray{T,1}) where {T} = begin
    let U = typejoin(map(typeof, s)...)
        return PS{U}(U[s...])
    end
end
_pset_construct(PS::Type, s::AbstractSet{T}) where {T} = begin
    let U = typejoin([typeof(u) for u in s]...)
        return PS{U}(U[s...])
    end
end
PVectorSet(s::AbstractArray{T,1})      where {T} = _pset_construct(PVectorSet, s)
PIdVectorSet(s::AbstractArray{T,1})    where {T} = _pset_construct(PIdVectorSet, s)
PEqualVectorSet(s::AbstractArray{T,1}) where {T} = _pset_construct(PEqualVectorSet, s)
PHashSet(s::AbstractArray{T,1})      where {T} = _pset_construct(PHashSet, s)
PIdHashSet(s::AbstractArray{T,1})    where {T} = _pset_construct(PIdHashSet, s)
PEqualHashSet(s::AbstractArray{T,1}) where {T} = _pset_construct(PEqualHashSet, s)
PSet(s::AbstractArray{T,1})      where {T} = _pset_construct(PHashSet, s)
PIdSet(s::AbstractArray{T,1})    where {T} = _pset_construct(PIdHashSet, s)
PEqualSet(s::AbstractArray{T,1}) where {T} = _pset_construct(PEqualHashSet, s)

PVectorSet(s::AbstractSet{T})      where {T} = _pset_construct(PVectorSet, s)
PIdVectorSet(s::AbstractSet{T})    where {T} = _pset_construct(PIdVectorSet, s)
PEqualVectorSet(s::AbstractSet{T}) where {T} = _pset_construct(PEqualVectorSet, s)
PHashSet(s::AbstractSet{T})      where {T} = _pset_construct(PHashSet, s)
PIdHashSet(s::AbstractSet{T})    where {T} = _pset_construct(PIdHashSet, s)
PEqualHashSet(s::AbstractSet{T}) where {T} = _pset_construct(PEqualHashSet, s)
PSet(s::AbstractSet{T})      where {T} = _pset_construct(PHashSet, s)
PIdSet(s::AbstractSet{T})    where {T} = _pset_construct(PIdHashSet, s)
PEqualSet(s::AbstractSet{T}) where {T} = _pset_construct(PEqualHashSet, s)
