################################################################################
# PSet
# A persistent set type that uses PVec and PMap types.

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
equalfn(::Type{PIdSet}) = (===)
equalfn(::Type{PEqualSet}) = isequal
equalfn(::Type{PEquivSet}) = isequiv
hashfn(::Type{PIdSet}) = objectid
hashfn(::Type{PEqualSet}) = hashcode
hashfn(::Type{PEquivSet}) = equivhash
equalfn(::PIdSet) = (===)
equalfn(::PEqualSet) = isequal
equalfn(::PEquivSet) = isequiv
hashfn(::PIdSet) = objectid
hashfn(::PEqualSet) = hashcode
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
    TT = typeof(t)
    SS = typeof(s)
    if TT <: PIdSet
        (SS <: Base.IdSet || SS <: PIdSet)  || return false
    elseif TT <: PEqualSet
        (SS <: Set || SS <: PEqualSet) || return false
    elseif TT <: PEquivSet
        (SS <: PEquivSet) || return false
    end
    for x in t
        in(x, s, eqfn) || return false
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
    PIdVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
end
struct PEqualVectorSet{T} <: PEqualSet{T}
    _items::PVec{T}
    PEqualVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
end
struct PVectorSet{T} <: PEquivSet{T}
    _items::PVec{T}
    PVectorSet{T}(u::PVec{T}) where {T} = new{T}(u)
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
Base.in(x::S, u::PVECSETS{T}) where {T, S<:T} = in(u, x, equalfn(u))
"""
    conj(set, x)
Yields the set that is formed of the given set conjoined with
the object x.
"""
conj(u::PVECSETS{T}, x::S) where {T,S<:T} = let U = typeof(u)
    (x in u) && return u
    return U(push(convert(PVec{T}, u), x))
end
"""
    disj(set, x)
Yields the set that is formed of the given set disjoined from
the object x.
"""
disj(u::PVECSETS{T}, x::S) where {T, S<:T} = begin
    let k, eq = equalfn(u), v = convert(PVec{T}, u), U = typeof(u)
        for (k,y) in enumerate(v)
            eq(y, x) && return U(PVec{T}(vcat(v[1:k-1], v[k+1:end])))
        end
        return u
    end
end
_pset_init(u0::PSet{T}, arr::AbstractArray{S, 1}) where {T,S<:T} = let u = u0
    for a in arr
        (a in u0) && continue
        u = conj(u, a)
    end
    return u
end
_pset_init(u0::PSet{T}, s::AbstractSet{S}) where {T,S<:T} = let u = u0
    for a in s
        (a in u0) && continue
        u = conj(u, a)
    end
    return u
end
PVectorSet{T}() where {T} = PVectorSet{T}(PVec{T}(T[]))
PVectorSet{T}(s::PVectorSet{S}) where {T,S<:T} = s
PVectorSet{T}(s::AbstractSet{S}) where {T,S<:T} = _pset_init(
    PVectorSet{T}(PVec{T}(T[])),
    s)
PVectorSet{T}(arr::AbstractArray{S, 1}) where {T,S<:T} = _pset_init(
    PVectorSet{T}(PVec{T}(T[])),
    arr)
PIdVectorSet{T}() where {T} = PIdVectorSet{T}(PVec{T}(T[]))
PIdVectorSet{T}(s::PIdVectorSet{S}) where {T,S<:T} = s
PIdVectorSet{T}(s::AbstractSet{S}) where {T,S<:T} = _pset_init(
    PIdVectorSet{T}(PVec{T}(T[])),
    s)
PIdVectorSet{T}(arr::AbstractArray{S, 1}) where {T,S<:T} = _pset_init(
    PIdVectorSet{T}(PVec{T}(T[])),
    arr)
PEqualVectorSet{T}() where {T} = PEqualVectorSet{T}(PVec{T}(T[]))
PEqualVectorSet{T}(s::PEqualVectorSet{S}) where {T,S<:T} = s
PEqualVectorSet{T}(s::AbstractSet{S}) where {T,S<:T} = _pset_init(
    PEqualVectorSet{T}(PVec{T}(T[])),
    s)
PEqualVectorSet{T}(arr::AbstractArray{S, 1}) where {T,S<:T} = _pset_init(
    PEqualVectorSet{T}(PVec{T}(T[])),
    arr)

# PHashSet is a a set that uses persistent maps to keep track of items in and
# not in the set. The expectation is that the sets that make up the map values
# are PVectorSet type sets.
struct PHashSet{T} <: PEquivSet{T}
    _table::PMap{PSet{T}}
end
struct PIdHashSet{T} <: PIdSet{T}
    _table::PMap{PSet{T}}
end
struct PEqualHashSet{T} <: PEqualSet{T}
    _table::PMap{PSet{T}}
end
const PHASHSETS = (Union{PHashSet{T}, PIdHashSet{T}, PEqualHashSet{T}} where {T})
Base.convert(::Type{PMap{PSet{T}}}, u::PHashSet{T}) where {T} = u._table
Base.convert(::Type{PMap{PSet{T}}}, u::PIdHashSet{T}) where {T} = u._table
Base.convert(::Type{PMap{PSet{T}}}, u::PEqualHashSet{T}) where {T} = u._table
Base.length(u::PHashSet)       = length(u._table)
Base.length(u::PIdHashSet)     = length(u._table)
Base.length(u::PEqualHashSet)  = length(u._table)
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
conj(u::PHASHSETS{T}, x::S) where {T, S<:T} = begin
    let U = typeof(u), t = convert(PMap{PSet{T}}, u), eq = equalfn(u),
        h = UInt(hashfn(u)(x)), v = get(t, h, t)
        if v === t
            return U(assoc(t, h, PIdVectorSet{T}(PVec1{T}(x))))
        elseif in(x, v, eq)
            return u
        else
            return U(assoc(t, h, push(v, x)))
        end
    end
end
disj(u::PHASHSETS{T}, x::S) where {T, S<:T} = begin
    let U = typeof(u), t = convert(PMap{PSet{T}}, u), eq = equalfn(u),
        h = UInt(hashfn(u)(x)), v = get(t, h, t)
        if v === t
            return u
        elseif in(x, v, eq)
            return U(length(v) == 1
                     ? dissoc(t, h)
                     : assoc(t, h, dissoc(v, h)))
        else
            return u
        end
    end
end

PHashSet{T}() where {T} = PHashSet(PMap{PSet{T}}())
PHashSet{T}(s::PHashSet{T}) where {T} = s
PHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(
    PHashSet(PMap{PSet{T}}()),
    s)
PHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(
    PHashSet(PMap{PSet{T}}()),
    arr)
PIdHashSet{T}() where {T} = PIdHashSet(PMap{PSet{T}}())
PIdHashSet{T}(s::PIdHashSet{T}) where {T} = s
PIdHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(
    PIdHashSet(PMap{PSet{T}}()),
    s)
PIdHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(
    PIdHashSet(PMap{PSet{T}}()),
    arr)
PEqualHashSet{T}() where {T} = PEqualHashSet(PMap{PSet{T}}())
PEqualHashSet{T}(s::PEqualHashSet{T}) where {T} = s
PEqualHashSet{T}(s::AbstractSet{S}) where {T, S<:T} = _pset_init(
    PEqualHashSet(PMap{PSet{T}}()),
    s)
PEqualHashSet{T}(arr::AbstractArray{S, 1}) where {T, S<:T} = _pset_init(
    PEqualHashSet(PMap{PSet{T}}()),
    arr)

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



