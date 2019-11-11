################################################################################
# PDict
# A persistent dictionary type that uses PVec and PMap types.

# Before we declare our own dicts, let's use a clever trick to make a mutable
# dict type based on equivalence.

"""
    EquivDict{T}
The EquivDict type is like the Dict and IdDict types, except that it works on
the isequiv() and equivhash() functions instead of the isequal() and hash()
functions or the === and objectid() functions.
"""
struct EquivDict{K,V} <: AbstractDict{K,V}
    dict::Dict{EquivRef{K},EquivRef{V}}
end
EquivDict{K,V}() where {K,V} = EquivDict{K,V}(Dict{EquivRef{K},EquivRef{V}}())
EquivDict() = EquivDict{Any,Any}()
EquivDict{K,V}(arr::AbstractArray{S,1}) where {K,V,S} = begin
    let kvs = [EquivRef{K}(k) => EquivRef{V}(v) for (k,v) in arr]
        return EquivDict{K,V}(Dict{EquivRef{K},EquivRef{V}}(kvs...))
    end
end
EquivDict{K,V}(arr::AbstractDict{M,U}) where {K,V,M<:K,U<:V} = begin
    let kvs = [EquivRef{K}(k) => EquivRef{V}(v) for (k,v) in arr]
        return EquivDict{K,V}(Dict{EquivRef{K},EquivRef{V}}(kvs...))
    end
end
EquivDict(arr::Union{AbstractArray,AbstractDict}) = begin
    let K = [], V = []
        for (k,v) in arr:
            push!(K, k)
            push!(V, v)
        end
        (K, V) = (typejoin(K...), typejoin(V...))
        return EquivDict{K,V}(arr)
    end
end
Base.length(s::EquivDict) = length(s.dict)
Base.in(x::Pair{KK,VV}, s::EquivDict{K,V}) where {K,V,KK<:K,VV<:V} = begin
    return in(EquivRef{K}(x[1]) => EquivRef{V}(x[2]), s.dict)
end
Base.in(x::Pair{KK,VV}, s::EquivDict{K,V}, f::Function) where {K,V,KK<:K,VV<:V} = begin
    return in(EquivRef{K}(x[1]) => EquivRef{V}(x[2]), s.dict, f)
end
Base.iterate(s::EquivDict{K,V}) where {K,V} = let u = iterate(s.dict)
    (u === nothing) && return nothing
    let (k,v) = u[1]
        return (k.object => v.object, u[2])
    end
end
Base.iterate(s::EquivDict{K,V}, it) where {K,V} = let u = iterate(s.dict, it)
    (u === nothing) && return nothing
    let (k,v) = u[1]
        return (k.object => v.object, u[2])
    end
end
Base.push!(s::EquivDict{K,V}, x::Pair{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
    push!(s.dict, EquivRef{K}(x[1]) => EquivRef{V}(x[2]))
    return s
end
Base.setindex!(s::EquivDict{K,V}, v::VV, k::KK) where {K,V,KK<:K,VV<:V} = begin
    setindex!(s.dict, EquivRef{V}(v), EquivRef{K}(k))
    return s
end
Base.delete!(s::EquivDict{K,V}, k::KK) where {K,V,KK<:K} = begin
    delete!(s.set, EquivRef{K}(k))
    return s
end
Base.show(io::IO, m::MIME"text/plain", s::EquivDict{K,V}) where {K,V} = begin
    print(io, "EquivDict")
    (K === Any && V === Any) || print(io, "{$(repr(m, K)),$(repr(m, V))}")
    print(io, "(")
    let ii = 1, st = iterate(s)
        st === nothing && break
        ii >= 20 && break
        print(io, (ii == 1? "" : ", "))
        (k,v) = st[1]
        print(io, repr(m, k.object => v.object))
        st = iterate(s, st[2])
        ii += 1
    end
    (st === nothing) || print(io, " ...")
    print(io, ")")
end
equalfn(::Type{EquivDict}) = isequiv
equalfn(::Type{EquivDict{K,V}}) where {K,V} = isequiv
hashfn(::Type{EquivDict}) = equivhash
hashfn(::Type{EquivDict{K,V}}) where {K,V} = equivhash
        

# Okay, now declare the persistent types:
abstract type PDict{K,V} <: AbstractDict{K,V} end
"""
    PIdDict{K,V}
A persistent version of IdDict{K,V}.
"""
abstract type PIdDict{K,V} <: PDict{K,V} end
"""
    PEqualDict{K,V}
A persistent form of Dict{K,V}.
"""
abstract type PEqualDict{K,V} <: PDict{K,V} end
"""
    PEquivDict{K,V}
A persistent dict based on object equivalence rather than equality or identity
(see isequiv and equivhash). This is the default kind of PDict.
"""
abstract type PEquivDict{K,V} <: PDict{K,V} end

# setup our equality functions functions
equalfn(::Type{T}) where {T <: PIdDict} = (===)
equalfn(::Type{T}) where {T <: PEqualDict} = isequal
equalfn(::Type{T}) where {T <: PEquivDict} = isequiv
hashfn(::Type{T}) where {T <: PIdDict} = objectid
hashfn(::Type{T}) where {T <: PEqualDict} = hash
hashfn(::Type{T}) where {T <: PEquivDict} = equivhash
equalfn(::PIdDict) = (===)
equalfn(::PEqualDict) = isequal
equalfn(::PEquivDict) = isequiv
hashfn(::PIdDict) = objectid
hashfn(::PEqualDict) = hash
hashfn(::PEquivDict) = equivhash
_pdict_show(io::IO, u::PDict{K,V}, name::String) where {K,V} = begin
    print(io, name)
    (K === Any && V === Any) || print(io, "{$(repr(m, K)),$(repr(m, V))}")
    print(io, "(")
    let ii = 1, st = iterate(s)
        st === nothing && break
        ii >= 20 && break
        (ii > 1) && print(io, ", ")
        print(io, repr(m, st[1]))
        st = iterate(s, st[2])
        ii += 1
    end
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", u::PDict{K,V}) where {K,V} = _pdict_show(io, u, "PDict")
Base.show(io::IO, ::MIME"text/plain", u::PIdDict{K,V}) where {K,V} = _pdict_show(io, u, "PIdDict")
Base.show(io::IO, ::MIME"text/plain", u::PEqualDict{K,V}) where {K,V} = _pdict_show(io, u, "PEqualDict")
Base.show(io::IO, ::MIME"text/plain", u::PEquivDict{K,V}) where {K,V} = _pdict_show(io, u, "PDict")

_iseq(t::PDict{K,V}, s::AbstractDict{KK,VV}, eqfn::Function) where {K,V,KK,VV} = begin
    (length(t) == length(s)) || return false
    # we want the general types to match...
    (equalfn(typeof(t)) === equalfn(typeof(s))) || return false
    for x in s
        in(x, t, eqfn) || return false
    end
    return true
end
Base.isequal(t::PDict{K,V}, s::AbstractDict{KK,VV}) where {K,V,KK,VV} = _iseq(t, s, isequal)
Base.isequal(s::AbstractDict{KK,VV}, t::PDict{K,V}) where {K,V,KK,VV} = _iseq(t, s, isequal)
isequiv(t::PDict{K,V}, s::AbstractDict{KK,VV}) where {K,V,KK,VV} = _iseq(t, s, isequiv)
isequiv(s::AbstractDict{KK,VV}, t::PDict{K,V}) where {K,V,KK,VV} = _iseq(t, s, isequiv)

# PVectorDict objects don't attenpt to differentiate objects; they just act in
# a naive linear fashion; this makes them good if, for example, you need to
# keep track of objects that have already had a hash collision.
struct PIdVectorDict{K,V} <: PIdDict{K,V}
    _items::PVec{Pair{K,V}}
    # When given a PVec or IdSet, PVectorSet-types don't have to check for correctness!
    PIdVectorDict{K,V}(u::IdDict{KK,VV}) where {K,V,KK<:K,VV<:V} = new{K,V}(PVec{Pair{K,V}}(Pair{K,V}[u...]))
    PIdVectorDict{K,V}(u::PIdDict{KK,VV}) where {K,V,KK<:K,VV<:V} = new{K,V}(PVec{Pair{K,V}}(Pair{K,V}[u...]))
    # Otherwise, they do:
    PIdVectorDict{T}(u::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === (===)
                return new{K,V}(PVec{Pair{K,V}}(Pair{K,V}[u...]))
            else
                return new{K,V}(PVec{Pair{K,V}}([IdDict{Pair{K,V}}(u)...]))
            end
        end
    end
    PIdVectorDict{K,V}(u::AbstractArray) where {K,V} = begin
        new{K,V}(PVec{T}(T[IdDict{K,V}(Pair{K,V}[u...])...]))
    end
    # However, it is too expensive to check every possible iteration of a dict for
    # correctness (e.g. as we build up or edit a dict), so we do not check when
    # we are passed a PVec that is of the appropriate type:
    PIdVectorDict{K,V}(u::PVec{Pair{K,V}}) where {T} = new{K,V}(u)
end
PIdVectorDict{K,V}(u::PIdVectorDict{K,V}) where {K,V} = u
struct PEqualVectorDict{K,V} <: PEqualDict{K,V}
    _items::PVec{Pair{K,V}}
    PEqualVectorDict{K,V}(u::Dict{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
        return new{T}(PVec{K,V}(Pair{K,V}[u...]))
    end
    PEqualVectorDict{K,V}(u::PEqualDict{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
        return new{T}(PVec{K,V}(Pair{K,V}[u...]))
    end
    PEqualVectorDict{K,V}(u::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === isequal
                return new{K,V}(PVec{Pair{K,V}}(Pair{K,V}[u...]))
            else
                return new{K,V}(PVec{Pair{K,V}}([Dict{K,V}(u)...]))
            end
        end
    end
    PEqualVectorDict{K,V}(u::AbstractArray) where {K,V} = begin
        return new{K,V}(PVec{Pair{K,V}}(T[Dict{K,V}(u)...]))
    end
    PEqualVectorDict{K,V}(u::PVec{Pair{K,V}}) where {K,V} = new{K,V}(u)
end
PEqualVectorDict{K,V}(u::PEqualVectorDict{K,V}) where {K,V} = u
struct PVectorDict{K,V} <: PEquivDict{K,V}
    _items::PVec{Pair{K,V}}
    PVectorDict{K,V}(u::PVec{Pair{K,V}}) where {K,V} = new{K,V}(u)
    PVectorDict{K,V}(u::EquivDict{KK,VV}) where {K,V,KK,VV} = new{K,V}(PVec{K,V}(Pair{K,V}[u...]))
    PVectorDict{K,V}(u::PDict{KK,VV}) where {K,V,KK,VV} = new{K,V}(PVec{K,V}(Pair{K,V}[u...]))
    PVectorDict{K,V}(u::AbstractDict{KK,VV}) where {K,V,KK,VV} = begin
        let U = typeof(u), eqfn = equalfn(U), hfn = hashfn(U)
            if eqfn === isequiv
                return new{K,V}(PVec{K,V}(Pair{K,V}[u...]))
            else
                return new{K,V}(PVec{K,V}([EquivDict{K,V}(u)...]))
            end
        end
    end
    PVectorDict{K,V}(u::AbstractArray) where {K,V} = new{K,V}(PVec{K,V}(Pair{K,V}[EquivDict{K,V}(u)...]))
end
const PVECDICTS = (Union{PIdVectorDict{K,V}, PEqualVectorDict{K,V}, PVectorDict{K,V}} where {K,V})
Base.convert(::Type{PVec{Pair{K,V}}}, u::PVectorDict{K,V}) where {K,V} = u._items
Base.convert(::Type{PVec{Pair{K,V}}}, u::PIdVectorDict{K,V}) where {K,V} = u._items
Base.convert(::Type{PVec{Pair{K,V}}}, u::PEqualVectorDict{K,V}) where {K,V} = u._items
Base.length(u::PVectorDict) = length(u._items)
Base.length(u::PIdVectorDict) = length(u._items)
Base.length(u::PEqualVectorDict) = length(u._items)
Base.iterate(u::PVectorDict) = iterate(u._items)
Base.iterate(u::PIdVectorDict) = iterate(u._items)
Base.iterate(u::PEqualVectorDict) = iterate(u._items)
Base.iterate(u::PVectorDict, x) = iterate(u._items, x)
Base.iterate(u::PIdVectorDict, x) = iterate(u._items, x)
Base.iterate(u::PEqualVectorDict, x) = iterate(u._items, x)
Base.in(x::Pair{KK,VV}, u::PVECDICTS{K,V}, eqfn::Function) where {K,V,KK<:K,VV<:V} = begin
    for y in convert(PVec{Pair{K,V}}, u)
        eqfn(y, x) && return true
    end
    return false
end
Base.in(x::Pair{KK,VV}, u::PVECDICTS{K,V}) where {K,V,KK<:K,VV<:V} = in(x, u, equalfn(u))
push(u::PVECDICTS{K,V}, x::Pair{KK,VV}) where {K,V,KK<:K,VV<:V} = let U = typeof(u)
    (x in u) && return u
    return U(push(convert(PVec{Pair{K,V}}, u), x))
end
delete(u::PVECDICTS{K,V}, k::KK) where {K,V,KK<:K} = begin
    let eq = equalfn(u), v = convert(PVec{Pair{K,V}}, u), U = typeof(u)
        for (ii,ky,vy) in enumerate(v)
            eq(ky, k) && return U(PVec{Pair{K,V}}(vcat(v[1:ii-1], v[ii+1:end])))
        end
        return u
    end
end

# PHashDict is a a dict that uses persistent maps to keep track of keys in and
# not in the dict. The expectation is that the dicts that make up the map values
# are PVectordict type sets.
struct PHashDict{K,V} <: PEquivDict{K,V}
    _n::Integer
    _table::PMap{Pair{K,PSet{V}}}
end
struct PIdHashDict{K,V} <: PIdDict{K,V}
    _n::Integer
    _table::PMap{Pair{K,PSet{V}}}
end
struct PEqualHashDict{K,V} <: PEqualDict{K,V}
    _n::Integer
    _table::PMap{Pair{K,PSet{V}}}
end
_pvecdict_type(::Type{PHashDict}) = PVectorDict
_pvecdict_type(::Type{PIdHashDict}) = PIdVectorDict
_pvecdict_type(::Type{PEqualHashDict}) = PEqualVectorDict
_pvecdict_type(::Type{PHashDict{K,V}}) where {K,V} = PVectorDict{K,V}
_pvecdict_type(::Type{PIdHashDict{K,V}}) where {K,V} = PIdVectorDict{K,V}
_pvecdict_type(::Type{PEqualHashDict{K,V}}) where {K,V} = PEqualVectorDict{K,V}
_pvecdict_type(::Type{PVectorDict}) = PVectorDict
_pvecdict_type(::Type{PIdVectorDict}) = PIdVectorDict
_pvecdict_type(::Type{PEqualVectorDict}) = PEqualVectorDict
_pvecdict_type(::Type{PVectorDict{K,V}}) where {K,V} = PVectorDict{K,V}
_pvecdict_type(::Type{PIdVectorDict{K,V}}) where {K,V} = PIdVectorDict{K,V}
_pvecdict_type(::Type{PEqualVectorDict{K,V}}) where {K,V} = PEqualVectorDict{K,V}
_pvecdict_type(::Type{EquivDict}) = PVectorDict
_pvecdict_type(::Type{IdDict}) = PIdVectorDict
_pvecdict_type(::Type{Dict}) = PEqualVectorDict
_pvecdict_type(::Type{EquivDict{K,V}}) where {K,V} = PVectorDict{K,V}
_pvecdict_type(::Type{IdDict{K,V}}) where {K,V} = PIdVectorDict{K,V}
_pvecdict_type(::Type{Dict{K,V}}) where {K,V} = PEqualVectorDict{K,V}
_pvecdict_type(::Type{AbstractDict}) = PEqualVectorDict
_pvecdict_type(::Type{AbstractDict{K,V}}) where {K,V} = PEqualVectorDict{K,V}
_phashdict_type(::Type{PVectorDict}) = PHashDict
_phashdict_type(::Type{PIdVectorDict}) = PIdHashDict
_phashdict_type(::Type{PEqualVectorDict}) = PEqualHashDict
_phashdict_type(::Type{PVectorDict{K,V}}) where {K,V} = PHashDict{K,V}
_phashdict_type(::Type{PIdVectorDict{K,V}}) where {K,V} = PIdHashDict{K,V}
_phashdict_type(::Type{PEqualVectorDict{K,V}}) where {K,V} = PEqualHashDict{K,V}
_phashdict_type(::Type{PVectorDict}) = PHashDict
_phashdict_type(::Type{PIdVectorDict}) = PIdHashDict
_phashdict_type(::Type{PEqualVectorDict}) = PEqualHashDict
_phashdict_type(::Type{PVectorDict{K,V}}) where {K,V} = PHashDict{K,V}
_phashdict_type(::Type{PIdVectorDict{K,V}}) where {K,V} = PIdHashDict{K,V}
_phashdict_type(::Type{PEqualVectorDict{K,V}}) where {K,V} = PEqualHashDict{K,V}
_phashdict_type(::Type{PHashDict}) = PHashDict
_phashdict_type(::Type{PIdHashDict}) = PIdHashDict
_phashdict_type(::Type{PEqualHashDict}) = PEqualHashDict
_phashdict_type(::Type{PHashDict{K,V}}) where {K,V} = PHashDict{K,V}
_phashdict_type(::Type{PIdHashDict{K,V}}) where {K,V} = PIdHashDict{K,V}
_phashdict_type(::Type{PEqualHashDict{K,V}}) where {K,V} = PEqualHashDict{K,V}
_phashdict_type(::Type{EquivDict}) = PVectorDict
_phashdict_type(::Type{IdDict}) = PIdHashDict
_phashdict_type(::Type{Dict}) = PEqualHashDict
_phashdict_type(::Type{EquivDict{K,V}}) where {K,V} = PVectorDict{K,V}
_phashdict_type(::Type{IdDict{K,V}}) where {K,V} = PIdHashDict{K,V}
_phashdict_type(::Type{Dict{K,V}}) where {K,V} = PEqualHashDict{K,V}
_phashdict_type(::Type{AbstractDict}) = PEqualHashDict
_phashdict_type(::Type{AbstractDict{K,V}}) where {K,V} = PEqualHashDict{K,V}
# For hash-type sets:
_pdict_init(P::DataType, arr::AbstractArray)= begin
    let P = (isa(P, UnionAll) ? P{K,V} : P), K = P.parameters[1], V = P.parameters[2], KK,VV
        D = _pvecdict_type(P), u = TMap{PDict{K,V}}(), eqfn = equalfn(P), hfn = hashfn(P), h, x,
        for (k,v) in arr
            h = hfn(k)
            x = get(u, h, nothing)
            if (x === nothing)
                u[h] = D([Pair{K,V}(k,v)])
            else
                u[h] = push(x, Pair{K,V}(k,v))
            end
        end
        return P(length(u), freeze(u))
    end
end
_pdict_init(P::DataType, s::AbstractDict{KK,VV}) where {KK,VV} = _pdict_init(P, Array(s))
# PHash type methods:
const PHASHDICTS = (Union{PHashDict{K,V}, PIdHashDict{K,V}, PEqualHashDict{K,V}} where {K,V})
Base.convert(::Type{PMap{PDict{K,V}}}, u::PHashDict{K,V}) where {K,V} = u._table
Base.convert(::Type{PMap{PDict{K,V}}}, u::PIdHashDict{K,V}) where {K,V} = u._table
Base.convert(::Type{PMap{PDict{K,V}}}, u::PEqualHashDict{K,V}) where {K,V} = u._table
Base.length(u::PHashDict)       = u._n
Base.length(u::PIdHashDict)     = u._n
Base.length(u::PEqualHashDict)  = u._n
Base.iterate(u::PHASHDICTS{K,V}) where {K,V} = begin
    let t = convert(PMap{PDict{K,V}}, u), x = iterate(t)
        if x === nothing
            return x
        else
            let ((kk,subd),addr1) = x, (kv,addr2) = iterate(subd)
                return (kv, (t, subd, addr1, addr2))
            end
        end
    end
end
Base.iterate(u::PHASHDICTS{K,V}, addr::Tuple) where {K,V} = beginlet (t,subd,ad1,ad2) = addr, x
    x = iterate(subd, ad2)
    if x === nothing
        x = iterate(t, ad1)
        (x === nothing) && return nothing
        let ((kk,subd),addr1) = x, (kv,addr2) = iterate(subd)
            return (kv, (t, subd, addr1, addr2))
        end
    else
        let (kv,addr2) = x
            return (kv, (t, k, vs, ad1, addr2))
        end
    end
end
Base.in(x::Pair{KK,VV}, u::PHASHDICTS{K,V}, eqfn::Function) where {K,V,KK<:K,VV<:V} = begin
    let h = UInt(hashfn(u)(x)), t = convert(PMap{PDict{K,V}}, u), ys = get(t, h, nothing)
        (ys === nothing) && return false
        return in(x, ys, eqfn)
    end
end
Base.in(x::S, u::PHASHDICTS{T}) where {T, S<:T} = in(x, u, equalfn(u))
push(u::PHASHDICTS{K,V}, p::Pair{KK,VV}) where {K,V,KK<:K,VV<:V} = begin
    let U = typeof(u), t = convert(PMap{PDict{K,V}}, u), eqfn = equalfn(u),
        (k,v) = x, h = UInt(hashfn(u)(k)), vv = get(t, h, nothing)
        (vv === nothing) && return U(u._n + 1, assoc(t, h, _pvecdict_type(U)([p])))
        if haskey(vv, k)
            (get(vv, k, nothing) === v) && return u
            return U(u._n, assoc(t, h, assoc(vv, k, v)))
        else
            let vvv = push(vv, p)
                (vvv === v) && return u
                return U(u._n + 1, assoc(t, h, vvv))
            end
        end
    end
end
delete(u::PHASHDICTS{K,V}, k::KK) where {K,V,KK<:K} = begin
    let U = typeof(u), t = convert(PMap{PDict{K,V}}, u), eqfn = equalfn(u),
        h = UInt(hashfn(u)(x)), v = get(t, h, nothing)
        (v === nothing || !haskey(v, k)) && return u
        return U(u._n - 1,
                 length(v) == 1
                 ? dissoc(t, h)
                 : assoc(t, h, delete(v, k)))
    end
end
# Constructors:
PHashDict{K,V}() where {K,V} = PHashDict{K,V}(0, PMap{Pair{K,PSet{V}}}())
PHashDict{K,V}(s::PHashDict{K,V}) where {K,V} = s
PHashDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = _pset_init(PHashDict{K,V}, s)
PHashDict{K,V}(arr::AbstractArray{S, 1}) where {K,V,KK<:K,VV<:V} = _pset_init(PHashDict{K,V}, arr)
PIdHashDict{K,V}() where {K,V} = PIdHashDict{K,V}(0, PMap{Pair{K,PSet{V}}}())
PIdHashDict{K,V}(s::PIdHashDict{K,V}) where {K,V} = s
PIdHashDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = _pset_init(PIdHashDict{K,V}, s)
PIdHashDict{K,V}(arr::AbstractArray{S, 1}) where {K,V,KK<:K,VV<:V} = _pset_init(PIdHashDict{K,V}, arr)
PEqualHashDict{K,V}() where {K,V} = PEqualHashSet(0, PMap{Pair{K,PSet{V}}}())
PEqualHashDict{K,V}(s::PEqualHashDict{K,V}) where {K,V} = s
PEqualHashDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = _pset_init(PEqualHashDict{K,V}, s)
PEqualHashDict{K,V}(arr::AbstractArray{S, 1}) where {K,V,KK<:K,VV<:V} = _pset_init(PEqualHashDict{K,V}, arr)

PDict{K,V}() where {K,V} = PHashDict{K,V}()
PDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = PHashDict{K,V}(s)
PDict{K,V}(s::AbstractArray{S,1}) where {K,V,KK<:K,VV<:V} = PHashDict{K,V}(s)
PIdDict{K,V}() where {K,V} = PIdHashDict{K,V}()
PIdDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = PIdHashDict{K,V}(s)
PIdDict{K,V}(s::AbstractArray{S,1}) where {K,V,KK<:K,VV<:V} = PIdHashDict{K,V}(s)
PEqualDict{K,V}() where {K,V} = PEqualHashDict{K,V}()
PEqualDict{K,V}(s::AbstractDict{KK,VV}) where {K,V,KK<:K,VV<:V} = PEqualHashDict{K,V}(s)
PEqualDict{K,V}(s::AbstractArray{S,1}) where {K,V,KK<:K,VV<:V} = PEqualHashDict{K,V}(s)
_pdict_construct(PS::Type, s::AbstractArray{T,1}) where {T} = begin
    let K = typejoin([typeof(k) for (k,v) in s]...), V = typejoin([typeof(v) for (k,v) in s]...)
        return PS{K,V}(s)
    end
end
_pdict_construct(PS::Type, s::AbstractDict{K,V}) where {K,V} = begin
    return PS{K,V}(Pair{K,V}[s...])
end
PVectorDict(s::AbstractArray{T,1})      where {T} = _pdict_construct(PVectorDict, s)
PIdVectorDict(s::AbstractArray{T,1})    where {T} = _pdict_construct(PIdVectorDict, s)
PEqualVectorDict(s::AbstractArray{T,1}) where {T} = _pdict_construct(PEqualVectorDict, s)
PHashDict(s::AbstractArray{T,1})      where {T} = _pdict_construct(PHashDict, s)
PIdHashDict(s::AbstractArray{T,1})    where {T} = _pdict_construct(PIdHashDict, s)
PEqualHashDict(s::AbstractArray{T,1}) where {T} = _pdict_construct(PEqualHashDict, s)
PDict(s::AbstractArray{T,1})      where {T} = _pdict_construct(PHashDict, s)
PIdDict(s::AbstractArray{T,1})    where {T} = _pdict_construct(PIdHashDict, s)
PEqualDict(s::AbstractArray{T,1}) where {T} = _pdict_construct(PEqualHashDict, s)

PVectorDict(s::AbstractDict{K,V})      where {K,V} = _pdict_construct(PVectorDict, s)
PIdVectorDict(s::AbstractDict{K,V})    where {K,V} = _pdict_construct(PIdVectorDict, s)
PEqualVectorDict(s::AbstractDict{K,V}) where {K,V} = _pdict_construct(PEqualVectorDict, s)
PHashDict(s::AbstractDict{K,V})      where {K,V} = _pdict_construct(PHashDict, s)
PIdHashDict(s::AbstractDict{K,V})    where {K,V} = _pdict_construct(PIdHashDict, s)
PEqualHashDict(s::AbstractDict{K,V}) where {K,V} = _pdict_construct(PEqualHashDict, s)
PDict(s::AbstractDict{K,V})      where {K,V} = _pdict_construct(PHashDict, s)
PIdDict(s::AbstractDict{K,V})    where {K,V} = _pdict_construct(PIdHashDict, s)
PEqualDict(s::AbstractDict{K,V}) where {K,V} = _pdict_construct(PEqualHashDict, s)
