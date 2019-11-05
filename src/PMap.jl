################################################################################
# PMap
# A persistent map type that maps integers to values; it works for hash's,
# objectid's, or for making sparse vectors.

const PMAP_KEY_T = typejoin(typeof(hash(nothing)), typeof(objectid(nothing)))
_keysplit(k::PMAP_KEY_T) = (Int(k & 0b11111) + 1, k >> 5)
abstract type PMap{T} <: AbstractDict{PMAP_KEY_T, T} end

# make sure pmaps have a hash and equality operator that is appropriate
_iseq(t::PMap{T}, s::AbstractDict{K,S}, eqfn::Function) where {T,K<:Integer,S} = begin
    (length(t) == length(s)) || return false
    for kv in s
        v = get(t, kv, t)
        (v === t) && return false
        eqfn(v, kv[2]) || return false
    end
    return true
end
isequiv(t::PMap{T}, s::AbstractDict{K,S}) where {T,K<:Integer,S} = _iseq(t, s, isequiv)
isequiv(s::AbstractDict{K,S}, t::PMap{T}) where {T,K<:Integer,S} = _iseq(t, s, isequiv)
Base.isequal(t::PMap{T}, s::AbstractDict{K,S}) where {T,K<:Integer,S} = _iseq(t, s, isequal)
Base.isequal(s::AbstractDict{K,S}, t::PMap{T}) where {T,K<:Integer,S} = _iseq(t, s, isequal)
equivhash(t::PMap{T}) where {T} = let h = ~hash(T)
    for (k,v) in t
        h += hasheq(v) * k
    end
    return h
end
# We need a type for PMap entries--they can be empty entries (not yet assigned),
# singleton entries (a key and a value), or an antire sub-map:
abstract type PMapNode{T} end
struct PMapEmptyNode{T} <: PMapNode{T} end
struct PMapSubmapNode{T} <: PMapNode{T}
    _map::PMap{T}
end
struct PMapSingleNode{T} <: PMapNode{T}
    _key::PMAP_KEY_T
    _val::T
end
PMapSingleNode{T}(p::Pair{PMAP_KEY_T,S}) where {T,S<:T} = PMapSingleNode{T}(p[1], p[2])

# How we print pmaps:
_print_pmap(io::IO, u::PMap{T}, head::String) where {T} = begin
    print(io, "$(head)[")
    n = length(u)
    if n < 50
        for (ii,kv) in enumerate(u)
            if   ii > 1 print(io, ", $(repr(kv))")
            else        print(io, "$(repr(kv))")
            end
        end
    else
        for (ii,kv) in enumerate(u)
            if     ii > 50 break
            elseif ii > 1  print(io, ", $(repr(kv))")
            else           print(io, "$(repr(kv))")
            end
        end
        print(io, " ...")
    end
    print(io, "]")
end
Base.show(io::IO, ::MIME"text/plain", pv::PMap{T}) where {T} =
    _print_pmap(io, pv, "PMap")

struct PMap0{T} <: PMap{T} end
Base.length(::PMap0) = 0
Base.get(::PMap0, k::PMAP_KEY_T, df) = df
Base.iterate(::PMap0) = nothing
Base.iterate(::PMap0, x) = nothing
dissoc(m::PMap0, k::PMAP_KEY_T) = m
assoc(m::PMap0{T}, k::PMAP_KEY_T, v::S) where {T, S<:T} = begin
    let ee = PMapEmptyNode{T}(), ar = PMapNode{T}[ee for x in 1:32], (k1,k2) = _keysplit(k)
        ar[k1] = PMapSingleNode{T}(k, v)
        return PMap32{T}(1, PVec32{PMapNode{T}}(ar...))
    end
end

struct PMap32{T} <: PMap{T}
    _n::Integer
    _data::PVec32{PMapNode{T}}
end
Base.length(m::PMap32) = m._n
Base.get(m::PMap32{T}, k::PMAP_KEY_T, df) where {T} = begin
    let (k1,k2) = _keysplit(k), u = m._data[k1]
        if isa(u, PMapEmptyNode{T})
            return df
        elseif isa(u, PMapSubmapNode{T})
            return get(u._map, k2, df)
        elseif k != u._key
            return df
        else
            return u._val
        end
    end
end
Base.in(m::PMap32{T}, kv::Pair{PMAP_KEY_T,T}) where {T} = begin
    let k = kv[1], v = kv[2], u = get(m, k, m)
        (u === m) && return false
        return u == v
    end
end
Base.iterate(m::PMap32{T}) where {T} = begin
    for (k,v) in enumerate(m._data)
        if isa(v, PMapEmptyNode{T})
            continue
        elseif isa(v, PMapSubmapNode{T})
            let (x,st) = iterate(v._map), kk = (x[1] << 5) | (k - 1)
                return (Pair{PMAP_KEY_T,T}(kk, x[2]), (k, st...))
            end
        else
            return (Pair{PMAP_KEY_T,T}(v._key, v._val), (k+1,))
        end
    end
end
Base.iterate(m::PMap32{T}, st::Tuple) where {T} = begin
    if length(st) == 0
        return iterate(m)
    else
        let k1 = st[1]
            for k in k1:32
                let u = m._data[k]
                    if isa(u, PMapEmptyNode{T})
                        continue
                    elseif isa(u, PMapSubmapNode{T})
                        let it = (k == k1 ? iterate(u._map, st[2:end]) : iterate(u._map))
                            if it === nothing
                                continue
                            else
                                let kv = it[1], kk = (kv[1] << 5) | (k - 1)
                                    return (Pair{PMAP_KEY_T,T}(kk, kv[2]), (k, it[2]...))
                                end
                            end
                        end
                    else
                        return (Pair{PMAP_KEY_T,T}(u._key, u._val), (k+1,))
                    end
                end
            end
        end
    end
    return nothing
end
assoc(m::PMap32{T}, k::PMAP_KEY_T, v::S) where {T, S<:T} = begin
    let (k1,k2) = _keysplit(k), u = m._data[k1], node, n = m._n
        if isa(u, PMapEmptyNode{T})
            n += 1
            node = PMapSingleNode{T}(k, v)
        elseif isa(u, PMapSubmapNode{T})
            let nu = length(u._map), uu = assoc(u._map, k2, v)
                (uu === u._map) && return m
                n += length(uu) - nu
                node = PMapSubmapNode{T}(uu)
            end
        elseif u._key != k
            n += 1
            node = assoc(assoc(PMap0{T}(), u._key >> 5, u._val), k2, v)
            node = PMapSubmapNode{T}(node)
        elseif u._val === v
            return m
        else
            node = PMapSingleNode{T}(k, v)
        end
        return PMap32{T}(n, assoc(m._data, k1, node))
    end
end
dissoc(m::PMap32{T}, k::PMAP_KEY_T) where {T} = begin
    let (k1,k2) = _keysplit(k), u = m._data[k1], node
        if isa(u, PMapEmptyNode{T})
            return m
        elseif isa(u, PMapSubmapNode{T})
            let uu = dissoc(u._map, k2)
                if uu === u._map
                    return m
                elseif length(uu) == 0
                    node = PMapEmptyNode{T}()
                elseif length(uu) == 1
                    let kv = first(uu), kk = (kv[1] << 5) | (k1 - 1)
                        node = PMapSingleNode{T}(kk, kv[2])
                    end
                else
                    node = PMapSubmapNode{T}(uu)
                end
            end
        elseif u._key == k
            node = PMapEmptyNode{T}()
        else
            return m
        end
        return PMap32{T}(m._n - 1, assoc(m._data, k1, node))
    end
end

# Constructors
PMap(kvs::Vararg{Union{Tuple{PMAP_KEY_T,S},Pair{PMAP_KEY_T,S}}}) where {S} = begin
    let m0 = PMap0{S}(), m = m0
        for kv in kvs
            m = assoc(m, kv[1], kv[2])
        end
        return m
    end
end
PMap{T}(kvs::Vararg{Union{Tuple{PMAP_KEY_T,S},Pair{PMAP_KEY_T,S}}}) where {T,S<:T} = begin
    let m0 = PMap0{T}(), m = m0
        for kv in kvs
            m = assoc(m, kv[1], kv[2])
        end
        return m
    end
end
