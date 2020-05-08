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
equalfn(::Type{EquivSet}) = isequiv
equalfn(::Type{EquivSet{T}}) where {T} = isequiv
hashfn(::Type{EquivSet}) = equivhash
hashfn(::Type{EquivSet{T}}) where {T} = equivhash

# We should define an isequiv function for abstract sets; this should work
# equally well on PSets or normal sets.
isequiv(s::SS, t::ST) where {S, T, SS <: AbstractSet{S}, ST <: AbstractSet{T}} = begin
    (ismuttype(ST) || ismuttype(ST)) && return (s === t)
    (length(t) == length(s)) || return false
    for ss in s
        tt = get(t, ss, t)
        (tt === t) && return false
        isequiv(ss, tt) || return false
    end
    return true
end


################################################################################
# PSet
# We construct PSets in an unfortunaetly complex way in order to reduce code
# duplication. The three kinds of persistent sets are very similar, just using
# different equality tests, so we can keep most of the code the same.
macro _pset_code(name::Symbol, eqfn, hashfn)
    let eq = gensym(), h = gensym(), _name = gensym(), q
        q = quote
            let $eq = $eqfn, $h = $hashfn
                struct $name{T} <: AbstractSet{T}
                    _n::Int
                    _tree::PTree{PList{T}}
                end
                # Document the equal/hash types.
                equalfn(::Type{T}) where {T <: $name} = $eq
                hashfn(::Type{T}) where {T <: $name} = $h
                equalfn(::$name) = $eq
                hashfn(::$name) = $h
                # Construction.
                $name{T}() where {T} = $name{T}(0, PTree{PList{T}}())
                $name{T}(itr) where {T} = reduce(push, itr, init=$name{T}())
                $name() = $name{Any}()
                $name(itr) = $_name(itr, Base.IteratorEltype(itr))
                $_name(itr, ::Base.HasEltype) = $name{eltype(itr)}(itr)
                $_name(itr, ::Base.EltypeUnknown) = begin
                    T = typejoin(map(typeof, itr)...)
                    return $name{T}(itr)
                end
                # Base methods.
                Base.empty(s::$name{T}, ::Type{U}=T) where {T,U} = $name{U}()
                Base.length(s::$name) = s._n
                Base.iterate(u::$name{T}) where {T} = begin
                    (u._n == 0) && return nothing
                    (lst,titer) = iterate(u._tree)
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    return (el,(lst,liter,titer))
                end
                Base.iterate(u::$name{T}, tup) where {T} = begin
                    (lst,liter,titer) = tup
                    if liter !== nothing
                        (el,liter) = iterate(lst, liter)
                        return (el,(lst,liter,titer))
                    end
                    q = iterate(u._tree, titer)
                    (q === nothing) && return q
                    (lst,titer) = q
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    return (el,(lst,liter,titer))
                end
                Base.in(x::S, u::$name{T}, eqfn::Function) where {T, S} = begin
                    hh = $h(x)
                    uu = get(u._tree, hh, PList{T}())
                    return in(x, uu, eqfn)
                end
                Base.in(x::S, u::$name{T}) where {S,T} = in(x, u, $eq)
                # And the Air methods.
                Air.push(u::$name{T}, x::S) where {T,S} = begin
                    hh = $h(x)
                    uu = get(u._tree, hh, PList{T}())
                    in(x, u, $eq) && return u
                    uu = pushfirst(uu, x)
                    tr = setindex(u._tree, uu, hh)
                    return $name{T}(u._n + 1, tr)
                end
                Air.delete(u::$name{T}, x::S) where {T,S} = begin
                    hh = $h(x)
                    uu = get(u._tree, hh, PList{T}())
                    (length(uu) === 0) && return u
                    vv = delete(uu, x, $eq)
                    (uu === vv) && return u
                    if length(uu) == 0
                        return $name{T}(u._n - 1, delete(u._tree, hh))
                    else
                        return $name{T}(u._n - 1, setindex(u._tree, hh, vv))
                    end
                end
            end
        end
        return esc(q)
    end
end
    
@_pset_code PSet isequiv equivhash
#@_pset_code PIdSet (===) objectid
#@_pset_code PEqualSet isequal hash
