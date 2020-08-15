################################################################################
# PSet
# A persistent set type that uses PVec and PMap types.

import Base.IdSet


# ==============================================================================
# AbstractPSet
"""
    AbstractPSet{T}

AbstractPSet is an abstract type extended by persistent set types such as PSet,
PIdSet, and PEqualSet, as well as the weighted persistent set types.
"""
abstract type AbstractPSet{T} <: AbstractSet{T} end

# ==============================================================================
# PSet
# We construct PSets in an unfortunaetly complex way in order to reduce code
# duplication. The three kinds of persistent sets are very similar, just using
# different equality tests, so we can keep most of the code the same.
macro _pset_code(name::Symbol, eqfn, hashfn)
    let eq = gensym(), h = gensym(), _name = gensym(), q,
        n = gensym("[private] length"), tree = gensym("[private] tree")
        return quote
            begin
                local $eq = $eqfn
                local $h = $hashfn
                struct $name{T} <: AbstractPSet{T}
                    $n::Int
                    $tree::PTree{PList{T}}
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
                $name(itr) = $name(itr, Base.IteratorEltype(itr))
                $name(itr, ::Base.HasEltype) = $name{eltype(itr)}(itr)
                $name(itr, ::Base.EltypeUnknown) = begin
                    T = typejoin(map(typeof, itr)...)
                    return $name{T}(itr)
                end
                # Base methods.
                Base.empty(s::$name{T}, ::Type{U}=T) where {T,U} = $name{U}()
                Base.length(s::$name) = s.$n
                Base.IteratorSize(::Type{$name{T}}) where {T} = Base.HasLength()
                Base.IteratorEltype(::Type{$name{T}}) where {T} = Base.HasEltype()
                Base.eltype(::Type{$name{T}}) where {T,N} = T
                Base.eltype(u::$name{T}) where {T,N} = T
                Base.iterate(u::$name{T}) where {T} = begin
                    (u.$n == 0) && return nothing
                    (lst,titer) = iterate(u.$tree)
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
                    q = iterate(u.$tree, titer)
                    (q === nothing) && return q
                    (lst,titer) = q
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    return (el,(lst,liter,titer))
                end
                Base.in(x::S, u::$name{T}, eqfn::Function) where {T, S} = begin
                    hh = $h(x)
                    uu = get(u.$tree, hh, PList{T}())::PList{T}
                    return in(x, uu, eqfn)
                end
                Base.in(x::S, u::$name{T}) where {S,T} = in(x, u, $eq)
                # And the Air methods.
                Air.push(u::$name{T}, x::S) where {T,S} = begin
                    hh = $h(x)
                    uu = get(u.$tree, hh, PList{T}())::PList{T}
                    in(x, u, $eq) && return u
                    uu = pushfirst(uu, x)
                    tr = setindex(u.$tree, uu, hh)
                    return $name{T}(u.$n + 1, tr)
                end
                Air.delete(u::$name{T}, x::S) where {T,S} = begin
                    hh = $h(x)
                    uu = get(u.$tree, hh, PList{T}())::PList{T}
                    (length(uu) === 0) && return u
                    vv = delete(uu, x, $eq)
                    (uu === vv) && return u
                    if length(vv) == 0
                        return $name{T}(u.$n - 1, delete(u.$tree, hh))
                    else
                        return $name{T}(u.$n - 1, setindex(u.$tree, vv, hh))
                    end
                end
            end
        end |> esc
    end
end
    
@_pset_code PSet isequiv equivhash
@_pset_code PIdSet (===) objectid
@_pset_code PEqualSet isequal hash

mutability(::Type{PSet}) = Immutable()
mutability(::Type{PSet{T}}) where {T} = Immutable()
mutability(::Type{PIdSet}) = Immutable()
mutability(::Type{PIdSet{T}}) where {T} = Immutable()
mutability(::Type{PEqualSet}) = Immutable()
mutability(::Type{PEqualSet{T}}) where {T} = Immutable()

# We should define an isequiv function for abstract sets; this should work
# equally well on PSets or normal sets.
isequiv(s::SS, t::ST) where {S,T,SS <: AbstractPSet{S},ST <: AbstractPSet{T}} = begin
    (length(t) == length(s)) || return false
    for k in s
        (k in t) || return false
    end
    (equalfn(SS) === equalfn(ST)) && return true
    for k in t
        (k in s) || return false
    end
    return true
end
equivhash(u::SS) where {K,SS<:AbstractPSet{K}} = begin
    h = hash(length(u))
    for k in u
        h += equivhash(k)
    end
    return h
end
