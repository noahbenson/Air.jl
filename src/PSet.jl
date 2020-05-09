################################################################################
# PSet
# A persistent set type that uses PVec and PMap types.

import Base.IdSet

# We should define an isequiv function for abstract sets; this should work
# equally well on PSets or normal sets.
isequiv(s::SS, t::ST) where {S, T, SS <: AbstractSet{S}, ST <: AbstractSet{T}} = begin
    (ismuttype(ST) || ismuttype(ST)) && return (s === t)
    (length(t) == length(s)) || return false
    for ss in s
        (ss in t) || return false
    end
    # If s and t use the same notion of equivalence, this is enough
    (equalfn(SS) === equalfn(ST)) && return true
    # Otherwise, we need to check the other direction also
    for tt in t
        (tt in s) || return false
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
                    if length(vv) == 0
                        return $name{T}(u._n - 1, delete(u._tree, hh))
                    else
                        return $name{T}(u._n - 1, setindex(u._tree, vv, hh))
                    end
                end
            end
        end
        return esc(q)
    end
end
    
@_pset_code PSet isequiv equivhash
@_pset_code PIdSet (===) objectid
@_pset_code PEqualSet isequal hash
