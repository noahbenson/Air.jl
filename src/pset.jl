################################################################################
# PSet
# A persistent set type that uses PVec and PMap types.

import Base.IdSet


# ==============================================================================
# PSet
# We construct PSets in an unfortunaetly complex way in order to reduce code
# duplication. The three kinds of persistent sets are very similar, just using
# different equality tests, so we can keep most of the code the same.
macro _pset_code(name::Symbol, eqfn, hashfn, linset)
    eq = gensym()
    h = gensym()
    return quote
        begin
            local $eq = $eqfn
            local $h = $hashfn
            struct $name{T} <: AbstractPSet{T}
                count::Int
                root::PTree{$linset{T}}
            end
            # Document the equal/hash types.
            equalfn(::Type{T}) where {T <: $name} = $eq
            hashfn(::Type{T}) where {T <: $name} = $h
            equalfn(::$name) = $eq
            hashfn(::$name) = $h
            # Construction.
            $name{T}() where {T} = $name{T}(0, PTree{$linset{T}}())
            $name{T}(itr) where {T} = reduce(push, itr, init=$name{T}())
            $name() = $name{Any}()
            $name(itr) = $name(itr, Base.IteratorEltype(itr))
            $name(itr, ::Base.HasEltype) = $name{eltype(itr)}(itr)
            $name(itr, ::Base.EltypeUnknown) = begin
                T = typejoin(map(typeof, itr)...)
                return $name{T}(itr)
            end
            # Base methods.
            Base.empty(s::$name{T}) where {T} = $name{T}()
            Base.empty(s::$name{T}, ::Type{S}) where {T,S} = $name{S}()
            Base.length(s::$name) = getfield(s, :count)
            Base.IteratorSize(::Type{$name{T}}) where {T} = Base.HasLength()
            Base.IteratorEltype(::Type{$name{T}}) where {T} = Base.HasEltype()
            Base.eltype(::Type{$name{T}}) where {T,N} = T
            Base.eltype(u::$name{T}) where {T,N} = T
            Base.iterate(u::$name{T}) where {T} = begin
                (getfield(u, :count) == 0) && return nothing
                (lst,titer) = iterate(getfield(u, :root))
                lst = lst[2]
                (el,liter) = iterate(lst)
                return (el,(lst,liter,titer))
            end
            Base.iterate(u::$name{T}, tup) where {T} = begin
                (lst,liter,titer) = tup
                if liter !== nothing
                    r = iterate(lst, liter)
                    if r !== nothing
                        (el,liter) = r
                        return (el,(lst,liter,titer))
                    end
                end
                q = iterate(getfield(u, :root), titer)
                (q === nothing) && return q
                (lst,titer) = q
                lst = lst[2]
                (el,liter) = iterate(lst)
                return (el,(lst,liter,titer))
            end
            Base.in(x::S, u::$name{T}) where {T,S} = begin
                hh = $h(x)
                uu = get(getfield(u, :root), hh, nothing)
                return uu === nothing ? false : in(x, uu::$linset{T})
            end
            # And the Air methods.
            Air.push(u::$name{T}, x::S) where {T,S} = begin
                hh = $h(x)
                root = getfield(u, :root)
                uu = get(root, hh, nothing)
                if uu === nothing
                    uu = $linset{T}(T[x], Val{:safe}())
                else
                    in(x, uu) && return u
                    uu = push(uu, x)
                end
                tr = setindex(root, uu, hh)
                return $name{T}(getfield(u, :count) + 1, tr)
            end
            Air.delete(u::$name{T}, x::S) where {T,S} = begin
                hh = $h(x)
                uu = get(getfield(u, :root), hh, nothing)
                (uu === nothing) && return u
                vv = delete(uu, x)
                (uu === vv) && return u
                if length(vv) == 0
                    return $name{T}(getfield(u, :count) - 1, delete(getfield(u, :root), hh))
                else
                    return $name{T}(getfield(u, :count) - 1, setindex(getfield(u, :root), vv, hh))
                end
            end
        end
    end |> esc
end
    
@_pset_code PSet isequal hash PLinearSet
@_pset_code PIdSet (===) objectid PIdLinearSet
@_pset_code PEquivSet isequiv equivhash PEquivLinearSet

mutability(::Type{PSet}) = Immutable()
mutability(::Type{PSet{T}}) where {T} = Immutable()
mutability(::Type{PIdSet}) = Immutable()
mutability(::Type{PIdSet{T}}) where {T} = Immutable()
mutability(::Type{PEquivSet}) = Immutable()
mutability(::Type{PEquivSet{T}}) where {T} = Immutable()

