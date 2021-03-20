################################################################################
# pset.jl
#
# Implementation of a persistent set type that uses the PTree type internally.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

import Base.IdSet

# #PLinearSet ==================================================================
# We construct these objects in an unfortunaetly complex way in order to reduce
# code duplication. The two kinds of persistent sets are very similar, just
# using different equality tests, so we can keep most of the code the same.
# See PDict.jl for information about this macro-based construction, used here.
macro _plinset_code(name::Symbol, eq, h)
    return quote
        struct $name{T} <: AbstractPSet{T}
            elements::Union{Nothing, Vector{T}}
            function $name{T}() where {T}
                return new{T}(nothing)
            end
            function $name{T}(u::Vector{T}, ::Val{:safe}) where {T}
                return length(u) == 0 ? new{T}(nothing) : new{T}(u)
            end
            function $name{T}(u::Vector{T}, ::Val{:unique}) where {T}
                return length(u) == 0 ? new{T}(nothing) : new{T}(copy(u))
            end
            function $name{T}(u::Vector{T}) where {T}
                (length(u) == 0) && return new{T}(nothing)
                v = Vector{T}(undef, length(u))
                kk = 0
                for ii in 1:length(u)
                    ui = u[ii]
                    found = false
                    for jj in 1:k
                        if v[ii] == ui
                            found = true
                            break
                        end
                    end
                    found && continue
                    k += 1
                    v[k] = ui
                end
                resize!(v, k)
                return new{T}(v)
            end
        end
        function $name{T}(d::$name{T}) where {T}
            return d
        end
        function $name{T}(t0, ts...) where {T}
            return $name{T}(T[t0, ts...])
        end
        function $name{T}(itr) where {T}
            return $name{T}(T[itr...])
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = $eq
        hashfn(::Type{T}) where {T <: $name} = $h
        equalfn(::$name) = $eq
        hashfn(::$name) = $h
        # Generic construction.
        $name() = $name{Any}()
        $name(::Tuple{}) = $name{Any}()
        $name(itr) = $name(itr, Base.IteratorEltype(itr))
        $name(itr, ::Base.HasEltype) = $name{eltype(itr)}(itr)
        $name(itr, ::Base.EltypeUnknown) = begin
            T = typejoin(map(typeof, itr)...)
            return $name{T}(itr)
        end
        # Base methods.
        Base.empty(s::$name{T}, ::Type{S}=T) where {T,S} = $name{S}()
        Base.length(s::$name) = begin
            ks = getfield(s, :elements)
            return (ks === nothing ? 0 : length(ks))
        end
        Base.IteratorEltype(::Type{$name{T}}) where {T} = Base.HasEltype()
        Base.eltype(::Type{$name{T}}) where {T} = T
        Base.eltype(::$name{T}) where {T} = T
        Base.IteratorSize(::Type{$name{T}}) where {T} = Base.HasLength()
        Base.iterate(u::$name{T}) where {T} = begin
            els = getfield(u, :elements)
            (els === nothing) && return nothing
            return iterate(els)
        end
        Base.iterate(u::$name{T}, state) where {T} = begin
            els = getfield(u, :elements)
            return iterate(els, state)
        end
        Base.in(x::S, u::$name{T}) where {S,T} = begin
            els = getfield(u, :elements)
            (els === nothing) && return false
            els::Vector{T}
            for el in els
                $eq(x,el) && return true
            end
            return false
        end
        # And the Air methods.
        Air.push(u::$name{T}, x::S) where {T,S} = begin
            els = getfield(u, :elements)
            (els === nothing) && return $name{T}(x)
            for el in els
                $eq(el,x) && return u
            end
            return $name{T}(push(els, x), Val{:safe}())
        end
        Air.delete(u::$name{T}, x::S) where {T,S} = begin
            els = getfield(u, :elements)
            (els === nothing) && return u
            els::Vector{T}
            n = length(els)
            if n == 1
                for ii in 1:n
                    if $eq(els[ii], x)
                        return $name{T}()
                    end
                end
            else
                for ii in 1:n
                    if $eq(els[ii], x)
                        return $name{T}(delete(elss, ii), Val{:safe}())
                    end
                end
            end
            return u
        end
        Base.propertynames(u::$name) = (:count, :elements)
        Base.getproperty(u::$name{T}, s::Symbol) where {T} = begin
            if s == :count
                return length(u)
            elseif s == :elements
                els = getfield(u, :elements)
                return els === nothing ? () : tuple(els...)
            else
                throw(ArgumentError("no such property $s of type $(typeof(u))"))
            end
        end
    end |> esc
end
# Macro-write the code fo these two types.
@_plinset_code PLinearSet isequal hash
@_plinset_code PIdLinearSet (===) objectid
# Document the types.
@doc """
    PLinearSet{T}

A persistent set type that uses a simple array representation internally;
accordingly `PLinearSet` does not peform persistent operations such as
`setindex`, `push`, and `delete` efficiently. It is intended to be used by Air
internally as a set that stores items whose hash values are identical.

See also: [`PIdLinearSet`](@ref), [`PSet`](@ref).
""" PLinearSet

@doc """
    PIdLinearSet{T}

A persistent `IdSet` type that uses a simple array representation internally;
accordingly `PIdLinearSet` does not peform persistent operations such as
`setindex`, `push`, and `delete` efficiently. It is intended to be used by Air
internally as a set that stores items whose hash values are identical.

See also: [`PLinearSet`](@ref), [`PSet`](@ref).
""" PIdLinearSet

# #PSet ========================================================================
# The PSet type employs both PTree and PLinearSet to create an unordered set.
macro _pset_code(name::Symbol, eq, h, linset)
    return quote
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
    end |> esc
end
# Declare the types.    
@_pset_code PSet isequal hash PLinearSet
@_pset_code PIdSet (===) objectid PIdLinearSet
# The show function is best written in terms of an AbstractPSet.
Base.show(io::IO, s::AbstractPSet{T}) where {T} = begin
    S = typeof(s)
    print(io, get(io, :typeinfo, Any) == s ? String(nameof(S)) : S)
    print(io, "(")
    isempty(s) || Base.show_vector(io, s)
    print(io, ")")
end
# Document them.
@doc """
    PSet{T}

A persistent set type that uses an array hash-mapped trie structure to allow for
efficient updating of the set using the persistent operations defined by Air
such as `setindex`, `push`, and `delete`.

See also: `Set`, [`PIdSet`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> PSet()
PSet{Any}()

julia> :a in PSet([:a, :b, :c])
true

julia> :d in PSet([:a, :b, :c])
false

julia> :d in push(PSet(), :d)
true
```
""" PSet
@doc """
    PIdSet{T}

A persistent identity-based set type that uses an array hash-mapped trie
structure to allow for efficient updating of the set using the persistent
operations defined by Air such as `push`, and `delete`.

See also: `Set`, [`PIdSet`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> PIdSet()
PIdSet{Any}()

julia> :a in PIdSet([:a, :b, :c])
true

julia> :d in PIdSet([:a, :b, :c])
false

julia> :d in push(PIdSet(), :d)
true
```
""" PIdSet

# Export these types.
export PSet, PIdSet

