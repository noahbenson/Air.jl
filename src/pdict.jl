################################################################################
# PDict
#
# Implementation of a persistent dictionary type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson


# #PLinearDict #################################################################
# We construct the types in this file in an unfortunaetly complex way in order
# to reduce code duplication. The three kinds of persistent sets are very
# similar, just using different equality tests, so we can keep most of the code
# the same. This method is used for PDict, PLinearDict, PSet, and PLinearSet.
macro _plindict_code(name::Symbol, eq, h)
    return quote
        struct $name{K,V} <: AbstractPDict{K,V}
            keys::Union{Nothing, Vector{K}}
            values::Union{Nothing, Vector{V}}
        end
        function $name{K,V}() where {K,V}
            return $name{K,V}(nothing, nothing)
        end
        function $name{K,V}(kv::Pair{J,U}) where {K,V,J,U}
            return $name{K,V}(K[kv[1]], V[kv[2]])
        end
        function $name{K,V}(kv::Tuple{J,U}) where {K,V,J,U}
            return $name{K,V}(K[kv[1]], V[kv[2]])
        end
        function $name{K,V}(d::$name{K,V}) where {K,V}
            return d
        end
        function $name{K,V}(d::AbstractDict) where {K,V}
            return reduce(push, d, init=$name{K,V}())
        end
        function $name{K,V}(kv::Pair, ps::Pair...) where {K,V}
            return reduce(push, ps, init=$name{K,V}(kv))
        end
        function $name{K,V}(kv::Tuple, ps::Tuple...) where {K,V}
            return reduce(push, ps, init=$name{K,V}(kv))
        end
        function $name{K,V}(itr) where {K,V}
            return reduce(push, itr, init=$name{K,V}())
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = $eq
        hashfn(::Type{T}) where {T <: $name} = $h
        equalfn(::$name) = $eq
        hashfn(::$name) = $h
        # Generic construction.
        $name() = $name{Any,Any}()
        $name(::Tuple{}) = $name{Any,Any}()
        $name(p::Pair{K,V}) where {K,V} = $name{K,V}(p)
        $name(d::AbstractDict{K,V}) where {K,V} = $name{K,V}(d)
        $name(ps::Pair...) = begin
            ps = _to_pairs(ps)
            T = Base.eltype(ps)
            K = T.parameters[1]
            V = T.parameters[2]
            return $name{K,V}(ps...)
        end
        $name(itr) = begin
            ps = _to_pairs(itr)
            T = Base.eltype(ps)
            K = T.parameters[1]
            V = T.parameters[2]
            return $name{K,V}(ps...)
        end
        # Base methods.
        Base.empty(s::$name{K,V}, ::Type{J}=K, ::Type{U}=V) where {K,V,J,U} = $name{J,U}()
        Base.length(s::$name) = begin
            ks = getfield(s, :keys)
            return (ks === nothing ? 0 : length(ks))
        end
        Base.IteratorEltype(::Type{$name{K,V}}) where {K,V} = Base.HasEltype()
        Base.eltype(::Type{$name{K,V}}) where {K,V} = Pair{K,V}
        Base.eltype(::$name{K,V}) where {K,V} = Pair{K,V}
        Base.IteratorSize(::Type{$name{K,V}}) where {K,V} = Base.HasLength()
        Base.iterate(u::$name{K,V}, ii::Int) where {K,V} = begin
            (ii > length(u)) && return nothing
            ks = getfield(u, :keys)::Vector{K}
            vs = getfield(u, :values)::Vector{V}
            return (Pair{K,V}(ks[ii], vs[ii]), #(@inbounds ks[ii], @inbounds vs[kk]),
                    ii + 1)
        end
        Base.iterate(u::$name{K,V}) where {K,V} = iterate(u, 1)
        Base.get(u::$name{K,V}, kk, df) where {K,V} = begin
            ks = getfield(u, :keys)
            (ks === nothing) && return df
            vs = getfield(u, :values)::Vector{V}
            for ii in 1:length(ks)
                k = ks[ii]
                $eq(k,kk) && return vs[ii]
            end
            return df
        end
        Base.in(kv::Pair, u::$name{K,V}, eqfn::F) where {K,V,F<:Function} = begin
            ks = getfield(u, :keys)
            (ks === nothing) && return false
            vs = getfield(u, :values)::Vector{V}
            for ii in 1:length(ks)
                k = ks[ii]
                $eq(k,kv[1]) && return eqfn(kv[2], vs[ii])
            end
            return false
        end
        Base.in(x::Pair, u::$name{K,V}) where {K,V} = in(x, u, (==))
        # And the Air methods.
        Air.push(u::$name{K,V}, kv::Pair) where {K,V} = begin
            ks = getfield(u, :keys)
            (ks === nothing) && return $name{K,V}(kv)
            vs = getfield(u, :values)::Vector{V}
            for ii in 1:length(ks)
                k = ks[ii]
                if $eq(k, kv[1])
                    if $eq(vs[ii], kv[2])
                        return u
                    else
                        return $name{K,V}(ks, setindex(vs, kv[2], ii))
                    end
                end
            end
            return $name{K,V}(push(ks, kv[1]), push(vs, kv[2]))
        end
        Air.setindex(u::$name{K,V}, v, k) where {K,V} = push(u, k => v)
        Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
            ks = getfield(u, :keys)
            (ks === nothing) && return u
            vs = getfield(u, :values)::Vector{V}
            n = length(ks)
            if n == 1
                if $eq(ks[1], x)
                    return $name{K,V}(nothing, nothing)
                end
            else
                for ii in 1:n
                    if $eq(ks[ii], x)
                        return $name{K,V}(delete(ks, ii), delete(vs, ii))
                    end
                end
            end
            return u
        end
        Base.propertynames(u::$name) = (:count, :keys, :values)
        Base.getproperty(u::$name{K,V}, s::Symbol) where {K,V} = begin
            if s == :count
                return length(u)
            elseif s == :keys
                ks = getfield(u, :keys)
                return ks === nothing ? () : tuple(ks...)
            elseif s == :values
                vs = getfield(u, :values)
                return vs === nothing ? () : tuple(vs...)
            elseif s == :pairs
                ks = getfield(u, :keys)
                (ks === nothing) && return ()
                vs = getfield(u, :values)::Vector{V}
                return Tuple([k => v for (k,v) in zip(ks,vs)])
            else
                throw(ArgumentError("no such property $s of type $(typeof(u))"))
            end
        end
    end |> esc
end
# Declare the linear types.
@_plindict_code PLinearDict isequal hash
@_plindict_code PIdLinearDict (===) objectid
# And document them.
@doc """
    PLinearDict{K,V}

A persistent dict type that uses a simple array representation internally;
accordingly `PLinearDict` does not peform persistent operations such as
`setindex`, `push`, and `delete` efficiently. It is intended to be used by Air
internally as a sub-dict that stores items whose hash values are identical.

See also: [`PIdLinearDict`](@ref), [`PDict`](@ref).
""" PLinearDict
@doc """
    PIdLinearDict{K,V}

An identity-based persistent dict type that uses a simple array representation
internally; accordingly `PIdLinearDict` does not peform persistent operations
such as `setindex`, `push`, and `delete` efficiently. It is intended to be used
by Air internally as a sub-dict that stores items whose hash values are
identical.

See also: [`PLinearDict`](@ref), [`PIdDict`](@ref).
""" PIdLinearDict

################################################################################
# PDict
macro _pdict_code(name::Symbol, eq, h, lindict)
    return quote
        struct $name{K,V} <: AbstractPDict{K,V}
            count::Int
            root::PTree{$lindict{K,V}}
        end
        function $name{K,V}() where {K,V}
            return $name{K,V}(0, PTree{$lindict{K,V}}())
        end
        function $name{K,V}(d::$name{K,V}) where {K,V}
            return d
        end
        function $name{K,V}(d::AbstractDict) where {K,V}
            return reduce(push, d, init=$name{K,V}())
        end
        function $name{K,V}(ps::Pair...) where {K,V}
            return reduce(push, ps, init=$name{K,V}())
        end
        function $name{K,V}(ps::Tuple...) where {K,V}
            return reduce(push, ps, init=$name{K,V}())
        end
        function $name{K,V}(itr) where {K,V}
            return reduce(push, itr, init=$name{K,V}())
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = $eq
        hashfn(::Type{T}) where {T <: $name} = $h
        equalfn(::$name) = $eq
        hashfn(::$name) = $h
        # Generic construction.
        $name() = $name{Any,Any}()
        $name(::Tuple{}) = $name{Any,Any}()
        $name(p::Pair{K,V}) where {K,V} = $name{K,V}(p)
        $name(d::AbstractDict{K,V}) where {K,V} = $name{K,V}(d)
        $name(ps::Pair...) = begin
            ps = _to_pairs(ps)
            T = Base.eltype(ps)
            K = T.parameters[1]
            V = T.parameters[2]
            return $name{K,V}(ps...)
        end
        $name(itr) = begin
            ps = _to_pairs(itr)
            T = Base.eltype(ps)
            K = T.parameters[1]
            V = T.parameters[2]
            return $name{K,V}(ps...)
        end
        # Base methods.
        Base.empty(s::$name{K,V}, ::Type{J}=K, ::Type{U}=V) where {K,V,J,U} = $name{J,U}()
        Base.length(s::$name) = getfield(s, :count)
        Base.IteratorEltype(::Type{$name{K,V}}) where {K,V} = Base.HasEltype()
        Base.eltype(::Type{$name{K,V}}) where {K,V} = Pair{K,V}
        Base.eltype(::$name{K,V}) where {K,V} = Pair{K,V}
        Base.IteratorSize(::Type{$name{K,V}}) where {K,V} = Base.HasLength()
        Base.iterate(u::$name{K,V}) where {K,V} = begin
            x = itervals(getfield(u, :root))
            (x === nothing) && return nothing
            (rootel,rootii) = x
            (el,reliter) = iterate(rootel)
            return (el,(rootii,reliter,length(rootel)))
        end
        Base.iterate(u::$name{K,V}, tup::Tuple{HASH_T,Int,Int}) where {K,V} = begin
            (rootii, reliter, rellen) = tup
            root = getfield(u, :root)
            if reliter < rellen
                # There are more in that node:
                rootel = get(root, rootii, nothing)::$lindict{K,V}
                (el,reliter) = iterate(rootel, reliter)
                return (el,(rootii, reliter, rellen))
            end
            x = itervals(root, rootii)
            (x === nothing) && return nothing
            (rootel,rootii) = x
            (el,reliter) = iterate(rootel)
            return (el,(rootii,reliter,length(rootel)))
        end
        Base.get(u::$name{K,V}, k, df) where {K,V} = begin
            ld = get(getfield(u, :root), $h(k), nothing)
            return (ld === nothing ? df : get(ld, k, df))
        end
        Base.in(kv::Pair, u::$name{K,V}, eqfn::F) where {K,V,F<:Function} = begin
            ld = get(getfield(u, :root), $h(kv[1]), nothing)
            return (ld === nothing ? false : in(kv, ld, eqfn))
        end
        Base.in(x::Pair, u::$name{K,V}) where {K,V} = in(x, u, (==))
        # And the Air methods.
        Air.push(u::$name{K,V}, x::Pair) where {K,V} =
            setindex(u, x.second, x.first)
        Air.setindex(u::$name{K,V}, v, k) where {K,V} = begin
            hh = $h(k)
            root = getfield(u, :root)
            ld0 = get(root, hh, nothing)
            if ld0 === nothing
                return $name{K,V}(
                    getfield(u, :count) + 1,
                    setindex(root, $lindict{K,V}(Pair{K,V}(k,v)), hh))
            else
                ld1 = setindex(ld0, v, k)
                (ld1 === ld0) && return u
                n = getfield(u, :count) + length(ld1) - length(ld0)
                return $name{K,V}(n, setindex(root, ld1, hh))
            end
        end
        Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
            hh = $h(x)
            root = getfield(u, :root)
            ld0 = get(root, hh, nothing)
            (ld0 === nothing) && return u
            ld1 = delete(ld0, x)
            (ld0 === ld1) && return u
            if length(ld1) == 0
                return $name{K,V}(getfield(u, :count) - 1, delete(root, hh))
            else
                return $name{K,V}(getfield(u, :count) - 1, setindex(root, ld1, hh))
            end
        end
    end |> esc
end
# Declare the types.    
(@_pdict_code PIdDict (===) objectid PIdLinearDict)
(@_pdict_code PDict isequal hash PLinearDict)
# Document the types.
@doc """
    PDict{K,V}

A dictionay type roughly equivalent to `Dict{K,V}` but that stores data using a
persistent hash array mapped trie system that allows for efficient persistent
operations. These operations are represented by non-mutating versions of the
standard dictionary functions, such as `setindex` in place of `setindex!`,
`push` in place of `push!`, and `pop` in place of `pop!`. These operations are
performed in `O(log n)` time, and minimal data duplication is performed in 
update operations.

See also: [`PIdDict`](@ref), [`Dict`](@ref), [`push`](@ref), [`pop`](@ref),
[`setindex`](@ref), [`delete`](@ref), [`insert`](@ref), [`getpair`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PDict{Any, ?Any}"
julia> PDict()
PDict{Any,Any}()
```

```jldoctest; filter=r"PDict{Symbol, ?Real} with 3 entries:"
julia> PDict(:a => 1, :b => 2, :c => 12.8)
PDict{Symbol,Real} with 3 entries:
  :c => 12.8
  :a => 1
  :b => 2
```

```jldoctest; filter=r"PDict{Symbol, ?Float64} with [34] entries:"
julia> d = PDict{Symbol,Float64}(:a => 1, :b => 2, :c => 12.8)
PDict{Symbol,Float64} with 3 entries:
  :c => 12.8
  :a => 1.0
  :b => 2.0

julia> :b in keys(d)
true

julia> d[:c]
12.8

julia> push(d, :d => 0.1)
PDict{Symbol,Float64} with 4 entries:
  :d => 0.1
  :c => 12.8
  :a => 1.0
  :b => 2.0
```
""" PDict
@doc """
    PIdDict{K,V}

A dictionay type roughly equivalent to `IdDict{K,V}` but that stores data using
a persistent hash array mapped trie system that allows for efficient persistent
operations. These operations are represented by non-mutating versions of the
standard dictionary functions, such as `setindex` in place of `setindex!`,
`push` in place of `push!`, and `pop` in place of `pop!`. These operations are
performed in `O(log n)` time, and minimal data duplication is performed in 
update operations.

See also: [`PDict`](@ref), [`IdDict`](@ref), [`push`](@ref), [`pop`](@ref),
[`setindex`](@ref), [`delete`](@ref), [`insert`](@ref), [`getpair`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PIdDict{Any, ?Any}\\(\\)"
julia> PIdDict()
PIdDict{Any,Any}()

```jldoctest; filter=r"PIdDict{Symbol, ?Real} with 3 entries:"
julia> PIdDict(:a => 1, :b => 2, :c => 12.8)
PIdDict{Symbol,Real} with 3 entries:
  :c => 12.8
  :a => 1
  :b => 2
```

```jldoctest; filter=r"PIdDict{Symbol, ?Float64} with [34] entries:"
julia> d = PIdDict{Symbol,Float64}(:a => 1, :b => 2, :c => 12.8)
PIdDict{Symbol,Float64} with 3 entries:
  :c => 12.8
  :a => 1.0
  :b => 2.0

julia> :b in keys(d)
true

julia> d[:c]
12.8

julia> push(d, :d => 0.1)
PIdDict{Symbol,Float64} with 4 entries:
  :d => 0.1
  :c => 12.8
  :a => 1.0
  :b => 2.0
""" PIdDict

# Export the relevant symbols.
export PDict, PIdDict
