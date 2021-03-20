################################################################################
# lazydict.jl
#
# A lazy persistent dictionary type that uses the PDict types.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

# This function is used for processing sequences of arguments when constructing
# a lazy dictionary: you want to typejoin across pairs/tuples, and you want to
# figure out Pair{K,V} for which all the given values are either of type V or
# of Delay{V}.
function _lazypairs_typejoin(args)
    n = length(args)
    ks = Vector{Type}(undef, n)
    vs = Vector{Type}(undef, n)
    top = 0
    for (k,v) in args
        top += 1
        ks[top] = typeof(k)
        vT = typeof(v)
        if vT <: Delay
            vT = vT.parameters[1]
        end
        vs[top] = vT
    end
    return (typejoin(ks...), typejoin(vs...))
end

# Like with PDicts, we construct LazyDicts in an unfortunaetly complex way in
# order to reduce code duplication. The three kinds of persistent sets are very
# similar, just using different equality tests, so we can keep most of the code
# the same.
macro _ldict_code(name::Symbol, eq, h, dtype)
    return quote
        struct $name{K,V} <: AbstractPDict{K,V}
            dict::$dtype{K,Delay{V}}
            function $name{K,V}() where {K,V}
                return new{K,V}($dtype{K,Delay{V}}())
            end
            function $name{K,V}(itr) where {K,V}
                return reduce(push, itr, init=new{K,V}($dtype{K,Delay{V}}()))
            end
            function $name{K,V}(d::$name{K,V}) where {K,V}
                return d
            end
            function $name{K,V}(ps::Union{Pair,Tuple}...) where {K,V}
                return reduce(push, ps, init=new{K,V}($dtype{K,Delay{V}}()))
            end
            function $name{K,V}(dict::$dtype{K,Delay{V}}) where {K,V}
                return new{K,V}(dict)
            end
        end
        # Document the equal/hash types.
        equalfn(::Type{T}) where {T <: $name} = $eq
        hashfn(::Type{T}) where {T <: $name} = $h
        equalfn(::$name) = $eq
        hashfn(::$name) = $h
        # Generic construction.
        $name() = $name{Any,Any}()
        $name(itr) = begin
            (K,V) = _lazypairs_typejoin(itr)
            return $name{K,V}(itr)
        end
        $name(d::AbstractDict{K,V}) where {K,V} = $name{K,V}(d)
        $name(::Tuple{}) = $name{Any,Any}()
        $name(p::Pair{K,V}) where {K,V} = $name{K,V}(p)
        $name(ps::Union{Pair,Tuple}...) = begin
            (K,V) = _lazypairs_typejoin(ps)
            return $name{K,V}(ps...)
        end
        # Base methods.
        Base.empty(s::$name{K,V}, ::Type{J}=K, ::Type{U}=V) where {K,V,J,U} = $name{J,U}()
        Base.length(s::$name) = length(getfield(s, :dict))
        Base.iterate(u::$name{K,V}) where {K,V} = begin
            (length(u) == 0) && return nothing
            it = iterate(getfield(u, :dict))
            (it === nothing) && return nothing
            ((k,v),state) = it
            return (k => v[], state)
        end
        Base.iterate(u::$name{K,V}, state) where {K,V} = begin
            it = iterate(getfield(u, :dict), state)
            (it === nothing) && return nothing
            ((k,v),state) = it
            return (k => v[], state)
        end
        Base.get(u::$name{K,V}, k, df) where {K,V} = begin
            uu = get(getfield(u, :dict), k, nothing)
            (uu === nothing) && return df
            uu::Delay{V}
            return uu[]
        end
        Base.haskey(u::$name{K,V}, k) where {K,V} = begin
            uu = get(getfield(u, :dict), k, nothing)
            return (uu !== nothing)
        end
        Base.keys(u::$name{K,V}, k) where {K,V} = keys(getfield(u, :dict))
        Base.in(kv::Pair, u::$name{K,V}, eqfn::F) where {K,V,F<:Function} = begin
            dict = getfield(u, :dict)
            v = get(dict, kv.first, nothing)
            (v === nothing) && return false
            return eqfn(kv.second, v[])
        end
        Base.in(x::Pair, u::$name{K,V}) where {K,V} = in(x, u, $eq)
        # And the Air methods.
        Air.push(u::$name{K,V}, x::Pair) where {K,V} = begin
            (k,v) = x
            dict = getfield(u, :dict)
            currval = get(dict, k, nothing)
            if currval !== nothing
                (currval === x) && return u
                (Base.isready(currval) && currval[] === x) && return u
            end
            isa(v, Delay{V}) || (v = Delay{V}(v))
            return $name{K,V}(setindex(dict, v, k))
        end
        Air.push(u::$name{K,V}, x::Pair{J,Delay{V}}) where {K,J,V} = begin
            (k,v) = x
            dict = getfield(u, :dict)
            currval = get(dict, k, nothing)
            if currval !== nothing
                (currval === x) && return u
                (Base.isready(currval) && currval[] === x) && return u
            end
            return $name{K,V}(setindex(dict, v, k))
        end
        Air.push(u::$name{K,V}, x::Pair{J,V}) where {K,J,V} = begin
            (k,v) = x
            dict = getfield(u, :dict)
            currval = get(dict, k, nothing)
            if currval !== nothing && Base.isready(currval) && currval[] === x
                return u
            else
                return $name{K,V}(setindex(dict, Delay{V}(v), k))
            end
        end
        Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
            dict = getfield(u, :dict)
            newd = delete(dict, x)
            (newd === dict) && return u
            return $name{K,V}(newd)
        end
        Base.isready(u::$name{K,V}, k) where {K,V} = begin
            dict = getfield(u, :dict)
            v = get(dict, k, nothing)
            (v === nothing) && throw(KeyError(k))
            return !isa(v._val, _DelayPending)
        end
    end |> esc
end
# Generate the types.    
@_ldict_code LazyDict isequal hash PDict
@_ldict_code LazyIdDict (===) objectid PIdDict
# Document the types.
@doc """
    LazyDict{K,V}

A dictionay type equivalent to `PDict{K,V}` in every way except that for any
value in the dict that is a `Delay` object, the dictionary hides the delay and
always returns the delay's value. This allows any value that has not yet been
requested to be lazily unevaluated.

See also: [`PDict`](@ref), `Base.Dict`, [`LazyIdDict`](@ref),
[`Delay`](@ref), [`@delay`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"LazyDict{Any, ?Any}()"
julia> LazyDict()
LazyDict{Any,Any}()
```

```jldoctest; filter=r"LazyDict{Symbol, ?Real} with 3 entries:"
julia> LazyDict(:a => 1, :b => 2, :c => 12.8)
LazyDict{Symbol,Real} with 3 entries:
  :c => 12.8
  :a => 1
  :b => 2
```

```jldoctest; filter=r"LazyDict{Symbol, ?Float64} with 3 entries:"
julia> d = LazyDict{Symbol,Float64}(:a => 1, :b => 2, :c => 12.8)
LazyDict{Symbol,Float64} with 3 entries:
  :c => 12.8
  :a => 1.0
  :b => 2.0

julia> d2 = push(d, :d => Delay{Real}(() -> (println("Running..."); 0.5))); haskey(d2, :d)
true

julia> d2[:d]
Running...
0.5

julia> d2[:d]
0.5
```
""" LazyDict
@doc """
    LazyIdDict{K,V}

A dictionay type equivalent to `PIdDict{K,V}` in every way except that for any
value in the dict that is a `Delay` object, the dictionary hides the delay and
always returns the delay's value. This allows any value that has not yet been
requested to be lazily unevaluated.

See also: [`PIdDict`](@ref), `Base.IdDict`, [`LazyDict`](@ref),
[`Delay`](@ref), [`@delay`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"LazyIdDict{Any, ?Any}()"
julia> LazyIdDict()
LazyIdDict{Any,Any}()
```

```jldoctest; filter=r"LazyIdDict{Symbol, ?Real} with 3 entries:"
julia> LazyIdDict(:a => 1, :b => 2, :c => 12.8)
LazyIdDict{Symbol,Real} with 3 entries:
  :b => 2
  :a => 1
  :c => 12.8
```

```jldoctest; filter=r"LazyIdDict{Symbol, ?Float64} with 3 entries:"
julia> d = LazyIdDict{Symbol,Float64}(:a => 1, :b => 2, :c => 12.8)
LazyIdDict{Symbol,Float64} with 3 entries:
  :c => 12.8
  :a => 1.0
  :b => 2.0

julia> d2 = push(d, :d => Delay{Real}(() -> (println("Running..."); 0.5))); haskey(d2, :d)
true

julia> d2[:d]
Running...
0.5

julia> d2[:d]
0.5
""" LazyIdDict

# Export the relevant symbols.
export LazyDict, LazyIdDict
