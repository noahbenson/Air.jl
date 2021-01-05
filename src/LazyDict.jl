################################################################################
# LazyDict
# A lazy persistent dictionary type that uses the PDict types.

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
    return (typejoin(ks), typejoin(vs))
end
# Line with PDicts, we construct LazyDicts in an unfortunaetly complex way in
# order to reduce code duplication. The three kinds of persistent sets are very
# similar, just using different equality tests, so we can keep most of the code
# the same.
const LazyDictLeaf{K,V} = PList{Pair{K,Union{V,Delay{V}}}} where {K,V}
macro _ldict_code(name::Symbol, eqfn, hashfn)
    eq = gensym()
    h = gensym()
    return quote
        begin
            local $eq = $eqfn
            local $h = $hashfn
            struct $name{K,V} <: AbstractPDict{K,V}
                dict::PDict{K,Delay{V}}
                function $name{K,V}() where {K,V}
                    return new{K,V}(PDict{K,Delay{V}}())
                end
                function $name{K,V}(itr) where {K,V}
                    return reduce(push, itr, init=new{K,V}(PDict{K,Delay{V}}()))
                end
                function $name{K,V}(d::$name{K,V}) where {K,V}
                    return d
                end
                function $name{K,V}(ps::Union{Pair,Tuple}...) where {K,V}
                    return reduce(push, ps, init=new{K,V}(PDict{K,Delay{V}}()))
                end
                function $name{K,V}(dict::PDict{K,Delay{V}}) where {K,V}
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
                (K,V,ps) = _lazypairs_typejoin(itr)
                return $name{K,V}(ps)
            end
            $name(d::AbstractDict{K,V}) where {K,V} = $name{K,V}(d)
            $name(::Tuple{}) = $name{Any,Any}()
            $name(p::Pair{K,V}) where {K,V} = $name{K,V}(p)
            $name(ps::Union{Pair,Tuple}...) = begin
                (K,V,ps) = _lazypairs_typejoin(itr)
                return $name{K,V}(ps)
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
            Base.in(kv::Pair, u::$name{K,V}, eqfn::F) where {K,V,F<:Function} = begin
                v = get(u, kv.first, nothing)
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
                !isa(v, Delay{V}) && (v = Delay{V}(v))
                return $name{K,V}(setindex(dict, v, k))
            end
            Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
                dict = getfield(u, :dict)
                newd = delete(dict, x)
                (newd === dict) && return u
                return $name{K,V}(newd)
            end
        end
    end |> esc
end
    
@_ldict_code LazyDict isequiv equivhash
@_ldict_code LazyIdDict (===) objectid
@_ldict_code LazyEqualDict isequal hash

mutability(::Type{LazyDict}) = Immutable
mutability(::Type{LazyDict{K,V}}) where {K,V} = Immutable
mutability(::Type{LazyIdDict}) = Immutable
mutability(::Type{LazyIdDict{K,V}}) where {K,V} = Immutable
mutability(::Type{LazyEqualDict}) = Immutable
mutability(::Type{LazyEqualDict{K,V}}) where {K,V} = Immutable
