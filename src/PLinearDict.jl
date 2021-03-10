################################################################################
# PLinearDict
# A persistent dictionary type that uses the PTree type to store a list of
# key-value pairs.


################################################################################
# PLinearDict
# See PDict.jl for information about the construction here 
macro _plindict_code(name::Symbol, eqfn, hashfn)
    eq = gensym()
    h = gensym()
    return quote
        begin
            local $eq = $eqfn
            local $h = $hashfn
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
        end
    end |> esc
end
    
@_plindict_code PLinearDict isequal hash
@_plindict_code PIdLinearDict (===) objectid
@_plindict_code PEquivLinearDict isequiv equivhash

mutability(::Type{PLinearDict}) = Immutable()
mutability(::Type{PLinearDict{K,V}}) where {K,V} = Immutable()
mutability(::Type{PIdLinearDict}) = Immutable()
mutability(::Type{PIdLinearDict{K,V}}) where {K,V} = Immutable()
mutability(::Type{PEquivLinearDict}) = Immutable()
mutability(::Type{PEquivLinearDict{K,V}}) where {K,V} = Immutable()
