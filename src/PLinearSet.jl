################################################################################
# PLinearSet
# A persistent dictionary type that uses the PTree type to store a list of
# key-value pairs.


################################################################################
# PLinearSet
# See PDict.jl for information about the construction here 
macro _plinset_code(name::Symbol, eqfn, hashfn)
    eq = gensym()
    h = gensym()
    return quote
        begin
            local $eq = $eqfn
            local $h = $hashfn
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
        end
    end |> esc
end
    
@_plinset_code PLinearSet isequiv equivhash
@_plinset_code PIdLinearSet (===) objectid
@_plinset_code PEqualLinearSet isequal hash

mutability(::Type{PLinearSet}) = Immutable()
mutability(::Type{PLinearSet{T}}) where {T} = Immutable()
mutability(::Type{PIdLinearSet}) = Immutable()
mutability(::Type{PIdLinearSet{T}}) where {T} = Immutable()
mutability(::Type{PEqualLinearSet}) = Immutable()
mutability(::Type{PEqualLinearSet{T}}) where {T} = Immutable()
