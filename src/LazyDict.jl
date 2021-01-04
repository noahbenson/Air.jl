################################################################################
# LazyDict
# A lazy persistent dictionary type that uses the PDict types.

# Line with PDicts, we construct LazyDicts in an unfortunaetly complex way in
# order to reduce code duplication. The three kinds of persistent sets are very
# similar, just using different equality tests, so we can keep most of the code
# the same.
const LazyDictLeaf{K,V} = PList{Pair{K,Union{V,Delay{V}}}} where {K,V}
macro _ldict_code(name::Symbol, eqfn, hashfn)
    let eq = gensym(), h = gensym(), _name = gensym(), q
        q = quote
            let $eq = $eqfn, $h = $hashfn
                struct $name{K,V} <: AbstractPDict{K,V}
                    _n::Int
                    _tree::PTree{LazyDictLeaf{K,V}}
                    function $name{K,V}() where {K,V}
                        return new(0, PTree{LazyDictLeaf{K,V}}())
                    end
                    function $name{K,V}(d::$name{K,V}) where {K,V}
                        return d
                    end
                    function $name{K,V}(d::AbstractDict) where {K,V}
                        return reduce(push, d, init=new(0, PTree{LazyDictLeaf{K,V}}()))
                    end
                    function $name{K,V}(ps::Pair...) where {K,V}
                        return reduce(push, ps, init=new(0, PTree{LazyDictLeaf{K,V}}()))
                    end
                    function $name{K,V}(ps::Tuple...) where {K,V}
                        return reduce(push, ps, init=new(0, PTree{LazyDictLeaf{K,V}}()))
                    end
                    function $name{K,V}(itr) where {K,V}
                        return reduce(push, itr, init=new(0, PTree{LazyDictLeaf{K,V}}()))
                    end
                    function $name{K,V}(n::Int, tree::PTree{LazyDictLeaf{K,V}}) where {K,V}
                        return new(n, tree)
                    end
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
                $name(d::AbstractDict{EquivRef{K},V}) where {K,V} = $name{K,V}(d)
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
                Base.length(s::$name) = s._n
                Base.iterate(u::$name{K,V}) where {K,V} = begin
                    (u._n == 0) && return nothing
                    (lst,titer) = iterate(u._tree)
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    v = el.second
                    el = Pair{K,V}(el.first,  isa(v, Delay) ? v[] : v)
                    return (el,(lst,liter,titer))
                end
                Base.iterate(u::$name{K,V}, tup) where {K,V} = begin
                    (lst,liter,titer) = tup
                    if liter !== nothing
                        (el,liter) = iterate(lst, liter)
                        v = el.second
                        el = Pair{K,V}(el.first, isa(v, Delay) ? v[] : v)
                        return (el,(lst,liter,titer))
                    end
                    q = iterate(u._tree, titer)
                    (q === nothing) && return q
                    (lst,titer) = q
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    v = el.second
                    el = Pair{K,V}(el.first, isa(v, Delay) ? v[] : v)
                    return (el,(lst,liter,titer))
                end
                Base.get(u::$name{K,V}, k, df) where {K,V} = begin
                    hh = $h(k)
                    uu = get(u._tree, hh, u)
                    (uu === u) && return df
                    for (kk,vv) in uu
                        $eq(k, kk) || continue
                        isa(vv, Delay) && (vv = vv[])
                        return vv
                    end
                    return df
                end
                Base.in(kv::Pair, u::$name{K,V}, eqfn::Function) where {K,V} = begin
                    v = get(u, kv.first, u)
                    (v === u) && return false
                    isa(v, Delay) && (v = v[])
                    return eqfn(kv.second, v)
                end
                Base.in(x::Pair, u::$name{K,V}) where {K,V} = in(x, u, $eq)
                # And the Air methods.
                Air.push(u::$name{K,V}, x::Pair) where {K,V} = begin
                    hh = $h(x.first)
                    uu = get(u._tree, hh, LazyDictLeaf{K,V}())
                    in(x, u, $eq) && return u
                    uu = pushfirst(uu, x)
                    tr = setindex(u._tree, uu, hh)
                    return $name{K,V}(u._n + 1, tr)
                end
                Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
                    hh = $h(x)
                    uu = get(u._tree, hh, LazyDictLeaf{K,V}())
                    (length(uu) === 0) && return u
                    vv = delete(uu, x, (k,kv) -> $eq(k, kv.first))
                    (uu === vv) && return u
                    if length(vv) == 0
                        return $name{K,V}(u._n - 1, delete(u._tree, hh))
                    else
                        return $name{K,V}(u._n - 1, setindex(u._tree, vv, hh))
                    end
                end
            end
        end
        return esc(q)
    end
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
