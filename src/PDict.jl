################################################################################
# PDict
# A persistent dictionary type that uses the PTree type.

# We should define an isequiv function for abstract dicts; this should work
# equally well on PSets or normal sets.
isequiv(s::DS, t::DT) where {KT,VT,KS,VS, DS <: AbstractDict{KS,VS}, DT <: AbstractDict{KT,VT}} = begin
    (ismuttype(DT) || ismuttype(DS)) && return (s === t)
    (length(t) == length(s)) || return false
    for (k,v) in s
        tt = get(t, k, t)
        (tt === t) && return false
        isequiv(tt, v) || return false
    end
    (equalfn(DS) === equalfn(DT)) && return true
    for (k,v) in t
        ss = get(s, k, s)
        (ss === s) && return false
        isequiv(ss, v) || return false
    end
    return true
end

################################################################################
# PDict
# We construct PDicts in an unfortunaetly complex way in order to reduce code
# duplication. The three kinds of persistent sets are very similar, just using
# different equality tests, so we can keep most of the code the same.
const PDictLeaf{K,V} = PList{Pair{K,V}} where {K,V}
macro _pdict_code(name::Symbol, eqfn, hashfn)
    let eq = gensym(), h = gensym(), _name = gensym(), q,
        n = gensym("[private] length"), tree = gensym("[private] tree"),
        empt = gensym("[private] emptylist")
        q = quote
            begin
                local $eq = $eqfn
                local $h = $hashfn
                struct $name{K,V} <: AbstractDict{K,V}
                    $n::Int
                    $tree::PTree{PDictLeaf{K,V}}
                    $empt::PDictLeaf{K,V}
                    function $name{K,V}() where {K,V}
                        return new(0, PTree{PDictLeaf{K,V}}(), PDictLeaf{K,V}())
                    end
                    function $name{K,V}(d::$name{K,V}) where {K,V}
                        return d
                    end
                    function $name{K,V}(d::AbstractDict) where {K,V}
                        return reduce(push, d, init=new(0, PTree{PDictLeaf{K,V}}(),
                                                        PDictLeaf{K,V}()))
                    end
                    function $name{K,V}(ps::Pair...) where {K,V}
                        return reduce(push, ps, init=new(0, PTree{PDictLeaf{K,V}}(),
                                                         PDictLeaf{K,V}()))
                    end
                    function $name{K,V}(ps::Tuple...) where {K,V}
                        return reduce(push, ps, init=new(0, PTree{PDictLeaf{K,V}}(),
                                                         PDictLeaf{K,V}()))
                    end
                    function $name{K,V}(itr) where {K,V}
                        return reduce(push, itr, init=new(0, PTree{PDictLeaf{K,V}}(),
                                                          PDictLeaf{K,V}()))
                    end
                    function $name{K,V}(n::Int, tree::PTree{PDictLeaf{K,V}}) where {K,V}
                        return new(n, tree, PDictLeaf{K,V}())
                    end
                    function $name{K,V}(n::Int, tree::PTree{PDictLeaf{K,V}},
                                        l::PDictLeaf{K,V}) where {K,V}
                        return new(n, tree, l)
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
                Base.length(s::$name) = s.$n
                Base.iterate(u::$name{K,V}) where {K,V} = begin
                    (u.$n == 0) && return nothing
                    (lst,titer) = iterate(u.$tree)
                    lst = lst[2]
                    (el,liter) = iterate(lst)
                    return (el,(lst,liter,titer))
                end
                Base.iterate(u::$name{K,V}, tup) where {K,V} = begin
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
                Base.get(u::$name{K,V}, k, df) where {K,V} = begin
                    hh = $h(k)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) == 0) && return df
                    for (kk,vv) in uu
                        $eq(k, kk) && return vv
                    end
                    return df
                end
                Base.in(kv::Pair, u::$name{K,V}, eqfn::Function) where {K,V} = begin
                    k = kv.first
                    v = kv.second
                    hh = $h(k)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) == 0) && return false
                    for (kk,vv) in uu
                        eqfn(k, kk) && return eqfn(v, vv)
                    end
                    return false
                end
                Base.in(kv::Pair, u::$name{K,V}, ::typeof(===)) where {K,V} = begin
                    k = kv.first
                    v = kv.second
                    hh = $h(k)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) == 0) && return false
                    for (kk,vv) in uu
                        (k === kk) && return (v === vv)
                    end
                    return false
                end
                Base.in(kv::Pair, u::$name{K,V}, ::typeof(isequiv)) where {K,V} = begin
                    k = kv.first
                    v = kv.second
                    hh = $h(k)
                    #uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    #(length(uu) == 0) && return false
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) == 0) && return false
                    for (kk,vv) in uu
                        isequiv(k, kk) && return isequiv(v, vv)
                    end
                    return false
                end
                Base.in(kv::Pair, u::$name{K,V}, ::typeof(isequal)) where {K,V} = begin
                    k = kv.first
                    v = kv.second
                    hh = $h(k)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) == 0) && return false
                    for (kk,vv) in uu
                        isequal(k, kk) && return isequal(v, vv)
                    end
                    return false
                end
                Base.in(x::Pair, u::$name{K,V}) where {K,V} = in(x, u, $eq)
                # And the Air methods.
                Air.push(u::$name{K,V}, x::Pair) where {K,V} = begin
                    hh = $h(x.first)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    in(x, u, $eq) && return u
                    uu = pushfirst(uu, x)
                    tr = setindex(u.$tree, uu, hh)
                    return $name{K,V}(u.$n + 1, tr, u.$empt)
                end
                Air.setindex(u::$name{K,V}, v, k) where {K,V} = push(u, k => v)
                Air.delete(u::$name{K,V}, x::J) where {K,V,J} = begin
                    hh = $h(x)
                    uu = get(u.$tree, hh, u.$empt)::PDictLeaf{K,V}
                    (length(uu) === 0) && return u
                    vv = delete(uu, x, (k,kv) -> $eq(k, kv.first))
                    (uu === vv) && return u
                    if length(vv) == 0
                        return $name{K,V}(u.$n - 1, delete(u.$tree, hh), u.$empt)
                    else
                        return $name{K,V}(u.$n - 1, setindex(u.$tree, vv, hh), u.$empt)
                    end
                end
            end
        end
        return esc(q)
    end
end
    
@_pdict_code PDict isequiv equivhash
@_pdict_code PIdDict (===) objectid
@_pdict_code PEqualDict isequal hash

mutability(::Type{PDict}) = Immutable
mutability(::Type{PDict{K,V}}) where {K,V} = Immutable
mutability(::Type{PIdDict}) = Immutable
mutability(::Type{PIdDict{K,V}}) where {K,V} = Immutable
mutability(::Type{PEqualDict}) = Immutable
mutability(::Type{PEqualDict{K,V}}) where {K,V} = Immutable
