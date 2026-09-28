################################################################################
# PWSet.jl
#
# The Persistent weighted collection PWSet which is built on top of the
# persistent heap type in PHeap.jl transient heap type.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020 Noah C. Benson

# ==============================================================================
# AbstractPWSet
"""
    AbstractPWSet{T,W}

AbstractPWSet is an abstract type implemented by all persistent weighted set
types.

See also: [`PWSet`](@ref), [`AbstractPSet`](@ref).
"""
abstract type AbstractPWSet{T,W<:Number} <: AbstractPSet{T} end

# ==============================================================================
# PWSet
macro _pwset_code(name::Symbol, dicttype::Symbol)
    return esc(
        quote
            struct $name{T,W} <: AbstractPWSet{T,W}
                heap::PHeap{T,W,typeof(>),$dicttype{T,Int}}
            end
            function $name{T,W}(itr::AbstractArray) where {T,W<:Number}
                return $name{T,W}(itr...)
            end
            function $name{T,W}(itr::AbstractSet) where {T,W<:Number}
                return $name{T,W}(itr...)
            end
            function $name{T,W}(itr::AbstractDict) where {T,W<:Number}
                return $name{T,W}(itr...)
            end
            function $name{T,W}() where {T,W<:Number}
                return $name{T,W}(PHeap{T,W,typeof(>),$dicttype{T,Int}}(>))
            end
            # The copy constructor needs the same `W<:Number` bound as the
            # iteration constructors above. Without it, `$name{T,W}(d)` is
            # ambiguous with `$name{T,W}(itr::AbstractSet)`: the latter is the
            # more specific method in `W` and the former in its argument.
            # The bound costs nothing, since the struct cannot be instantiated
            # with a non-`Number` weight in the first place.
            function $name{T,W}(d::$name{T,W}) where {T,W<:Number}
                return d
            end
            function $name{T,W}(ps::Union{Tuple,Pair}...) where {T,W<:Number}
                return reduce(push, ps, init=($name{T,W}()))
            end
            function $name{T}(itr) where {T}
                return $name{T}(itr...)
            end
            function $name{T}(::Tuple{}) where {T}
                return $name{T,Float64}(PHeap{T,Float64,typeof(>),$dicttype{T,Int}}(>))
            end
            function $name{T}() where {T}
                return $name{T,Float64}(PHeap{T,Float64,typeof(>),$dicttype{T,Int}}(>))
            end
            function $name{T}(d::$name{T,W}) where {T,W}
                return d
            end
            function $name{T}(ps::Union{Tuple,Pair}...) where {T}
                return $name{T,Float64}(ps...)
            end
            function $name(itr)
                return $name(itr...)
            end
            function $name(d::$name{T,W}) where {T,W}
                return d
            end
            function $name(ps::Union{Tuple,Pair}...)
                return $name{Any,Float64}(ps...)
            end
            function $name(::Tuple{})
                return $name{Any,Float64}()
            end
            function $name()
                return $name{Any,Float64}(
                    PHeap{Any,Float64,typeof(>),$dicttype{Any,Int}}(>)
                )
            end
            # Document the equal/hash types.
            equalfn(::Type{T}) where {T<:$name} = equalfn($(dicttype))
            hashfn(::Type{T}) where {T<:$name} = hashfn($(dicttype))
            equalfn(::$name) = equalfn($(dicttype))
            hashfn(::$name) = hashfn($(dicttype))
            # Generic construction.
            # Base methods.
            Base.length(s::$name) = length(s.heap)
            Base.iterate(u::$name{T,W}) where {T,W} = iterate(u.heap)
            Base.iterate(u::$name{T,W}, x) where {T,W} = iterate(u.heap, x)
            Base.in(s::S, u::$name{T,W}, eqfn::Function) where {S,T,W} =
                in(s, u.heap, eqfn)
            Base.in(s::S, u::$name{T,W}) where {S,T,W} = in(s, u.heap)
            # And the Air methods.
            Air.push(u::$name{T,W}, x::Tuple{S,X}) where {T,W,S,X<:Number} = begin
                heap = push(u.heap, x)
                (heap === u.heap) && return u
                return $name{T,W}(heap)
            end
            Air.push(u::$name{T,W}, x::Pair{S,X}) where {T,W,S,X<:Number} =
                push(u, (x[1], x[2]))
            Air.pop(u::$name{T,W}) where {T,W} = begin
                (length(u) == 0) && throw(ArgumentError("$name must be non-empty"))
                return $name{T,W}(pop(u.heap))
            end
            Base.first(u::$name{T,W}) where {T,W} = begin
                (length(u) == 0) && throw(ArgumentError("$name must be non-empty"))
                return first(u.heap)
            end
            Air.delete(u::$name{T,W}, s::S) where {T,W,S} = begin
                heap = delete(u.heap, s)
                (heap === u.heap) && return u
                return $name{T,W}(heap)
            end
            Air.getweight(u::$name{T,W}, x::J) where {T,W,J} = getweight(u.heap, x)
            Air.setweight(u::$name{T,W}, x::J, w::X) where {T,W,J,X<:Number} = begin
                heap = setweight(u.heap, x, w)
                (heap === u.heap) && return u
                return $name{T,W}(heap)
            end
            Random.rand(u::$name{T,W}) where {T,W} = begin
                (length(u) > 0) || throw(ArgumentError("PWSet must be non-empty"))
                return rand(u.heap)
            end
        end,
    )
end
# Declare the types.    
@_pwset_code PWSet PDict
@_pwset_code PWIdSet PIdDict
# Document the types.
@doc """
    PWSet{K,V}

A persistent set with weighted elements. As such, a `PWSet` supports the
typical operations of a `PSet` as well as the following:
* The `first` function yields the element with the highest weight.
* The `pop` function yields a copy of the set without the element that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PSet`](@ref), [`getweight`](@ref), [`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PWSet{Symbol, ?Float64} with 3 elements:"
julia> PWSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWSet{Symbol,Float64} with 3 elements:
  :c
  :b
  :a
```
""" PWSet
@doc """
    PWIdSet{T}

A persistent set with weighted elements. As such, a `PWIdSet` supports the
typical operations of a `PIdSet` as well as the following:
* The `first` function yields the element with the highest weight.
* The `pop` function yields a copy of the set without the element that has
  the highest weight.
* Iteration occurs in the order of greatest to least weight.
* The weights can be changed with the `getweight` and `setweight` functions;
  `setweight` yields a duplicate dictionary with updated weights.

See also: [`PIdSet`](@ref), [`PWSet`](@ref), [`getweight`](@ref),
[`setweight`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"PWIdSet{Symbol, ?Float64} with 3 elements:"
julia> PWIdSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWIdSet{Symbol,Float64} with 3 elements:
  :c
  :b
  :a
```
""" PWIdSet
# A few functions needed here instead of inside the macro.
Base.empty(s::PWSet{T,W}, ::Type{S}=T, ::Type{X}=W) where {T,W,S,X} = PWSet{S,X}()
Base.empty(s::PWIdSet{T,W}, ::Type{S}=T, ::Type{X}=W) where {T,W,S,X} = PWIdSet{S,X}()
# As in `pwdict.jl`: a weighted set is never equal to an unweighted one, and
# these fallbacks take `AbstractSet` rather than `Any` so that they do not
# overlap every other `isequal` method and turn it ambiguous.
Base.isequal(s::ST, t::AbstractSet) where {ST<:AbstractPWSet} = false
Base.isequal(t::AbstractSet, s::ST) where {ST<:AbstractPWSet} = false
# The mixed cases; see the analogous note in `pwdict.jl`. Without these,
# `isequal(::PSet, ::PWSet)` — one argument from each side — is ambiguous.
Base.isequal(s::AbstractPSet, t::AbstractPWSet) = false
Base.isequal(s::AbstractPWSet, t::AbstractPSet) = false
# As in `pwdict.jl`: the abstract `AbstractPWSet` in the signature is what makes
# this the most specific method for two weighted sets, rather than a
# parameterised form that `api.jl`'s set `isequal` overlaps.
function Base.isequal(s::AbstractPWSet, t::AbstractPWSet)
    (length(t) == length(s)) || return false
    for x in s
        (x in t) || return false
        (getweight(s, x) == getweight(t, x)) || return false
    end
    (equalfn(typeof(s)) === equalfn(typeof(t))) && return true
    for x in t
        (x in s) || return false
    end
    return true
end
Base.hash(u::PWS) where {T,W,PWS<:AbstractPWSet{T,W}} = begin
    h = length(u)
    for s in u.heap
        w = getweight(u.heap, s)
        h += hash(s) ⊻ hash(w)
    end
    return h
end

# Export the relevant symbols.
export PWSet, PWIdSet

# The weighted counterparts of `filter` and `replace` (`api.jl` holds the
# unweighted ones; these live here because `AbstractPWSet` is defined in this
# file, which is included after `api.jl`). Each element's weight follows it: a
# filtered or replaced element keeps the weight it had.
function Base.filter(f, s::AbstractPWSet)
    out = empty(s)
    # `pset_view` yields each element with its weight in one pass, where iterating
    # the set itself costs a pop per element (see the note above it).
    for (x, w) in pset_view(s)
        f(x) && (out = push(out, x => w))
    end
    return out
end
function Base.replace(s::AbstractPWSet, pairs::Pair...)
    alt = _altlookup(pairs...)
    out = empty(s)
    for (x, w) in pset_view(s)
        out = push(out, get(alt, x, x) => w)
    end
    return out
end

# ==============================================================================
# Combining weights

# When an operation merges two weighted sets, an element they share arrives with a
# weight from each, and `weight` says what to do about it:
#
#   * `:first`, `:last`, `:min`, `:max`, `:sum` and `:mean` are folds over the
#     weights in argument order, so an element found in one argument folds to that
#     argument's weight — which is the rule for everything but a caller's own
#     function;
#   * `:median` needs the weights themselves;
#   * a caller's own `f(el, weights)` gets the element and the vector of weights it
#     was found with, and is called for *every* element of the result, not only the
#     ones several arguments shared.
#
# The folds are what make this cheap: because the weights arrive as an iterator,
# `:sum` and friends never build the vector that only `:median` and a caller's
# function need. `Statistics` is not a dependency of this package, so the two
# statistics are spelled out here rather than imported.
const _WEIGHT_RULES = (:first, :last, :sum, :mean, :min, :max, :median)

# The rule as it is dispatched on: a `Val` for one of the names, the function
# itself for a caller's own.
function _weight_rule(w)
    (w isa Symbol) || return w
    (w in _WEIGHT_RULES) || throw(
        ArgumentError("weight: $w is not one of $(_WEIGHT_RULES), nor a function")
    )
    return Val(w)
end

"""
    _wcombine(rule, el, weights)

Yields the weight an element takes, given the weights it was found with. See the
note above for what each `rule` means.
"""
_wcombine(::Val{:first}, el, ws) = first(ws)
_wcombine(::Val{:last}, el, ws) = last(ws)
_wcombine(::Val{:sum}, el, ws) = sum(ws)
_wcombine(::Val{:min}, el, ws) = minimum(ws)
_wcombine(::Val{:max}, el, ws) = maximum(ws)
_wcombine(::Val{:mean}, el, ws) = begin
    (s, n) = (zero(first(ws)), 0)
    for w in ws
        s += w
        n += 1
    end
    return s / n
end
_wcombine(::Val{:median}, el, ws) = _median(collect(ws))
_wcombine(f::Function, el, ws) = f(el, collect(ws))
# The median of a set of weights: the middle one of an odd count, the mean of the
# two middle ones of an even count.
function _median(v::AbstractVector)
    W = eltype(v)
    (isempty(v)) && throw(ArgumentError("median: no weights"))
    sort!(v)
    n = length(v)
    return isodd(n) ? v[(n + 1) ÷ 2] : (v[n ÷ 2] + v[n ÷ 2 + 1]) / 2
end

# The type the combined weight will have, so that the result's type is known
# before the elements are visited. The folds preserve the weights' own type; the
# two statistics divide, and a caller's function is asked of the compiler.
_wtype(::Union{Val{:first},Val{:last},Val{:sum},Val{:min},Val{:max}}, ::Type{E}, ::Type{W}) where {E,W} = W
_wtype(::Union{Val{:mean},Val{:median}}, ::Type{E}, ::Type{W}) where {E,W} = typeof(zero(W) / 1)
_wtype(f::Function, ::Type{E}, ::Type{W}) where {E,W} = Base.promote_op(f, E, Vector{W})

# The elements of `s` that `keep` accepts, each with the weight it takes from the
# arguments that hold it. `keep` is given the element and the argument list.
function _pwset_combine(s, args, keep, elems, rule)
    E = _pwset_eltype(args)
    W = _pwset_wtype(args)
    Wr = _wtype(rule, E, W)
    out = empty(s, E, Wr)
    for x in elems
        keep(x, args) || continue
        w = _wcombine(rule, x, (getweight(t, x) for t in args if x in t))
        out = push(out, x => w)
    end
    return out
end
# The element and weight types of the result: the arguments' promoted.
_pwset_eltype(args) = promote_type(map(eltype, args)...)
_pwset_wtype(args) = promote_type(map(Air._pwset_w, args)...)
_pwset_w(::AbstractPWSet{T,W}) where {T,W} = W

"""
    union(s::AbstractPWSet, others::AbstractPWSet...; weight=:first)

Yields the persistent weighted set of every element of the given sets. An element
that several of them hold takes its weight from `weight`; see the note in this
file for what that may be. An element held by only one keeps that one's weight.
"""
function Base.union(s::AbstractPWSet, ss::AbstractPWSet...; weight = :first)
    args = (s, ss...)
    return _pwset_combine(s, args, (x, a) -> true, _pwiter_all(args), _weight_rule(weight))
end

# The elements of every argument, without repeats: a `Set` of the elements, which
# is what keeps a shared element from being visited once per argument.
_pwiter_all(args) = begin
    # Via the views: iterating a `PWSet` costs a pop per element, and all this
    # needs is the elements.
    ks = Set{_pwset_eltype(args)}()
    for t in args, (x, _) in pset_view(t)
        push!(ks, x)
    end
    ks
end

"""
    intersect(s::AbstractPWSet, t::AbstractPWSet; weight=:first)

Yields the persistent weighted set of the elements both hold, each taking its
weight from `weight`.
"""
function Base.intersect(s::AbstractPWSet, t::AbstractPWSet; weight = :first)
    args = (s, t)
    return _pwset_combine(s, args, (x, a) -> x in t, s, _weight_rule(weight))
end

"""
    setdiff(s::AbstractPWSet, t::AbstractPWSet; weight=:first)

Yields the persistent weighted set of the elements of `s` that are not in `t`.
Each is held by only one of the arguments, so `weight` does not come into it — a
caller's own function is still called, as its contract says.
"""
function Base.setdiff(s::AbstractPWSet, t::AbstractPWSet; weight = :first)
    args = (s, t)
    return _pwset_combine(s, args, (x, a) -> !(x in t), s, _weight_rule(weight))
end

"""
    symdiff(s::AbstractPWSet, t::AbstractPWSet; weight=:first)

Yields the persistent weighted set of the elements in exactly one of the two. As
with `setdiff`, each is held by only one argument.
"""
function Base.symdiff(s::AbstractPWSet, t::AbstractPWSet; weight = :first)
    args = (s, t)
    return _pwset_combine(
        s, args, (x, a) -> (x in t) ⊻ (x in s),
        union(_pwiter_all((s,)), _pwiter_all((t,))), _weight_rule(weight),
    )
end

# ==============================================================================
# A weighted set as a set of pairs

# A weighted set is a heap whose nodes are `(value, weight, total subtree weight)`,
# and the heap's own traversal — which `PWSet` and `PWDict` inherit — yields those
# values in weight order, at a sort's cost: O(n log n) comparisons, measured at
# 16.5 KB and 54 us for a thousand elements. The weight is already stored beside
# the value, so a view that walks the node vector directly is O(1) to build and a
# single pass to traverse: 32 bytes and 11 us for the same thousand. That is the
# trade the views make — the storage's order rather than the weights'.

"""
    pset_view(s::AbstractPWSet) -> AbstractPSet{Tuple{T,W}}

Yields a view of a weighted set as a set of `(element, weight)` pairs. Building
one is O(1), and traversing it is one pass over the collection's storage, where
traversing the set itself yields the same elements ordered by weight at a sort's
cost. For a thousand elements the view is 32 bytes and 11 us against the ordered
traversal's 16.5 KB and 54 us.

The view iterates in the order of the underlying storage, not the
weights-sorted order that `PHeap`'s own traversal gives. That is what it trades
for the difference in cost.

It is read-only, having no storage of its own: `push`, `delete` and `empty` raise.
"""
pset_view(s::P) where {T,W,P<:AbstractPWSet{T,W}} = PWSetView{T,W,P}(s)
# The type `pset_view` returns; part of `Air`'s internal implementation details.
struct PWSetView{T,W,P<:AbstractPWSet{T,W}} <: AbstractPSet{Tuple{T,W}}
    p::P
end
Base.length(v::PWSetView) = length(v.p)
Base.eltype(::Type{<:PWSetView{T,W}}) where {T,W} = Tuple{T,W}
Base.eltype(v::PWSetView{T,W}) where {T,W} = Tuple{T,W}
Base.IteratorEltype(::Type{<:PWSetView}) = Base.HasEltype()
Base.IteratorSize(::Type{<:PWSetView}) = Base.HasLength()
function Base.iterate(v::PWSetView, st...)
    nxt = iterate(getfield(getfield(v.p, :heap), :_heap), st...)
    (nxt === nothing) && return nothing
    ((k, w, _), s) = nxt
    return ((k, w), s)
end
function Base.in(x::Tuple, v::PWSetView)
    (length(x) == 2) || return false
    (k, w) = x
    return (k in v.p) && (getweight(v.p, k) == w)
end
_readonly(v) = error(
    "$(typeof(v)) is a read-only view; its argument is the collection itself"
)
push(v::PWSetView, x) = _readonly(v)
delete(v::PWSetView, x) = _readonly(v)
Base.empty(v::PWSetView) = _readonly(v)

"""
    pdict_view(d::AbstractPWDict) -> AbstractPDict{K,Tuple{V,W}}

Yields a view of a weighted dictionary as a dictionary of `key => (value,
weight)` pairs, with the same guarantees as [`pset_view`](@ref): O(1) to build,
one pass to traverse, read-only, and in the storage's order rather than the
weights-sorted one.

One thing the view cannot avoid: the weight is stored beside the key in the heap,
but the value lives in the dictionary, so every element costs a dictionary lookup
on top of the pass. That is a property of the split rather than of the view.
Iterating the dictionary and looking each weight up instead is worse rather than
better — measured at 428 KB and 468 us per thousand against 32 bytes and 57 us
here — because the dictionary's own iteration and the heap index's lookup are
both paid per element either way, and the index lookup is the dearer of the two.
"""
pdict_view(d::P) where {K,V,W,P<:AbstractPWDict{K,V,W}} = PWDictView{K,V,W,P}(d)
# The type `pdict_view` returns; part of `Air`'s internal implementation details.
struct PWDictView{K,V,W,P<:AbstractPWDict{K,V,W}} <: AbstractPDict{K,Tuple{V,W}}
    d::P
end
Base.length(v::PWDictView) = length(v.d)
Base.eltype(::Type{<:PWDictView{K,V,W}}) where {K,V,W} = Pair{K,Tuple{V,W}}
Base.eltype(v::PWDictView{K,V,W}) where {K,V,W} = Pair{K,Tuple{V,W}}
Base.IteratorEltype(::Type{<:PWDictView}) = Base.HasEltype()
Base.IteratorSize(::Type{<:PWDictView}) = Base.HasLength()
function Base.iterate(v::PWDictView{K,V,W}, st...) where {K,V,W}
    nxt = iterate(getfield(getfield(v.d, :heap), :_heap), st...)
    (nxt === nothing) && return nothing
    ((k, w, _), s) = nxt
    return (k => (getfield(v.d, :dict)[k], w), s)
end
Base.get(v::PWDictView, k, df) = haskey(v.d, k) ? (v.d[k], getweight(v.d, k)) : df
Base.haskey(v::PWDictView, k) = haskey(v.d, k)
# The *value* of the view's dictionary is the tuple, as for any dictionary; it is
# only the iteration that yields `key => value` pairs, as for any dictionary.
function Base.getindex(v::PWDictView, k)
    return (v.d[k], getweight(v.d, k))
end
push(v::PWDictView, x::Pair) = _readonly(v)
delete(v::PWDictView, x) = _readonly(v)
Base.empty(v::PWDictView) = _readonly(v)
export pset_view, pdict_view
