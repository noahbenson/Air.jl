################################################################################
# api.jl
#
# The core API of the Air library. Abstract types are defined here as well as
# a number of functions that are common across many of the persistent
# collections. Functions that are defined for core Julia types (like Tuple and
# Array) are defined here as well.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

import Base.setindex
import Base.delete!

################################################################################
# Abstract Types
"""
    AbstractPDict{K,V}

`AbstractPDict` is a subtype of `AbstractDict` that is extended only by
persistent dictionary types such as `PDict` and `LazyDict`.

See also [`PDict`](@ref), [`PWDict`](@ref), [`LazyDict`](@ref),
[`PIdDict`](@ref), [`LazyIdDict`](@ref), [`PWIdDict`](@ref).
"""
abstract type AbstractPDict{K,V} <: AbstractDict{K,V} end
"""
    AbstractPSet{T}

`AbstractPSet` is an abstract type extended by persistent set types such as
 `PSet` and `PWSet.

See also: [`PSet`](@ref), [`PIdSet`](@ref), [`PWSet`](@ref), [`PWIdSet`](@ref).
"""
abstract type AbstractPSet{T} <: AbstractSet{T} end
"""
    AbstractPArray{T,N}

The Abstract persistent array type represents any `Array`-like type that is
persistent. The Air library provides the reified type `PArray{T,N}`.

See also: [`PArray`](@ref), [`PVector`](@ref).
"""
abstract type AbstractPArray{T,N} <: AbstractArray{T,N} end
export AbstractPDict, AbstractPSet, AbstractPArray
"""
    ReentrantRef{T}

A `ReentrantRef` is a type of ref object that is thread-safe. There are a few
strategies for this, each of which are encoded in a different object type that
inheits from `ReentrantRef`. These types are as follows.

`Var`: `Var` objects act like `Ref` objects except that changes to them are
exclusively task-local and must be performed in specific scoped expressions.
However, the state of all `Var` objects for the current task can also be
saved and restored at a later point.

`Actor`: `Actor` objects obey the actor pattern; you can call
`send(fn, actor, args...)` where `fn` is a function that is, in another thread
at some point, called as `fn(actor[], args...)`. The new value of the actor
after running `fn` is the return value of the call.

`Volatile`: `Volatile` objects are `Ref`s that can be changed by any thread but
that must be changed only within a synchronized transaction that ensures that
all reads and writes of volatiles, as well as reads from and sends to actors,
are atomic: either they all happen successfully, or none of them do.
"""
abstract type ReentrantRef{T} <: Ref{T} end
export ReentrantRef
"""
    TransactionalRef{T}

A `TransactionalRef` is a reentrant reference that additional participates in
transactions. Transactional refs include `Actor`s and `Volatile`s.
"""
abstract type TransactionalRef{T} <: ReentrantRef{T} end
export TransactionalRef


# #equalfn and #hashfn #########################################################
"""
    equalfn(x)

If `x` is an object (such as a `PSet` or `Dict`) that has an opinion about
equality, `equalfn(x)` returns the function that it uses.

See also: [`hashfn`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> equalfn(Dict) == isequal
true

julia> equalfn(IdDict) == (===)
true

julia> equalfn(PSet) == isequal
true

julia> equalfn(PArray)
ERROR: ArgumentError: no equalfn for type PArray
```
"""
equalfn(x::T) where {T} = equalfn(T)
equalfn(x::Type{T}) where {T} = throw(ArgumentError("no equalfn for type $T"))
equalfn(x::Type{Dict}) = isequal
equalfn(x::Type{Set}) = isequal
equalfn(x::Type{IdDict}) = (===)
equalfn(x::Type{Base.IdSet}) = (===)
equalfn(x::Type{Dict{K,V}}) where {K,V} = isequal
equalfn(x::Type{Set{T}}) where {T} = isequal
equalfn(x::Type{IdDict{K,V}}) where {K,V} = (===)
equalfn(x::Type{Base.IdSet{T}}) where {T} = (===)
"""
    hashfn(x)

If `x` is an object (such as a `PSet` or `Dict`) that has an opinion about how
it hashes objects, `hashfn(x)` returns the function that it uses. It is
sufficient in almost all circumstances to define `equalfn(T)`; the `hashfn`
should always match the `equalfn` regardless.

See also: [`equalfn`](@ref)

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> hashfn(Dict) == hash
true

julia> hashfn(IdDict) == objectid
true

julia> hashfn(PSet) == hash
true

julia> hashfn(PArray)
ERROR: ArgumentError: no equalfn for type PArray
```
"""
hashfn(x::T) where {T} = hashfn(equalfn(T))
hashfn(::Type{T}) where {T} = hashfn(equalfn(T))
# equalfn and hashfn can also be used on their respective functions
equalfn(::typeof(objectid))  = (===)
equalfn(::typeof(hash))      = isequal
hashfn(::typeof(===))     = objectid
hashfn(::typeof(isequal)) = hash
export equalfn, hashfn

# #setindex ====================================================================
"""
    setindex(coll, val, index) 

Yields a copy of the given collection `coll` with the value `val` set at the
given index. This is a persistent version of `setindex!` and works for most
collection types including `Array`s, `Dict`s, `IdDict`s, and Air's persistent
versions of these.

Note that setindex() *always* returns a copy of the argument `coll` or fails.
For Air's persistent collections these operations are efficient, but for the
mutable counterparts, they are typically `O(n)`.

See also: `Base.setindex!`, [`push`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> setindex((:a,:b,:e,:d), :c, 3)
(:a, :b, :c, :d)
```

```jldoctest; filter=r"4-element (Array{Symbol,1}|Vector{Symbol}):"
julia> u = [:a,:b,:e,:d]; v = setindex(u, :c, 3)
4-element Array{Symbol,1}:
 :a
 :b
 :c
 :d

julia> u == v
false

julia> u[3]
:e
```

```jldoctest; filter=[r"Dict{Any, ?Any} with 1 entry:", r"Dict{Any, ?Any}()"]
julia> d1 = Dict(); setindex(d1, 10, :a)
Dict{Any,Any} with 1 entry:
  :a => 10

julia> d1
Dict{Any,Any}()
```
"""
setindex(u::AbstractArray{T,N}, x::S, I...) where {T,N,S} =
    Base.setindex!(copy(u), x, I...)
setindex(d::AbstractDict{K,V}, v::U, k::J) where {K,V,U,J} =
    Base.setindex!(copy(d), v, k)
export setindex

# #push ========================================================================
"""
    push(coll, val)
    push(coll1, val1, val2...)

Yields a copy of the given collection `coll` with the given value `val`
appended. This function is essentially a persistent equivalent of the `push!`
function that never modifies the object `coll`. For persistent collections in 
the Air library, this operation is efficient, but for most mutable objects, it
is `O(n)`.

See also: `Base.push!`, [`pop`](@ref), `Base.setindex!`.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"5-element (Array{Symbol,1}|Vector{Symbol}):"
julia> push((:a,:b,:c,:d), :e)
(:a, :b, :c, :d, :e)
```

```jldoctest; filter=r"5-element (Array{Symbol,1}|Vector{Symbol}):"
julia> u = [:a,:b,:c,:d]; push(u, :e)
5-element Array{Symbol,1}:
 :a
 :b
 :c
 :d
 :e
```

```jldoctest
julia> s = Set([:a, :b, :c]); push(s, :d)
Set{Symbol} with 4 elements:
  :a
  :b
  :d
  :c
```

```jldoctest; filter=r"Dict{Symbol, ?Int64} with 4 entries:"
julia> d = Dict(:a => 1, :b => 2); push(d, :c => 3, :d => 4)
Dict{Symbol,Int64} with 4 entries:
  :a => 1
  :b => 2
  :d => 4
  :c => 3
```
"""
function push end
push(A, a, b, c...) = reduce(push, c, init=push(push(A, a), b))
push(u::NTuple{N,T}, val::S) where {T,N,S} = begin
    return NTuple{N+1,T}(T[u..., val])
end
push(u::AbstractVector{T}, val::S) where {T,S} = T[u..., val]
push(d::AbstractDict{K,V}, kv::Pair{J,U}) where {K,V,J,U} = begin
    return setindex(d, kv[2], kv[1])
end
push(s::AbstractSet{T}, u::S) where {T,S} = push!(copy(s), u)
export push

# #pop =========================================================================
"""
    pop(coll)

Yields a tuple `(last, most)` where `last` is the last element of the given
collection and most is a duplicate tuple of all but the last element of tup.
This is basically a persistent equivalent to the `pop!` function that never
modifies the given `coll`. For persistent collections in Air, these operations
are efficient, but for most mutable types, it is `O(n)`.

    pop(coll, key[, default])

Similar to `pop!`, pops the specific `key` from the given collection `coll` and
yields `(val, rest)` where `val` is the value associated with `key` in `coll`
and `rest` is the remainder of the collection without `key`. If the `key` is not
in `coll`, then `default` is yielded or an error is thrown.

See also: `Base.pop!`, [`push`](@ref), [`popat`](@ref),
[`delete`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> pop((:a,:b,:c,:d))
(:d, (:a, :b, :c))

julia> u = [:a,:b,:c,:d]; pop(u)
(:d, [:a, :b, :c])

julia> length(u)
4

julia> s = Set([:a, :b, :c]); pop(s, :c)
(:c, Set([:a, :b]))

julia> :c in s
true
```

```jldoctest; filter=r"\\(3, Dict\\(:a => 1, ?:b => 2\\)\\)"
julia> d = Dict(:a => 1, :b => 2, :c => 3); pop(d, :c)
(3, Dict(:a => 1,:b => 2))

julia> d[:c]
3
```
"""
function pop end
pop(u::NTuple{N,T}) where {N,T} = (u[end], u[1:end-1])
pop(u::Tuple{}) = throw(ArgumentError("n-tuple must be non-empty"))
pop(u::AbstractArray{T,1}) where {T} = (u[end], copy(u[1:end-1]))
pop(d::AbstractDict{K,V}) where {K,V} = begin
    rem = copy(d)
    v = pop!(rem)
    return (v, rem)
end
pop(d::AbstractDict{K,V}, k::J) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k)
    return (v, rem)
end
pop(d::AbstractPDict{K,V}) where {K,V} = begin
    (length(d) == 0) && throw(
        ArgumentError("cannot pop from an empty collection"))
    kv = iterate(d)[1]
    return (kv, delete(d, kv[1]))
end
pop(d::AbstractPDict{K,V}, k::J) where {K,V,J} = begin
    kv = getpair(d, k)
    (kv === missing) && throw(KeyError(k))
    return (kv, delete(d, k))
end
pop(d::AbstractDict{K,V}, k::J, dv) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k, dv)
    return (v, rem)
end
pop(s::AbstractSet{T}) where {T,S} = begin
    rem = copy(s)
    u = pop!(rem)
    return (u, rem)
end
pop(s::AbstractSet{T}, u::S) where {T,S} = begin
    rem = copy(s)
    u = pop!(rem, u)
    return (u, rem)
end
pop(s::AbstractPSet{T}) where {T,S} = begin
    (length(d) == 0) && throw(
        ArgumentError("cannot pop from an empty collection"))
    x = iterate(d)[1]
    return (x, delete(d, x))
end
pop(s::AbstractPSet{T}, u::S) where {T,S} = begin
    rem = copy(s)
    u = pop!(rem, u)
    return (u, rem)
end
export pop

# #pushfirst ===================================================================
"""
    pushfirst(coll, val)
    pushfirst(coll1, val1, val2...)

Yields a copy of the given collection `coll` with the given value `val`
prepended. This function is essentially a persistent equivalent of the
`pushfirst!` function that never modifies the object `coll`. For persistent
collections in the Air library, this operation is efficient, but for most
mutable objects, it is `O(n)`.

See also: `Base.pushfirst!`, [`popfirst`](@ref), `Base.setindex!`,
[`push`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest; filter=r"5-element (Array{Symbol,1}|Vector{Symbol}):"
julia> pushfirst((:a,:b,:c,:d), :e)
(:e, :a, :b, :c, :d)

julia> u = [:a,:b,:c,:d]; pushfirst(u, :e)
5-element Array{Symbol,1}:
 :e
 :a
 :b
 :c
 :d

julia> length(u)
4
```
"""
function pushfirst end
pushfirst(u::NTuple{N,T}, val::S) where {T,N,S} = begin
    return NTuple{N+1,T}(T[val, u...])
end
pushfirst(u::AbstractVector{T}, val::S) where {T,S} = T[val, u...]
pushfirst(A, a, b...) = push(push(A, a), b...)
export pushfirst

# #popfirst ====================================================================
"""
    popfirst(coll)

Yields a tuple `(first, rest)` where `first` is the first element of the given
collection and rest is a duplicate tuple of all but the first element of `tup`.
This is basically a persistent equivalent to the `popfirst!` function that never
modifies the given `coll`. For persistent collections in Air, these operations
are efficient, but for most mutable types, it is `O(n)`.

See also: `Base.popfirst!`, [`pushfirst`](@ref), [`popat`](@ref),
[`delete`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> popfirst((:a,:b,:c,:d))
(:a, (:b, :c, :d))

julia> u = [:a,:b,:c,:d]; popfirst(u)
(:a, [:b, :c, :d])

julia> length(u)
4
```
"""
function popfirst end
popfirst(u::NTuple{N,T}) where {N,T} = (u[1], u[2:end])
popfirst(u::Tuple{}) = throw(ArgumentError("n-tuple must be non-empty"))
popfirst(u::AbstractArray{T,1}) where {T} = (u[1], copy(u[2:end]))
export popfirst

# #popat =======================================================================
"""
    popat(coll, k)

Yields a tuple `(el, rest)` where `el` is the `k`th element of the given
collection `coll` and `rest` is a duplicate of all but the `k`th element of
`coll`. This is basically a persistent equivalent to the `popat!` function that
never modifies the given `coll`. For persistent collections in Air, these
operations are efficient, but for most mutable types, it is `O(n)`.

See also: `Base.popat!`, [`pop`](@ref), [`popfirst`](@ref),
[`delete`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> popat((:a,:b,:c,:d), 2)
(:b, (:a, :c, :d))

julia> u = [:a,:b,:c,:d]; popat(u, 3)
(:c, [:a, :b, :d])

julia> length(u)
4
```
"""
function popat end
popat(u::NTuple{N,T}, k::K) where {N,T,K<:Integer} =
    (u[k], (u[1:k-1]..., u[k+1:end]...))
popat(u::Tuple{}) = throw(ArgumentError("n-tuple must be non-empty"))
popat(u::AbstractArray{T,1}, k::K) where {T,K<:Integer} = begin
    c = copy(u)
    el = popat!(c, k)
    return (el, c)
end
export popat

# #insert ======================================================================
"""
    insert(coll, idx, val)

Yields a copy of the given collection `coll` with the given value `cal` inserted
at the given index. Roughly equivalent to `insert!(copy(arr), idx, va)`: the
`insert` function never modifies its arguments and always yields a copy. For the
persistent collections defined in Air, this operation is efficient, but for most
mutable objects, this is `O(n)`.

See also: `Base.insert!`, [`push`](@ref), [`delete`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> insert((:a,:b,:c,:d), 2, :x)
(:a, :x, :b, :c, :d)
```

```jldoctest; filter=r"5-element (Array{Symbol,1}|Vector{Symbol}):"
julia> u = [:a,:b,:c,:d]; insert(u, 1, :x)
5-element Array{Symbol,1}:
 :x
 :a
 :b
 :c
 :d

julia> length(u)
4
```

```jldoctest; filter=r"5-element Bit.+:"
julia> b = BitArray([0,0,0,0]); insert(b, 2, 1)
5-element BitArray{1}:
 0
 1
 0
 0
 0

julia> b[2]
false
```
"""
function insert end
insert(a::AbstractVector{T}, k::K, val::S) where {T,K<:Integer,S} =
    insert!(copy(a), k, val)
insert(a::Vector{T}, k::K, val::S) where {T,K<:Integer,S} = begin
    n = length(a) + 1
    (k > n) && throw(
        ArgumentError("insert: $k is out of range for Vector of length $(n-1)"))
    out = Vector{T}(undef, n)
    (k > 1) && copyto!(out, 1, a, 1, k - 1)
    @inbounds out[k] = val
    (k < n) && copyto!(out, k + 1, a, k, n - k)
    return out
end
insert(a::Tuple, idx::II, val::T) where {II<:Integer,T} =
    (a[1:idx]..., val, a[idx:end]...)
insert(a::Tuple{}, idx::II, val::T) where {II<:Integer,T} = begin
    if idx == 1
        return (val,)
    else
        throw(ArgumentError("invalid index $idx for tuple of size 0"))
    end
end
# For tuples, we want to generate versions of this for NTuples up to size 64;
# beyond that we can use a generic function.
macro _tuple_insert_gencode(N::Int)
    # We'll want to refer to the tuple elements:
    els = [:(tup[$k]) for k in 1:N]
    # Build up the if-elseif-else expression, starting with the else:
    ifexpr = :(throw(ArgumentError(
        "insert: $k if out of range for a Tuple of length $(length(tup))")))
    ifexpr = Expr(:elseif, :(k == $(N+1)), :(push(tup, el)), ifexpr) 
    for k in N:-1:1
        ifexpr = Expr(k == 1 ? :if : :elseif,
                      :(k == $k),
                      :(($(els[1:k-1]...), el, $(els[k:end]...))),
                      ifexpr)
    end
    return quote
        insert(tup::NTuple{$N,T}, k::K, el::S) where {T,K<:Integer,S} = $ifexpr
    end |> esc
end
# Now generate functions for up to 64:
(@_tuple_insert_gencode  1)
(@_tuple_insert_gencode  2)
(@_tuple_insert_gencode  3)
(@_tuple_insert_gencode  4)
(@_tuple_insert_gencode  5)
(@_tuple_insert_gencode  6)
(@_tuple_insert_gencode  7)
(@_tuple_insert_gencode  8)
(@_tuple_insert_gencode  9)
(@_tuple_insert_gencode 10)
(@_tuple_insert_gencode 11)
(@_tuple_insert_gencode 12)
(@_tuple_insert_gencode 13)
(@_tuple_insert_gencode 14)
(@_tuple_insert_gencode 15)
(@_tuple_insert_gencode 16)
(@_tuple_insert_gencode 17)
(@_tuple_insert_gencode 18)
(@_tuple_insert_gencode 19)
(@_tuple_insert_gencode 20)
(@_tuple_insert_gencode 21)
(@_tuple_insert_gencode 22)
(@_tuple_insert_gencode 23)
(@_tuple_insert_gencode 24)
(@_tuple_insert_gencode 25)
(@_tuple_insert_gencode 26)
(@_tuple_insert_gencode 27)
(@_tuple_insert_gencode 28)
(@_tuple_insert_gencode 29)
(@_tuple_insert_gencode 30)
(@_tuple_insert_gencode 31)
(@_tuple_insert_gencode 32)
(@_tuple_insert_gencode 33)
(@_tuple_insert_gencode 34)
(@_tuple_insert_gencode 35)
(@_tuple_insert_gencode 36)
(@_tuple_insert_gencode 37)
(@_tuple_insert_gencode 38)
(@_tuple_insert_gencode 39)
(@_tuple_insert_gencode 40)
(@_tuple_insert_gencode 41)
(@_tuple_insert_gencode 42)
(@_tuple_insert_gencode 43)
(@_tuple_insert_gencode 44)
(@_tuple_insert_gencode 45)
(@_tuple_insert_gencode 46)
(@_tuple_insert_gencode 47)
(@_tuple_insert_gencode 48)
(@_tuple_insert_gencode 49)
(@_tuple_insert_gencode 50)
(@_tuple_insert_gencode 51)
(@_tuple_insert_gencode 52)
(@_tuple_insert_gencode 53)
(@_tuple_insert_gencode 54)
(@_tuple_insert_gencode 55)
(@_tuple_insert_gencode 56)
(@_tuple_insert_gencode 57)
(@_tuple_insert_gencode 58)
(@_tuple_insert_gencode 59)
(@_tuple_insert_gencode 60)
(@_tuple_insert_gencode 61)
(@_tuple_insert_gencode 62)
(@_tuple_insert_gencode 63)
(@_tuple_insert_gencode 64)
export insert

# #delete ======================================================================
"""
    delete(coll, k)

Yields a copy of the collection `coll` with the item associated with the given
index or key `k` deleted. This function works on dictionaries, vectors, tuples,
and sets, and it always yields a copy of the collection without modifying its
aguments. For the persistent collections defined in Air, this is efficient, but
for most mutable objects, this is `O(n)`.

See also: `Base.delete!`, [`insert`](@ref), [`pop`](@ref).

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> delete((:a,:b,:c,:d), 2)
(:a, :c, :d)
```

```jldoctest; filter=r"3-element (Array{Symbol,1}|Vector{Symbol}):"
julia> u = [:a,:b,:c,:d]; delete(u, 1)
3-element Array{Symbol,1}:
 :b
 :c
 :d
```

```jldoctest
julia> s1 = Set([:a,:b,:c,:d]); s2 = delete(s1, :d); (:d in s1, :d in s2)
(true, false)
```

```jldoctest; filter=r"4-element Bit.+:"
julia> b = BitArray([0,0,0,1,0]); delete(b, 4)
4-element BitArray{1}:
 0
 0
 0
 0

julia> b[2]
false
```
"""
delete(d::AbstractDict{K,V}, k::J) where {K,V,J} = delete!(copy(d), k)
delete(d::AbstractSet{V}, k::U) where {V,U} = delete!(copy(d), k)
delete(u::AbstractVector{T}, k::K) where {T,K<:Integer} = begin
    n = length(u)
    V = typeof(u)
    (n == 0 || k < 0 || k > n) && throw(
        ArgumentError("delete: $k is out of range for Vector of length $n"))
    out = V(undef, n - 1)
    (k > 1) && copyto!(out, 1, u, 1, k - 1)
    (k < n) && copyto!(out, k, u, k + 1, n - k)
    return out
end
delete(a::Tuple, idx::K) where {K<:Integer,T} = (a[1:idx-1]..., a[idx+1:end]...)
delete(a::Tuple{}, idx::II) where {II<:Integer} =
    throw(ArgumentError("delete: $idx is out of range for a Tuple{}"))
# For tuples, we want to generate versions of this for NTuples up to size 64;
# beyond that we can use a generic function.
macro _tuple_delete_gencode(N::Int)
    # We'll want to refer to the tuple elements:
    els = [:(tup[$k]) for k in 1:N]
    # Build up the if-elseif-else expression, starting with the else:
    ifexpr = :(throw(ArgumentError(
        "delete: $k if out of range for a Tuple of length $(length(tup))")))
    for k in N:-1:1
        ifexpr = Expr(k == 1 ? :if : :elseif,
                      :(k == $k),
                      :(($(els[1:k-1]...), $(els[k+1:end]...))),
                      ifexpr)
    end
    return quote
        delete(tup::NTuple{$N,T}, k::K) where {T,K<:Integer} = $ifexpr
    end |> esc
end
# Now generate functions for up to 64:
(@_tuple_delete_gencode  1)
(@_tuple_delete_gencode  2)
(@_tuple_delete_gencode  3)
(@_tuple_delete_gencode  4)
(@_tuple_delete_gencode  5)
(@_tuple_delete_gencode  6)
(@_tuple_delete_gencode  7)
(@_tuple_delete_gencode  8)
(@_tuple_delete_gencode  9)
(@_tuple_delete_gencode 10)
(@_tuple_delete_gencode 11)
(@_tuple_delete_gencode 12)
(@_tuple_delete_gencode 13)
(@_tuple_delete_gencode 14)
(@_tuple_delete_gencode 15)
(@_tuple_delete_gencode 16)
(@_tuple_delete_gencode 17)
(@_tuple_delete_gencode 18)
(@_tuple_delete_gencode 19)
(@_tuple_delete_gencode 20)
(@_tuple_delete_gencode 21)
(@_tuple_delete_gencode 22)
(@_tuple_delete_gencode 23)
(@_tuple_delete_gencode 24)
(@_tuple_delete_gencode 25)
(@_tuple_delete_gencode 26)
(@_tuple_delete_gencode 27)
(@_tuple_delete_gencode 28)
(@_tuple_delete_gencode 29)
(@_tuple_delete_gencode 30)
(@_tuple_delete_gencode 31)
(@_tuple_delete_gencode 32)
(@_tuple_delete_gencode 33)
(@_tuple_delete_gencode 34)
(@_tuple_delete_gencode 35)
(@_tuple_delete_gencode 36)
(@_tuple_delete_gencode 37)
(@_tuple_delete_gencode 38)
(@_tuple_delete_gencode 39)
(@_tuple_delete_gencode 40)
(@_tuple_delete_gencode 41)
(@_tuple_delete_gencode 42)
(@_tuple_delete_gencode 43)
(@_tuple_delete_gencode 44)
(@_tuple_delete_gencode 45)
(@_tuple_delete_gencode 46)
(@_tuple_delete_gencode 47)
(@_tuple_delete_gencode 48)
(@_tuple_delete_gencode 49)
(@_tuple_delete_gencode 50)
(@_tuple_delete_gencode 51)
(@_tuple_delete_gencode 52)
(@_tuple_delete_gencode 53)
(@_tuple_delete_gencode 54)
(@_tuple_delete_gencode 55)
(@_tuple_delete_gencode 56)
(@_tuple_delete_gencode 57)
(@_tuple_delete_gencode 58)
(@_tuple_delete_gencode 59)
(@_tuple_delete_gencode 60)
(@_tuple_delete_gencode 61)
(@_tuple_delete_gencode 62)
(@_tuple_delete_gencode 63)
(@_tuple_delete_gencode 64)
export delete

# #getpair =====================================================================
"""
    getpair(d, k)

If the key `k` is found in the dictionary `d`, yields the pair (`k => d[k]`);
otherwise yields `missing`.

See also: `Base.get`, `Base.getindex`.

# Examples

```@meta
DocTestSetup = quote
    using Air
end
```

```jldoctest
julia> d = Dict(:a => 1, :b => 2, :c => 3); getpair(d, :a)
:a => 1

julia> getpair(d, :x)
missing
```
"""
function getpair end
let _private_symbol = gensym("private_getpair_c61d3fee7deae4d3d2237de9044247b2")
    global getpair(d::AbstractDict, k) = begin
        v = get(d, k, _private_symbol)
        return (v === _private_symbol) ? missing : (k => v)
    end
end
export getpair

# #isequal #####################################################################
function _equalfn_isequal(
    l::DL, r::DR
) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}
}
    # Identical dictionaries are always equal.
    (l === r) && return true
    # Equality functions must match.
    (equalfn(DL) === equalfn(DR)) || return false
    # Dictionaries must be the same length to be equal.
    (length(l) == length(r)) || return false
    # All pairs from one must be in the other.
    for pair in l
        in(pair, r, isequal) || return false
    end
    return true
end
Base.isequal(l::AbstractDict, r::AbstractPDict) = _equalfn_isequal(l, r)
Base.isequal(l::AbstractPDict, r::AbstractDict) = _equalfn_isequal(l, r)
Base.isequal(l::AbstractPDict, r::AbstractPDict) = _equalfn_isequal(l, r)
function _equalfn_isequal(
    l::SL, r::SR
) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}
}
    # Sets must have the same equality base to be equal.
    (equalfn(SL) === equalfn(SR)) || return false
    # Identical sets are always equal.
    (l === r) && return true
    # Sets must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all elements must be found in both.
    for el in l
        in(el, r) || return false
    end
    return true
end
Base.isequal(l::AbstractSet, r::AbstractPSet) = _equalfn_isequal(l, r)
Base.isequal(l::AbstractPSet, r::AbstractSet) = _equalfn_isequal(l, r)
Base.isequal(l::AbstractPSet, r::AbstractPSet) = _equalfn_isequal(l, r)
function _equalfn_eq(
    l::DL, r::DR
) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}
}
    # Dictionaries must have the same equality base to be equal.
    (equalfn(DL) === equalfn(DR)) || return false
    # Identical dictionaries are always equal.
    (l === r) && return true
    # Dictionaries must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all pairs must be found in both.
    anymissing = false
    for pair in l
        isin = in(pair, r)
        if ismissing(isin)
            anymissing = true
        else
            isin || return false
        end
    end
    return anymissing ? missing : true
end
Base.:(==)(l::AbstractDict, r::AbstractPDict) = _equalfn_eq(l, r)
Base.:(==)(l::AbstractPDict, r::AbstractDict) = _equalfn_eq(l, r)
Base.:(==)(l::AbstractPDict, r::AbstractPDict) = _equalfn_eq(l, r)
function _equalfn_eq(
    l::SL, r::SR
) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}
}
    # Sets must have the same equality base to be equal.
    (equalfn(SL) === equalfn(SR)) || return false
    # Identical sets are always equal.
    (l === r) && return true
    # Sets must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all elements must be found in both.
    for el in l
        in(el, r) || return false
    end
    return true
end
Base.:(==)(l::AbstractSet, r::AbstractPSet) = _equalfn_eq(l, r)
Base.:(==)(l::AbstractPSet, r::AbstractSet) = _equalfn_eq(l, r)
Base.:(==)(l::AbstractPSet, r::AbstractPSet) = _equalfn_eq(l, r)

