# Persistent Weighted Sets

`Air` supports, in additon to the standard persistent set type
[`PSet`](pset.md), a weighted set `PWSet` (and `PWIdSet`). These types are
broadly similar to `PSet` and `PIdSet` except that their entries additionally
have non-negative numeric weights. These sets act as (1) persistent sets, (2)
priority queues, and (3) discrete random variables.

## Examples

Elements are given with their weights as pairs, and are ordered from highest to
lowest weight:

```julia
julia> using Air, Random

julia> s = PWSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWSet{Symbol, Float64}([:c, :b, :a])

julia> collect(s)
3-element Vector{Symbol}:
 :c
 :b
 :a
```

A weighted set supports the usual set operations:

```julia
julia> :b in s
true

julia> length(s)
3

julia> delete(s, :b)
PWSet{Symbol, Float64}([:c, :a])
```

It also supports the weighted extras:

```julia
julia> first(s)          # the element with the highest weight
:c

julia> pop(s)            # a copy without that element
PWSet{Symbol, Float64}([:b, :a])

julia> getweight(s, :b)
0.2

julia> setweight(s, :a, 9.9)   # a copy with one weight changed
PWSet{Symbol, Float64}([:a, :c, :b])

julia> rand(s)           # draws an element with probability proportional to weight
:c
```

See also: [`PWDict`](@ref), [`PSet`](@ref), [Persistent Heaps](pheap.md),
[`getweight`](@ref), [`setweight`](@ref).
