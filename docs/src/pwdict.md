# Persistent Weighted Dictionaries

`Air` supports, in additon to the standard persistent dictionary type
[`PDict`](pdict.md), a weighted dictionary `PWDict` (and `PWIdDict`). These types
are broadly similar to `PDict` and `PIdDict` except that their entries
additionally have non-negative numeric weights. These dictoinaries act as (1)
dictionaries, (2) priority queues, and (3) discrete random variables.

## Examples

Weights are required whenever an entry is inserted, so they are given as the
second half of a nested pair, `(key => value) => weight`:

```julia
julia> using Air, Random

julia> d = PWDict{Symbol,Int,Float64}();

julia> d = push(d, (:a => 1) => 1.0)
PWDict(:a => 1)

julia> d = push(d, (:b => 2) => 5.0)
PWDict(:b => 2, :a => 1)

julia> d = push(d, (:c => 3) => 3.0)
PWDict(:b => 2, :c => 3, :a => 1)
```

A weighted dictionary behaves as an ordinary dictionary:

```julia
julia> d[:a]
1

julia> d[:missing]
ERROR: KeyError: key :missing not found

julia> length(d)
3
```

It also behaves as a priority queue, ordered by weight:

```julia
julia> first(d)          # the key-value pair with the highest weight
:b => 2

julia> pop(d)            # a copy without that entry
PWDict(:c => 3, :a => 1)

julia> d                 # the original is unchanged
PWDict(:b => 2, :c => 3, :a => 1)

julia> collect(d)        # iterates from highest to lowest weight
3-element Vector{Pair{Symbol, Int64}}:
 :b => 2
 :c => 3
 :a => 1
```

Weights can be inspected and changed, and they determine the distribution used
by `rand`:

```julia
julia> getweight(d, :b)
5.0

julia> d = setweight(d, :a, 10.0)
PWDict(:a => 1, :b => 2, :c => 3)

julia> first(d)
:a => 1

julia> rand(d)           # :a is now ten times as likely as :c
:a => 1
```

See also: [`PWSet`](@ref), [`PDict`](@ref), [Persistent Heaps](pheap.md),
[`getweight`](@ref), [`setweight`](@ref).
