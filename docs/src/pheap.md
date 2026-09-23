# Persistent Heaps

`Air` provides the internal `PHeap{T,W,F,D}` type, a persistent (immutable)
priority queue in which each element carries a non-negative numeric weight. A
persistent binary heap keeps a cached total weight for every subtree, which
allows three operations to be efficient:

1. **Priority queue**: `first` yields the element with the highest weight, and
   `pop` yields a copy of the heap without it.
2. **Ordered iteration**: a heap iterates from highest to lowest weight.
3. **Weighted sampling**: `rand` draws an element with probability proportional
   to its weight, using the cached subtree totals rather than scanning.

You will rarely construct a `PHeap` directly. It is the backing store for the
[weighted dictionaries](pwdict.md) and [weighted sets](pwset.md), which expose
this behaviour through the usual dictionary and set interfaces.

## Examples

A weighted set is the simplest way to see the ordering and sampling behaviour:

```julia
julia> using Air, Random

julia> s = PWSet{Symbol}(:a => 0.1, :b => 0.2, :c => 0.3)
PWSet{Symbol, Float64}([:c, :b, :a])

julia> first(s)          # the highest-weighted element
:c

julia> getweight(s, :b)
0.2

julia> pop(s)            # a copy without the highest-weighted element
PWSet{Symbol, Float64}([:b, :a])

julia> s                 # the original is unchanged
PWSet{Symbol, Float64}([:c, :b, :a])

julia> setweight(s, :a, 1.0)   # a copy with one weight changed
PWSet{Symbol, Float64}([:a, :c, :b])

julia> rand(s)           # weighted sampling: :c is three times as likely as :a
:c
```

`PWDict` adds the same weights to a dictionary. Note that weights are required
whenever an entry is inserted, since there is no sensible default:

```julia
julia> d = push(PWDict{Symbol,Int,Float64}(), (:a => 1) => 1.0)
PWDict(:a => 1)

julia> d = push(d, (:b => 2) => 5.0)
PWDict(:b => 2, :a => 1)

julia> first(d)          # the key-value pair with the highest weight
:b => 2

julia> d[:a]
1
```

The ordering function is a parameter of the heap. It defaults to `>`, so high
weights come first; passing `<` reverses that, while weighted sampling is
unaffected.

The docstrings for `PHeap`, [`getweight`](@ref) and [`setweight`](@ref) are in
the full [API reference](API.md).

See also: [`PWDict`](@ref), [`PWSet`](@ref), [`PDict`](@ref), [`PSet`](@ref).
