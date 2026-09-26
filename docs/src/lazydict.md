# Persistent Lazy Dictionaries

One of the core pieces of `Air` is the persistent dictionary type,
[`PDict`](pdict.md). A nearly identical type is the `LazyDict`. Lazy
dictionaries differ from normal persistent dictionaries in that a
`LazyDict{K,V}` may be given a key and values of types `K` and `Delay{V}`,
respectively, and in this case, only calculates and yields the actual value (of
type `V`) if and when requested. Otherwise, the `LazyDict` and `LazyIdDict` types
are broadly similar to the `PDict` and `PIdDict` types.

## Examples

Values are wrapped in a [`Delay`](util.md), which computes its value only when
it is first read:

```julia
julia> using Air

julia> d = LazyDict{Symbol,Int}(:a => Delay{Int}(() -> (println("computing :a..."); 7)));

julia> d[:a]             # the first read computes the value
computing :a...
7

julia> d[:a]             # subsequent reads reuse it
7
```

Because the values are delays, a `LazyDict` is a natural way to describe a
mapping of keys to expensive computations that may never be needed. Entries that
are never read are never computed:

```julia
julia> d = LazyDict{Symbol,Int}(
           :a => Delay{Int}(() -> (println("computing :a..."); 1)),
           :b => Delay{Int}(() -> (println("computing :b..."); 2)),
       );

julia> d[:a]             # :b is never computed
computing :a...
1
```

[`@delay`](@ref) provides a concise way to write the values, and a trailing type
annotation gives the delay its element type:

```julia
julia> d = LazyDict{Symbol,Float64}(:x => (@delay 1 / 4)::Float64);

julia> d[:x]
0.25
```

Aside from the delaying of values, `LazyDict` behaves like [`PDict`](@ref), and
`LazyIdDict` like [`PIdDict`](@ref).

See also: [`PDict`](@ref), [`PIdDict`](@ref), [`Delay`](@ref), [`@delay`](@ref).
