# Persistent Arrays

One of the core pieces of `Air` is the persistent array type. This type is
broadly similar to Julia's native array type; the `PArray{T,N}` type imitates
the `Array{T,N}` type. This page walks through several simple examples of the
usage of this type.

## Examples

```julia
julia> using Air

# Persistent arrays are most commonly made from other arrays.
julia> v = PArray([1.0, 4.4, 2.9])
3-element PArray{Float64,1}:
 1.0
 4.4
 2.9

# They can also be made using pfill(), pones(), and pzeros(), which are all
# similar to their non-persistent counterparts fill(), ones(), and zeros().
julia> m = pfill(5, (2,3))
2×3 PArray{Int64,2}:
 5  5  5
 5  5  5

# The operations push(), pushfirst(), pop(), and popfirst() are all similar to
# their mutable equivalents (push!(), pushfirst!(), pop!(), and popfirst!()),
# and all are efficient with PVector (PArray{*,1}) objects, allocating only
# the path they change.
julia> v = push(v, -1.8)
4-element PArray{Float64,1}:
  1.0
  4.4
  2.9
 -1.8

julia> v = pushfirst(v, -4.3)
5-element PArray{Float64,1}:
 -4.3
  1.0
  4.4
  2.9
 -1.8

julia> pop(v)
4-element PArray{Float64,1}:
 -4.3
  1.0
  4.4
  2.9
  
# Updates can be performed with setindex(), which is like setindex!().
julia> setindex(m, -100, 2, 1)
2×3 PArray{Int64,2}:
    5  5  5
 -100  5  5
```


## Operations that keep the array persistent

A `PArray` is a dense/sparse hybrid. Every array has a default value, and any
position that is not stored explicitly reads as that default. Operations on a
`PArray` preserve both halves of that: a broadcast maps the *default* as well as
the stored entries, and a result equal to the new default is simply not stored.
So a sparse array stays sparse, and a batch of values that happen to equal the
default stores nothing at all.

That is true of broadcasting, `map`, `filter`,
`reverse`, indexing with a vector of positions, `vcat` and
`hcat`, the elementwise arithmetic operators, and `reshape` — which
reuses the same tree with a different shape and so costs nothing.

It is also true of the operations that reorder a vector's values — `sort`,
`unique`, `circshift`, `deleteat` and `splice` — and of `repeat` and
`cat`, which build a new array out of one or more arguments.

Note what these preserve. A `PArray`'s positions are positions and not merely the
entries that happen to be stored, so `sort` keeps the array's *length*, and the
positions that were unset stay unset rather than vanishing. And the default
travels with the entries into the result, so neither `repeat` nor `cat` stores an
entry that only restates it: repeating a four-element sparse vector that holds
two entries gives two stored entries per copy rather than four. Where the
arguments disagree about the default the first one's wins, and the entries that
disagree with it become stored, since they can no longer be the default.

`LinearAlgebra`'s `triu` and `rotl90` belong here too, and `triu` is worth a
note. It *zeroes* what lies below the `k`-th diagonal rather than setting it to
the array's default, and those two coincide only when the default is zero. In
that case the positions below the diagonal already read as zero and `triu` only
has to drop the stored entries there; in the other case it must store an explicit
zero at each of them, so the operation costs the shape rather than the entries.
That asymmetry is the operation's rather than the representation's, and no
representation avoids it. `rotl90` both rotates and transposes, so the result's
shape is the argument's reversed; like everything here it moves the stored
entries and carries the default, so a sparse matrix stays sparse.

Two things are
deliberately otherwise: `similar` returns a mutable `Array`, because that is what
it promises its caller and an immutable array cannot be it; and matrix
multiplication is left to Base, because `A * B` is not elementwise and a mutable
`Matrix` is the right answer for it.

## Transients

[`transient`](@ref) yields a mutable view of an array for batch updates. It
supports `push!` and `pop!`, and `setindex!` in place, with the interface of an
`Array`; [`persistent!`](@ref) hands back a `PArray` in O(1) when the batch is
done. A transient is an `AbstractArray` but deliberately *not* an
`AbstractPArray`, which is the hierarchy for persistent collections.

The gain is measured rather than assumed, and it is narrower than it looks. For
appends a transient is a clear win, since consecutive appends reuse the same
rightmost node. For updates to existing entries the crossover is around ten
updates per batch, and a single update is *slower* through a transient than
directly. The per-type docstrings — [`TArray`](@ref), [`TDict`](@ref),
[`TSet`](@ref) and their identity-keyed variants — carry the tables.
