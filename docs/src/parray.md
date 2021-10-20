# Persistent Array Examples

One of the core pieces of `Air` is the persistent array type. This type is
broadly similar to Julia's native array type; the `PArray{T,N}` type imitates
the `Array{T,N}` type. This page walks through several simple examples of the
usage of this type.

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
# and all are very efficient with PVector (PArray{*,1}) objects.
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

