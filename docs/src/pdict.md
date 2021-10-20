# Persistent Dictionary Examples

One of the core pieces of `Air` is the persistent dictionary types. These types
are broadly similar to Julia's native dictionary types; the `PDict{K,V}` type
imitates the `Dict{K,V}` type, and the `PIdDict{K,V}` type imitates
`IdDict{K,V}`. This page walks through several simple examples of the usage of
these types.

```julia
julia> using Air

# Creating a persistent dictionary is just like using Dict().
julia> cube = PDict{Symbol,Float64}()
PDict{Symbol,Float64} with 0 entries

# Instead of push!() and setindex!(), PDict objects efficiently produce
# persistent duplicates of themselves with any requested updates using
# functions like push() and setindex().
julia> cube = push(cube, :height => 0.2)
PDict{Symbol,Float64} with 1 entry:
  :height => 0.2

julia> cube = push(cube, :width => 9.4, :depth => 1.1)
PDict{Symbol,Float64} with 3 entries:
  :depth  => 1.1
  :height => 0.2
  :width  => 9.4

# The delete!() function is replaced by delete().
julia> delete(cube, :height)
PDict{Symbol,Float64} with 2 entries:
  :depth => 1.1
  :width => 9.4

# This does not modify the original object.
julia> cube
PDict{Symbol,Float64} with 3 entries:
  :depth  => 1.1
  :height => 0.2
  :width  => 9.4

# Lookup operatoins are nearly as fast as with native Dict objects.
julia> cube[:height]
0.2

# The behavioir of PDict is generally similar to Dict.
julia> cube[:name]
ERROR: KeyError: key :name not found
```
