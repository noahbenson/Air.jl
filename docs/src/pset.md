# Persistent Sets

The persistent set types are a core feature of `Air`. These types are broadly
similar to Julia's native set types; the `PSet{T}` type imitates the `Set{T}`
type, and the `PIdSet{T}` type imitates `Core.IdSet{T}`. This page walks through
several simple examples of the usage of these types.

```julia
julia> using Air

# Creating a persistent set is just like using Set().
julia> names = PSet(:Emily, :Lena, :Stephanie)
PSet{Symbol} with 3 elements:
  :Emily
  :Lena
  :Stephanie

# Instead of push!(), PSet objects efficiently produce persistent duplicates of
# themselves with any requested updates using functions like push().
julia> names = push(names, :Sheri)
PSet{Symbol} with 4 elements:
  :Emily
  :Lena
  :Stephanie
  :Sheri

julia> names = push(names, :Caitlin, :Joseph)
PSet{Symbol} with 6 elements:
  :Joseph
  :Emily
  :Lena
  :Stephanie
  :Caitlin
  :Sheri

# The delete!() function is replaced by delete().
julia> delete(names, :Emily)
PSet{Symbol} with 5 elements:
  :Joseph
  :Lena
  :Stephanie
  :Caitlin
  :Sheri

# This does not modify the original object.
julia> names
PSet{Symbol} with 6 elements:
  :Joseph
  :Emily
  :Lena
  :Stephanie
  :Caitlin
  :Sheri

# Lookup operatoins are nearly as fast as with native Set objects.
julia> in(:Joseph, names)
true
```
