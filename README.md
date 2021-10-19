<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://noahbenson.github.io/Air.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://noahbenson.github.io/Air.jl/dev)
[![Build Status](https://travis-ci.com/noahbenson/Air.jl.svg?branch=master)](https://travis-ci.com/noahbenson/Air.jl)
[![Codecov](https://codecov.io/gh/noahbenson/Air.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/noahbenson/Air.jl)


# Air ##########################################################################

An **A**tomic **I**mmutable **R**esource for Julia that is as light as air.


## Author ######################################################################

Noah C. Benson &lt;<nben@uw.edu>&gt;


## Description #################################################################

Air is an [immutable data structure](https://en.wikipedia.org/wiki/Persistent_data_structure)
and [software transactional memory](https://en.wikipedia.org/wiki/Software_transactional_memory)
tool.  It provides several persistent data structures including dictionaries,
sets, N-dimensional arrays, heaps, weighted dictionaries, and weighted sets,
and it includes a transaction system that allows one to write atomic operations
over such data.

Air is a Julia library that takes advantage of Julia's builtin immutable
paradigms to provide this set of functional utilities. Air is currently 
under development but includes substantial testing and is generally stable.
Inspiration for Air's design is derived largely from paradigms in
[Clojure](https://en.wikipedia.org/wiki/Clojure).

### Examples

The following code blocks demonstrate some simple examples of how one can use
the `Air` library.


#### `PDict`: Persistent Dictionaries

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

#### `PArray`: Persistent Arrays

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

#### `Actor`, `Volatile`, and `tx()`: Thread-safe transactions

```julia
julia> using Air

# Volatile objects are like Ref objects. They are designed to be always safe to
# read but to only be settable inside a transaction.
julia> vol = Volatile(PDict(:a => 1, :b => 2, :c => 3))
Volatile{PDict{Symbol,Int64}}(@CjK4yBeEKUm: PDict(:c => 3,:a => 1,:b => 2))

# Like with Ref objects, Volatiles can be access with getindex(), inside or
# outside of a transaction.
julia> vol[]
PDict{Symbol,Int64} with 3 entries:
  :c => 3
  :a => 1
  :b => 2

# However, outside of a transaction, it is an error to change a volatile.
julia> vol[] = push(vol[], :d => 4)
ERROR: cannot set volatile outside of transaction

# Inside of a transaction, it is guaranteed that no other thread will change any
# volatile that is read or changed during the transaction. This is the only safe
# time to edit volatiles.
julia> @tx vol[] = push(vol[], :d => 4)
PDict{Symbol,Int64} with 4 entries:
  :d => 4
  :c => 3
  :a => 1
  :b => 2

julia> vol[]
PDict{Symbol,Int64} with 4 entries:
  :d => 4
  :c => 3
  :a => 1
  :b => 2

# The tx() function is an alternate way to write transactions.
julia> tx() do
         vol[] = push(vol[], :e => 5)
         nothing
       end

julia> vol[]
PDict{Symbol,Int64} with 5 entries:
  :d => 4
  :e => 5
  :c => 3
  :a => 1
  :b => 2

# Actors are for managing side-effects of transactions or for managing data that
# needs to be updated through sequential single-threaded operations.
julia> notes = Actor{PVector{String}}(pfill("", 0))
Actor{PArray{String,1}}(@1j6NGQzgacH: String[])

# It is safe to read actors at any time; this similar to reading Ref values.
julia> notes[]
0-element PArray{String,1}

# We can define a function that sends a log to the message. The send function
# requires that the first argument be a function, to be run in the actor's
# thread, whose one parameter is the current value of actor and whose return
# value is the updated value for the actor.
julia> logmsg(s) = begin
         send(val -> pushfirst(val, s), notes)
         nothing
       end
logmsg (generic function with 1 method)

# We can test out the log function. Note that the pushfirst action in the
# logmsg function above is run in a separate thread (the actor's thread).
julia> logmsg("Log initialized.")

julia> notes[]
1-element PArray{String,1}:
 "Log initialized."

# Actors can be run and read in transactions, and they are guaranteed to be
# unchanged for the duration of a transaction, and no other thread will send
# it any functions during the transaction.
julia> @tx logmsg("Next log message.")

julia> notes[]
2-element PArray{String,1}:
 "Next log message."
 "Log initialized."

# If a transaction fails for any reason, no sends occur. Actors are only sent
# functions when a transacton succeeds.
julia> tx() do
         logmsg("Failure message.")
         error()
       end
ERROR: 

julia> notes[]
2-element PArray{String,1}:
 "Next log message."
 "Log initialized."
```

### Plans

Note that many of the core components for Air already have working
implementations. Several more are currently undergoing testing. In particular,
the existing persistent data structures are fairly well tested and have a
performance comparable to their clojure counterparts. Additionally, initial
tests of the thread-safe transaction system using the `@tx` transaction block
macro, `Actor`s and `Volatile`s are promising (though more testing is neeeded).

* Completed plans
  * Persistent data structures:
    * `PArray`, a persistent array type that mimics Julia's native `Array`
    * `PDict`, a persistent dictionary type that mimics Julia's native `Dict`
    * `PSet`, a persistent set type.
    * `PWSet`, a weighted persistent set type. Each element of a `PWSet` has a
      weight, and the sete itself acts both as a piority queue (`first(pwset)`
      always yields the element with the highest weight in `O(1)` time) and as
      a distribution (`rand(pwset)` yields a random element of the set where the
      probability of an element being chosen is equal to its weight divided by
      the total weight of the elements in the set.
    * `PWDict`, a persistent weighted dictionary type is like the `PWSet` type,
      but insted of elements with weights, the dictionary contains key-valye
      pairs with weights.
    * `LazyDict`, a persistent lazy dictionary type.
* Plans with incomplete testing:
  * Multi-threading utilities, inspired by Clojure
    * Thread-safe `Delay` type.
    * Thread-local `Var` type.
    * A `Volatile` type that operates with the `@tx` macro to ensure
      that all updates to references within a synchronized block are performed
      atomically.
    * An `Actor` type for sending asynchronous jobs to independent threads which
      also respects the atomic requirements of transactional `tx` blocks.
* Plans that are implemented but require testing
  * Multi-threading utilities, inspired by Clojure
    * Thread-safe `Promise` types.
* Plans that are not yet implemented
  * **Persistent Record Types.** Tool for defining data structures that
    support lazy reification of members.
  * **Inheritable Method Parameters.** A common problem in many languages is
    that an API function that calls another lower-level function often needs to
    support some subset of the argument preprocessing that occurs in the
    lower-level function. The preprocessing and documentation of the shared
    arguments shouldn't need to be duplicated in multiple places; rather a
    simple macro should make it easy for common parameters to occur in both
    places and be documented identically in both.
  * Better query/build/update API tools for the data structures?


## License

MIT License

Copyright (c) 2019-2021 Noah C. Benson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
