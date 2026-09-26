# Utilities

Air provides a handful of small, independent utilities that do not depend on the
collection types or the transaction system.

## Lazy values: `Delay` and `@delay`

A `Delay{T}` computes its value the first time it is requested, and returns that
same value on every later request. Unlike a plain closure, a `Delay` is
thread-safe: if several tasks ask for the value at once, the underlying function
still runs exactly once.

```julia
julia> using Air

julia> d = Delay{Int}(() -> (println("computing..."); 7));

julia> isready(d);      # nothing has been computed yet

julia> d[]
computing...
7

julia> d[]              # the function is not run again
7
```

[`@delay`](@ref) builds a `Delay` from an expression. The expression may be a
closure, whose arguments are closed over, or a plain expression, which behaves
as `() -> expression`:

```julia
julia> d = @delay 10;

julia> d[]
10

julia> n = 3;

julia> d2 = @delay (n + 1);

julia> n = 100;         # the value was captured, not re-read

julia> d2[]
4
```

An optional trailing type annotation on the expression sets the element type of
the resulting `Delay`:

```julia
julia> d = @delay (2 + 3)::Int;

julia> typeof(d)
Delay{Int64}
```

`Delay` objects are used by [`LazyDict`](@ref) to represent values that have not
been computed yet.

## Memoized functions: `@memoize`

[`@memoize`](@ref) declares a function whose results are cached, keyed by the
arguments. The function body is evaluated at most once per distinct set of
arguments, and never concurrently with itself.

```julia
julia> @memoize square(n::Int) = (println("computing..."); n * n)
square (generic function with 1 method)

julia> square(4)
computing...
16

julia> square(4)
16
```

A trailing type annotation names the type cached in the memoization dictionary:

```julia
julia> @memoize half(n::Int) = (n / 2)::Float64;
```

Arguments are compared by equality, so passing mutable arguments that are later
modified leads to undefined behaviour.

## One-shot values: `Promise`

A `Promise{T}` is a container for a single value that is produced at some later
point, typically by another task. Reading a `Promise` blocks until the value is
available.

```julia
julia> p = Promise{Int}();

julia> isready(p)
false

julia> Threads.@spawn put!(p, 7);

julia> take(p)
7
```

## Locking several objects: `lockall`

`lockall` takes a collection of locks, acquires all of them (in a consistent
order to avoid deadlock), runs the given function, and releases them:

```julia
julia> r1, r2 = ReentrantLock(), ReentrantLock();

julia> lockall(() -> :done, r1, r2)
:done
```

## API

```@docs
@delay
@memoize
take
```

The remaining docstrings for this page — `Delay`, `Promise`, and `lockall` —
are in the full [API reference](API.md).
