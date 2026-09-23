# Task-Local Variables

A `Var{T}` is a reference-like object whose value is **task-local**: each task
sees its own binding, and a task that has never bound the `Var` sees its initial
value. This is the mechanism behind dynamically-scoped configuration, and it is
the same mechanism Air uses internally to track the current transaction.

`Var`s differ from `Ref`s in an important way: a `Var` cannot be assigned
directly. Its binding is established for the duration of a scope with
[`withvars`](@ref), [`setvars`](@ref), or the functions returned by
[`wrapwithvars`](@ref) and [`wrapsetvars`](@ref).

## Declaring a Var

Use [`@var`](@ref) to declare a global binding. An optional type annotation on
the right-hand side gives the `Var` its element type:

```julia
julia> using Air

julia> @var depth = 0
Var{Int64}(@JfLHavRdgnZ: 0; init=0)

julia> @var name = :root::Symbol
Var{Symbol}(@cFgU9Kqe8: :root; init=:root)

julia> depth[]
0
```

`Var`s can equally be constructed directly:

```julia
julia> v = Var{Int}(0)
Var{Int64}(@JfLHavRdgnZ: 0; init=0)
```

## Binding a Var

`withvars` runs a function with new bindings in place, and restores the previous
bindings afterwards:

```julia
julia> v = Var{Int}(0);

julia> withvars(() -> v[], v => 42)
42

julia> v[]          # the binding was scoped to the call
0
```

The `do`-block form is usually the most readable:

```julia
julia> withvars(v => 10) do
           v[] + 1
       end
11
```

`setvars` takes a dictionary instead of individual pairs, and `vars()` returns
the bindings visible in the current task:

```julia
julia> setvars(() -> v[], IdDict{Var,Any}(v => 7))
7

julia> withvars(v => 3) do
           vars()[v]
       end
3
```

Bindings nest, with the innermost one winning:

```julia
julia> withvars(v => 1) do
           withvars(v => 2) do
               v[]
           end
       end
2
```

## Task locality

Bindings are per task, so a task spawned inside a `withvars` block does **not**
inherit them:

```julia
julia> withvars(v => 42) do
           (v[], fetch(Threads.@spawn v[]))
       end
(42, 0)
```

## Wrapping a function

`wrapwithvars` returns a function that installs the given bindings whenever it
is called, which is convenient for callbacks and for spawning tasks that should
see particular bindings:

```julia
julia> v = Var{Int}(0);

julia> f = wrapwithvars(v => 5) do x
           v[] + x
       end;

julia> f(1)
6

julia> v[]          # unchanged outside the wrapped call
0
```

## API

The docstrings for `Var`, [`@var`](@ref), `vars`, `withvars`, `setvars`,
`wrapwithvars`, and `wrapsetvars` are in the full [API reference](API.md).
