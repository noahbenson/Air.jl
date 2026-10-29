# Utilities

Air provides a handful of small, independent utilities that do not depend on the
collection types. Most also do not depend on the transaction system, with one
exception: a `Promise` may be delivered inside a transaction, in which case the
delivery is staged until the transaction commits (see
[Promises and transactions](stm.md#Promises-and-transactions)).

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
available. When the values form a sequence rather than a single result, use a
[`Stream`](@ref) instead; see [Streams](streams.md).

```julia
julia> p = Promise{Int}();

julia> isready(p)
false

julia> Threads.@spawn put!(p, 7);

julia> take(p)
7
```

A promise settles exactly once. It may be settled successfully with `put!`, or
as failed with [`fail!`](@ref); after that, its result cannot change. `isready`
is true after either settlement, [`isfailed`](@ref) tells the two apart, and
[`geterror`](@ref) reads a failure's payload without throwing.

```julia
julia> p = Promise{Int}();

julia> fail!(p, "no value");

julia> isready(p), isfailed(p)
(true, true)

julia> try
           take(p)
       catch e
           sprint(showerror, e)
       end
"promise failed: \"no value\""

julia> geterror(p)
"no value"
```

A promise may be given an `ondelivery` and/or an `onfailure` callback when it is
constructed. It is called with the settlement's payload, exactly once, after the
settlement has been made; its return value is ignored, and an exception
propagates while the promise stays settled. See
[Promises and transactions](stm.md#Promises-and-transactions) for what a
settlement inside a transaction means.

```julia
julia> p = Promise{Int}(ondelivery = x -> println("got ", x),
                        onfailure = e -> println("failed: ", e));

julia> put!(p, 7);
got 7

julia> take(p)
7
```

Promises compare by identity: a promise is `isequal` only to itself, and its
hash never waits for a settlement, so an undelivered promise is a perfectly good
dictionary key.

## Cancellation and timeouts

`take` can bound its wait with `timeout`, in seconds. If the promise is still
pending when the timeout expires, `take` throws [`ReaderTimeout`](@ref) instead
of waiting forever; a `Stream`'s `take` and `put!` throw the same exception.

```julia
julia> p = Promise{Int}();

julia> try
           take(p; timeout = 0.1)
       catch e
           sprint(showerror, e)
       end
"reader timed out after 0.1 seconds"
```

Every task also has a [`CancelToken`](@ref): it is created on first use and is
available as [`canceltoken`](@ref). Blocking Air calls observe the calling
task's token by default, so a task is cancellable without any setup — another
task only needs the token, typically recorded when the task is created.
[`cancel!`](@ref) cancels a token and wakes every task waiting under it with a
[`Cancelled`](@ref) exception, carrying the reason given to `cancel!`.
Cancellation is one-way and idempotent; the first reason given is kept.

```julia
julia> tok = CancelToken();

julia> cancel!(tok, "shutting down");

julia> try
           take(Promise{Int}(); token = tok)
       catch e
           sprint(showerror, e)
       end
"cancelled: shutting down"
```

Cancellation wins over settlement: a cancelled token makes `take` throw even
when the promise is already settled, because a task whose work was cancelled
should not go on to use the value. To wait regardless of cancellation — in the
cleanup of a cancelled task, for example — pass `token = nothing`. Inside a
transaction, a `take` on a *promise* that has to wait is a deadlock, since the
settlement it waits for cannot happen until the body returns (a `Stream` may be
written by another task while the body waits; see
[Streams and transactions](streams.md#Streams-and-transactions)). A `Cancelled`
thrown inside a transaction body aborts that transaction and is not retried. The
same holds for `ReaderTimeout`.

[`checkcancel`](@ref) observes cancellation without blocking (`iscancelled` is
its non-throwing form), and [`setcancaltoken!`](@ref) replaces or clears a
task's token. Air's internal `Air.sleep_interruptible(seconds; token)` is the
cancellable counterpart of `sleep`; it is not exported, because it would shadow
`Base.sleep` for anyone who does `using Air`.

Cancellation is observed by blocking waits, not everywhere: a transaction that
is retrying, and a `Delay`'s `[]`, do not check a token, so a cancelled task
should call `checkcancel` between transactions.

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
fail!
isfailed
geterror
```

The remaining docstrings for this page — `Delay`, `Promise`, and `lockall` —
are in the full [API reference](API.md).
