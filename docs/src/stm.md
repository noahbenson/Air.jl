# Transactions (Software Transactional Memory)

Air includes a software transactional memory (STM) system, in the spirit of
Clojure's `ref`s and `dosync`. It lets you take a set of reads and writes that
must happen together and have them either all take effect or none of them do,
without taking locks by hand.

Three object types participate:

* [`Volatile`](@ref) — a reference-like object whose value may only be changed
  inside a transaction. This is the analogue of a Clojure `ref`.
* [`Actor`](@ref) — a worker with its own queue and value, to which functions
  are *sent* rather than called. This is the analogue of a Clojure agent.
* [`Var`](@ref) — a task-local binding; see the [Task-Local Variables](var.md)
  page. `Var`s are not transactional, but they are how the current transaction
  is tracked internally.

## Volatiles and `tx`

A `Volatile` behaves like a `Ref` except that it may only be written inside a
transaction. The transaction body is run with [`tx`](@ref) (or the
[`@tx`](@ref) macro):

```julia
julia> using Air

julia> v = Volatile{Int}(0)
Volatile{Int64}(@Hx4dz8Y0A3w)

julia> @tx v[] = 5
5

julia> v[]
5

julia> v[] = 6          # writing outside a transaction is an error
ERROR: cannot set volatile outside of transaction
```

Because the body may be re-run, it should be free of side effects that are not
themselves transactional:

```julia
julia> @tx (v[] = v[] + 1)
6
```

## Atomicity and retries

A transaction records the values it reads and the values it writes. When it
commits, it checks that none of the values it read have changed in the
meantime. If any have, the transaction is discarded and run again, so the body
of a transaction may execute more than once. Once a transaction commits, all of
its writes become visible at once.

If the body throws, nothing is written:

```julia
julia> v = Volatile{Int}(0);

julia> try
           @tx begin
               v[] = 99
               error("abort")
           end
       catch
       end

julia> v[]
0
```

A transaction may deliberately restart itself by throwing `TxRetryException`.
Nested transactions are folded into the enclosing one, so calling `tx` inside a
transaction simply runs the body in the current transaction.

## Filters and finalizers

A `Volatile` may carry two functions that constrain its value. A *filter* runs
whenever the volatile is set directly (`v[] = x`), and a *finalizer* runs
immediately before the transaction attempts to commit. Both replace the value
that would otherwise be stored, and both may throw to abort the transaction.

```julia
julia> v = Volatile{Int}(0);

julia> tx() do
           setfinalize!(v, x -> x * 2)   # committed values are doubled
       end

julia> @tx v[] = 3;

julia> v[]
6
```

Use [`getfilter`](@ref), [`getfinalize`](@ref), [`setfilter!`](@ref) and
[`setfinalize!`](@ref) to inspect and change them. Note that changing a filter
does not re-apply it to the current value.

## Actors

An [`Actor`](@ref) owns a value and a queue of pending functions. Sending a
function to an actor schedules it to run on the actor's own task; its return
value becomes the actor's new value. Reads and sends made inside a transaction
are only committed if the transaction commits.

```julia
julia> a = Actor{Symbol}(:start)
Actor{Symbol}(@97gWy32lwvV)

julia> send(x -> Symbol(string(x), "_done"), a);

julia> a[]
:start_done
```

Sending is asynchronous, so reading the actor immediately after `send` may still
yield the old value. A function that throws puts the actor into an error state:
subsequent reads and sends raise the captured exception until the actor is
reset.

```julia
julia> b = Actor{Symbol}(:start);

julia> send(x -> error("boom"), b);

julia> geterror(b) isa Air.ActorException
true

julia> reset(b, :ok)     # restart with a new value
Air.ActorException{Symbol}(ErrorException)

julia> b[]
:ok
```

## IO

[`TxIO`](@ref) wraps an `IO` object so that output is performed on a separate
task and writes made inside a transaction only happen if the transaction
commits. [`airout`](@ref) is a `TxIO`-like wrapper around `stdout`.

## API

The docstrings for every type and function discussed above are in the full
[API reference](API.md).

See also: [`Var`](@ref), [`@var`](@ref).
