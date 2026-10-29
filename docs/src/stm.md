# Transactions (Software Transactional Memory)

Air includes a software transactional memory (STM) system, in the spirit of
Clojure's `ref`s and `dosync`. It lets you take a set of reads and writes that
must happen together and have them either all take effect or none of them do,
without taking locks by hand.

Four object types participate:

* [`Volatile`](@ref) — a reference-like object whose value may only be changed
  inside a transaction. This is the analogue of a Clojure `ref`.
* [`Actor`](@ref) — a worker with its own queue and value, to which functions
  are *sent* rather than called. This is the analogue of a Clojure agent.
* [`Promise`](@ref) — a one-shot value container whose delivery may be staged
  in a transaction; see [Promises and transactions](#Promises-and-transactions).
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

A sent function is called in the latest world, so it may be one that was defined
after the actor was created — a closure typed at the REPL, for example, which an
actor's long-lived task could otherwise not call.

An actor may be given two hooks, at construction or later:

* a **post-processor**, called after each scheduled function with the old value
  and the function's return value. The value it returns is stored instead, which
  is useful when the scheduled function's job is a side effect — an `IO` stream,
  for example, whose `println` returns `nothing`;
* an **error handler**, called with the actor, the value the failing function
  was given, and the exception. It observes the failure; the actor still enters
  its error state.

```julia
julia> a = Actor{Int}(0; postprocess = (old, r) -> old + r);

julia> send(_ -> 3, a);

julia> a[]                      # after the actor has processed the message
3

julia> send(_ -> 4, a);

julia> a[]
7
```

Use [`getpostprocess`](@ref), [`setpostprocess!`](@ref),
[`geterrorhandler`](@ref) and [`seterrorhandler!`](@ref) to inspect and change
them. A setter affects the next message the actor dequeues; setting a hook is
immediate configuration and is not rolled back by a transaction.

A `send` made from a volatile finalizer or a promise delivery callback runs in
the transaction's finalize phase, where no transaction is active, but the
message still belongs to the transaction that is finalizing: it is staged and
queued only if that transaction commits, once rather than once per attempt.

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

## Promises and transactions

A [`Promise`](@ref) is a one-shot container: `put!` delivers a value, `fail!`
settles it as failed, and `take` (or `p[]`) waits for either settlement. A
promise may also be given an `ondelivery` and/or an `onfailure` callback when it
is constructed.

Inside a transaction, `put!` and `fail!` do not settle: they stage the
settlement, and it happens only if the transaction commits. If the body throws,
or the transaction exhausts its retries, nothing is settled. Inside the
transaction that staged it, `take` (which throws `PromiseFailed` for a staged
failure), `isready`, `isfailed` and `geterror` see the staged settlement, so a
`take` after a `put!` does not deadlock waiting for a settlement that cannot
happen until the body returns.

The promise's callback runs in the transaction's finalize phase — the
`ondelivery` callback for a staged delivery, the `onfailure` callback for a
staged failure — after the body and after the volatile finalizers, but before
the commit locks or validates anything. It runs on **every attempt**, so under
contention it may run more than once, and it may run for an attempt that never
settles. It is a notification, not a commit hook: its return value is ignored,
and an exception it throws aborts the transaction. Outside a transaction, the
settlement happens immediately and the callback runs in the calling task once the
promise's lock is released; the promise then stays settled even if the callback
throws.

A promise settles only once. A second settlement — staged twice in one
transaction, of either kind, or attempted after another task has settled the
promise — raises `ArgumentError`. That is a permanent failure rather than a
retryable conflict, so the transaction fails at once and commits nothing.

Callbacks run with no transaction active, `currtx() === nothing`, so a volatile
read inside one sees the last committed value rather than the transaction's
staged write. A `send` made there is nonetheless staged in the transaction that
is finalizing, so a callback whose transaction loses its race queues nothing.
Every other side effect a callback performs is immediate and is not rolled back
if the transaction fails — in particular, a `put!` or `fail!` there settles its
promise at once. Since callbacks run on every attempt, a callback that settles
another promise must be careful: doing so on the first attempt makes the second
attempt settle it twice, which raises `ArgumentError` and aborts the
transaction. A callback may start a transaction of its own, but it should not
block, and it must not wait on the promise it is being called for.

## IO

[`TxIO`](@ref) wraps an `IO` object so that output is performed on a separate
task and writes made inside a transaction only happen if the transaction
commits. [`airout`](@ref) is a `TxIO`-like wrapper around `stdout`.

## API

The docstrings for every type and function discussed above are in the full
[API reference](API.md).

See also: [`Var`](@ref), [`@var`](@ref).
