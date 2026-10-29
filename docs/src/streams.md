# Streams

A `Stream{T}` is a FIFO of values that may be written and read repeatedly. It
is the reusable counterpart of a [`Promise`](@ref): where a promise carries one
value to one or more readers, a stream carries a sequence of values from any
number of writers to any number of readers, and it may be closed when the
sequence ends.

```julia
julia> using Air

julia> s = Stream{Int}();

julia> put!(s, 1); put!(s, 2);

julia> take(s)
1

julia> take(s)
2
```

`put!` appends and [`take`](@ref) removes the oldest value. Both block when they
cannot proceed — a read on an empty stream, a write on a full one — and both
accept the `timeout` and `token` keywords of `take(::Promise)`. Streams compare
by identity, so one is a perfectly good dictionary key.

## Bounding a stream

`maxsize` bounds how many values may be waiting. With `maxsize = 2`, two values
may be pending and the third `put!` waits for a `take`; this is backpressure,
so a producer cannot run away from a consumer. The default, `nothing`, leaves
the stream unbounded, and then `put!` never waits.

```julia
julia> s = Stream{Int}(maxsize = 1);

julia> put!(s, 1);

julia> try
           put!(s, 2; timeout = 0.01)
       catch e
           sprint(showerror, e)
       end
"reader timed out after 0.01 seconds"

julia> take(s)
1
```

[`maxsize`](@ref) reports the bound, and `length`, `isempty` and `isready`
inspect the stream without blocking. Here `isready` means "a `take` would return or throw
rather than wait", so a closed stream is ready even when it is empty, just as a
failed promise is.

[`snapshot`](@ref) copies the values that are waiting, oldest first, without
consuming them. Its result is an independent `Vector`, so it is the form to log,
keep, or serialize; a stream describes itself further with `isopen` and
`maxsize`.

```julia
julia> s = Stream{Int}(maxsize = 3); put!(s, 1); put!(s, 2);

julia> snapshot(s)
2-element Vector{Int64}:
 1
 2

julia> length(s)        # nothing was taken
2
```

## Closing a stream

`close` ends a stream. A later `put!` throws
[`StreamClosed`](@ref), and a `take` yields the values that remain and then
throws it too, so a reader can drain whatever was written before the end. Every
waiting task is woken, and closing twice is not an error.

```julia
julia> s = Stream{Int}();

julia> put!(s, 7); close(s);

julia> take(s)
7

julia> try
           take(s)
       catch e
           sprint(showerror, e)
       end
"stream is closed"
```

## Cancellation and timeouts

A blocking `take` or `put!` observes the calling task's [`CancelToken`](@ref)
exactly as `take(::Promise)` does: by default it uses [`canceltoken`](@ref), a
[`Cancelled`](@ref) exception wins even over an available value, and
`token = nothing` opts out. `timeout`, in seconds, bounds either wait and throws
[`ReaderTimeout`](@ref). A value or a free slot that arrives as the timer
expires wins over the timer, so a stream at the edge of its timeout does not
lose a race it actually won.

## Streams and transactions

Streams take part in transactions. Inside a transaction body, `put!`, `take` and
`close` are **staged**: they happen only if the transaction commits, and the
body reads its own writes. A staged `take` yields the same values a
non-transactional `take` would, in the same order, and it does not consume
anything if the attempt aborts.

```julia
julia> s = Stream{Int}();

julia> tx(() -> put!(s, 1));

julia> take(s)
1

julia> try
           tx() do
               put!(s, 2)
               error("abort")
           end
       catch
       end

julia> isempty(s)   # the staged put did not happen
true
```

A concurrent change is a conflict. If another task — one outside a transaction,
or one whose transaction commits first — writes to or reads from the stream
while this transaction is running, the commit's version check fails and the
body runs again from the start, with its staged operations discarded. This is
the same optimistic rule the STM applies to [`Volatile`](@ref)s: a transaction
always commits against the values it first saw. A stream that a transaction
only inspects with `length`, `isempty`, `isready`, `isopen` or `snapshot` is not
part of that validation; those are advisory snapshots, not validated reads, so a
body should not branch on them and then assume the branch still holds at commit.
(Inside a transaction, `snapshot` shows the body's own staged puts and takes.)

Two consequences are worth stating plainly:

* **Blocking inside a transaction is allowed.** Nothing is locked until the
  commit phase, so waiting for another task to put a value does not hold up
  anyone's lock. It does, however, make the body long-lived: use `timeout` and
  `token` to bound the wait.
* **A stream created in the current transaction cannot be waited on.** No other
  task can have a reference to it until the transaction commits, so a `take`
  that would block, or a `put!` on a full newborn stream, raises
  [`StreamDeadlock`](@ref) instead of hanging forever. Using such a stream
  without blocking — including taking a value this same body put — is fine.

Streams created outside a transaction are unaffected by the second rule, even
when a transaction uses them, and the marker does not outlive the transaction
that created the stream.

Like every other Air primitive, a stream operation made from a promise callback
or a volatile finalizer runs with no transaction in effect, so it takes effect
immediately and is not rolled back; see
[Promises and transactions](stm.md#Promises-and-transactions).

## Working with several tasks

Any number of tasks may read and write one stream. Every waiting reader is
woken by a write and every waiting writer by a read; a task that is woken but
loses the race for the value or the slot simply waits again, so no wake-up is
lost. Among several waiting readers or writers the winner is unspecified —
streams are not fair and impose no priority.

Consumption through transactions is optimistic, so many transactional
consumers on one hot stream conflict and retry. For a high rate of consumption,
prefer one non-transactional consumer per stream.

## API

The docstrings for `Stream`, [`StreamClosed`](@ref) and
[`StreamDeadlock`](@ref) are in the full [API reference](API.md), along with
those for `take`, `put!`, `close`, [`snapshot`](@ref), `isopen`, [`maxsize`](@ref),
`isempty`, `length` and `isready`. All of those but `take`, `maxsize` and
`snapshot` extend `Base`; the `Promise` methods that `take` and `put!` share
with streams are documented there too.
