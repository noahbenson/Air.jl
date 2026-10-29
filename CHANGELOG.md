# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Under `0.x`, SemVer puts breaking changes in the minor position and compatible
fixes in the patch position, which is how the versions below are chosen.

## [0.3.0] - 2026-10-28

### Added

- **Actor hooks.** An `Actor` may now be given a post-processing function,
  called after each scheduled function with `(old value, function result)` and
  whose return value is stored instead, and an error handler, called with
  `(actor, value, exception)` when a scheduled function or the post-processor
  throws. Both are constructor keywords and have `get`/`set!` accessors:
  `getpostprocess`/`setpostprocess!` and `geterrorhandler`/`seterrorhandler!`.
  The `Actor` docstring has promised both since the original release without
  implementing them.
- **Promise delivery callbacks.** `Promise{T}` accepts an `ondelivery` keyword:
  a function called with the delivered value, exactly once, whose return value
  is ignored.
- **Transactional promise delivery.** A `put!` inside a transaction is now
  staged and performed only if the transaction commits, with the promise's lock
  joining the commit's sorted lock set. `take` and `isready` see a delivery
  staged by the current transaction. A second delivery raises `ArgumentError`
  and is not retried.
- **Promise failure.** `fail!(promise, err)` settles a promise as failed; a
  later `take` (or `p[]`) throws `PromiseFailed(err)` and `isready` is true.
  `isfailed` tells a failure from a delivery, and `geterror(promise)` reads the
  payload without throwing. A failed promise also accepts an `onfailure`
  keyword callback, with the same contract as `ondelivery`. Inside a
  transaction, a failure is staged like a delivery and performed only if the
  transaction commits; `take`, `isready`, `isfailed` and `geterror` see the
  staged settlement.
- **Streams.** `Stream{T}` is a reusable FIFO: `put!` appends, `take` removes
  the oldest value, and both block when they cannot proceed, with the same
  `timeout` and `token` keywords as `take(::Promise)`. A stream may be bounded
  with `maxsize`, so a full stream makes `put!` wait and a producer cannot run
  away from its consumer, and it may be closed with `close`: a later `put!`
  throws `StreamClosed`, while a `take` drains what remains and then throws it.
  Inside a transaction, `put!`, `take` and `close` are staged and take effect
  only if the transaction commits, and a concurrent change makes the attempt
  retry rather than commit against a stale view; a stream created in the
  current transaction raises `StreamDeadlock` instead of waiting, since no other
  task can reach it yet. `snapshot` copies the values that are waiting, oldest
  first, as an independent vector — the form to log, keep or serialize — and
  `isopen`, `maxsize`, `length`, `isempty` and `isready` inspect a stream
  without blocking. Streams compare by identity, like promises.
- **Timeouts and cancellation.** `take` accepts a `timeout` keyword (in
  seconds), throwing the new `ReaderTimeout` if the promise is still pending
  when it expires, and a `token` keyword. Every task has a `CancelToken`
  (`canceltoken()`, replaceable with `setcancaltoken!`); `cancel!` wakes every
  task waiting under it with the new `Cancelled` exception. Blocking Air calls
  observe the calling task's token by default, and cancellation wins even over
  an already settled promise, so `token = nothing` opts out. `checkcancel` and
  `iscancelled` observe cancellation without blocking. `Air.sleep_interruptible`
  is the cancellable counterpart of `sleep`; it is not exported. Consulting the
  token first makes `take` on an already settled promise cost about 74 ns rather
  than 4 ns; it still allocates nothing, and a `take` that has to wait is
  dominated by the wait itself.

### Changed

- `PromiseTimeout` is now `ReaderTimeout`, and serves both promise and stream
  waits. It was added in this same unreleased version, so no released API
  changes.
- The commit's sorted lock set gained a fourth group, streams, after promises;
  a stream that a transaction staged is locked and validated like a volatile.
- `Promise{T}` is now a mutable struct with a persistent mutex, a condition
  bound to it, a stored settlement, and lazily allocated callback fields and
  waiter set; the private field layout and the implicit
  `Promise{T}(::Threads.Condition)` constructor are gone. All of it is
  documented as internal.
- `isready(p::Promise)` is now true for a failed promise as well as a delivered
  one, since both are settlements.
- `Promise`s now compare by identity: `isequal(a, b)` is `a === b`, and `hash`
  no longer waits for a settlement. The old value-comparing forms called `take`,
  so an undelivered promise used as a dictionary key hung the calling task.
- A transaction started from a volatile finalizer or a promise delivery (or
  failure) callback now gets a fresh `Transaction` rather than the pooled one,
  which the finalizing transaction is still using.
- A `send` made from a volatile finalizer or a promise settlement callback is now
  staged in the transaction that is finalizing and queued only if it commits.
  It used to be queued immediately, escaping a transaction that then failed —
  once per attempt.

### Fixed

- Corrected stale `Actor` and `ActorException` docstrings: `restart` is `reset`,
  the exception's fields are `error`, `value`, `fn` and `queue`, and the
  post-processing example used a three-argument `send` that never existed.
- A task's pooled `Transaction` is now cleared when an attempt ends rather than
  when the next one begins, so the last transaction's reads, writes and staged
  messages — and through them whatever values were stored in them — are no
  longer kept alive until that task happens to run another transaction. The
  clear still happens exactly once per transaction, so nothing gets slower.
- An actor's scheduled function, post-processor and error handler are now called
  in the latest world. A message is run by the actor's long-lived task, which
  could not call a function defined after that task was started — a closure
  written at the REPL or in a test body, typically — and raised `MethodError`
  ("the applicable method may be too new") instead. The cost is about 2 ns per
  message on an actor that already costs ~165 ns, with no extra allocation.
- A promise's `ondelivery` and `onfailure` callbacks are now called in the
  latest world, for the same reason as an actor's functions: the task that
  settles a promise may be long-lived, and the callback may have been written
  after that task was started, which raised `MethodError` ("the applicable
  method may be too new") and failed the settling task instead of running the
  callback.

### Performance

- A stream operation outside a transaction — one mutex acquisition, a push or a
  pop, and a version bump — costs about 110 ns, the same as `take` on a settled
  promise. A `take` allocates nothing; a `put!` pays only the live vector's
  amortized growth. Staging the same operation in a transaction costs about
  550 ns and 0.7 KB, which is what being able to roll it back costs.

### Documentation

- Documented the actor hooks, transactional promise settlement, promise failure,
  and timeouts and cancellation in the manual, including what a settlement
  callback may safely do during a transaction and which waits observe
  cancellation.
- Added a [Streams](https://noahbenson.github.io/Air.jl/stable/streams/) page
  covering bounds, closing, cancellation, and the transactional rules, and a
  "Streams and transactions" section to the STM page. The utilities page now
  says that the in-transaction deadlock warning for `take` applies to promises,
  not to streams.

## [0.2.0] - 2026-10-02

The first release since 0.1.1. Both the persistent collections and the
transactional layer were reworked, the manual was written, and the test suite,
benchmarks and CI were rebuilt. It also folds in a month of documentation work
from late 2021 that was committed but never released.

### Changed (breaking)

- **Julia 1.10 is now the minimum**, up from 1.5.
- **`Documenter` is no longer a dependency.** It was only ever needed to build
  the manual, which now has an environment of its own.
- **`setindex` is `Air`'s own function.** It was an extension of
  `Base.setindex` — method piracy, and the source of the package's two type
  piracies — so it no longer carries `Base`'s meanings for types `Air` does not
  own.
- Removed the exports `Source`, `AbstractSourceKernel` and `receive`, none of
  which was ever defined, and deleted `src/metadata.jl`, which was never
  included.

### Added

- **Transients** — `TArray`, `TDict`, `TSet` and their identity-keyed variants: a
  mutable view for batch updates, with `persistent!` returning a persistent
  collection in O(1). Measured, the gain is narrower than it looks: appends win
  clearly, but the crossover for updating existing entries is around ten updates
  per batch, and a single update is *slower* through a transient.
- **Array generators** — `psparse`, `prand`, `prandn`, `pdiagm`, `blockdiag` and
  `permutedims`.
- **Weighted views** — `pset_view` and `pdict_view`, which traverse a weighted
  collection in a single pass rather than one pop per element.
- **A `weight` argument for the weighted set operations** — `:first`, `:last`,
  `:sum`, `:mean`, `:min`, `:max`, `:median`, or a caller's `f(element, weights)`.
- **Operations that keep the argument's kind** rather than answering with a
  mutable `Array`: broadcasting, `map`, `filter`, `reverse`, indexing with a
  vector of positions, `vcat`, `hcat`, the elementwise arithmetic, `reshape`,
  `sort`, `unique`, `circshift`, `deleteat`, `splice`, `repeat`, `cat`, `triu`
  and `rotl90`.
- `filter` and `replace` for dictionaries and sets, and `union`, `intersect`,
  `setdiff` and `symdiff` for `PSet` — the last of which previously dropped a
  `PWSet`'s weights silently, which was the worst of the wrong-kind operations.
- `PHeap` declares its element type, so `collect` and comprehensions over a heap
  build a `Vector{T}` instead of a `Vector{Any}`.
- A benchmark suite (`bench/`), type-stability baselines, and a CI matrix over
  Julia 1.10/1.11/nightly and Linux/macOS/Windows, including thread legs above
  four — nothing above four threads had ever run.

### Fixed

- **Seventeen latent defects** across the collections and the transactional
  layer, none of which was reachable from a test.
- **The persistent heap's cached subtree totals**, which drifted on deletion and
  skewed weighted sampling.
- **Twenty method ambiguities**, reported by Aqua.
- **Two type piracies**, both in `setindex`.
- An allocation on every dictionary lookup; an unset-position read on a
  default-less `PArray` that returned whatever was in memory; and a type
  instability on every actor read, which is now 1.8x faster and allocation-free.
- Tree iteration, which re-descended from the root at every step.

### Performance

- A transaction's overhead fell by 61%, and an empty transaction now allocates
  nothing.
- Traversing a `PHeap` — and a `PWSet` or `PWDict` through it — went from 60 MB
  and 10.2 ms per thousand elements to 16.5 KB and 54 µs, in the same order.
- A failed transaction attempt now yields before retrying, rather than
  immediately re-locking the volatiles that just invalidated it. Measured
  neutral in the intended shape (many volatiles, short transactions touching few
  of them) and a fifth to two-fifths fewer wasted attempts under contention.

### Documentation

- A complete manual: a page per collection type, plus the transactional layer,
  task-local variables and utilities, and an API reference covering every
  exported symbol. The build is warning-free, and `checkdocs = :exports` keeps
  it that way.

## [0.1.1] - 2021-10-19

The last release before the rework. Its tree is
`019001d20ce1fb82bd989f60e2368e772eb9e6ca`, which is what the General registry
records for it.
