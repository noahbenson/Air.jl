# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Under `0.x`, SemVer puts breaking changes in the minor position and compatible
fixes in the patch position, which is how the versions below are chosen.

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
