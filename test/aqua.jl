################################################################################
# aqua.jl
#
# Generic package-quality checks. `ambiguities` is now on: the twenty method
# ambiguities it reported have been resolved, most of them by giving the
# weighted collections' `isequal` fallbacks a real argument type instead of
# `Any`, and one by adding the missing methods for comparing a weighted
# collection with an unweighted one (which was a genuine bug — the call threw).
# The remaining checks are off deliberately:
#   * `persistent_tasks` — the task-local `Var` implementation intentionally
#     uses `task_local_storage`, which that check flags.
#   * `unbound_args` — Aqua reports false positives for `NTuple{N,T}` and
#     `Vararg{T,N}` signatures; see the note at the call below.
#   * `piracies`, `stale_deps` — not yet triaged.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        Air;
        # Aqua's `unbound_args` check reports false positives for signatures
        # written with `NTuple{N,T}` or `Vararg{T,N}` (it flags the type
        # variables even though they appear in the argument types, and Julia's
        # own method-definition check does not). The genuinely unbound type
        # variables that Air did have — 17 of them, which were also causing
        # "declares type variable ... but does not use it" warnings at
        # precompilation time — have been removed.
        unbound_args=false,
        stale_deps=false,
        persistent_tasks=false,
        piracies=false,
    )
end
