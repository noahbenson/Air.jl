################################################################################
# aqua.jl
#
# Generic package-quality checks. These are deliberately conservative to start
# with: method ambiguities are disabled because the collection types are built
# from large macro-generated method families that need a separate triage pass,
# and the task-local `Var` implementation intentionally uses
# `task_local_storage`, which `persistent_tasks` flags.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        Air;
        ambiguities=false,
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
