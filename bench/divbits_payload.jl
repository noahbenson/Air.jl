################################################################################
# divbits_payload.jl
#
# One run of the trie-geometry sweep driven by bench/divbits.sh. Reports whether
# the current `PTREE_NODE_SHIFT`/`PTREE_TWIG_SHIFT` constants are feasible, checks
# correctness against `Base`, and times the collection paths that depend on the
# trie's shape. Prints a single machine-readable line so the driver can collect
# the runs into a table.
#
# The constants are compile-time (they determine the occupancy-bit width and the
# node-id encoding), so each combination needs its own Julia process; this file
# is that process's payload.
#
# # Results
#
# The sweep (bench/divbits.sh, 2026-09-24, Julia 1.11, arm64 macOS) over
# node/twig shifts of 4, 5, and 6:
#
#   node twig levels  ok  construct(us)  push(ns)  failing checks
#     4    4     16  yes          428        81
#     4    5     16   NO          402        74  vector-build, vector-pop
#     4    6     16   NO          409        62  vector-build, vector-pop
#     5    4     13   NO          393        74  vector-build, vector-pop
#     5    5     13  yes          396        62
#     5    6     13  yes          395        62
#     6    4     11   NO          413        66  vector-build, vector-pop
#     6    5     11  yes          381        59
#     6    6     11  yes          383        61   (the current default)
#
# Two conclusions. First, there is nothing to win here: the working geometries
# are within a few percent of each other, which is inside this machine's noise,
# so the current 6/6 stands. pcollections' finding that 5 bits beat 6 does not
# carry over — their AMT varies a single `AMT_DIVBITS` with a narrower last
# level, whereas Air's geometry is a node/twig split with a derived root shift.
#
# Second, and more important: four geometries are *silently wrong*, and only on
# the `PArray` path — they corrupt vector construction and `pop` while dicts and
# sets stay correct, and they fail no assert on load. The documented constraint
# (twigs >= 3 and shifts <= 7) does not capture this. Anyone retuning these
# constants gets wrong data rather than an error, which is worth fixing whether
# or not the constants are ever changed again.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Air, BenchmarkTools, Printf

const A = Air

# The timed work. These live at module scope because `@benchmark` runs the
# expression it is given in the enclosing module, not in the caller's frame.
t_construct() = begin
    d = PDict{Symbol,Int}()
    for i in 1:1000
        d = push(d, Symbol("k", i) => i)
    end
    d
end
t_push() = begin
    v = PVector{Int}()
    for i in 1:1000
        v = push(v, i)
    end
    v
end

function main()
    # A node's depth is stored in the low `PTREE_TWIG_SHIFT` bits of its id, so
    # the tree must not need more distinct depths than that field can hold.
    if A.PTREE_LEVELS > (1 << A.PTREE_TWIG_SHIFT)
        @printf(
            "SHIFTS node=%d twig=%d feasible=false levels=%d depthfield=%d\n",
            A.PTREE_NODE_SHIFT,
            A.PTREE_TWIG_SHIFT,
            A.PTREE_LEVELS,
            1 << A.PTREE_TWIG_SHIFT,
        )
        return
    end

    # Correctness: a moderate dict with deletions, a set, a vector, and a
    # matrix, each against `Base`. A geometry that mis-encodes keys or depths
    # tends to show up here rather than in the timings, so each check is
    # labelled: knowing *which* one breaks is what identifies the constraint the
    # constants have to satisfy.
    n = 2000
    fails = String[]
    check(label, cond) = cond || push!(fails, label)

    d = PDict{Symbol,Int}()
    for i in 1:n
        d = push(d, Symbol("k", i) => i)
    end
    check("dict-build", length(d) == n)
    check("dict-read", all(d[Symbol("k", i)] == i for i in 1:n))
    for i in 1:2:n
        d = delete(d, Symbol("k", i))
    end
    check("dict-delete-len", length(d) == n ÷ 2)
    check("dict-delete-read", all(d[Symbol("k", i)] == i for i in 2:2:n))
    check("dict-delete-gone", all(!haskey(d, Symbol("k", i)) for i in 1:2:n))

    s = PSet(Symbol("k", i) for i in 1:n)
    check("set-build", length(s) == n)
    check("set-in", all(Symbol("k", i) in s for i in 1:n))

    v = PVector(collect(1:n))
    check("vector-build", all(v[i] == i for i in 1:n))
    v = push(v, -1)
    check("vector-push", v[n + 1] == -1)
    check("vector-pop", pop(v)[n] == n)

    m = PMatrix(0.0, (40, 50))
    m = setindex(m, 3.5, 7, 9)
    check("matrix-setindex", m[7, 9] == 3.5)

    ok = isempty(fails)

    t_construct()
    t_push()
    tt = @benchmark t_construct() samples = 100 seconds = 1
    tp = @benchmark t_push() samples = 100 seconds = 1

    @printf(
        "SHIFTS node=%d twig=%d feasible=true ok=%s fails=%s levels=%d construct=%.0f push=%.0f\n",
        A.PTREE_NODE_SHIFT,
        A.PTREE_TWIG_SHIFT,
        ok,
        join(fails, ","),
        A.PTREE_LEVELS,
        minimum(tt).time / 1000,
        minimum(tp).time / 1000,
    )
end

main()
