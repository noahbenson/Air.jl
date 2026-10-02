################################################################################
# stm_contention.jl
#
# What contention costs a transaction, and what the retry policy does about it.
#
# This exists because the earlier answer to "is the retry policy right?" was
# wrong, and the way it was wrong is easy to repeat. The measurements here are
# built to avoid the three traps that produced it:
#
#  * **Latency is not throughput.** A per-op figure computed as wall time over
#    total operations across T tasks divides by T, which cancels exactly the
#    contention a caller feels. At T=4 the same run gives 1,117 ns/op that way
#    and 4,267 ns/op per task. Both are real; only the second is what a caller
#    waits for. Every number below is per-task latency, with the aggregate
#    printed beside it so the difference is visible rather than assumed.
#
#  * **A contention measurement needs a matched control.** Running more tasks
#    slows a task down for reasons that have nothing to do with the STM, so each
#    workload is also run with every task on its own volatile under the same
#    load. The private column is the floor: if it moves, the measurement is
#    measuring the scheduler.
#
#  * **Conflicts come from sharing, not from task count.** The intended shape is
#    many volatiles with each transaction touching a few, so the sweep varies the
#    volatile count and the per-transaction working set, not just the task count.
#    One volatile written by every task is the worst case, not a typical one, and
#    is included only to bound the other end.
#
# The retry policy: a failed attempt yields before retrying, because the retrying
# transaction would otherwise re-lock the volatiles that just invalidated it. It
# is a scheduler hand-off rather than a sleep — `sleep`'s resolution is a
# millisecond, three orders of magnitude coarser than a transaction. Attempts per
# commit is the metric to read; it is far less noisy than latency and it is what
# the policy actually changes.
#
# Run:  julia --project=bench -t 8 bench/stm_contention.jl
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Air, Random, Printf

const M = 3000                     # transactions per task

# One task's worth of work: short transactions, each touching `k` volatiles drawn
# at random from `vols`. The atomic counter records how many times the body ran,
# which is the attempt count — one per retry, one per commit.
function _worker(vols, seed, k, m, cnt)
    rng = MersenneTwister(seed)
    idx = Vector{Int}(undef, k)
    return @elapsed for _ in 1:m
        for j in 1:k
            idx[j] = rand(rng, 1:length(vols))
        end
        tx() do
            Threads.atomic_add!(cnt, 1)
            for j in 1:k
                x = vols[idx[j]][]
                vols[idx[j]][] = x + 1
            end
        end
    end
end

# Run T tasks over N volatiles, each transaction touching k of them, and report
# per-task latency, attempts per commit, and the aggregate figure that hides it.
# `private` gives every task its own volatile, which is the control: same task
# count, same work, no sharing, so nothing can conflict.
function _cell(N, k, T, m; private = false)
    vols = [Volatile{Int}(0) for _ in 1:(private ? T : N)]
    counts = [Threads.Atomic{Int}(0) for _ in 1:T]
    lat = zeros(T)
    wall = @elapsed begin
        tasks = [Threads.@spawn begin
            vs = private ? [vols[i]] : vols
            lat[$i] = _worker(vs, $i, k, m, counts[$i])
        end for i in 1:T]
        foreach(wait, tasks)
    end
    tot = sum(c[] for c in counts)
    (sum(lat) / T / m * 1e9, tot / (T * m), wall / (T * m) * 1e9)
end

function main()
    # Warm every shape before any of it is measured: an unwarmed first cell pays
    # compilation on each task and reports a latency several times too large.
    for (_, N, k, T) in (("", 1000, 3, 8), ("", 100, 3, 8), ("", 1, 1, 8))
        _cell(N, k, T, 300)
    end

    @printf("threads: %d\n\n", Threads.nthreads())
    @printf("%-26s %14s %12s %14s\n",
            "workload", "per-task ns/op", "attempts/op", "agg ns/op")
    println("-"^70)
    for (nm, N, k, T) in (("many volatiles, small ws", 10_000, 1, 8),
                          ("", 1_000, 1, 8),
                          ("", 1_000, 3, 8),
                          ("", 100, 3, 8),
                          ("one volatile, all writers", 1, 1, 8),
                          ("", 2, 2, 8))
        r = _cell(N, k, T, M)
        @printf("%-26s %14.0f %12.3f %14.0f\n", nm, r[1], r[2], r[3])
    end

    println()
    @printf("%-26s %14s %12s %14s\n",
            "shared vs private control", "per-task ns/op", "attempts/op", "agg ns/op")
    println("-"^70)
    for K in (0, 1, 2, 3), priv in (false, true)
        # the measuring task plus K background writers, all on one volatile when
        # shared; each on its own when private
        r = _cell(1, 1, K + 1, M; private = priv)
        @printf("K=%-2d %-20s %14.0f %12.3f %14.0f\n",
                K, priv ? "private (control)" : "shared", r[1], r[2], r[3])
    end
end

main()
