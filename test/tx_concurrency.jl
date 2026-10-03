################################################################################
# tx_concurrency.jl
#
# The STM under real contention: several threads touching the same volatiles.
#
# `test/TX.jl`'s `tx1` already runs up to eleven workers over a shared array and
# asserts an exact invariant — the root of a summed tree — so the core is
# covered. What is added here are the properties a user relies on, stated
# directly rather than through a larger workload, and exercised at thread counts
# above the four the CI matrix has always used:
#
#  * every increment lands exactly once, so there are no lost updates;
#  * a transaction sees a consistent snapshot across several volatiles, so there
#    are no torn reads — this is what the commit's read-set validation is for;
#  * a transaction that aborts writes nothing, while other transactions commit.
#
# Every loop is bounded and every assertion is exact, so a failure names itself
# and a regression cannot hide behind a loose threshold.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

@testset "STM under contention" begin
    # At least two tasks, so the tests are meaningful even on a single-threaded
    # run; otherwise one task per thread, so the counts rise with the machine and
    # the CI leg that runs eight threads exercises more than the four-thread one.
    T = max(2, Threads.nthreads())

    @testset "no lost updates" begin
        # T tasks each increment one volatile M times. Every increment must land
        # exactly once, so the total is exact: a lost update shows up as a short
        # count, and a duplicated one as a long count.
        M = 200
        v = Volatile{Int}(0)
        tasks = [Threads.@spawn begin
            for _ in 1:M
                tx() do
                    v[] = v[] + 1
                end
            end
        end for _ in 1:T]
        foreach(wait, tasks)
        @test v[] == T * M
    end

    @testset "a consistent snapshot across volatiles" begin
        # Two volatiles that are only ever changed together, in one transaction.
        # A reader that can see one without the other has seen a torn snapshot.
        # This is the property the commit's read-set validation exists to provide,
        # and a read-only transaction must fail validation and retry like any
        # other.
        M = 200
        a = Volatile{Int}(0)
        b = Volatile{Int}(0)
        torn = Threads.Atomic{Int}(0)
        nw = max(1, T ÷ 2)
        writers = [Threads.@spawn begin
            for _ in 1:M
                tx() do
                    a[] = a[] + 1
                    b[] = b[] + 1
                end
            end
        end for _ in 1:nw]
        readers = [Threads.@spawn begin
            for _ in 1:M
                # What matters is what a *committed* transaction observed. A
                # reader may legitimately see a torn pair inside an attempt that
                # is then invalidated and retried — that is what optimistic
                # concurrency permits — so counting a violation from inside the
                # body would count speculative attempts. `tx`'s return value is
                # produced only by an attempt that commits, so counting here
                # counts exactly the observations we care about.
                bad = tx() do
                    a[] != b[]
                end
                bad && Threads.atomic_add!(torn, 1)
            end
        end for _ in 1:(T - nw)]
        foreach(wait, writers)
        foreach(wait, readers)
        @test torn[] == 0
        # and the writers' totals agree with each other, which they only can if
        # every pair of increments committed together
        @test a[] == b[] == nw * M
    end

    @testset "an aborting transaction writes nothing, under contention" begin
        # Half the tasks commit one increment each; the other half stage a large
        # write and then abort. The aborted writes must leave no trace, so the
        # total is exactly what the committing half contributed.
        M = 200
        v = Volatile{Int}(0)
        ncommit = max(1, T ÷ 2)
        nabor = T - ncommit
        committers = [Threads.@spawn begin
            for _ in 1:M
                tx() do
                    v[] = v[] + 1
                end
            end
        end for _ in 1:ncommit]
        aborters = [Threads.@spawn begin
            for _ in 1:M
                try
                    tx() do
                        v[] = v[] + 1000
                        error("abort")
                    end
                catch
                end
            end
        end for _ in 1:nabor]
        foreach(wait, committers)
        foreach(wait, aborters)
        @test v[] == ncommit * M
    end
end
