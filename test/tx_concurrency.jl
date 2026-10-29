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

    # Bounded polling for the promise tests below, which watch the waiter
    # registry to know that every task is parked: a missed wake-up fails the
    # test instead of hanging it.
    waituntil(f; timeout=30.0) = begin
        t0 = time()
        while !f() && (time() - t0) < timeout
            sleep(0.01)
        end
        f()
    end
    nwaiters(d) = (ws = getfield(d, :waiters); ws === nothing ? 0 : length(ws))

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

    @testset "exactly one transaction delivers a shared promise" begin
        # A promise is one-shot, so among T transactions racing to deliver the
        # same one exactly one must commit and the rest must fail with an
        # `ArgumentError` rather than retry.
        p = Promise{Int}()
        results = Vector{Any}(undef, T)
        tasks = [Threads.@spawn begin
            try
                tx() do
                    put!(p, i)
                end
                results[i] = :ok
            catch e
                results[i] = e
            end
        end for i in 1:T]
        foreach(wait, tasks)
        @test count(==(:ok), results) == 1
        @test count(r -> r isa ArgumentError, results) == T - 1
        @test take(p) in 1:T
    end

    @testset "a two-promise transaction stages in either order" begin
        # Each transaction stages the same two promises, alternating the order.
        # Without the commit's sorted promise-lock group two such transactions
        # could each hold one lock and wait for the other; with it, one commits
        # and the other fails cleanly.
        p1 = Promise{Int}()
        p2 = Promise{Int}()
        results = Vector{Any}(undef, T)
        tasks = [Threads.@spawn begin
            try
                tx() do
                    if isodd(i)
                        put!(p1, i)
                        put!(p2, i)
                    else
                        put!(p2, i)
                        put!(p1, i)
                    end
                end
                results[i] = :ok
            catch e
                results[i] = e
            end
        end for i in 1:T]
        foreach(wait, tasks)
        @test count(==(:ok), results) == 1
        @test count(r -> r isa ArgumentError, results) == T - 1
        @test isready(p1) && isready(p2)
        @test take(p1) == take(p2)
    end

    @testset "deliveries survive volatile contention" begin
        # Each iteration commits a volatile increment and a delivery in the same
        # transaction, so a retry has to discard and restage the delivery. The
        # volatile total is exact and every promise is delivered with its own
        # tag. Callback counts are deliberately not asserted: a retry runs the
        # callback for an attempt that never delivers.
        M = 100
        v = Volatile{Int}(0)
        ps = [[Promise{Int}() for _ in 1:M] for _ in 1:T]
        tasks = [Threads.@spawn begin
            for j in 1:M
                tx() do
                    v[] = v[] + 1
                    put!(ps[i][j], j)
                end
            end
        end for i in 1:T]
        foreach(wait, tasks)
        @test v[] == T * M
        @test all(isready(ps[i][j]) for i in 1:T for j in 1:M)
        @test all(take(ps[i][j]) == j for i in 1:T for j in 1:M)
    end

    @testset "a finalize-phase send happens once per commit" begin
        # Each iteration increments the volatile and delivers a promise whose
        # callback sends to one shared actor. A losing attempt re-runs the
        # callback and restages its send, so the actor must receive exactly one
        # message per committed transaction — one per increment, not one per
        # attempt. The sentinel send at the end is queued behind every committed
        # message, so taking it proves the actor has drained.
        M = 100
        v = Volatile{Int}(0)
        count = Threads.Atomic{Int}(0)
        a = Actor{Int}(0)
        # The sentinel closure is defined before anything is sent. The actor's
        # task is still running when the backlog is drained, and a task cannot
        # call a method defined after it was created: defining the sentinel
        # afterwards makes the drain a `MethodError` inside the actor.
        drained = Promise{Nothing}()
        sentinel = y -> (put!(drained, nothing); y)
        tasks = [Threads.@spawn begin
            for j in 1:M
                p = Promise{Int}(ondelivery = _ -> begin
                    send(y -> (Threads.atomic_add!(count, 1); y), a)
                end)
                tx() do
                    v[] = v[] + 1
                    put!(p, j)
                end
            end
        end for _ in 1:T]
        foreach(wait, tasks)
        send(sentinel, a)
        take(drained)
        @test v[] == T * M
        @test count[] == T * M
    end

    @testset "a commit wakes every parked taker" begin
        # Park T takers on one promise, commit a delivery in a transaction, and
        # every one of them must receive the value and unregister.
        p = Promise{Int}()
        takers = [Threads.@spawn take(p; token = nothing) for _ in 1:T]
        @test waituntil(() -> nwaiters(p) == T)
        tx() do
            put!(p, 9)
        end
        @test all(fetch(t) == 9 for t in takers)
        @test getfield(p, :waiters) === nothing
    end

    @testset "a committed failure wakes every parked taker" begin
        # The same for the failure settlement: every parked taker throws rather
        # than continuing to wait.
        p = Promise{Int}()
        takers = [Threads.@spawn (try
                     take(p; token = nothing)
                 catch e
                     e
                 end) for _ in 1:T]
        @test waituntil(() -> nwaiters(p) == T)
        tx() do
            fail!(p, :boom)
        end
        es = fetch.(takers)
        @test all(e isa PromiseFailed for e in es)
        @test all(e.err === :boom for e in es)
        @test getfield(p, :waiters) === nothing
    end

    @testset "cancellation and a commit wake only their own waiters" begin
        # Half the takers wait under their own token, half under the promise
        # alone. Cancelling releases exactly the first half, while the promise is
        # still pending; failing it then releases the second.
        p = Promise{Int}()
        n = max(2, T ÷ 2)
        toks = [CancelToken() for _ in 1:n]
        cancelled = [Threads.@spawn (try
                         take(p; token = toks[i])
                     catch e
                         e
                     end) for i in 1:n]
        takers = [Threads.@spawn (try
                     take(p; token = nothing)
                 catch e
                     e
                 end) for _ in 1:n]
        @test waituntil(() -> nwaiters(p) == 2n)
        foreach(cancel!, toks)
        @test all(e isa Cancelled for e in fetch.(cancelled))
        @test waituntil(() -> nwaiters(p) == n)
        @test !isready(p)
        tx() do
            fail!(p, :boom)
        end
        @test all(e isa PromiseFailed for e in fetch.(takers))
        @test getfield(p, :waiters) === nothing
    end
end
