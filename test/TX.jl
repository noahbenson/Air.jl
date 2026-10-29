# Tests for the Transaction system.
# Author: Noah C. Benson <n@nben.net>

@testset "TX" begin
    @testset "tx1" begin
        vals = Int[]
        vols = Vector{Volatile{Int}}[]
        log = Actor{typeof(stdout)}(stdout)
        worker() = begin
            global vals
            global vols
            while vols[end][1][] < 0
                # Pick a random layer and cell to calculate
                l = rand((2, 3, 4, 5))
                layer = vols[l]
                k = rand(1:length(layer))
                # synconize for the rest.
                tx() do
                    #println("$l $k $(current_tx() === nothing) $(getfield(_current_tx, :values)))")
                    res = layer[k]
                    if res[] < 0
                        up = vols[l - 1]
                        a = up[k * 2 - 1]
                        b = up[k * 2]
                        if a[] >= 0 && b[] >= 0
                            res[] = a[] + b[]
                            send(log) do log
                                println(
                                    "Layer $l, item $k complete $(objectid(current_task())).",
                                )
                                log
                            end
                        end
                    end
                end
            end
            send(log) do log
                println("Ending Woker $(objectid(current_task())).")
                log
            end
        end
        tx1_run(n::Int) = begin
            global vals
            global vols
            vals = rand(0:1000, 16)
            vols = Vector{Volatile{Int}}[
                [Volatile{Int}(u) for u in vals],
                [Volatile{Int}(-1) for _ in 1:8],
                [Volatile{Int}(-1) for _ in 1:4],
                [Volatile{Int}(-1) for _ in 1:2],
                [Volatile{Int}(-1)],
            ]
            # Spawn n threads.
            threads = [(Threads.@spawn worker()) for _ in 1:n]
            for th in threads
                Threads.wait(th)
            end
            println("  ---->  $(vols[end][1][]) == $(sum(vals))")
            return (vols[end][1][] == sum(vals))
        end
        @test tx1_run(1)
        @test tx1_run(2)
        @test tx1_run(3)
        @test tx1_run(4)
        @test tx1_run(5)
        @test tx1_run(6)
        @test tx1_run(7)
        @test tx1_run(8)
        @test tx1_run(9)
        @test tx1_run(10)
        @test tx1_run(11)
    end

    # The `tx1` testset above only ever exercises successful transactions. The
    # tests below cover the failure modes that make a transaction system
    # transactional: abort/rollback, conflict detection, and retry.
    @testset "STM semantics" begin
        @testset "commit and out-of-transaction writes" begin
            v = Volatile{Int}(0)
            @tx v[] = 1
            @test v[] == 1
            # A Volatile may only be written inside a transaction.
            @test_throws ErrorException (v[] = 2)
            @test v[] == 1
        end

        @testset "aborting a transaction discards its writes" begin
            v = Volatile{Int}(0)
            w = Volatile{Int}(0)
            @test_throws ErrorException tx() do
                v[] = 5
                w[] = 6
                error("abort")
            end
            @test v[] == 0
            @test w[] == 0
        end

        @testset "nested transactions fold into the outer one" begin
            v = Volatile{Int}(0)
            @tx begin
                v[] = 1
                tx(() -> (v[] = 2))
            end
            @test v[] == 2
        end

        @testset "an explicit TxRetryException restarts the transaction" begin
            v = Volatile{Int}(0)
            n = Ref(0)
            tx() do
                n[] += 1
                # Throw once, then succeed on the retry.
                (n[] == 1) && throw(TxRetryException())
                v[] = 7
            end
            @test n[] == 2
            @test v[] == 7
        end

        @testset "the pooled transaction releases the attempt's data" begin
            # `tx_clear!` runs on the way out of an attempt as well as at the
            # top of the next one. Without the trailing clear the task's pooled
            # `Transaction` would keep the last attempt's reads, writes and
            # staged messages — and through them whatever values they hold —
            # reachable until the task happened to run another transaction.
            v = Volatile{Vector{Int}}(Int[])
            w = Volatile{Int}(0)
            p = Promise{Int}()
            tx() do
                w[] = length(v[])
                put!(p, 2)
            end
            @test take(p) == 2
            t = Air._transaction()
            @test t.state === :running
            @test isempty(t.rvolatiles)
            @test isempty(t.wvolatiles)
            @test isempty(t.actors)
            @test isempty(getfield(t, :promises))
            # An error escaping the body releases the staged data too.
            @test_throws ErrorException tx() do
                w[] = length(v[])
                error("abort")
            end
            t = Air._transaction()
            @test t.state === :running
            @test isempty(t.rvolatiles)
            @test isempty(t.wvolatiles)
        end

        @testset "a conflicting write forces a retry" begin
            # Task A reads `v` and then blocks. The main task commits a write to
            # `v`, invalidating A's read set. When A resumes it must detect the
            # conflict, retry, observe the new value, and commit a write derived
            # from it. Ordering is established with channels rather than sleeps,
            # so the test does not depend on scheduler timing.
            v = Volatile{Int}(0)
            started = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            attempts = Threads.Atomic{Int}(0)
            t = Threads.@spawn tx() do
                # `atomic_add!` returns the previous value, so `n` is this
                # attempt's 1-based number. The handshake happens only on the
                # first attempt: on a retry the main task is no longer waiting
                # to rendezvous, so repeating it would deadlock.
                n = Threads.atomic_add!(attempts, 1) + 1
                x = v[]
                if n == 1
                    put!(started, nothing)
                    take!(resume)
                end
                v[] = x + 1
            end
            take!(started)
            tx() do
                v[] = 100
            end
            put!(resume, nothing)
            wait(t)
            # Exactly two: the handshake makes the first attempt fail for certain
            # and nothing writes `v` again after the main task's commit, so the
            # second must succeed. Asserting the exact count rather than `>= 2`
            # pins the retry policy — an extra attempt here is a regression.
            @test attempts[] == 2
            @test v[] == 101
        end

        @testset "Actor" begin
            # The actor runs its functions on an internally spawned task, so we
            # poll for the observable outcome with a bounded timeout rather than
            # sleeping a fixed amount.
            waituntil(f; timeout=30.0) = begin
                t0 = time()
                while !f() && (time() - t0) < timeout
                    sleep(0.01)
                end
                f()
            end

            a = Actor{Symbol}(:start)
            @test a[] == :start
            send(a) do _x
                :finished
            end
            @test waituntil(() -> a[] === :finished)
            @test geterror(a) === nothing

            # An error inside the actor is captured: the actor enters an error
            # state and reading it rethrows the captured exception. `reset`
            # restores it to a usable state with a new value.
            b = Actor{Symbol}(:start)
            send(b) do _x
                error("actor failure")
            end
            @test waituntil(() -> geterror(b) !== nothing)
            @test geterror(b) isa Air.ActorException
            @test_throws Air.ActorException b[]
            reset(b, :ok)
            @test b[] == :ok
            @test geterror(b) === nothing

            # Post-processing: the hook is given (old value, raw result), and
            # the value it returns is what the actor stores. The next message
            # therefore sees the post-processed value.
            seen = Channel{Tuple{Int,Int}}(4)
            pp = Actor{Int}(0; postprocess = (old, r) -> begin
                put!(seen, (old, r))
                old + r
            end)
            send(_ -> 3, pp)
            @test take!(seen) == (0, 3)
            send(_ -> 4, pp)
            @test take!(seen) == (3, 4)
            @test waituntil(() -> pp[] == 7)
            @test getpostprocess(pp) !== nothing
            @test geterrorhandler(pp) === nothing

            # Hooks default to nothing, and a setter affects the next message
            # dequeued and can be cleared again.
            plain = Actor{Int}(0)
            @test getpostprocess(plain) === nothing
            @test geterrorhandler(plain) === nothing
            hits = Channel{Int}(4)
            setpostprocess!(plain, (old, r) -> (put!(hits, old); r + 100))
            send(_ -> 1, plain)
            @test take!(hits) == 0
            @test waituntil(() -> plain[] == 101)
            setpostprocess!(plain, nothing)
            @test getpostprocess(plain) === nothing
            send(_ -> 5, plain)
            @test waituntil(() -> plain[] == 5)

            # An error handler is given (actor, value, exception); the actor
            # still enters its error state, and `reset` restores it without
            # clearing the hook.
            errs = Channel{Any}(1)
            eh = Actor{Int}(0; errorhandler = (act, val, e) -> put!(errs, (act, val, e)))
            send(_ -> error("boom"), eh)
            @test waituntil(() -> geterror(eh) !== nothing)
            act0, val0, e0 = take!(errs)
            @test act0 === eh
            @test val0 == 0
            @test e0 isa ErrorException
            @test geterror(eh) isa Air.ActorException
            reset(eh, 3)
            @test eh[] == 3
            @test geterror(eh) === nothing
            @test geterrorhandler(eh) !== nothing

            # A post-processor that throws is an actor error: the handler sees
            # the post-processor's exception, while `.fn` remains the scheduled
            # function and `.value` the value it was given.
            errs2 = Channel{Any}(1)
            ppf = Actor{Int}(0;
                postprocess = (old, r) -> error("hook boom"),
                errorhandler = (act, val, e) -> put!(errs2, (val, e)),
            )
            f = _ -> 1
            send(f, ppf)
            @test waituntil(() -> geterror(ppf) !== nothing)
            val1, e1 = take!(errs2)
            @test val1 == 0
            @test e1 isa ErrorException
            @test geterror(ppf).fn === f
            @test geterror(ppf).value == 0

            # The hooks are private fields, not properties.
            @test propertynames(pp) == (:state, :value, :error)
            @test_throws ErrorException pp.postprocess

            # Post-processing serialises concurrent sends: one hook call per
            # message means N concurrent increments total exactly N.
            N = 8
            counter = Actor{Int}(0; postprocess = (old, _) -> old + 1)
            for _ in 1:N
                Threads.@spawn send(_ -> 0, counter)
            end
            @test waituntil(() -> counter[] == N)
            @test counter[] == N

            # A message defined after the actor's task was created still runs.
            # The task is kept alive by a blocking message, so `latef` is
            # compiled in a newer world than the one that task is running in; a
            # direct call would fail with "the applicable method may be too
            # new", exactly as a closure typed at the REPL during a backlog
            # would. This is why the actor calls its functions with
            # `invokelatest`.
            wa = Actor{Int}(0)
            started = Channel{Nothing}(1)
            gate = Channel{Nothing}(1)
            send(_ -> (put!(started, nothing); take!(gate); 0), wa)
            @test take!(started) === nothing
            done = Promise{Nothing}()
            latef = @eval x -> (put!($done, nothing); x + 1)
            send(latef, wa)
            put!(gate, nothing)
            @test waituntil(() -> isready(done))
            @test geterror(wa) === nothing
            @test wa[] == 1
        end

        @testset "Transactional promises" begin
            # A delivery is staged: it happens only if the transaction commits.
            # The callback runs in the finalize phase, and `take`/`isready` see
            # the staged value from inside the transaction (read-your-writes).
            log = Symbol[]
            p = Promise{Int}(ondelivery = _ -> push!(log, :delivery))
            res = tx() do
                put!(p, 4)
                @test take(p) == 4
                @test isready(p)
                :committed
            end
            @test res === :committed
            @test isready(p)
            @test take(p) == 4
            @test log == [:delivery]

            # An aborted body delivers nothing and runs no callback.
            log2 = Symbol[]
            q = Promise{Int}(ondelivery = _ -> push!(log2, :delivery))
            @test_throws ErrorException tx() do
                put!(q, 1)
                error("abort")
            end
            @test !isready(q)
            @test isempty(log2)

            # A throwing callback aborts the transaction and rolls back its
            # volatile writes.
            v = Volatile{Int}(0)
            r = Promise{Int}(ondelivery = _ -> error("callback boom"))
            @test_throws ErrorException tx() do
                v[] = 9
                put!(r, 1)
            end
            @test v[] == 0
            @test !isready(r)

            # Volatile finalizers run before delivery callbacks.
            order = Symbol[]
            v2 = Volatile{Int}(0, nothing, x -> (push!(order, :finalize); x))
            s = Promise{Int}(ondelivery = _ -> push!(order, :delivery))
            tx() do
                v2[] = 5
                put!(s, 1)
            end
            @test v2[] == 5
            @test order == [:finalize, :delivery]

            # Two deliveries in one transaction are rejected immediately, and
            # so is a delivery of a promise that has already been delivered.
            t = Promise{Int}()
            @test_throws ArgumentError tx() do
                put!(t, 1)
                put!(t, 2)
            end
            @test !isready(t)
            u = Promise{Int}()
            put!(u, 1)
            @test_throws ArgumentError tx() do
                put!(u, 2)
            end

            # A callback may start its own transaction. The transaction that is
            # finalizing must not hand out its pooled `Transaction` object, or
            # the nested one would clear it mid-commit.
            w = Volatile{Int}(0)
            z = Volatile{Int}(0)
            y = Promise{Int}(ondelivery = _ -> tx(() -> (z[] = 7)))
            tx() do
                w[] = 3
                put!(y, 1)
            end
            @test w[] == 3
            @test z[] == 7
            @test take(y) == 1
        end

        @testset "a delivery callback runs on every attempt" begin
            # The first attempt is invalidated while its callback is blocked in
            # the finalize phase; the retry runs the callback again, but the
            # promise is delivered only once.
            v = Volatile{Int}(0)
            calls = Threads.Atomic{Int}(0)
            entered = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            p = Promise{Int}(ondelivery = _ -> begin
                n = Threads.atomic_add!(calls, 1) + 1
                if n == 1
                    put!(entered, nothing)
                    take!(resume)
                end
            end)
            t = Threads.@spawn begin
                take!(entered)
                tx() do
                    v[] = 1
                end
                put!(resume, nothing)
            end
            tx() do
                x = v[]
                put!(p, x)
            end
            wait(t)
            @test calls[] == 2
            @test v[] == 1
            @test take(p) == 1
        end

        @testset "a losing delivery does not retry and commits nothing" begin
            # Another task delivers `lost` while this transaction's callback is
            # blocked. The conflict is permanent — a one-shot promise cannot be
            # delivered twice — so the commit throws rather than retrying, and
            # the volatile write is not committed.
            v = Volatile{Int}(0)
            lost = Promise{Int}()
            entered = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            gate = Promise{Int}(ondelivery = _ -> begin
                put!(entered, nothing)
                take!(resume)
            end)
            t = Threads.@spawn begin
                take!(entered)
                put!(lost, 7)
                put!(resume, nothing)
            end
            @test_throws ArgumentError tx() do
                v[] = 5
                put!(gate, 1)
                put!(lost, 2)
            end
            wait(t)
            @test v[] == 0
            @test take(lost) == 7
            @test !isready(gate)
        end

        @testset "transactional promise failures" begin
            # A failure is staged exactly like a delivery: the body sees it
            # through the read-your-writes functions, and it happens only if the
            # transaction commits.
            log = Symbol[]
            p = Promise{Int}(onfailure = _ -> push!(log, :failed))
            res = tx() do
                fail!(p, :boom)
                @test isready(p)
                @test isfailed(p)
                @test geterror(p) === :boom
                @test_throws PromiseFailed take(p)
                :committed
            end
            @test res === :committed
            @test isready(p)
            @test isfailed(p)
            @test geterror(p) === :boom
            @test_throws PromiseFailed take(p)
            @test log == [:failed]

            # Aborting leaves the promise pending and runs no failure callback.
            log2 = Symbol[]
            q = Promise{Int}(onfailure = _ -> push!(log2, :failed))
            @test_throws ErrorException tx() do
                fail!(q, :nope)
                error("abort")
            end
            @test !isready(q)
            @test !isfailed(q)
            @test geterror(q) === nothing
            @test isempty(log2)

            # A throwing failure callback aborts the transaction like a throwing
            # delivery callback does.
            v = Volatile{Int}(0)
            r = Promise{Int}(onfailure = _ -> error("failure callback boom"))
            @test_throws ErrorException tx() do
                v[] = 4
                fail!(r, :x)
            end
            @test v[] == 0
            @test !isready(r)

            # A promise settles once, whatever the kind of settlement.
            m = Promise{Int}()
            @test_throws ArgumentError tx() do
                put!(m, 1)
                fail!(m, 2)
            end
            @test !isready(m)
            n = Promise{Int}()
            @test_throws ArgumentError tx() do
                fail!(n, 1)
                put!(n, 2)
            end
            @test !isready(n)
            already = Promise{Int}()
            fail!(already, :settled)
            @test_throws ArgumentError tx() do
                fail!(already, :again)
            end
            @test_throws ArgumentError tx() do
                put!(already, 1)
            end
            @test geterror(already) === :settled

            # `Promise{Any}` is the case where the delivered and failed slots
            # have the same type; staging it must not be ambiguous.
            ap = Promise()
            tx() do
                fail!(ap, :any)
            end
            @test geterror(ap) === :any
            @test_throws PromiseFailed take(ap)
        end

        @testset "a failure callback runs on every attempt" begin
            # Same contract as a delivery callback: it runs in the finalize
            # phase of every attempt, but the promise is failed only once.
            v = Volatile{Int}(0)
            calls = Threads.Atomic{Int}(0)
            entered = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            p = Promise{Int}(onfailure = _ -> begin
                n = Threads.atomic_add!(calls, 1) + 1
                if n == 1
                    put!(entered, nothing)
                    take!(resume)
                end
            end)
            t = Threads.@spawn begin
                take!(entered)
                tx() do
                    v[] = 1
                end
                put!(resume, nothing)
            end
            tx() do
                fail!(p, v[])
            end
            wait(t)
            @test calls[] == 2
            @test v[] == 1
            @test geterror(p) === 1
        end

        @testset "a losing failure does not retry and commits nothing" begin
            # The failure conflict is permanent too: the promise was settled by
            # another task during finalize, so the commit throws and the volatile
            # write is not committed.
            v = Volatile{Int}(0)
            lost = Promise{Int}()
            entered = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            gate = Promise{Int}(onfailure = _ -> begin
                put!(entered, nothing)
                take!(resume)
            end)
            t = Threads.@spawn begin
                take!(entered)
                fail!(lost, 7)
                put!(resume, nothing)
            end
            @test_throws ArgumentError tx() do
                v[] = 5
                fail!(gate, 1)
                fail!(lost, 2)
            end
            wait(t)
            @test v[] == 0
            @test geterror(lost) === 7
            @test !isready(gate)
        end

        @testset "actor sends in the finalize phase" begin
            # A send from a volatile finalizer or a promise delivery callback
            # runs with no transaction active (`currtx() === nothing`), but it
            # is staged in the transaction that is finalizing: it is queued if
            # and only if that transaction commits, and once rather than once
            # per attempt. The actor's task runs asynchronously, so each case
            # flushes the queue with a sentinel and waits for it; messages are
            # processed in order, so when the sentinel's promise is delivered
            # every earlier message has been processed.
            # The sentinel is one closure, defined before any actor in this
            # testset exists. A task cannot call a method defined after it was
            # created, and an actor's task may still be draining its queue when
            # the sentinel is sent.
            drainref = Ref{Promise{Nothing}}()
            drainmsg = y -> begin
                put!(drainref[], nothing)
                y
            end
            drain!(a) = begin
                drainref[] = Promise{Nothing}()
                send(drainmsg, a)
                take(drainref[])
            end

            # A volatile finalizer's send commits with the transaction.
            a1 = Actor{Int}(0)
            hits1 = Threads.Atomic{Int}(0)
            v1 = Volatile{Int}(0, nothing, x -> begin
                send(y -> (Threads.atomic_add!(hits1, 1); y), a1)
                x
            end)
            tx() do
                v1[] = 5
            end
            @test v1[] == 5
            drain!(a1)
            @test hits1[] == 1

            # It is discarded if the transaction cannot commit. Another task
            # delivers `lost` while the finalizer is parked, which is the
            # permanent conflict; the volatile write and the staged send both
            # vanish.
            a2 = Actor{Int}(0)
            hits2 = Threads.Atomic{Int}(0)
            entered = Channel{Nothing}(1)
            resume = Channel{Nothing}(1)
            lost = Promise{Int}()
            v2 = Volatile{Int}(0, nothing, x -> begin
                send(y -> (Threads.atomic_add!(hits2, 1); y), a2)
                put!(entered, nothing)
                take!(resume)
                x
            end)
            t = Threads.@spawn begin
                take!(entered)
                put!(lost, 7)
                put!(resume, nothing)
            end
            @test_throws ArgumentError tx() do
                v2[] = 1
                put!(lost, 2)
            end
            wait(t)
            @test v2[] == 0
            @test take(lost) == 7
            drain!(a2)
            @test hits2[] == 0

            # A failed attempt's send is discarded and restaged, so a retried
            # transaction sends once: the finalizer runs twice and only the
            # committed attempt queues a message.
            a3 = Actor{Int}(0)
            hits3 = Threads.Atomic{Int}(0)
            nfin = Threads.Atomic{Int}(0)
            entered3 = Channel{Nothing}(1)
            resume3 = Channel{Nothing}(1)
            gate3 = Volatile{Int}(0)
            v3 = Volatile{Int}(0, nothing, x -> begin
                send(y -> (Threads.atomic_add!(hits3, 1); y), a3)
                n = Threads.atomic_add!(nfin, 1) + 1
                if n == 1
                    put!(entered3, nothing)
                    take!(resume3)
                end
                x
            end)
            t3 = Threads.@spawn begin
                take!(entered3)
                tx() do
                    gate3[] = gate3[] + 10
                end
                put!(resume3, nothing)
            end
            tx() do
                gate3[]                  # read, so the other write invalidates it
                v3[] = 1
            end
            wait(t3)
            @test nfin[] == 2
            @test v3[] == 1
            drain!(a3)
            @test hits3[] == 1

            # The same for a delivery callback.
            a4 = Actor{Int}(0)
            hits4 = Threads.Atomic{Int}(0)
            p4 = Promise{Int}(ondelivery = _ -> begin
                send(y -> (Threads.atomic_add!(hits4, 1); y), a4)
            end)
            tx() do
                put!(p4, 1)
            end
            drain!(a4)
            @test hits4[] == 1

            # ... including a callback whose transaction loses its race.
            a5 = Actor{Int}(0)
            hits5 = Threads.Atomic{Int}(0)
            entered5 = Channel{Nothing}(1)
            resume5 = Channel{Nothing}(1)
            lost5 = Promise{Int}()
            p5 = Promise{Int}(ondelivery = _ -> begin
                send(y -> (Threads.atomic_add!(hits5, 1); y), a5)
                put!(entered5, nothing)
                take!(resume5)
            end)
            t5 = Threads.@spawn begin
                take!(entered5)
                put!(lost5, 7)
                put!(resume5, nothing)
            end
            @test_throws ArgumentError tx() do
                put!(p5, 1)
                put!(lost5, 2)
            end
            wait(t5)
            @test !isready(p5)
            drain!(a5)
            @test hits5[] == 0

            # Sends from the body and from the finalize phase reach the actor in
            # that order: both are enqueued into the same staged queue.
            order = Symbol[]
            a6 = Actor{Int}(0)
            v6 = Volatile{Int}(0, nothing, x -> begin
                send(y -> (push!(order, :finalize); y), a6)
                x
            end)
            tx() do
                send(y -> (push!(order, :body); y), a6)
                v6[] = 1
            end
            drain!(a6)
            @test order == [:body, :finalize]

            # A nested transaction started by a callback is an independent
            # transaction, so its send commits even when the outer one loses.
            a7 = Actor{Int}(0)
            hits7 = Threads.Atomic{Int}(0)
            entered7 = Channel{Nothing}(1)
            resume7 = Channel{Nothing}(1)
            lost7 = Promise{Int}()
            p7 = Promise{Int}(ondelivery = _ -> begin
                tx() do
                    send(y -> (Threads.atomic_add!(hits7, 1); y), a7)
                end
                put!(entered7, nothing)
                take!(resume7)
            end)
            t7 = Threads.@spawn begin
                take!(entered7)
                put!(lost7, 7)
                put!(resume7, nothing)
            end
            @test_throws ArgumentError tx() do
                put!(p7, 1)
                put!(lost7, 2)
            end
            wait(t7)
            drain!(a7)
            @test hits7[] == 1

            # Sending to an actor in an error state throws, and from a finalizer
            # that exception aborts the transaction.
            bad = Actor{Int}(0)
            send(_ -> error("boom"), bad)
            t0 = time()
            while (geterror(bad) === nothing) && ((time() - t0) < 30.0)
                sleep(0.01)
            end
            @test geterror(bad) isa Air.ActorException
            v8 = Volatile{Int}(0, nothing, x -> begin
                send(y -> y, bad)
                x
            end)
            @test_throws Air.ActorException tx() do
                v8[] = 3
            end
            @test v8[] == 0
        end
    end

    @testset "TxIO and airout" begin
        # `TxIO` writes through an `Actor`, so concurrent writers are serialised —
        # each `write` is one message. And because a send inside a transaction is
        # itself transactional (the message is staged and released only at commit),
        # output written inside a transaction appears only if it succeeds. That is
        # the contract its docstring states, and the abort case below is the part
        # of it nothing else covers.
        buf = IOBuffer()
        io = TxIO(buf)
        content() = String(take!(copy(buf)))
        # The writes below are asynchronous — they are performed by the TxIO's
        # actor task — so each has to be waited for. (`waituntil` is defined
        # inside the neighbouring testset, so this one defines its own.)
        waitio(f; timeout = 10.0) = begin
            t0 = time()
            while !f() && (time() - t0) < timeout
                sleep(0.005)
            end
            f()
        end
        @test iswritable(io)
        @test !isreadable(io)
        @test eof(io) == eof(buf)
        @test_throws ErrorException TxIO(IOBuffer(read=true, write=false))

        # outside a transaction, a write goes through
        write(io, "direct")
        @test waitio(() -> content() == "direct")
        # inside a transaction that commits
        tx() do
            write(io, "-committed")
        end
        @test waitio(() -> content() == "direct-committed")
        # inside a transaction that aborts, it must never be sent. The check is
        # made deterministic rather than by waiting and hoping: a marker is
        # written *after* the aborted attempt and waited for. The actor processes
        # its queue in order, so once the marker has appeared, anything the
        # aborted transaction had sent would have appeared too.
        try
            tx() do
                write(io, "-ABORTED")
                error("abort")
            end
        catch
        end
        write(io, "-marker")
        @test waitio(() -> content() == "direct-committed-marker")
        @test !occursin("ABORTED", content())

        # the rest of the IO interface
        print(io, "-print")
        println(io, "-line")
        @test waitio(() -> content() == "direct-committed-marker-print-line\n")

        # concurrent writers do not interleave, since one write is one message
        if Threads.nthreads() >= 2
            n = 20
            @sync for i in 1:n
                Threads.@spawn write(io, "|$(i)")
            end
            @test waitio(() -> all(occursin("|$(i)", content()) for i in 1:n))
        end

        # `airout` is the `stdout` counterpart
        @test iswritable(airout)
        @test !isreadable(airout)
        @test Air.AirOut() isa Air.AirOut
    end
end
