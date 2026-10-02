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
