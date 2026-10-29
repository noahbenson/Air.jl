################################################################################
# cancel.jl
#
# Tests for cancellation tokens and the blocking calls that observe them.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

@testset "cancel" begin
    # Bounded polling: a missed wake-up fails the test instead of hanging it.
    waituntil(f; timeout=30.0) = begin
        t0 = time()
        while !f() && (time() - t0) < timeout
            sleep(0.01)
        end
        f()
    end

    @testset "timeouts" begin
        # A timeout expires with the promise still pending, and leaves no
        # waiter behind.
        p = Promise{Int}()
        e = try
            take(p; timeout = 0.05)
        catch err
            err
        end
        @test e isa ReaderTimeout
        @test e.timeout == 0.05
        @test !isready(p)
        @test p.waiters === nothing
        @test_throws ReaderTimeout take(p; timeout = 0.05)

        # A delivery that beats the timeout wins, and a settled promise never
        # waits, even with a zero timeout.
        q = Promise{Int}()
        Threads.@spawn (sleep(0.02); put!(q, 7))
        @test take(q; timeout = 5.0) == 7
        @test take(q; timeout = 0) == 7

        # A failure that beats the timeout wins as well.
        pf = Promise{Int}()
        Threads.@spawn (sleep(0.02); fail!(pf, :bad))
        @test_throws PromiseFailed take(pf; timeout = 5.0)

        # Nor does an already failed promise wait.
        f = Promise{Int}()
        fail!(f, :bad)
        @test_throws PromiseFailed take(f; timeout = 0)

        # A negative timeout is a programming error.
        @test_throws ArgumentError take(Promise{Int}(), timeout = -1)
    end

    @testset "cancellation" begin
        # A parked waiter wakes with `Cancelled`, and unregisters from both the
        # token and the promise.
        tok = CancelToken()
        @test !iscancelled(tok)
        @test checkcancel(tok) === nothing
        @test checkcancel(nothing) === nothing
        p = Promise{Int}()
        t = Threads.@spawn take(p; token = tok)
        @test waituntil(() -> length(tok.waiters) == 1)
        cancel!(tok, :stop)
        @test iscancelled(tok)
        @test waituntil(() -> istaskdone(t))
        e = t.result
        @test e isa Cancelled
        @test e.token === tok
        @test e.reason === :stop
        @test occursin("cancelled", sprint(showerror, e))
        @test occursin("stop", sprint(showerror, e))
        @test waituntil(() -> isempty(tok.waiters))
        @test p.waiters === nothing

        # Cancelling first is observed before any wait, and cancellation wins
        # even over an already settled promise; `token = nothing` opts out.
        tok2 = CancelToken()
        cancel!(tok2)
        @test_throws Cancelled take(Promise{Int}(), token = tok2)
        s = Promise{Int}()
        put!(s, 3)
        @test_throws Cancelled take(s; token = tok2)
        @test take(s; token = nothing) == 3

        # Cancelling twice changes nothing and is not an error; the first
        # reason is kept.
        @test cancel!(tok2, :again) === tok2
        again = try
            checkcancel(tok2)
        catch x
            x
        end
        @test again isa Cancelled
        @test again.reason === nothing
        @test sprint(showerror, again) == "cancelled (no reason given)"

        # Cancellation that races registration is still observed: a task
        # cancelled immediately after it is spawned never sleeps forever.
        raced = [begin
                     tk = CancelToken()
                     th = Threads.@spawn take(Promise{Int}(); token = tk)
                     cancel!(tk)
                     th
                 end for _ in 1:100]
        @test waituntil(() -> all(istaskdone, raced); timeout = 20.0)
        @test all(th.result isa Cancelled for th in raced)
        @test all(isempty(th.result.token.waiters) for th in raced)
    end

    @testset "cancellation in a transaction" begin
        # A `Cancelled` thrown by the body aborts the transaction and is not
        # retried: nothing it wrote is committed.
        v = Volatile{Int}(0)
        tok = CancelToken()
        cancel!(tok)
        @test_throws Cancelled tx() do
            v[] = 3
            take(Promise{Int}(); token = tok)
        end
        @test v[] == 0
        # Cancellation wins over a settlement staged by the same transaction,
        # so the promise is left pending.
        p = Promise{Int}()
        @test_throws Cancelled tx() do
            put!(p, 1)
            take(p; token = tok)
        end
        @test !isready(p)
    end

    @testset "task-local tokens" begin
        # A task's token is created lazily, stable, and may be replaced or
        # cleared; `take` observes it by default.
        setcancaltoken!(nothing)
        tk = canceltoken()
        @test tk isa CancelToken
        @test canceltoken() === tk
        @test setcancaltoken!(nothing) === nothing
        tk2 = canceltoken()
        @test tk2 !== tk
        try
            setcancaltoken!(tk2)
            cancel!(tk2)
            @test_throws Cancelled take(Promise{Int}())
            # `token = nothing` opts out of the task's cancelled token.
            s = Promise{Int}()
            th = Threads.@spawn take(s; token = nothing)
            @test waituntil(() -> s.waiters !== nothing)
            put!(s, 5)
            @test fetch(th) == 5
        finally
            setcancaltoken!(nothing)
        end
    end

    @testset "sleep_interruptible" begin
        # It returns normally when the sleep finishes...
        t0 = time()
        @test Air.sleep_interruptible(0.02) === nothing
        @test time() - t0 >= 0.015
        # ...and throws `Cancelled` as soon as the token is cancelled.
        tk = CancelToken()
        th = Threads.@spawn Air.sleep_interruptible(30; token = tk)
        @test waituntil(() -> length(tk.waiters) == 1)
        t0 = time()
        cancel!(tk)
        @test waituntil(() -> istaskdone(th); timeout = 5.0)
        @test th.result isa Cancelled
        @test time() - t0 < 5.0
        @test isempty(tk.waiters)
        # An already-cancelled token throws immediately.
        tkm = CancelToken()
        cancel!(tkm)
        @test_throws Cancelled Air.sleep_interruptible(30; token = tkm)
        # `token = nothing` sleeps uninterruptibly.
        t0 = time()
        @test Air.sleep_interruptible(0.01; token = nothing) === nothing
        @test time() - t0 >= 0.005
    end
end
