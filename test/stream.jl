################################################################################
# stream.jl
#
# Tests for Streams: the reusable, bounded, transactional FIFO.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

@testset "stream" begin
    # Bounded polling: a missed wake-up fails the test instead of hanging it.
    waituntil(f; timeout=30.0) = begin
        t0 = time()
        while !f() && (time() - t0) < timeout
            sleep(0.01)
        end
        f()
    end
    nwaiters(s) = (ws = getfield(s, :waiters); ws === nothing ? 0 : length(ws))

    @testset "basics" begin
        s = Stream{Int}()
        @test eltype(s) === Int
        @test eltype(Stream{Int}) === Int
        @test maxsize(s) === nothing
        @test isopen(s)
        @test isempty(s)
        @test length(s) == 0
        @test !isready(s)

        # FIFO order, and `put!` yields the converted value.
        @test put!(s, 1) === 1
        @test put!(s, 2.0) === 2
        @test length(s) == 2
        @test isready(s)
        @test take(s) === 1
        @test take(s) === 2
        @test isempty(s)
        @test_throws InexactError put!(s, 2.5)

        # `Stream()` is a `Stream{Any}`; `maxsize` must be nothing or >= 1.
        a = Stream()
        @test a isa Stream{Any}
        @test put!(a, :x) === :x
        @test take(a) === :x
        @test Stream{Int}() isa Stream{Int}
        @test maxsize(Stream{Int}(maxsize = 2)) == 2
        @test_throws ArgumentError Stream{Int}(maxsize = 0)
        @test_throws ArgumentError Stream(maxsize = -3)
        @test_throws TypeError Stream{Int}(maxsize = 2.5)

        # Streams are handles: identity equality and hashing, so one is a fine
        # dictionary key.
        b = Stream{Int}(maxsize = 2)
        @test isequal(s, s)
        @test !isequal(s, b)
        @test hash(s) == hash(s)
        d = Dict(s => :a, b => :b)
        @test d[s] === :a
        @test d[b] === :b

        # `snapshot` copies what is queued, oldest first, without consuming it,
        # and the copy is independent of the stream.
        q = Stream{Int}(maxsize = 3)
        @test snapshot(q) == Int[]
        @test snapshot(q) isa Vector{Int}
        put!(q, 1)
        put!(q, 2)
        snap = snapshot(q)
        @test snap == [1, 2]
        @test length(snap) == length(q)
        push!(snap, 99)
        empty!(snap)
        @test snapshot(q) == [1, 2]
        @test take(q) == 1
        @test snapshot(q) == [2]
        @test snapshot(q) isa Vector{Int}

        # `show` never blocks or enumerates.
        @test occursin("Stream{Int64}", sprint(show, s))
        @test occursin("open", sprint(show, s))
        @test occursin("maxsize 2", sprint(show, b))
        close(b)
        @test occursin("closed", sprint(show, b))
        @test sprint(show, MIME"text/plain"(), s) == sprint(show, s)
    end

    @testset "blocking" begin
        # A taker parks on an empty stream and is woken by a put.
        s = Stream{Int}()
        t = Threads.@spawn take(s; token = nothing)
        @test waituntil(() -> nwaiters(s) == 1)
        put!(s, 5)
        @test fetch(t) == 5
        @test waituntil(() -> nwaiters(s) == 0)

        # A writer parks on a full stream and is woken by a take.
        b = Stream{Int}(maxsize = 1)
        put!(b, 1)
        w = Threads.@spawn put!(b, 2; token = nothing)
        @test waituntil(() -> nwaiters(b) == 1)
        @test take(b) == 1
        @test fetch(w) === 2
        @test waituntil(() -> nwaiters(b) == 0)
        @test take(b) == 2

        # Both sides time out, and leave no waiter behind.
        e = try
            take(Stream{Int}(); timeout = 0.02)
        catch err
            err
        end
        @test e isa ReaderTimeout
        @test e.timeout == 0.02
        @test sprint(showerror, e) == "reader timed out after 0.02 seconds"
        full = Stream{Int}(maxsize = 1)
        put!(full, 1)
        @test_throws ReaderTimeout put!(full, 2; timeout = 0.02)
        @test nwaiters(full) == 0
        @test_throws ArgumentError take(Stream{Int}(); timeout = -1)
        @test_throws ArgumentError put!(full, 2; timeout = -1)

        # A put that races the timer wins, and a full-but-ready stream never
        # waits even with a zero timeout.
        r = Stream{Int}(maxsize = 1)
        put!(r, 1)
        w2 = Threads.@spawn put!(r, 2; timeout = 5.0)
        @test waituntil(() -> nwaiters(r) == 1)
        take(r)
        @test fetch(w2) === 2
        z = Stream{Int}()
        put!(z, 9)
        @test take(z; timeout = 0) === 9
        @test put!(z, 10; timeout = 0) === 10

        # Every parked waiter is woken by a mutation, and those that lose the
        # race re-park rather than missing a later one.
        m = Stream{Int}()
        takers = [Threads.@spawn take(m; token = nothing) for _ in 1:8]
        @test waituntil(() -> nwaiters(m) == 8)
        for i in 1:8
            put!(m, i)
        end
        @test sort(fetch.(takers)) == 1:8
        @test waituntil(() -> nwaiters(m) == 0)
    end

    @testset "cancellation" begin
        # A parked reader and a parked writer both wake with `Cancelled`.
        tok = CancelToken()
        s = Stream{Int}()
        t = Threads.@spawn take(s; token = tok)
        @test waituntil(() -> length(tok.waiters) == 1)
        cancel!(tok, :stop)
        @test waituntil(() -> istaskdone(t))
        e = t.result
        @test e isa Cancelled
        @test e.reason === :stop
        @test waituntil(() -> nwaiters(s) == 0)
        @test isempty(tok.waiters)

        b = Stream{Int}(maxsize = 1)
        put!(b, 1)
        tok2 = CancelToken()
        w = Threads.@spawn put!(b, 2; token = tok2)
        @test waituntil(() -> length(tok2.waiters) == 1)
        cancel!(tok2)
        @test waituntil(() -> istaskdone(w))
        @test w.result isa Cancelled
        @test take(b) == 1

        # Cancellation wins over an available value, and `token = nothing` opts
        # out of the task's token.
        tok3 = CancelToken()
        cancel!(tok3)
        c = Stream{Int}()
        put!(c, 3)
        @test_throws Cancelled take(c; token = tok3)
        @test_throws Cancelled put!(c, 4; token = tok3)
        @test take(c; token = nothing) == 3
    end

    @testset "close" begin
        s = Stream{Int}()
        @test isopen(s)
        @test close(s) === nothing
        @test !isopen(s)
        # A closed stream is ready: `take` returns or throws without waiting.
        @test isready(s)
        close(s)
        @test !isopen(s)
        @test_throws StreamClosed put!(s, 1)
        @test_throws StreamClosed take(s)
        @test sprint(showerror, StreamClosed()) == "stream is closed"

        # What remains is drained, then `take` throws.
        d = Stream{Int}()
        put!(d, 1)
        put!(d, 2)
        close(d)
        @test length(d) == 2
        @test !isempty(d)
        @test take(d) == 1
        @test take(d) == 2
        @test_throws StreamClosed take(d)

        # Close wakes a parked reader and a parked writer.
        r = Stream{Int}()
        rt = Threads.@spawn (try
                     take(r)
                 catch e
                     e
                 end)
        @test waituntil(() -> nwaiters(r) == 1)
        close(r)
        @test waituntil(() -> istaskdone(rt))
        @test rt.result isa StreamClosed

        w = Stream{Int}(maxsize = 1)
        put!(w, 1)
        wt = Threads.@spawn (try
                     put!(w, 2)
                 catch e
                     e
                 end)
        @test waituntil(() -> nwaiters(w) == 1)
        close(w)
        @test waituntil(() -> istaskdone(wt))
        @test wt.result isa StreamClosed
        @test take(w) == 1
    end

    @testset "transactions" begin
        # Staged puts commit, and vanish on abort.
        s = Stream{Int}()
        tx(() -> put!(s, 1))
        @test take(s) == 1
        try
            tx() do
                put!(s, 2)
                error("abort")
            end
        catch
        end
        @test isempty(s)

        # Staged takes commit, and return their item to the stream on abort.
        put!(s, 3)
        @test tx(() -> take(s)) == 3
        @test isempty(s)
        put!(s, 4)
        try
            tx() do
                @test take(s) == 4
                error("abort")
            end
        catch
        end
        @test take(s) == 4

        # Successive staged takes see successive items, and the transaction
        # reads its own writes.
        put!(s, 5)
        put!(s, 6)
        pair = tx() do
            (take(s), take(s))
        end
        @test pair == (5, 6)
        @test isempty(s)
        ryw = tx() do
            put!(s, 7)
            x = take(s)
            (x, length(s), isempty(s))
        end
        @test ryw == (7, 0, true)
        @test isempty(s)

        # The accessors are staged-aware inside a transaction.
        t = Stream{Int}()
        put!(t, 1)
        view = tx() do
            a = length(t)
            take(t)
            b = length(t)
            c = isempty(t)
            put!(t, 2)
            d = length(t)
            (a, b, c, d, isready(t))
        end
        @test view == (1, 0, true, 1, true)
        @test length(t) == 1
        @test take(t) == 2

        # A snapshot is staged-aware in the same way, and agrees with `length`.
        sn = Stream{Int}()
        put!(sn, 1)
        put!(sn, 2)
        shot = tx() do
            take(sn)
            put!(sn, 3)
            @test length(snapshot(sn)) == length(sn)
            snapshot(sn)
        end
        @test shot == [2, 3]
        @test snapshot(sn) == [2, 3]
        @test take(sn) == 2

        # A staged close commits, and one that aborts leaves the stream open.
        c = Stream{Int}()
        put!(c, 1)
        tx(() -> close(c))
        @test !isopen(c)
        @test take(c) == 1
        @test_throws StreamClosed take(c)
        d = Stream{Int}()
        try
            tx() do
                close(d)
                error("abort")
            end
        catch
        end
        @test isopen(d)
        @test put!(d, 1) === 1

        # A put after a staged close throws, and `isopen` sees the staged close.
        e = Stream{Int}()
        @test_throws StreamClosed tx() do
            close(e)
            @test !isopen(e)
            put!(e, 1)
        end
        @test isopen(e)

        # A transaction that only reads a stream changes nothing, not even its
        # version.
        f = Stream{Int}()
        put!(f, 1)
        v = getfield(f, :version)
        tx() do
            length(f)
            isready(f)
            isopen(f)
            isempty(f)
        end
        @test getfield(f, :version) === v
        @test take(f) == 1

        # A concurrent non-transactional take is a conflict even when the
        # transaction has staged nothing yet: the attempt restarts and both
        # consumers receive distinct items.
        g = Stream{Int}()
        for i in 1:3
            put!(g, i)
        end
        attempts = Threads.Atomic{Int}(0)
        ready = Channel{Nothing}(1)
        go = Channel{Nothing}(1)
        cons = Threads.@spawn begin
            take!(ready)
            take!(go)
            take(g)
        end
        got = tx() do
            attempts[] += 1
            x = take(g)
            if attempts[] == 1
                put!(ready, nothing)
                put!(go, nothing)
                sleep(0.2)
            end
            x
        end
        @test attempts[] == 2
        @test got == 2
        @test fetch(cons) == 1
        @test take(g) == 3
        @test isempty(g)

        # A change that lands between two accesses in one body restarts the
        # attempt as well, so the second staged take cannot index the wrong
        # items.
        h = Stream{Int}()
        for i in 1:4
            put!(h, i)
        end
        n2 = Threads.Atomic{Int}(0)
        ready2 = Channel{Nothing}(1)
        go2 = Channel{Nothing}(1)
        cons2 = Threads.@spawn begin
            take!(ready2)
            take!(go2)
            take(h)
        end
        pair = tx() do
            n2[] += 1
            a = take(h)
            if n2[] == 1
                put!(ready2, nothing)
                put!(go2, nothing)
                sleep(0.2)
            end
            (a, take(h))
        end
        @test n2[] == 2
        @test pair == (2, 3)
        @test fetch(cons2) == 1
        @test take(h) == 4
    end

    @testset "StreamDeadlock" begin
        # A stream created in the transaction cannot be reached by anyone else,
        # so a would-block read or write raises instead of deadlocking.
        @test_throws StreamDeadlock tx() do
            take(Stream{Int}())
        end
        @test_throws StreamDeadlock tx() do
            b = Stream{Int}(maxsize = 1)
            put!(b, 1)
            put!(b, 2)
        end
        @test occursin("transaction", sprint(showerror, StreamDeadlock()))

        # Nothing is written when the guard fires.
        v = Volatile{Int}(0)
        @test_throws StreamDeadlock tx() do
            v[] = 3
            take(Stream{Int}())
        end
        @test v[] == 0

        # But a newborn stream works normally as long as it does not block, and
        # a staged close makes a take throw `StreamClosed`, not `StreamDeadlock`.
        newborn = tx() do
            s = Stream{Int}()
            put!(s, 1)
            (take(s), close(s), isopen(s))
        end
        @test newborn == (1, nothing, false)
        @test_throws StreamClosed tx() do
            s = Stream{Int}()
            close(s)
            take(s)
        end

        # The born marker is per attempt: a stream made in one transaction may
        # block in a later one, and one used after commit outside is fine.
        s = tx(() -> Stream{Int}())
        put!(s, 1)
        @test tx(() -> take(s)) == 1
        # A staged put made in the creating transaction commits with it.
        b = tx() do
            x = Stream{Int}()
            put!(x, 5)
            x
        end
        @test length(b) == 1
        @test take(b) == 5
        @test put!(b, 6) === 6
        @test take(b) == 6
    end

    @testset "stress" begin
        # Several producers and consumers deliver every item exactly once.
        s = Stream{Int}(maxsize = 4)
        nprod = 4
        ncons = 4
        each = 25
        total = nprod * each
        prod = [Threads.@spawn(begin
                    for i in 1:each
                        put!(s, (p - 1) * each + i)
                    end
                end) for p in 1:nprod]
        # The `let` keeps each consumer's buffer its own: a bare `got = Int[]`
        # here would be a boxed local of the enclosing scope, shared by every
        # consumer as soon as the test also mentions `got` outside.
        cons = [Threads.@spawn let got = Int[]
                    while length(got) < total ÷ ncons
                        push!(got, take(s))
                    end
                    got
                end for _ in 1:ncons]
        wait.(prod)
        chunks = fetch.(cons)
        @test sum(length, chunks) == total
        @test sort(reduce(vcat, chunks)) == 1:total
        @test isempty(s)
        @test nwaiters(s) == 0
    end
end
