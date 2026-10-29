################################################################################
# util.jl
#
# Tests for the Air utilities.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

@testset "util" begin
    @testset "Delay" begin
        k = 1
        d = Delay{Int}(() -> (k += 1; k))
        @test !isready(d)
        @test_throws ArgumentError (d[] = 20)
        @test d[] == 2
        @test d[] == 2
        @test isready(d)
        d2 = Delay{Int}(() -> 2)
        d3 = Delay{Int}(() -> (k += 1; k))
        @test hash(d) == hash(d2)
        @test hash(d) != hash(d3)
        # An already-computed delay is ready immediately.
        d4 = Delay{Int}(5)
        @test isready(d4)
        @test d4[] == 5
    end

    @testset "@delay" begin
        counter = Ref(0)
        d = @delay begin
            counter[] += 1
            counter[]
        end
        @test d[] == 1
        @test d[] == 1
        @test counter[] == 1
        # The macro used to require an `Expr`: a literal was rejected outright
        # ("no method matching @delay(::LineNumberNode, ::Module, ::Int64)"),
        # and a block was unwrapped by index without skipping `LineNumberNode`s,
        # which threw. All of these forms should work.
        @test (@delay 10)[] == 10
        @test (@delay (2 + 3))[] == 5
        @test (@delay () -> 11)[] == 11
        @test (@delay () -> 12.0::Float64)[] === 12.0
        @test (@delay begin
            20
        end)[] == 20
        @test (@delay begin
            x = 1
            x + 1
        end)[] == 2
    end

    @testset "@memoize" begin
        @testset "untyped arguments" begin
            # Regression: an argument without a type annotation arrives at the
            # macro as a plain Symbol, which `_memoize_fixarg` did not accept.
            @eval module MemoizeUntyped
                using Air
                @memoize memf(x) = x + 1
            end
            @test MemoizeUntyped.memf(1) == 2
            @test MemoizeUntyped.memf(1) == 2
        end

        @testset "type-tagged result" begin
            # Regression: the type tag was read from the expression head (the
            # Symbol `::`) rather than its type argument.
            @eval module MemoizeTagged
                using Air
                @memoize memg(x::Int) = (x + 1.0)::Float64
            end
            @test MemoizeTagged.memg(1) === 2.0
        end

        @testset "values are computed once" begin
            @eval module MemoizeOnce
                using Air
                const COUNT = Ref(0)
                @memoize memh(x::Int) = begin
                    COUNT[] += 1
                    x * 10
                end
            end
            @test MemoizeOnce.memh(3) == 30
            @test MemoizeOnce.memh(3) == 30
            @test MemoizeOnce.COUNT[] == 1
        end
    end

    @testset "Promise" begin
        p = Promise{Int}()
        @test !isready(p)
        Threads.@spawn put!(p, 7)
        @test take(p) == 7
        @test isready(p)

        # A delivery callback is fixed at construction and runs exactly once,
        # after the promise's lock has been released. Its return value is
        # ignored: it is a notification, not a transformation.
        seen = Ref{Any}(nothing)
        waslocked = Ref(true)
        q = Promise{Int}(ondelivery = x -> (seen[] = x; waslocked[] = islocked(q.mutex)))
        @test put!(q, 3) == 3
        @test seen[] == 3
        @test !waslocked[]
        @test take(q) == 3
        @test_throws ArgumentError put!(q, 4)
        @test seen[] == 3
        @test take(q) == 3

        # The callback sees the value after conversion, as `take` does.
        vals = Float64[]
        r = Promise{Float64}(ondelivery = x -> push!(vals, x))
        @test put!(r, 1) === 1.0
        @test vals == [1.0]
        @test take(r) === 1.0

        # A throwing callback propagates out of `put!`, but the promise has
        # already been delivered and cannot be rolled back.
        bad = Promise{Int}(ondelivery = _ -> error("callback boom"))
        @test_throws ErrorException put!(bad, 9)
        @test isready(bad)
        @test take(bad) == 9

        # Every taker receives the value.
        s = Promise{Int}()
        ts = [Threads.@spawn take(s) for _ in 1:4]
        put!(s, 11)
        @test all(fetch(t) == 11 for t in ts)

        # A callback runs with no lock held, so it may deliver another promise.
        inner = Promise{Int}()
        outer = Promise{Int}(ondelivery = x -> put!(inner, 2x))
        put!(outer, 5)
        @test take(inner) == 10

        # A promise may also be settled as failed; `take` then throws
        # `PromiseFailed` carrying the payload, `isready` is true, `isfailed`
        # tells the two settlements apart, and `geterror` reads the payload
        # without throwing.
        p = Promise{Int}()
        @test !isfailed(p)
        @test geterror(p) === nothing
        @test fail!(p, :boom) === :boom
        @test isready(p)
        @test isfailed(p)
        @test geterror(p) === :boom
        @test_throws PromiseFailed take(p)
        @test_throws PromiseFailed p[]
        err = try
            take(p)
        catch e
            e
        end
        @test err isa PromiseFailed
        @test err.err === :boom
        @test sprint(showerror, err) == "promise failed: :boom"
        @test_throws ArgumentError fail!(p, :again)
        @test_throws ArgumentError put!(p, 1)
        @test geterror(p) === :boom

        # The other order: a delivered promise cannot be failed.
        q = Promise{Int}()
        @test put!(q, 5) == 5
        @test_throws ArgumentError fail!(q, :late)
        @test !isfailed(q)
        @test geterror(q) === nothing
        @test take(q) == 5

        # `show` prints the failure without blocking.
        @test occursin("failed: :boom", sprint(show, MIME"text/plain"(), p))
        @test occursin("<...>", sprint(show, MIME"text/plain"(), Promise{Int}()))

        # The failure callback has the same contract as the delivery callback:
        # fixed at construction, run once, outside the lock, and its exception
        # propagates while the promise stays failed.
        seen = Ref{Any}(nothing)
        waslocked = Ref(true)
        f = Promise{Int}(onfailure = x -> (seen[] = x; waslocked[] = islocked(f.mutex)))
        @test fail!(f, 8) === 8
        @test seen[] == 8
        @test !waslocked[]
        @test_throws PromiseFailed take(f)
        @test seen[] == 8
        badfail = Promise{Int}(onfailure = _ -> error("failure callback boom"))
        @test_throws ErrorException fail!(badfail, 9)
        @test isfailed(badfail)
        @test geterror(badfail) === 9

        # Every taker parked on a promise wakes when it fails.
        g = Promise{Int}()
        ts = [Threads.@spawn (try
                    take(g)
                catch e
                    e
                end) for _ in 1:4]
        fail!(g, :shared)
        es = fetch.(ts)
        @test all(e isa PromiseFailed for e in es)
        @test all(e.err === :shared for e in es)

        # A settlement callback written after the settling task was created still
        # runs. A long-lived worker may settle a promise whose callback is newer
        # than the worker's world; without `invokelatest` that is a `MethodError`
        # saying the applicable method is too new. The bounded wait turns the
        # regression into a failure rather than a hung test.
        started = Channel{Nothing}(1)
        keep = Channel{Tuple{Promise{Int},Promise{Int}}}(1)
        ranok = Channel{Int}(1)
        ranbad = Channel{Symbol}(1)
        worker = Threads.@spawn begin
            put!(started, nothing)
            p, q = take!(keep)
            put!(p, 1)
            fail!(q, :bad)
        end
        take!(started)
        latecbs = @eval (x -> (put!($ranok, x); x),
                         e -> (put!($ranbad, e); e))
        okp = Promise{Int}(ondelivery = latecbs[1])
        badp = Promise{Int}(onfailure = latecbs[2])
        put!(keep, (okp, badp))
        @test Base.timedwait(() -> isready(ranok) && isready(ranbad), 10.0) === :ok
        @test take!(ranok) == 1
        @test take!(ranbad) === :bad
        try
            wait(worker)
        catch
        end
        @test !istaskfailed(worker)

        # Equality and hashing are by identity, so an undelivered promise can be
        # a dictionary key and delivering it does not move the key. The
        # value-based forms this replaces blocked on an undelivered promise; the
        # bounded wait turns that regression into a failure instead of a hung
        # suite.
        a = Promise{Int}()
        b = Promise{Int}()
        eq = Threads.@spawn begin
            d = Dict(a => :a, b => :b)
            before = (isequal(a, a), isequal(a, b), hash(a) == hash(a),
                      length(d), d[a], d[b])
            put!(a, 1)
            after = (d[a], length(d), take(a), isequal(a, a))
            (before..., after...)
        end
        t0 = time()
        while !istaskdone(eq) && (time() - t0) < 10.0
            sleep(0.01)
        end
        @test istaskdone(eq)
        if istaskdone(eq)
            @test fetch(eq) === (true, false, true, 2, :a, :b, :a, 2, 1, true)
        end
    end

    @testset "lockall" begin
        r1, r2, r3 = ReentrantLock(), ReentrantLock(), ReentrantLock()
        @test lockall(() -> :success, [r1, r2, r3]) === :success
        @test lockall(() -> :success, (r1, r2, r3)) === :success
        # The locks must be released once the body has run.
        @test all(!islocked(r) for r in (r1, r2, r3))
    end

    @testset "_to_pairs" begin
        # Regression: the empty-iterator path referenced an undefined `itr`.
        @test Air._to_pairs(Pair{Symbol,Int}[]) isa Vector{Pair{Symbol,Int}}
        @test isempty(Air._to_pairs(Pair{Symbol,Int}[]))
        @test Air._to_pairs([:a => 1, :b => 2]) == [:a => 1, :b => 2]
    end
end
