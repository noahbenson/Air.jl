################################################################################
# util.jl
#
# Tests for the Air utilities.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

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
