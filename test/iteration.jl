################################################################################
# iteration.jl
#
# Tests for iteration over the tree-backed collections.
#
# Iteration is driven by an explicit stack whose state is `(path, todo)`. That
# makes two things worth pinning down that the randomized comparisons elsewhere
# would not necessarily catch: that every element is yielded exactly once, and
# that a finished iterator stays finished rather than starting over (an earlier
# version of the stack could not tell "exhausted" from "not yet started").
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "iteration" begin
    # Sizes chosen to straddle the node fan-out (64 children per node), where
    # the stack gains and loses levels.
    SIZES = (0, 1, 2, 3, 63, 64, 65, 128, 1000)

    @testset "PDict" begin
        for n in SIZES
            d = PDict{Symbol,Int}()
            for i in 1:n
                d = push(d, Symbol("k", i) => i)
            end
            @test length(d) == n
            @test collect(keys(d)) == collect(keys(d))          # stable
            @test sort(collect(keys(d)); by=string) ==
                sort([Symbol("k", i) for i in 1:n]; by=string)
            @test Set(collect(d)) == Set(Symbol("k", i) => i for i in 1:n)
            @test sum(values(d); init=0) == div(n * (n + 1), 2)
            @test sort(collect(values(d))) == collect(1:n)
        end
        # elements are yielded exactly once
        d = PDict{Symbol,Int}(:a => 1, :b => 2, :c => 3)
        seen = Pair{Symbol,Int}[]
        x = iterate(d)
        while x !== nothing
            push!(seen, x[1])
            x = iterate(d, x[2])
        end
        @test length(seen) == 3
        @test Set(seen) == Set([:a => 1, :b => 2, :c => 3])
        # an exhausted iterator keeps yielding nothing rather than starting over
        d2 = PDict{Symbol,Int}(:a => 1, :b => 2)
        s1 = iterate(d2)
        s2 = iterate(d2, s1[2])
        @test s2 !== nothing
        @test iterate(d2, s2[2]) === nothing
    end

    @testset "PSet" begin
        for n in SIZES
            s = PSet(Symbol("k", i) for i in 1:n)
            @test length(s) == n
            got = collect(s)
            @test length(got) == n
            @test Set(got) == Set(Symbol("k", i) for i in 1:n)
        end
        s2 = PSet([:a, :b])
        y1 = iterate(s2)
        y2 = iterate(s2, y1[2])
        @test y2 !== nothing
        @test iterate(s2, y2[2]) === nothing
    end

    @testset "identity variants" begin
        # These come from the same macros, so exercise them too: keys are
        # compared by identity, which routes through different bucket types.
        for n in (0, 1, 65, 200)
            d = PIdDict{Int,Int}()
            for i in 1:n
                d = push(d, i => i)
            end
            @test length(d) == n
            @test sort(collect(keys(d))) == collect(1:n)

            s = PIdSet(i for i in 1:n)
            @test length(s) == n
            @test sort(collect(s)) == collect(1:n)
        end
    end

    @testset "lazy dictionaries iterate like dicts" begin
        d = LazyDict{Symbol,Int}(:a => 1, :b => 2)
        @test Set(collect(d)) == Set([:a => 1, :b => 2])
        @test d[:a] == 1
    end

    @testset "buckets with several entries" begin
        # Force hash collisions so that iteration has to walk more than one
        # element within a single bucket.
        d = PDict{Int,Int}()
        for i in 1:200
            d = push(d, i => i)
        end
        @test length(d) == 200
        @test sort(collect(keys(d))) == collect(1:200)
        # deleting from the middle leaves iteration consistent
        for k in (1, 7, 64, 65, 200)
            d = delete(d, k)
        end
        @test sort(collect(keys(d))) == setdiff(collect(1:200), [1, 7, 64, 65, 200])
    end
end
