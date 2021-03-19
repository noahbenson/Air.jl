# Tests for the LazyDict type.
# Author: Noah C. Benson <n@nben.net>

@testset "LazyDict" begin
    syms = Symbol[:a, :b, :c, :d, :e, :f, :g, :h, :i, :j]
    nums = Real[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    n = 100
    # First just test them like normal dicts (they must work the same);
    # then we check delays. We use these functions, where the order they
    # are called will produce different results.
    # (Note that it's not a good idea to modify values like this inside
    # anonymous functions.
    fnval = 0
    delayfn1 = () -> begin
        fnval += 1
        return fnval
    end
    delayfn2 = () -> begin
        fnval *= 2
        return fnval
    end
    @testset "LazyIdDict" begin
        compare_test(Air.LazyIdDict{Symbol,Real}(),
                     Base.IdDict{Symbol,Real}(),
                     syms, nums, n)
        compare_test(Air.LazyIdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     Base.IdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     syms, nums, n)
        d = Air.LazyIdDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                       :b => Delay{Int}(delayfn2))
        fnval = 0
        @test (d[:a], d[:b]) == (1, 2)
        @test (d[:b], d[:a]) == (2, 1)
        d = Air.LazyIdDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                       :b => Delay{Int}(delayfn2))
        fnval = 0
        @test (d[:b], d[:a]) == (0, 1)
        @test (d[:a], d[:b]) == (1, 0)
    end
    @testset "LazyDict" begin
        compare_test(Air.LazyDict{Symbol,Real}(),
                     Base.Dict{Symbol,Real}(),
                     syms, nums, n)
        compare_test(Air.LazyDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     Base.Dict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     syms, nums, n)
        d = Air.LazyDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                     :b => Delay{Int}(delayfn2))
        fnval = 0
        @test (d[:a], d[:b]) == (1, 2)
        @test (d[:b], d[:a]) == (2, 1)
        d = Air.LazyDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                     :b => Delay{Int}(delayfn2))
        fnval = 0
        @test (d[:b], d[:a]) == (0, 1)
        @test (d[:a], d[:b]) == (1, 0)
    end
end    
