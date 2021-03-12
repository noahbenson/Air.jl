# Tests for the PDict type.
# Author: Noah C. Benson <n@nben.net>

function compare_test(
    p::DIMM, s::DMUT, ks::Vector{K}, vs::Vector{V}, n::Int
) where {K,V, DIMM <: AbstractDict{K,V}, DMUT <: AbstractDict{K,V}}
    pd = [:set, :delete, :get]
    for i in 1:n
        q = rand(pd)
        if q == :set
            k = rand(ks)
            v = rand(vs)
            p = Air.push(p, k => v)
            push!(s, k => v)
        elseif q == :delete && length(s) > 0
            k = rand(keys(s))
            p = Air.delete(p, k)
            delete!(s, k)
        else
            k = rand(ks)
        end
        @test isequal(p, s)
        @test length(p) == length(s)
        @test get(p, k, nothing) == get(s, k, nothing)
        @test (==)((k => get(p, k, nothing)) in p,
                   (k => get(s, k, nothing)) in s)
    end
end
@testset "PDict" begin
    # Next, using the above function
    syms = Symbol[:a, :b, :c, :d, :e, :f, :g, :h, :i, :j]
    nums = Real[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    n = 100
    @testset "PIdDict" begin
        compare_test(Air.PIdDict{Symbol,Real}(),
                     Base.IdDict{Symbol,Real}(),
                     syms, nums, n)
        compare_test(Air.PIdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     Base.IdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     syms, nums, n)
    end
    @testset "PEquivDict" begin
        compare_test(Air.PEquivDict{Symbol,Real}(),
                     Air.EquivDict{Symbol,Real}(),
                     syms, nums, n)
        compare_test(Air.PEquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     Air.EquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     syms, nums, n)
    end
    @testset "PDict" begin
        compare_test(Air.PDict{Symbol,Real}(),
                     Base.Dict{Symbol,Real}(),
                     syms, nums, n)
        compare_test(Air.PDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     Base.Dict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                     syms, nums, n)
    end
end
