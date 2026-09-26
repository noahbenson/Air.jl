################################################################################
# arrayops.jl
#
# Tests for the array operations that answer a persistent array with a
# persistent array.
#
# The property under test is that an operation over a persistent collection does
# not quietly stop being persistent, and that a sparse operand stays sparse: each
# operation carries the operand's default along. Values are compared against the
# same operation on the equivalent `Array`.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "array operations keep the kind" begin
    u = PVector(collect(1:5))
    # a `PArray` is immutable, so sharing is a copy
    @test copy(u) === u
    @test copy(u) isa PVector{Int}
    @test map(x -> x + 1, u) isa PVector{Int}
    @test collect(map(x -> x + 1, u)) == collect(u) .+ 1

    # a sparse operand keeps its default through every one of these, so the
    # positions it does not store stay unstored
    s = setindex(setindex(PVector{Float64}(2.0, (6,)), 9.0, 2), 9.0, 5)
    @test nnz(s) == 2
    @test Air.defaultvalue(s) == 2.0
    m = map(x -> x + 1, s)
    @test m isa PVector{Float64}
    @test Air.defaultvalue(m) == 3.0
    @test nnz(m) == 2
    @test collect(m) == collect(s) .+ 1

    @testset "filter" begin
        f = filter(isodd, u)
        @test f isa PVector{Int}
        @test collect(f) == [1, 3, 5]
        @test collect(filter(iseven, u)) == filter(iseven, collect(u))
        # keeping an entry equal to the default does not store it
        keen = filter(x -> true, s)
        @test keen isa PVector{Float64}
        @test Air.defaultvalue(keen) == 2.0
        @test nnz(keen) == 2
        @test length(keen) == 6
        @test collect(keen) == collect(s)
        @test collect(filter(iseven, s)) == filter(iseven, collect(s))
        @test isempty(filter(x -> false, u))
    end

    @testset "reverse" begin
        r = reverse(u)
        @test r isa PVector{Int}
        @test collect(r) == reverse(collect(u))
        rs = reverse(s)
        @test rs isa PVector{Float64}
        @test Air.defaultvalue(rs) == 2.0
        @test nnz(rs) == 2
        @test length(rs) == 6
        @test collect(rs) == reverse(collect(s))
        @test reverse(reverse(s)) == s
    end

    @testset "indexing with a vector of positions" begin
        @test u[1:2] isa PVector{Int}
        @test collect(u[1:2]) == [1, 2]
        @test collect(u[[1, 3, 5]]) == [1, 3, 5]
        @test collect(u[5:-1:3]) == [5, 4, 3]
        @test length(u[2:4]) == 3
        @test collect(Int[]) == collect(u[Int[]])
        # a selected element equal to the default is not stored
        si = s[[1, 2, 6]]
        @test si isa PVector{Float64}
        @test Air.defaultvalue(si) == 2.0
        @test nnz(si) == 1
        @test collect(si) == collect(s)[[1, 2, 6]]
        # out-of-range positions are an error, as for an `Array`
        @test_throws BoundsError u[1:6]
        @test_throws BoundsError u[[0]]
    end

    @testset "a transient answers with a transient" begin
        t = transient(u)
        @test map(x -> x + 1, t) isa TVector{Int}
        @test copy(t) isa TVector{Int}
        @test collect(map(x -> x + 1, t)) == collect(u) .+ 1
        @test length(t) == 5                 # the transient is unchanged
    end
end
