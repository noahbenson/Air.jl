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
        @test iserady(d)
        d2 = Delay{Int}(() -> 2)
        d3 = Delay{Int}(() -> (k += 1; k))
        @test hash(d) == hash(d2)
        @test hash(d) != hash(d3)        
    end

end
