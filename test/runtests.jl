using Air
using Test

@testset "Air.jl" begin
    # Some simple tests of the basic functions for PVec's:
    u = Air.PVec(1:100)
    @test length(u) == 100
    @test length(Air.dissoc(u, 100)) == 99
    @test length(Air.assoc(u, 101, 0)) == 101

    uu = [(1:100)...]
    @test u == uu
    v = Air.assoc(u, 10, 0)
    vv = [(1:100)...]
    vv[10] = 0
    @test v == vv
    @test v[10] == 0
    @test v[1] == 1
    @test v[100] == 100
end
