using Air
using Test

@testset "Air.jl" begin
    u = Air.PVec(1:100)
    @test length(u) == 100
end
