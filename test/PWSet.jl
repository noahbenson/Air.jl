# Tests for the PWDict type.
# Author: Noah C. Benson <n@nben.net>

@testset "PWSet" begin
    els = [:b => 20.0, :c => 30.0, :a => 10.0, :d => 40.0]
    p = PWSet{Symbol,Float64}()
    for el in els
        p = push(p, el)
    end
    (lst,mst) = (first(p), pop(p))
    @test lst == :d
    (lst,mst) = (first(mst), pop(mst))
    @test lst == :c
    (lst,mst) = (first(mst), pop(mst))
    @test lst == :b
    (lst,mst) = (first(mst), pop(mst))
    @test lst == :a
    @test isempty(mst)
    # make many samples and make sure they resemble the distribution
    counts = Dict(:a => 0.0, :b => 0.0, :c => 0.0, :d => 0.0)
    for _ in 1:100000
        sym = rand(p)
        counts[sym] += 1/1000
    end
    @test  8.5 < counts[:a] < 12.5
    @test 18.5 < counts[:b] < 22.5
    @test 28.5 < counts[:c] < 32.5
    @test 38.5 < counts[:d] < 42.5
end
