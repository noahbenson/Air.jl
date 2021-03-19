# Tests for the PWDict type.
# Author: Noah C. Benson <n@nben.net>

@testset "PWDict" begin
    intels = [gensym() => (x,rand(Float64)) for x in 1:2000]
    els = [:b => (2,20.0), :c => (3,30.0), :a => (1,10.0), :d => (4,40.0)]
    p = PWDict{Symbol,Int,Float64}(els...)
    for el in els
        p = push(p, el)
    end
    (lst,mst) = (first(p), pop(p))
    @test lst == (:d => 4)
    (lst,mst) = (first(mst), pop(mst))
    @test lst == (:c => 3)
    (lst,mst) = (first(mst), pop(mst))
    @test lst == (:b => 2)
    (lst,mst) = (first(mst), pop(mst))
    @test lst == (:a => 1)
    @test isempty(mst)
    # make many samples and make sure they resemble the distribution
    counts = Dict(:a => 0.0, :b => 0.0, :c => 0.0, :d => 0.0)
    vcounts = [0.0,0.0,0.0,0.0]
    for _ in 1:100000
        (k,v) = rand(p)
        counts[k] += 1/1000
        vcounts[v] += 1/1000
    end
    @test  8.5 <  counts[:a] < 12.5
    @test 18.5 <  counts[:b] < 22.5
    @test 28.5 <  counts[:c] < 32.5
    @test 38.5 <  counts[:d] < 42.5
    @test  8.5 < vcounts[1] < 12.5
    @test 18.5 < vcounts[2] < 22.5
    @test 28.5 < vcounts[3] < 32.5
    @test 38.5 < vcounts[4] < 42.5
end
