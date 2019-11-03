using Air
using Test

@testset "Air.jl" begin
    @testset "PVec" begin
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
    @testset "PMap" begin
        m0 = Air.PMap{Symbol}()
        pairs = Dict{UInt64,Symbol}(100 => :a, 5 => :b, 25 => :c,
                                    1 => :d, 0b100001 => :e,
                                    0b1000001 => :f)
        (m, n) = (m0, 0)
        for kv in pairs
            @test length(m) == n
            n += 1
            m = assoc(m, kv[1], kv[2])
        end
        for kv in pairs
            @test haskey(m, kv[1])
            @test m[kv[1]] === kv[2]
        end
        for kv in m
            @test haskey(pairs, kv[1])
            @test pairs[kv[1]] === kv[2]
        end
        mm = m
        for k in UInt64[0b1000001, 0b100001, 1]
            @test length(mm) == n
            mm = dissoc(mm, k)
            n -= 1
            @test !haskey(mm, k)
        end
    end
end
