using Air
using Test
using Random.Random

@testset "Air.jl" begin

    # #PVec and #PMap, simple tests ############################################
    #@testset "PVec" begin
    #    # Some simple tests of the basic functions for PVec's:
    #    u = Air.PVec(1:100)
    #    @test length(u) == 100
    #    @test length(Air.dissoc(u, 100)) == 99
    #    @test length(Air.assoc(u, 101, 0)) == 101
    #
    #    uu = [(1:100)...]
    #    @test u == uu
    #    v = Air.assoc(u, 10, 0)
    #    vv = [(1:100)...]
    #    vv[10] = 0
    #    @test v == vv
    #    @test v[10] == 0
    #    @test v[1] == 1
    #    @test v[100] == 100
    #end
    #@testset "PMap" begin
    #    m0 = Air.PMap{Symbol}()
    #    pairs = Dict{UInt64,Symbol}(100 => :a, 5 => :b, 25 => :c,
    #                                1 => :d, 0b100001 => :e,
    #                                0b1000001 => :f)
    #    (m, n) = (m0, 0)
    #    for kv in pairs
    #        @test length(m) == n
    #        n += 1
    #        m = Air.assoc(m, kv[1], kv[2])
    #    end
    #    for kv in pairs
    #        @test haskey(m, kv[1])
    #        @test m[kv[1]] === kv[2]
    #    end
    #    for kv in m
    #        @test haskey(pairs, kv[1])
    #        @test pairs[kv[1]] === kv[2]
    #    end
    #    mm = m
    #    for k in UInt64[0b1000001, 0b100001, 1]
    #        @test length(mm) == n
    #        mm = Air.dissoc(mm, k)
    #        n -= 1
    #        @test !haskey(mm, k)
    #    end
    #end
    @testset "PArray" begin
        numops = 100
        @testset "1D" begin
            a = Real[]
            p = Air.PVector{Real}()
            ops = [:push, :pushfirst, :pop, :popfirst, :set, :get]
            for ii in 1:numops
                q = rand(ops)
                if q == :push
                    x = rand(Float64)
                    push!(a, x)
                    p = Air.push(p, x)
                elseif q == :pop && length(a) > 0
                    x = pop!(a)
                    @test x == p[end]
                    p = Air.pop(p)
                elseif q == :pushfirst
                    x = rand(Float64)
                    pushfirst!(a, x)
                    p = Air.pushfirst(p, x)
                elseif q == :popfirst && length(a) > 0
                    x = popfirst!(a)
                    @test x == p[1]
                    p = Air.popfirst(p)
                elseif q == :set && length(a) > 0
                    k = rand(1:length(a))
                    v = rand(Float64)
                    a[k] = v
                    p = Air.setindex(p, v, k)
                elseif q == :get && length(a) > 0
                    k = rand(1:length(a))
                    @test a[k] == p[k]
                end
                @test a == p
                @test size(a) == size(p)
            end
        end
        #@testset "2D" begin
        #end
        #@testset "3D" begin
        #end
    end

    # #PSet ####################################################################
    function compare_test(p::Air.PSet{T}, s::AbstractSet{T}, ks::AbstractArray{T,1},
                          n::Integer) where {T}
        let k, q, pd = [:push, :delete], ks = collect(ks)
            for i in 1:n
                k = rand(ks)
                q = rand(pd)
                if q == :push
                    p = Air.push(p, k)
                    push!(s, k)
                else
                    p = Air.delete(p, k)
                    delete!(s, k)
                end
                @test length(p) == length(s)
                @test (k in p) == (k in s)
                @test isequal(p, s)
            end
        end
    end
    if false
    let syms = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j], n = 100
        @testset "PIdSet" begin
            compare_test(Air.PIdSet{Symbol}(), Base.IdSet{Symbol}(), syms, n)
            compare_test(Air.PIdSet{Symbol}([:b, :d, :e]), Base.IdSet{Symbol}([:b, :d, :e]), syms, n)
        end
        @testset "PSet" begin
            compare_test(Air.PSet{Symbol}(), Air.EquivSet{Symbol}(), syms, n)
            compare_test(Air.PSet{Symbol}([:b, :d, :e]), Air.EquivSet{Symbol}([:b, :d, :e]), syms, n)
        end
        @testset "PEqualSet" begin
            compare_test(Air.PEqualSet{Symbol}(), Set{Symbol}(), syms, n)
            compare_test(Air.PEqualSet{Symbol}([:b, :d, :e]), Set{Symbol}([:b, :d, :e]), syms, n)
        end
    end
    end
end
