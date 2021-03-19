# Tests for the PSet types.
# Author: Noah C. Benson <n@nben.net>

function compare_test(p::SIMM, s::SMUT, ks::Vector{T}, n::Int) where {T,SIMM<:AbstractSet{T},SMUT<:AbstractSet{T}}
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
@testset "PSet" begin
    # First, do a standard test:
    numops = 100
    numits = 10
    for it in 1:numits
        mut = Set{Real}()
        imm = Air.PSet{Real}()
        ops = [:push, :pop, :get]
        for ii in 1:numops
            q = rand(ops)
            if q == :push
                x = rand(Float64)
                push!(mut, x)
                imm = Air.push(imm, x)
            elseif q == :pop && length(mut) > 0
                x = rand(mut)
                delete!(mut, x)
                imm = Air.delete(imm, x)
            elseif q == :get && length(mut) > 0
                x = rand(mut)
                @test in(x, imm)
            end
            @test length(mut) == length(imm)
            @test mut == imm
        end
    end
    # Next, using the above function
    let syms = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j], n = 100
        @testset "PIdSet" begin
            compare_test(Air.PIdSet{Symbol}(), Base.IdSet{Symbol}(), syms, n)
            idset = Base.IdSet{Symbol}()
            for k in [:b, :d, :e]
                push!(idset, k)
            end
            compare_test(Air.PIdSet{Symbol}([:b, :d, :e]), idset, syms, n)
        end
        @testset "PSet" begin
            compare_test(Air.PSet{Symbol}(), Set{Symbol}(), syms, n)
            compare_test(Air.PSet{Symbol}([:b, :d, :e]), Set{Symbol}([:b, :d, :e]), syms, n)
        end
    end
end
