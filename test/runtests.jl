using Air
using Test
using Random.Random

struct TestType
    a::Int
    b::Symbol
end

@testset "Air.jl" begin

    @testset "API" begin
        @testset "PArray" begin
            a = collect(reshape(1:100, (10,10)))
            b = convert(Array{Float64,2}, collect(reshape(1:100, (10,10))))
            p = PArray(a)
            q = PArray(b)
            @test a == b
            @test p == q
            @test a == p
            @test a == q
            @test b == p
            @test b == q
            @test !isequiv(a, b)
            @test !isequiv(a, p)
            @test !isequiv(a, q)
            @test !isequiv(b, p)
            @test !isequiv(b, q)
            @test isequiv(p, q)
            @test !(a === b)
            @test !(a === p)
            @test !(a === q)
            @test !(b === p)
            @test !(b === q)
            @test !(p === q)
            @test hash(a) == hash(b)
            @test hash(p) == hash(q)
            @test hash(a) == hash(p)
            @test hash(a) == hash(q)
            @test hash(b) == hash(p)
            @test hash(b) == hash(q)
            @test equivhash(a) != equivhash(b)
            @test equivhash(p) == equivhash(q)
            @test equivhash(a) != equivhash(p)
            @test equivhash(a) != equivhash(q)
            @test equivhash(b) != equivhash(p)
            @test equivhash(b) != equivhash(q)
        end
        @testset "PDict" begin
            intels = [gensym() => x for x in 1:2000]
            fltels = [k => Float64(v) for (k,v) in intels]
            a = Dict(intels...)
            b = Dict(fltels...)
            p = PDict(a)
            q = PDict(b)
            @test a == b
            @test p == q
            @test a == p
            @test a == q
            @test b == p
            @test b == q
            @test !isequiv(a, b)
            @test !isequiv(a, p)
            @test !isequiv(a, q)
            @test !isequiv(b, p)
            @test !isequiv(b, q)
            @test isequiv(p, q)
            @test !(a === b)
            @test !(a === p)
            @test !(a === q)
            @test !(b === p)
            @test !(b === q)
            @test !(p === q)
            @test hash(a) == hash(b)
            @test hash(p) == hash(q)
            @test hash(a) == hash(p)
            @test hash(a) == hash(q)
            @test hash(b) == hash(p)
            @test hash(b) == hash(q)
            @test equivhash(a) != equivhash(b)
            @test equivhash(p) == equivhash(q)
            @test equivhash(a) != equivhash(p)
            @test equivhash(a) != equivhash(q)
            @test equivhash(b) != equivhash(p)
            @test equivhash(b) != equivhash(q)
        end
        @testset "PSet" begin
            intels = collect(1:1000)
            fltels = [Float64(v) for v in intels]
            a = Set(intels)
            b = Set(fltels)
            p = PSet(a)
            q = PSet(b)
            @test a == b
            @test p == q
            @test a == p
            @test a == q
            @test b == p
            @test b == q
            @test !isequiv(a, b)
            @test !isequiv(a, p)
            @test !isequiv(a, q)
            @test !isequiv(b, p)
            @test !isequiv(b, q)
            @test isequiv(p, q)
            @test !(a === b)
            @test !(a === p)
            @test !(a === q)
            @test !(b === p)
            @test !(b === q)
            @test !(p === q)
            @test hash(a) == hash(b)
            @test hash(p) == hash(q)
            @test hash(a) == hash(p)
            @test hash(a) == hash(q)
            @test hash(b) == hash(p)
            @test hash(b) == hash(q)
            @test equivhash(a) != equivhash(b)
            @test equivhash(p) == equivhash(q)
            @test equivhash(a) != equivhash(p)
            @test equivhash(a) != equivhash(q)
            @test equivhash(b) != equivhash(p)
            @test equivhash(b) != equivhash(q)
        end
        @testset "PWDict" begin
            intels = [gensym() => (x,rand(Float64)) for x in 1:2000]
            fltels = [k => (Float64(v),w) for (k,(v,w)) in intels]
            p = PWDict(intels...)
            q = PWDict(fltels...)
            @test p == q
            @test isequiv(p, q)
            @test !(p === q)
            @test hash(p) == hash(q)
            @test equivhash(p) == equivhash(q)
        end
        @testset "PWSet" begin
            intels = [gensym() => rand(Float64) for x in 1:2000]
            p = PSet(intels)
            q = PSet(reverse(intels))
            @test p == q
            @test isequiv(p, q)
            @test !(p === q)
            @test hash(p) == hash(q)
            @test equivhash(p) == equivhash(q)
        end
    end
    
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
        @testset "2D" begin
            # These should be more fleshed out, but for now, we can do just a
            # few simple tests.
            a = Array(reshape(1:200, (20,10)))
            p = PArray(a)
            @test p == a
            @test !isequiv(p, a)
            @test p[5,8] == a[5,8]
            @test p[:,4] == a[:,4]
            @test p[9,:] == a[9,:]
            p1 = setindex(p, -5, 15, 8)
            @test p1 != a
            @test p1[:,4] == a[:,4]
            @test p1[9,:] == a[9,:]
            @test p1[:,8] != a[:,8]
            @test p1[15,:] != a[15,:]
        end
        @testset "3D" begin
            a = Array(reshape(1:1000, (20,10,5)))
            p = PArray(a)
            @test p == a
            @test !isequiv(p, a)
            @test p[5,8,2] == a[5,8,2]
            @test p[:,4,1] == a[:,4,1]
            @test p[9,:,3] == a[9,:,3]
            @test p[:,2,4] == a[:,2,4]
            p1 = setindex(p, -5, 15, 8, 5)
            @test p1 != a
            @test p1[9,:,:] == a[9,:,:]
            @test p1[:,3,:] == a[:,3,:]
            @test p1[:,:,4] == a[:,:,4]
            @test p1[15,:,:] != a[15,:,:]
            @test p1[:,8,:] != a[:,8,:]
            @test p1[:,:,5] != a[:,:,5]
        end
    end

    # #PSet ####################################################################
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

    # #PSet ####################################################################
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
            @test ((k => get(p, k, nothing)) in p) == ((k => get(s, k, nothing)) in s)
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
        @testset "PDict" begin
            compare_test(Air.PDict{Symbol,Real}(),
                         Air.EquivDict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.PDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Air.EquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
        end
        @testset "PEqualDict" begin
            compare_test(Air.PEqualDict{Symbol,Real}(),
                         Base.Dict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.PEqualDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Base.Dict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
        end
    end
    @testset "general" begin
        @testset "Air.iscoll" begin
            @test Air.iscoll([1,2,3])
            @test Air.iscoll(Set([:a, :b, :c]))
            @test Air.iscoll(Air.PDict(:a => 1, :b => 2))
            @test !Air.iscoll(TestType(10, :x))
            @test !Air.iscoll(10)
            @test !Air.iscoll(1.5)
            @test !Air.iscoll(:abc)
            @test !Air.iscoll("string")

            @test Air.iscolltype(Array{Int,2})
            @test Air.iscolltype(Set{Symbol})
            @test Air.iscolltype(Air.PDict{Symbol,Int})
            @test !Air.iscolltype(TestType)
            @test !Air.iscolltype(Int)
            @test !Air.iscolltype(Float64)
            @test !Air.iscolltype(Symbol)
            @test !Air.iscolltype(String)
        end
        @testset "collections" begin
            @test !Air.issingle([1,2,3])
            @test !Air.issingle(Set([:a, :b, :c]))
            @test !Air.issingle(Air.PDict(:a => 1, :b => 2))
            @test Air.issingle(TestType(10, :x))
            @test Air.issingle(10)
            @test Air.issingle(1.5)
            @test Air.issingle(:abc)
            @test Air.issingle("string")
            
            @test !Air.issingletype(Array{Int,2})
            @test !Air.issingletype(Set{Symbol})
            @test !Air.issingletype(Air.PDict{Symbol,Int})
            @test Air.issingletype(TestType)
            @test Air.issingletype(Int)
            @test Air.issingletype(Float64)
            @test Air.issingletype(Symbol)
            @test Air.issingletype(String)
        end
        
    end
end
