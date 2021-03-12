using Air
using Test
using Random.Random
import Base.IdSet
import Base.delete!
import Base.isready

struct TestType
    a::Int
    b::Symbol
end

make_idset(u::Vector{T}) where {T} = begin
    s = IdSet{T}()
    for x in u
        push!(s, x)
    end
    return s
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
            for (DT,PDT) in ((Dict,      PDict),
                             (EquivDict, PEquivDict),
                             (IdDict,    PIdDict))
                a = DT(intels...)
                b = DT(fltels...)
                p = PDT(a)
                q = PDT(b)
                @test a == b
                @test a == p
                @test a == q
                @test b == p
                @test b == q
                @test isequiv(p, q)
                @test !isequiv(a, b)
                @test !isequiv(a, p)
                @test !isequiv(a, q)
                @test !isequiv(b, p)
                @test !isequiv(b, q)
                @test !(a === b)
                @test !(a === p)
                @test !(a === q)
                @test !(b === p)
                @test !(b === q)
                @test !(p === q)
                @test hash(a) == hash(b)
                @test hash(p) == hash(q)
                if DT === Dict
                    @test hash(a) == hash(p)
                    @test hash(a) == hash(q)
                    @test hash(b) == hash(p)
                    @test hash(b) == hash(q)
                end
                @test equivhash(a) != equivhash(b)
                @test equivhash(p) == equivhash(q)
                @test equivhash(a) != equivhash(p)
                @test equivhash(a) != equivhash(q)
                @test equivhash(b) != equivhash(p)
                @test equivhash(b) != equivhash(q)
            end
        end
        @testset "PSet" begin
            intels = collect(1:1000)
            fltels = [Float64(v) for v in intels]
            for (ST,PST) in ((Set,        PSet),
                             (EquivSet,   PEquivSet),
                             (IdSet,      PIdSet))
                if ST === IdSet
                    a = make_idset(intels)
                    b = make_idset(fltels)
                else
                    a = ST(intels)
                    b = ST(fltels)
                end
                p = PST(a)
                q = PST(b)
                @test a == p
                @test b == q
                @test (a == b) == (p == q)
                @test !isequiv(a, b)
                @test !isequiv(a, p)
                @test !isequiv(a, q)
                @test !isequiv(b, p)
                @test !isequiv(b, q)
                @test isequiv(p, q) == (ST != IdSet)
                @test !(a === b)
                @test !(a === p)
                @test !(a === q)
                @test !(b === p)
                @test !(b === q)
                @test !(p === q)
                @test hash(a) == hash(b)
                @test hash(p) == hash(q)
                if ST === Set
                    @test hash(a) == hash(p)
                    @test hash(a) == hash(q)
                    @test hash(b) == hash(p)
                    @test hash(b) == hash(q)
                end
                @test equivhash(a) != equivhash(b)
                @test equivhash(p) == equivhash(q)
                @test equivhash(a) != equivhash(p)
                @test equivhash(a) != equivhash(q)
                @test equivhash(b) != equivhash(p)
                @test equivhash(b) != equivhash(q)
            end
        end
        @testset "PWDict" begin
            intels = [gensym() => (x,rand(Float64)) for x in 1:2000]
            fltels = [k => (Float64(v),w) for (k,(v,w)) in intels]
            p = PWDict{Symbol,Int,Float64}(intels...)
            q = PWDict{Symbol,Float64,Float64}(fltels...)
            @test p == q
            @test isequiv(p, q)
            @test !(p === q)
            @test hash(p) == hash(q)
            @test equivhash(p) == equivhash(q)
        end
        @testset "PWSet" begin
            intels = [gensym() => rand(Float64) for x in 1:2000]
            p = PWSet{Symbol,Float64}(intels)
            q = PWSet{Symbol,Float64}(reverse(intels))
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
            p1 = Air.setindex(p, -5, 15, 8)
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
            p1 = Air.setindex(p, -5, 15, 8, 5)
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
                idset = Base.IdSet{Symbol}()
                for k in [:b, :d, :e]
                    push!(idset, k)
                end
                compare_test(Air.PIdSet{Symbol}([:b, :d, :e]), idset, syms, n)
            end
            @testset "PEquivSet" begin
                compare_test(Air.PEquivSet{Symbol}(), Air.EquivSet{Symbol}(), syms, n)
                compare_test(Air.PEquivSet{Symbol}([:b, :d, :e]), Air.EquivSet{Symbol}([:b, :d, :e]), syms, n)
            end
            @testset "PSet" begin
                compare_test(Air.PSet{Symbol}(), Set{Symbol}(), syms, n)
                compare_test(Air.PSet{Symbol}([:b, :d, :e]), Set{Symbol}([:b, :d, :e]), syms, n)
            end
        end
    end

    # #PDict ###################################################################
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
            @test (==)((k => get(p, k, nothing)) in p,
                       (k => get(s, k, nothing)) in s)
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
        @testset "PEquivDict" begin
            compare_test(Air.PEquivDict{Symbol,Real}(),
                         Air.EquivDict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.PEquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Air.EquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
        end
        @testset "PDict" begin
            compare_test(Air.PDict{Symbol,Real}(),
                         Base.Dict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.PDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Base.Dict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
        end
    end
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
    @testset "LazyDict" begin
        syms = Symbol[:a, :b, :c, :d, :e, :f, :g, :h, :i, :j]
        nums = Real[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        n = 100
        # First just test them like normal dicts (they must work the same);
        # then we check delays. We use these functions, where the order they
        # are called will produce different results.
        # (Note that it's not a good idea to modify values like this inside
        # anonymous functions.
        fnval = 0
        delayfn1 = () -> begin
            fnval += 1
            return fnval
        end
        delayfn2 = () -> begin
            fnval *= 2
            return fnval
        end
        @testset "LazyIdDict" begin
            compare_test(Air.LazyIdDict{Symbol,Real}(),
                         Base.IdDict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.LazyIdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Base.IdDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
            d = Air.LazyIdDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                           :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:a], d[:b]) == (1, 2)
            @test (d[:b], d[:a]) == (2, 1)
            d = Air.LazyIdDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                           :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:b], d[:a]) == (0, 1)
            @test (d[:a], d[:b]) == (1, 0)
        end
        @testset "LazyEquivDict" begin
            compare_test(Air.LazyEquivDict{Symbol,Real}(),
                         Air.EquivDict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.LazyEquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Air.EquivDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
            d = Air.LazyEquivDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                              :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:a], d[:b]) == (1, 2)
            @test (d[:b], d[:a]) == (2, 1)
            d = Air.LazyEquivDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                              :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:b], d[:a]) == (0, 1)
            @test (d[:a], d[:b]) == (1, 0)
        end
        @testset "LazyDict" begin
            compare_test(Air.LazyDict{Symbol,Real}(),
                         Base.Dict{Symbol,Real}(),
                         syms, nums, n)
            compare_test(Air.LazyDict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         Base.Dict{Symbol,Real}(:b=>20, :d=>40, :e=>50),
                         syms, nums, n)
            d = Air.LazyDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                         :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:a], d[:b]) == (1, 2)
            @test (d[:b], d[:a]) == (2, 1)
            d = Air.LazyDict{Symbol,Int}(:a => Delay{Int}(delayfn1),
                                         :b => Delay{Int}(delayfn2))
            fnval = 0
            @test (d[:b], d[:a]) == (0, 1)
            @test (d[:a], d[:b]) == (1, 0)
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
