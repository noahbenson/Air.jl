# Tests for the code in the API.jl file.
# Author: Noah C. Benson <n@nben.net>

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

    make_idset(u::Vector{T}) where {T} = begin
        s = IdSet{T}()
        for x in u
            push!(s, x)
        end
        return s
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

    struct TestType
        a::Int
        b::Symbol
    end
    @testset "iscoll" begin
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

    @testset "issingle" begin
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
