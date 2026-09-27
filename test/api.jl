# Tests for the code in the API.jl file.
# Author: Noah C. Benson <n@nben.net>

@testset "API" begin
    @testset "setindex" begin
        d1 = Dict()
        d2 = setindex(d1, 10, :a)
        @test d1 != d2
        @test d2[:a] == 10
        @test length(d1) == 0
        d1 = IdDict()
        d2 = setindex(d1, 10, :a)
        @test d1 != d2
        @test d2[:a] == 10
        @test length(d1) == 0
        u1 = [:x, :y, :z]
        u2 = setindex(u1, :a, 1)
        @test u1 != u2
        @test u2[1] == :a
        @test u1[1] == :x
        u1 = view(u1, 1:2)
        u2 = setindex(u1, :a, 1)
        @test u2[1] == :a
        @test u1[1] == :x
        u1 = reshape(u1, (1, 2))
        u2 = setindex(u1, :a, 1, 1)
        @test u2[1, 1] == :a
        @test u1[1, 1] == :x
        t1 = (:x, :y, :z)
        t2 = setindex(t1, :a, 1)
        @test t2[1] == :a
        @test length(t2) == 3
        b1 = BitArray([1, 1, 1, 1, 1])
        b2 = setindex(b1, 0, 3)
        @test b1[3] == true
        @test b2[3] == false
    end

    @testset "PArray" begin
        a = collect(reshape(1:100, (10, 10)))
        b = convert(Array{Float64,2}, collect(reshape(1:100, (10, 10))))
        p = PArray(a)
        q = PArray(b)
        @test a == b
        @test p == q
        @test a == p
        @test a == q
        @test b == p
        @test b == q
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
    end

    @testset "PDict" begin
        intels = [gensym() => x for x in 1:2000]
        fltels = [k => Float64(v) for (k, v) in intels]
        for (DT, PDT) in ((Dict, PDict), (IdDict, PIdDict))
            a = DT(intels...)
            b = DT(fltels...)
            p = PDT(a)
            q = PDT(b)
            @test a == b
            @test a == p
            @test a == q
            @test b == p
            @test b == q
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
        for (ST, PST) in ((Set, PSet), (IdSet, PIdSet))
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
        end
    end

    @testset "PWDict" begin
        intels = [gensym() => (x, rand(Float64)) for x in 1:2000]
        fltels = [k => (Float64(v), w) for (k, (v, w)) in intels]
        p = PWDict{Symbol,Int,Float64}(intels...)
        q = PWDict{Symbol,Float64,Float64}(fltels...)
        @test p == q
        @test !(p === q)
        @test hash(p) == hash(q)
    end

    @testset "PWSet" begin
        intels = [gensym() => rand(Float64) for x in 1:2000]
        p = PWSet{Symbol,Float64}(intels)
        q = PWSet{Symbol,Float64}(reverse(intels))
        @test p == q
        @test !(p === q)
        @test hash(p) == hash(q)
    end

@testset "operations that keep the collection's kind" begin
    # `Base` does these for `AbstractDict`/`AbstractSet` by building a mutable
    # `Dict`/`Set`, which for a persistent collection either raises (its `empty`
    # is immutable) or returns the wrong kind.
    d = PDict(:a => 1, :b => 2, :c => 3)
    f = filter(kv -> isodd(kv.second), d)
    @test f isa PDict{Symbol,Int}
    @test length(f) == 2 && f[:a] == 1 && f[:c] == 3 && !haskey(f, :b)
    @test length(d) == 3                        # the source is untouched

    r = replace(d, 1 => 9)
    @test r isa PDict{Symbol,Int}
    @test r[:a] == 9 && r[:b] == 2 && length(r) == 3
    @test d[:a] == 1

    a = PSet([1, 2, 3])
    b = PSet([3, 4])
    @test filter(isodd, a) isa PSet{Int}
    @test collect(filter(isodd, a)) == [1, 3]
    @test replace(a, 1 => 9) isa PSet{Int}
    @test 9 in replace(a, 1 => 9)
    @test union(a, b) isa PSet{Int}
    @test sort(collect(union(a, b))) == [1, 2, 3, 4]
    @test sort(collect(intersect(a, b))) == [3]
    @test sort(collect(setdiff(a, b))) == [1, 2]
    @test sort(collect(symdiff(a, b))) == [1, 2, 4]
    @test union(a, b, PSet([9])) isa PSet{Int}
    # a mixed call with a plain `Set` is `Base`'s business, and gives `Base`'s
    # answer
    @test union(a, Set([9])) isa Set{Int}

    # a weighted set keeps each element's weight, including through a filter and
    # a replacement
    w = PWSet{Int,Float64}(1 => 1.0, 2 => 2.0, 3 => 3.0)
    fw = filter(isodd, w)
    @test fw isa PWSet{Int,Float64}
    @test sort(collect(fw)) == [1, 3]
    @test getweight(fw, 3) == 3.0
    @test getweight(fw, 1) == 1.0
    rw = replace(w, 1 => 9)
    @test rw isa PWSet{Int,Float64}
    @test 9 in rw && !(1 in rw)
    @test getweight(rw, 9) == 1.0
end
end
