# Tests for the PWDict type.
# Author: Noah C. Benson <n@nben.net>

@testset "PWSet" begin
    els = [:b => 20.0, :c => 30.0, :a => 10.0, :d => 40.0]
    p = PWSet{Symbol,Float64}()
    for el in els
        p = push(p, el)
    end
    (lst, mst) = (first(p), pop(p))
    @test lst == :d
    (lst, mst) = (first(mst), pop(mst))
    @test lst == :c
    (lst, mst) = (first(mst), pop(mst))
    @test lst == :b
    (lst, mst) = (first(mst), pop(mst))
    @test lst == :a
    @test isempty(mst)
    # make many samples and make sure they resemble the distribution
    counts = Dict(:a => 0.0, :b => 0.0, :c => 0.0, :d => 0.0)
    for _ in 1:100000
        sym = rand(p)
        counts[sym] += 1/1000
    end
    @test 8.5 < counts[:a] < 12.5
    @test 18.5 < counts[:b] < 22.5
    @test 28.5 < counts[:c] < 32.5
    @test 38.5 < counts[:d] < 42.5
end

# As for `PWDict`: a `PWSet` weighs each element, so it is never `isequal` to a
# `PSet` with the same elements. `isequal(::PSet, ::PWSet)` was ambiguous rather
# than false before, so both orders are covered below.
@testset "PWSet equality" begin
    a = PWSet{Symbol,Float64}(:x => 2.0, :y => 3.0)
    b = PWSet{Symbol,Float64}(:y => 3.0, :x => 2.0)
    @test isequal(a, b)
    @test isequal(b, a)
    @test isequal(a, a)
    # same elements, one weight different: not equal, in both directions
    c = PWSet{Symbol,Float64}(:x => 2.0, :y => 9.0)
    @test !isequal(a, c)
    @test !isequal(c, a)
    # same elements, no weights at all: not equal, in both directions
    # (note that `PSet` has no element-vararg constructor, unlike `PSet{T,W}`)
    d = PSet([:x, :y])
    @test !isequal(a, d)
    @test !isequal(d, a)
    # and neither equal to something that is not a set
    @test !isequal(a, nothing)
    @test !isequal(nothing, a)
    @test !isequal(a, missing)
    @test !isequal(missing, a)
    @test !isequal(a, :x)
end

@testset "weighted set operations" begin
    a = PWSet{Int,Float64}(1 => 1.0, 2 => 2.0, 3 => 3.0)
    b = PWSet{Int,Float64}(3 => 30.0, 4 => 4.0)

    # An element only one argument holds keeps that argument's weight, whatever
    # the rule is — including `:mean`, which would otherwise divide by one.
    for rule in (:first, :last, :sum, :mean, :min, :max, :median)
        u = union(a, b; weight = rule)
        @test u isa PWSet{Int,Float64}
        @test sort(collect(u)) == [1, 2, 3, 4]
        @test getweight(u, 1) == 1.0        # only in `a`
        @test getweight(u, 2) == 2.0
        @test getweight(u, 4) == 4.0        # only in `b`
    end

    # ... and the rule decides for the one they share
    @test getweight(union(a, b; weight = :first), 3) == 3.0
    @test getweight(union(a, b; weight = :last), 3) == 30.0
    @test getweight(union(a, b; weight = :sum), 3) == 33.0
    @test getweight(union(a, b; weight = :mean), 3) == 16.5
    @test getweight(union(a, b; weight = :min), 3) == 3.0
    @test getweight(union(a, b; weight = :max), 3) == 30.0
    @test getweight(union(a, b; weight = :median), 3) == 16.5
    @test union(a, b) == union(a, b; weight = :first)   # the default is `:first`

    # A caller's own function is given the element and the weights it was found
    # with, and is called for every element of the result, not only the shared
    # ones. The weights arrive in argument order.
    f = (el, ws) -> sum(ws) * 10
    @test getweight(union(a, b; weight = f), 3) == 330.0
    @test getweight(union(a, b; weight = f), 1) == 10.0
    @test getweight(intersect(a, b; weight = f), 3) == 330.0
    @test getweight(symdiff(a, b; weight = f), 4) == 40.0
    # the weights arrive as a `Vector` of the arguments' weight type, in argument
    # order — and the value a function returns must be positive, since a weighted
    # set is a heap and `PHeap` rejects weights that are not
    seen = Ref{Any}(nothing)
    union(a, b; weight = (el, ws) -> (seen[] = ws; 1.0))
    @test seen[] isa Vector{Float64}
    @test_throws ArgumentError union(a, b; weight = (el, ws) -> -1.0)

    # the other three
    @test sort(collect(intersect(a, b))) == [3]
    @test getweight(intersect(a, b), 3) == 3.0
    @test sort(collect(setdiff(a, b))) == [1, 2]
    @test getweight(setdiff(a, b), 1) == 1.0
    @test sort(collect(symdiff(a, b))) == [1, 2, 4]
    @test getweight(symdiff(a, b), 4) == 4.0
    @test intersect(a, b) isa PWSet{Int,Float64}
    @test setdiff(a, b) isa PWSet{Int,Float64}
    @test symdiff(a, b) isa PWSet{Int,Float64}

    # a name that is not one of the rules is an error rather than a default
    @test_throws ArgumentError union(a, b; weight = :nope)

    # and the identity-keyed kind behaves the same way
    ia = PWIdSet{Int,Float64}(1 => 1.0, 3 => 3.0)
    ib = PWIdSet{Int,Float64}(3 => 30.0)
    @test union(ia, ib; weight = :sum) isa PWIdSet{Int,Float64}
    @test getweight(union(ia, ib; weight = :sum), 3) == 33.0
    @test getweight(union(ia, ib; weight = :sum), 1) == 1.0
end
