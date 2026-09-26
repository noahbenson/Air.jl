################################################################################
# regressions.jl
#
# Regression tests for long-standing defects that were never exercised by the
# test suite. Most of these lived in code paths that no test reached, so they
# only ever surfaced as `UndefVarError`s for users who happened to hit them.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "regressions" begin
    @testset "PTree internals" begin
        # `Base.empty(::PTree)` referenced an undefined `r`.
        @test isempty(Base.empty(Air.PTree{Int}()))
        # `isequal(::PTree,::PTree)` had a `getfiield` typo.
        @test isequal(Air.PTree{Int}(), Air.PTree{Int}())
        # The `convert` method for 1-d arrays had its eltype and dimensionality
        # swapped (`AbstractArray{1,U}`), so it could never be called.
        @test hasmethod(convert, Tuple{Type{Air.PTree{Int}},Vector{Int}})
        @test isempty(convert(Air.PTree{Int}, Int[]))
        # `Base.in` for a PTree referenced undefined `k` and `df`.
        @test !in(Air.HASH_T(1) => 1, Air.PTree{Int}(), (===))
    end

    @testset "PLinearSet construction and deletion" begin
        # The de-duplicating constructor used an undefined counter `k` and
        # compared the wrong element.
        @test length(Air.PLinearSet([1, 2, 2, 3])) == 3
        @test Set(Air.PLinearSet([1, 2, 2, 3])) == Set([1, 2, 3])
        @test collect(Air.PLinearSet([1, 2, 2, 3])) == [1, 2, 3]
        @test length(Air.PLinearSet(Int[])) == 0
        # `Air.delete` on a set referenced an undefined `elss`.
        for S in (Air.PLinearSet, Air.PIdLinearSet)
            @test length(delete(S([1, 2, 3]), 2)) == 2
            @test 2 ∉ delete(S([1, 2, 3]), 2)
            @test length(delete(S([1]), 1)) == 0
        end
        # A plain `Vector`-backed set is reachable from the public API too.
        @test length(PSet([1, 2, 2, 3])) == 3
    end

    @testset "empty typed pair iterators" begin
        # `_to_pairs` referenced an undefined `itr` on the empty-iterator path,
        # which every `PDict`/`PSet` constructor call could reach.
        @test isempty(PDict(Pair{Symbol,Int}[]))
        @test PDict(Pair{Symbol,Int}[]) isa PDict{Symbol,Int}
        @test isempty(PSet(Pair{Symbol,Int}[]))
    end

    @testset "@memoize" begin
        # Plain (untyped) arguments reached `_memoize_fixarg` as a `Symbol`.
        @eval module RegMemoizePlain
            using Air
            @memoize f(x) = x + 1
        end
        @test RegMemoizePlain.f(1) == 2
        # The result type tag was read from the `::` expression head.
        @eval module RegMemoizeTagged
            using Air
            @memoize g(x::Int) = (x + 1.5)::Float64
        end
        @test RegMemoizeTagged.g(1) === 2.5
    end

    @testset "PMatrix(default, sz)" begin
        # Referenced an undefined `len`.
        m = PMatrix(0.0, (2, 3))
        @test size(m) == (2, 3)
        @test all(m[i, j] == 0.0 for i in 1:2, j in 1:3)
    end

    @testset "weighted collections default constructors" begin
        # `PWSet()`/`PWIdSet()` referred to a type variable `T` that was not in
        # scope; `PHeap()` and friends used a bare `PDict` where the struct
        # requires `AbstractPDict{T,Int}`.
        for W in (PWSet, PWIdSet)
            @test isempty(W())
            @test length(push(W(), (1, 2.0))) == 1
        end
        @test isempty(PWDict())
        @test isempty(Air.PHeap())
        @test isempty(Air.PHeap{Int}())
        @test isempty(Air.PHeap{Int,Float64}())
        @test isempty(Air.PHeap(>))
    end

    @testset "Transaction properties" begin
        # `getproperty(::Transaction, …)` used the property Symbol in place of
        # the transaction itself, and omitted the `::Symbol` annotation, making
        # it ambiguous with `Base.getproperty`.
        t = Air.Transaction()
        @test t.state === :running
        @test isempty(t.rvolatiles)
        @test isempty(t.wvolatiles)
        @test isempty(t.actors)
        @test_throws ErrorException t.nonexistent
    end

    @testset "Volatile filter and finalize" begin
        @testset "filter applies exactly once" begin
            v = Volatile{Int}(0)
            Air.tx() do
                setfilter!(v, x -> x + 1)
            end
            Air.tx() do
                v[] = 10
            end
            @test v[] == 11
        end

        @testset "finalize applies at commit" begin
            v = Volatile{Int}(0)
            Air.tx() do
                setfinalize!(v, x -> x * 2)
            end
            Air.tx() do
                v[] = 3
            end
            @test v[] == 6
        end

        @testset "setfilter! does not refilter the value" begin
            v = Volatile{Int}(5)
            Air.tx() do
                setfilter!(v, x -> x + 100)
            end
            @test v[] == 5
        end

        @testset "getfilter/getfinalize" begin
            v = Volatile{Int}(0)
            Air.tx() do
                setfilter!(v, x -> x + 1)
                setfinalize!(v, x -> x * 2)
            end
            Air.tx() do
                @test getfilter(v) !== nothing
                @test getfinalize(v) !== nothing
            end
        end
    end

    @testset "@var error handling" begin
        # These error paths used the non-existent `ArgumentException`.
        @test_throws LoadError @eval @var 3 = 1
        @test_throws LoadError @eval @var
    end

    @testset "pop on a persistent set" begin
        # `pop(::AbstractPSet)` referenced an undefined `d` instead of its
        # `s` argument.
        s = PSet([1, 2, 3])
        (x, rest) = pop(s)
        @test x in (1, 2, 3)
        @test length(rest) == 2
        @test x ∉ rest
        # `pop` is persistent: the original set is unchanged.
        @test length(s) == 3
        @test_throws ArgumentError pop(PSet{Int}())
    end

    @testset "construction from a size vector" begin
        # `PArray(default, size::Vector)` referenced an unbound `N`.
        a = PArray(0.0, [2, 3])
        @test size(a) == (2, 3)
        @test all(a[i, j] == 0.0 for i in 1:2, j in 1:3)
    end

    @testset "exported names are defined" begin
        for s in (:Var, :Volatile, :Actor, :tx, :send, :geterror, :reset, Symbol("@p"))
            @test isdefined(Air, s)
        end
        # `Source`, `AbstractSourceKernel`, and `receive` were exported but
        # never defined anywhere; they must no longer be claimed as API.
        for s in (:Source, :AbstractSourceKernel, :receive)
            @test !isdefined(Air, s)
        end
    end
end
