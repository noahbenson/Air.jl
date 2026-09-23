################################################################################
# pheap.jl
#
# Tests for the persistent heap that backs the weighted collections. `PHeap` is
# an internal type, so these exercise it directly via the `Air` module; the
# public `PWSet`/`PWDict` wrappers are covered in their own files.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "PHeap" begin
    @testset "construction" begin
        @test isempty(Air.PHeap())
        @test isempty(Air.PHeap{Int}())
        @test isempty(Air.PHeap{Int,Float64}())
        @test isempty(Air.PHeap(>))
        @test isempty(Air.PHeap{Int,Float64}(>))
    end

    @testset "push and first" begin
        h = Air.PHeap{Int,Float64}(>)
        h = push(h, (1, 1.0))
        h = push(h, (2, 5.0))
        h = push(h, (3, 3.0))
        @test length(h) == 3
        # `>` means the largest weight comes first.
        @test first(h) == 2
        @test collect(h) == [2, 3, 1]
        # Pushing an existing item updates its weight rather than duplicating.
        h = push(h, (1, 10.0))
        @test length(h) == 3
        @test first(h) == 1
    end

    @testset "custom comparison" begin
        h = Air.PHeap{Int,Float64}(<)
        h = push(h, (1, 1.0))
        h = push(h, (2, 5.0))
        h = push(h, (3, 3.0))
        # With `<`, the smallest weight comes first.
        @test first(h) == 1
        @test collect(h) == [1, 3, 2]
    end

    @testset "pop" begin
        h = Air.PHeap{Int,Float64}(>)
        h = push(h, (1, 1.0))
        h = push(h, (2, 5.0))
        h = push(h, (3, 3.0))
        h0 = h
        h = pop(h)
        @test length(h) == 2
        @test first(h) == 3
        # `pop` is persistent: the original heap is unaffected.
        @test length(h0) == 3
        @test first(h0) == 2
    end

    @testset "getweight/setweight" begin
        h = Air.PHeap{Int,Float64}(>)
        h = push(h, (1, 1.0))
        h = push(h, (2, 5.0))
        @test Air.getweight(h, 1) == 1.0
        @test Air.getweight(h, 2) == 5.0
        h = Air.setweight(h, 1, 9.0)
        @test Air.getweight(h, 1) == 9.0
        @test first(h) == 1
        @test_throws ErrorException Air.setweight(h, 99, 1.0)
    end

    @testset "weighted collections" begin
        s = PWSet{Int,Float64}()
        s = push(s, (1, 1.0))
        s = push(s, (2, 5.0))
        s = push(s, (3, 3.0))
        @test length(s) == 3
        @test first(s) == 2
        @test Air.getweight(s, 2) == 5.0
        s = Air.setweight(s, 1, 10.0)
        @test first(s) == 1

        d = PWDict{Symbol,Int,Float64}()
        d = push(d, (:a => 1) => 1.0)
        d = push(d, (:b => 2) => 5.0)
        @test length(d) == 2
        # `first` on a PWDict yields the key-value pair with the highest weight.
        @test first(d) == (:b => 2)
        @test d[:a] == 1
        @test Air.getweight(d, :b) == 5.0
    end

    @testset "pop preserves the subtree-weight totals" begin
        # KNOWN BUG (pre-existing, not introduced here): deleting from a
        # `PHeap` leaves the cached subtree-weight totals inconsistent with the
        # actual weights of the remaining nodes. The heap ordering itself is
        # correct, but `Random.rand` scales its draw by the root total, so
        # weighted sampling is skewed after any `pop`/`delete`. The `@test_broken`
        # below records this; it should start passing when the totals are fixed
        # in `_pheap_swap`/`_pheap_delete`.
        h = Air.PHeap{Int,Float64}(>)
        for it in ((1, 1.0), (2, 5.0), (3, 3.0), (4, 2.0), (5, 4.0))
            h = push(h, it)
        end
        h = pop(h)
        hp = Air.getfield(h, :_heap)
        root_tot = hp[1][3]
        weight_sum = sum((hp[i][2] for i in 1:length(hp)); init=0.0)
        @test_broken root_tot == weight_sum
    end
end
