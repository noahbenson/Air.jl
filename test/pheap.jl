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

    @testset "the cached subtree totals stay consistent" begin
        # `PHeap` caches, for every node, the total weight of its subtree, and
        # `Random.rand` scales its draw by the root's total. If those totals
        # drift, weighted sampling is silently skewed, so check the defining
        # relation — a node's total is its own weight plus its children's —
        # after every operation.
        badtotals(h) = begin
            hp = Air.getfield(h, :_heap)
            n = length(hp)
            nbad = 0
            for i in 1:n
                (_, w, tot) = hp[i]
                expected = w
                (2i <= n) && (expected += hp[2i][3])
                (2i + 1 <= n) && (expected += hp[2i + 1][3])
                tot == expected || (nbad += 1)
            end
            nbad
        end
        # `>` puts the largest weight first, `<` the smallest.
        for (cmp, isdesc) in ((>, true), (<, false))
            h = Air.PHeap{Int,Float64}(cmp)
            live = Int[]
            for step in 1:200
                r = mod(step, 4)
                if isempty(live) || r == 0 || r == 1
                    h = push(h, (step, float(mod(step * 7, 13) + 1)))
                    push!(live, step)
                elseif r == 2
                    h = Air.setweight(h, rand(live), float(mod(step * 5, 13) + 1))
                elseif r == 3 && mod(step, 8) == 3
                    k = first(h); h = pop(h); filter!(!=(k), live)
                else
                    k = rand(live); h = delete(h, k); filter!(!=(k), live)
                end
                @test badtotals(h) == 0
            end
            # the root total must equal the sum of the remaining weights
            @test Air.getfield(h, :_heap)[1][3] ==
                sum(Air.getweight(h, k) for k in live; init=0.0)
        end
    end

    @testset "weighted sampling is proportional to weight" begin
        # The user-visible consequence of drifting totals: draws are scaled by
        # the root total, so a wrong total skews them.
        h = Air.PHeap{Int,Float64}(>)
        h = push(h, (1, 1.0))
        h = push(h, (2, 9.0))
        h = push(h, (3, 5.0))
        h = delete(h, 3)                 # deleting used to be what broke them
        counts = Dict(1 => 0, 2 => 0)
        for _ in 1:20_000
            counts[rand(h)] += 1
        end
        ratio = counts[2] / max(counts[1], 1)
        @test 7.0 < ratio < 11.0         # expected 9.0
    end
end
