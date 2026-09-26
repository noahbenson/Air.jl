################################################################################
# ptree.jl
#
# Structural invariants of the persistent tree that backs every collection.
# `PTree` is internal, so these reach into it directly rather than through a
# collection's interface.
#
# The invariant worth testing here is the *minimal tree*: a branch always has at
# least two children, so a branch left holding a single child is replaced by
# that child. Without it a sequence of deletions leaves chains of one-child
# branches behind, and every later lookup and insertion pays an extra level per
# chain link. Nothing else observes it: a non-minimal tree still answers every
# query correctly, it is just slower, which is exactly the kind of thing that
# needs an explicit test rather than a behavioural one.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

"""
    _ptree_check(u)

Walks the `PTree` `u`, testing its structural invariants at every node, and
yields the number of leaves beneath it:

 * every branch has at least two children (the minimal-tree invariant);
 * a node's occupancy bitmap has exactly as many set bits as it has cells;
 * a twig's cell count matches its occupancy and its leaf count;
 * every child lies within its parent's leaf range; and
 * a node's cached leaf count is the sum over its children.
"""
function _ptree_check(u::Air.PTree{T}) where {T}
    id = Air.getfield(u, :id)
    bits = Air.getfield(u, :bits)
    n = Air.getfield(u, :numel)
    nocc = count_ones(bits)
    # No node reachable from a persistent collection may be marked as owned by a
    # transient: that flag is what licenses an in-place change, so a persistent
    # node carrying it could be mutated out from under its other holders.
    @test !Air.ptree_owned(id)
    # The canonical empty node holds no cells at all.
    if n == 0
        @test getfield(u, :cells) === nothing
        return 0
    end
    if Air.ptree_depth(id) == Air.PTREE_TWIG_DEPTH
        cells = Air.getfield(u, :cells)::Vector{T}
        @test length(cells) == nocc == n
        return n
    end
    cells = Air.getfield(u, :cells)::Vector{Air.PTree{T}}
    @test length(cells) == nocc
    @test nocc >= 2
    total = 0
    for c in cells
        cid = Air.getfield(c, :id)
        # the child lies beneath this node, and nowhere else
        @test Air.ptree_isbeneath(id, Air.ptree_minleaf(cid))
        total += _ptree_check(c)
    end
    @test total == n
    return n
end

@testset "PTree invariants" begin
    @testset "minimal tree after deletions" begin
        # Deleting keys is what can leave a branch with one child, so exercise
        # mixed insertion and deletion, including deleting down to a very small
        # tree and back up again.
        d = PDict{Symbol,Int}()
        for i in 1:512
            d = push(d, Symbol("k", i) => i)
            _ptree_check(Air.getfield(d, :root))
        end
        # remove every other key, checking at each step
        for i in 1:2:512
            d = delete(d, Symbol("k", i))
            _ptree_check(Air.getfield(d, :root))
        end
        @test length(d) == 256
        # and the survivors are still reachable, with the right values
        for i in 2:2:512
            @test d[Symbol("k", i)] == i
        end
        # drain it completely, then refill: both directions must stay minimal
        for i in 2:2:512
            d = delete(d, Symbol("k", i))
        end
        @test length(d) == 0
        @test Air.getfield(Air.getfield(d, :root), :numel) == 0
        for i in 1:64
            d = push(d, Symbol("r", i) => i)
            _ptree_check(Air.getfield(d, :root))
        end
        @test length(d) == 64
    end

    @testset "minimal tree under random operations" begin
        # A deterministic pseudo-random sequence of insert/delete, so that
        # arbitrary interleavings are covered without depending on the RNG
        # state left behind by other testsets.
        rng = Random.MersenneTwister(0xa17)
        d = PDict{Int,Int}()
        live = Set{Int}()
        for step in 1:2000
            k = rand(rng, 1:300)
            if k in live && rand(rng, Bool) && length(live) > 1
                d = delete(d, k)
                delete!(live, k)
            else
                d = push(d, k => k)
                push!(live, k)
            end
            _ptree_check(Air.getfield(d, :root))
        end
        @test length(d) == length(live)
        for k in live
            @test d[k] == k
        end
    end

    @testset "collections still agree with Base" begin
        # The structure changed, so re-check that behaviour is untouched across
        # the collection types that sit on top of the tree.
        rng = Random.MersenneTwister(0xbee)
        for _ in 1:20
            d = PDict{Symbol,Int}()
            ref = Dict{Symbol,Int}()
            for _ in 1:200
                k = Symbol("v", rand(rng, 1:60))
                if rand(rng) < 0.4 && !isempty(ref)
                    k = rand(rng, collect(keys(ref)))
                    d = delete(d, k)
                    delete!(ref, k)
                else
                    v = rand(rng, 1:1000)
                    d = push(d, k => v)
                    ref[k] = v
                end
            end
            @test length(d) == length(ref)
            @test all(d[k] == v for (k, v) in ref)
            @test sort(collect(keys(d))) == sort(collect(keys(ref)))

            s = PSet{Int}()
            refs = Set{Int}()
            for _ in 1:200
                x = rand(rng, 1:60)
                if rand(rng) < 0.4 && !isempty(refs)
                    x = rand(rng, collect(refs))
                    s = delete(s, x)
                    delete!(refs, x)
                else
                    s = push(s, x)
                    push!(refs, x)
                end
            end
            @test sort(collect(s)) == sort(collect(refs))
        end
    end
end
