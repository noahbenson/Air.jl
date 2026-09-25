################################################################################
# transient.jl
#
# Tests for the transient (mutable) collections.
#
# The behaviour to pin down is that a transient is a *view* for batch updates
# and never a licence to change a collection the caller still holds. Two things
# follow, and neither is visible from the results of ordinary operations:
#
#  * updates through a transient must leave the collection it was made from
#    alone, even after `persistent!` has handed back a collection that shares
#    its nodes; and
#  * a collection produced by `persistent!` must contain no node marked as owned
#    by a transient, since such a node could be changed in place by a later
#    transient. That is what `persistent!`'s walk over the owned nodes is for,
#    and it is asserted directly here.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

"""
    _owned_count(u)

The number of nodes in the tree `u` that a transient owns. Should be zero for
any tree reachable from a collection obtained from `persistent!`.
"""
function _owned_count(u::Air.PTree{T}) where {T}
    id = Air.getfield(u, :id)
    n = Air.ptree_owned(id) ? 1 : 0
    cells = Air.getfield(u, :cells)
    (cells === nothing || Air.ptree_depth(id) == Air.PTREE_TWIG_DEPTH) && return n
    for c in cells::Vector{Air.PTree{T}}
        n += _owned_count(c)
    end
    return n
end
_owned_count(v::PVector) = _owned_count(Air.getfield(v, :_tree))
_owned_count(d::Union{Air.AbstractPDict,Air.AbstractPSet}) =
    _owned_count(Air.getfield(d, :root))
@testset "transients" begin
    @testset "a transient equals the persistent equivalent" begin
        for n in (0, 1, 64, 65, 1000)
            base = PVector(collect(1:n))
            t = transient(base)
            for i in 1:n
                push!(t, i)
            end
            out = persistent!(t)
            @test length(out) == 2n
            @test all(out[i] == i for i in 1:n)
            @test all(out[n + i] == i for i in 1:n)
            # and the same as an equivalent persistent build
            ref = base
            for i in 1:n
                ref = push(ref, i)
            end
            @test length(out) == length(ref)
            @test all(out[i] == ref[i] for i in 1:length(ref))
        end
    end

    @testset "the source collection is left alone" begin
        base = PVector(collect(1:500))
        before = collect(base)
        t = transient(base)
        for i in 1:1000
            push!(t, -i)
        end
        for _ in 1:100
            pop!(t)
        end
        out = persistent!(t)
        @test collect(base) == before
        @test length(out) == 500 + 1000 - 100
        @test all(out[i] == i for i in 1:500)
    end

    @testset "persistent! leaves no node owned" begin
        t = transient(PVector{Int}())
        for i in 1:2000
            push!(t, i)
        end
        out = persistent!(t)
        @test _owned_count(out) == 0
        # ... including when the batch started from an existing collection, and
        # when nothing was pushed at all.
        t2 = transient(PVector(collect(1:100)))
        @test _owned_count(persistent!(t2)) == 0
        t3 = transient(PVector{Int}())
        push!(t3, 1)
        @test _owned_count(persistent!(t3)) == 0
    end

    @testset "a second transient cannot write through the first" begin
        # This is the case the ownership flag exists to prevent: the collection
        # `persistent!` returns shares nodes with the transient that made it, so
        # if those nodes stayed marked as owned, a later transient would treat
        # them as its own and change them in place.
        t1 = transient(PVector{Int}())
        for i in 1:1000
            push!(t1, i)
        end
        b = persistent!(t1)
        before = collect(b)
        t2 = transient(b)
        for i in 1:500
            push!(t2, 10_000 + i)      # appends land in b's last twig
        end
        @test collect(b) == before
        @test length(b) == 1000
        b2 = persistent!(t2)
        @test length(b2) == 1500
        @test all(b2[i] == i for i in 1:1000)
        @test all(b2[1000 + i] == 10_000 + i for i in 1:500)
    end

    @testset "the transient and persistent! are cheap in the size of the tree" begin
        # Both are meant to be O(1) in the collection's size: making a transient
        # shares structure, and persisting hands it back. If either started
        # copying, a batch update over an existing collection would lose more
        # than it gains.
        small = PVector(collect(1:10))
        large = PVector(collect(1:100_000))
        tsmall = transient(small)
        tlarge = transient(large)
        @test (@allocated transient(large)) <= (@allocated transient(small)) + 64
        @test (@allocated persistent!(tlarge)) <= (@allocated persistent!(tsmall)) + 64
    end

    @testset "TDict" begin
        for n in (0, 1, 100, 1000)
            base = PDict{Symbol,Int}()
            for i in 1:n
                base = push(base, Symbol("k", i) => i)
            end
            # a batch of updates, half overwriting and half new
            t = transient(base)
            @test t isa TDict{Symbol,Int}
            for i in 1:n
                t[Symbol("k", i)] = 1000 + i
            end
            for i in 1:n
                push!(t, Symbol("new", i) => i)
            end
            out = persistent!(t)
            @test out isa PDict{Symbol,Int}
            @test length(out) == 2n
            @test all(out[Symbol("k", i)] == 1000 + i for i in 1:n)
            @test all(out[Symbol("new", i)] == i for i in 1:n)
            # and it matches the persistent equivalent
            ref = base
            for i in 1:n
                ref = push(ref, Symbol("k", i) => 1000 + i)
                ref = push(ref, Symbol("new", i) => i)
            end
            @test length(out) == length(ref)
            @test all(out[k] == v for (k, v) in ref)
        end
    end

    @testset "TDict leaves the source alone, and no node owned" begin
        base = PDict{Symbol,Int}()
        for i in 1:500
            base = push(base, Symbol("k", i) => i)
        end
        t = transient(base)
        for i in 1:500
            t[Symbol("k", i)] = -i
            delete!(t, Symbol("k", i + 250))
        end
        out = persistent!(t)
        @test all(base[Symbol("k", i)] == i for i in 1:500)
        @test length(base) == 500
        @test _owned_count(out) == 0
        # a second transient over the result cannot write through to it
        before = Dict(out)
        t2 = transient(out)
        for i in 1:100
            t2[Symbol("t2", i)] = i
        end
        @test Dict(out) == before
    end

    @testset "transient delete matches persistent delete" begin
        function _mkdict(n)
            d = PDict{Symbol,Int}()
            for i in 1:n
                d = push(d, Symbol("k", i) => i)
            end
            return d
        end
        for n in (0, 1, 100, 1000)
            base = _mkdict(n)
            t = transient(base)
            for i in 1:(n ÷ 2)
                delete!(t, Symbol("k", i))
            end
            out = persistent!(t)
            ref = base
            for i in 1:(n ÷ 2)
                ref = delete(ref, Symbol("k", i))
            end
            @test length(out) == length(ref) == n - n ÷ 2
            @test all(out[k] == v for (k, v) in ref)
            @test all(!haskey(out, Symbol("k", i)) for i in 1:(n ÷ 2))
            @test all(base[Symbol("k", i)] == i for i in 1:n)  # source untouched
            @test _owned_count(out) == 0
            # The result must satisfy the same structural invariants a
            # persistently-built tree does, the minimal-tree collapse included —
            # `_ptree_check` walks it and asserts them.
            @test _ptree_check(Air.getfield(out, :root)) == n - n ÷ 2
        end
    end

    @testset "TDict accessors" begin
        t = transient(PDict{Symbol,Int}(:a => 1))
        @test length(t) == 1
        @test t[:a] == 1
        @test haskey(t, :a)
        @test !haskey(t, :b)
        @test get(t, :b, -1) == -1
        @test_throws KeyError t[:b]
        t[:b] = 2
        @test t[:b] == 2 && length(t) == 2
        t[:b] = 3                        # overwrite, not insert
        @test t[:b] == 3 && length(t) == 2
        delete!(t, :b)
        @test !haskey(t, :b) && length(t) == 1
        delete!(t, :nothere)             # deleting a missing key is a no-op
        @test length(t) == 1
    end

    @testset "TIdDict keeps keys by identity" begin
        t = transient(PIdDict{Any,Int}())
        @test t isa TIdDict{Any,Int}
        x, y = [1], [1]                  # equal but not identical
        t[x] = 10
        t[y] = 20
        @test length(t) == 2
        @test t[x] == 10 && t[y] == 20
        out = persistent!(t)
        @test out isa PIdDict{Any,Int}
        @test length(out) == 2
        @test out[x] == 10 && out[y] == 20
    end

    @testset "TSet" begin
        function _mkset(n)
            s = PSet{Symbol}()
            for i in 1:n
                s = push(s, Symbol("k", i))
            end
            return s
        end
        for n in (0, 1, 100, 1000)
            base = _mkset(n)
            t = transient(base)
            @test t isa TSet{Symbol}
            @test length(t) == n
            for i in 1:n
                delete!(t, Symbol("k", i))          # remove half...
                i <= n ÷ 2 && push!(t, Symbol("k", i))  # ...and put half back
            end
            for i in 1:n
                push!(t, Symbol("new", i))
            end
            out = persistent!(t)
            @test out isa PSet{Symbol}
            @test length(out) == n + n ÷ 2
            @test all(Symbol("new", i) in out for i in 1:n)
            @test all(Symbol("k", i) in out for i in 1:(n ÷ 2))
            @test all(!(Symbol("k", i) in out) for i in ((n ÷ 2) + 1):n)
            # matches the persistent equivalent
            ref = base
            for i in 1:n
                ref = delete(ref, Symbol("k", i))
                i <= n ÷ 2 && (ref = push(ref, Symbol("k", i)))
            end
            for i in 1:n
                ref = push(ref, Symbol("new", i))
            end
            @test length(out) == length(ref)
            @test all(x in ref for x in out)
            # the source is untouched, no node is left owned, and the result has
            # the same structure a persistent build would
            @test length(base) == n
            @test _owned_count(out) == 0
            @test _ptree_check(Air.getfield(out, :root)) == n + n ÷ 2
        end

        # duplicates are no-ops, and a second transient cannot write through
        t = transient(PSet{Symbol}())
        push!(t, :a)
        push!(t, :a)
        @test length(t) == 1
        b = persistent!(t)
        t2 = transient(b)
        for i in 1:100
            push!(t2, Symbol("x", i))
        end
        @test length(b) == 1 && :a in b
        @test length(persistent!(t2)) == 101
    end

    @testset "TIdSet keeps elements by identity" begin
        t = transient(PIdSet{Any}())
        @test t isa TIdSet{Any}
        x, y = [1], [1]                     # equal but not identical
        push!(t, x)
        push!(t, y)
        @test length(t) == 2
        @test x in t && y in t
        out = persistent!(t)
        @test out isa PIdSet{Any}
        @test length(out) == 2
        # deleting one of two equal-but-distinct elements through a transient
        t3 = transient(out)
        delete!(t3, x)
        out2 = persistent!(t3)
        @test length(out2) == 1
        @test !(x in out2) && y in out2
        @test length(out) == 2               # the original is untouched
    end

    @testset "pop!" begin
        t = transient(PVector(collect(1:200)))
        for _ in 1:190
            pop!(t)
        end
        out = persistent!(t)
        @test length(out) == 10
        @test all(out[i] == i for i in 1:10)
        t2 = transient(PVector{Int}())
        @test_throws ArgumentError pop!(t2)
    end

    @testset "the shape survives a round trip" begin
        # A `TArray` carries its source's shape, so an array that is not a vector
        # comes back as the same shape. Only the vector operations can change a
        # transient's length, so nothing else can invalidate the shape.
        m = PMatrix(0.0, (5, 7))
        m = setindex(m, 2.5, 3, 4)
        tm = transient(m)
        @test tm isa TArray{Float64,2}
        @test length(tm) == 35
        @test persistent!(tm) == m
        # and a vector's shape follows the length it was pushed to
        t = transient(PVector(collect(1:10)))
        @test t isa TVector{Int}
        for i in 1:90
            push!(t, i)
        end
        out = persistent!(t)
        @test size(out) == (100,)
        @test all(out[i] == i for i in 1:10)
    end

    @testset "appending equal to the default needs no entry" begin
        # PVector's own push skips the tree when the value equals the array's
        # default; the transient must do the same, or it would build a tree of
        # identical values.
        base = PVector{Int}(0, (5,))
        t = transient(base)
        push!(t, 0)
        out = persistent!(t)
        @test length(out) == 6
        @test all(out[i] == 0 for i in 1:6)
        @test isempty(Air.getfield(Air.getfield(out, :_tree), :cells) === nothing ? () :
                      Air.getfield(Air.getfield(out, :_tree), :cells))
    end
end
