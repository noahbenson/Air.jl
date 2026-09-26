################################################################################
# broadcast.jl
#
# Tests for broadcasting over the persistent arrays.
#
# The behaviour to pin down is that a broadcast maps both halves of a `PArray` —
# its explicit entries *and* its default — and that the result is a `PArray`
# rather than a mutable `Array`. Values are compared against broadcasting the
# equivalent `Array`, so they are checked against Base's own semantics rather
# than against a restatement of the implementation.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "broadcasting" begin
    @testset "the result is a PArray" begin
        u = PVector{Int}(0, (4,))
        r = u .+ 1
        @test r isa PVector{Int}
        @test length(r) == 4
        @test collect(r) == collect(u) .+ 1
        # ... while `similar` still hands out mutable scratch space
        @test similar(u) isa Vector{Int}
        @test similar(u, Float64) isa Vector{Float64}
        @test similar(u, (2, 2)) isa Matrix{Int}
        # a plain array operand also yields a PArray
        r2 = u .+ [1, 2, 3, 4]
        @test r2 isa PVector{Int}
        @test collect(r2) == collect(u) .+ [1, 2, 3, 4]
        # and so does a scalar, on either side of the operator
        @test (1 .+ u) isa PVector{Int}
        @test collect(1 .+ u) == 1 .+ collect(u)
    end

    @testset "the default is mapped as well as the entries" begin
        u = PVector{Float64}(1.5, (5,))
        @test nnz(u) == 0
        r = u .+ 1
        @test r isa PVector{Float64}
        @test Air.defaultvalue(r) == 2.5
        @test nnz(r) == 0                   # still nothing stored
        @test collect(r) == fill(2.5, 5)
        # an explicit entry moves with the default, so their offset is preserved
        v = setindex(u, 10.0, 3)
        @test nnz(v) == 1
        rv = v .+ 1
        @test Air.defaultvalue(rv) == 2.5
        @test collect(rv) == collect(v) .+ 1
    end

    @testset "an entry that lands on the new default is dropped" begin
        # `abs` maps the default (-1.0) and the stored entry (1.0) onto the same
        # value, so after the broadcast the entry has nothing left to store.
        u = PVector{Float64}(-1.0, (5,))
        v = setindex(u, 1.0, 3)
        @test nnz(v) == 1
        r = abs.(v)
        @test Air.defaultvalue(r) == 1.0
        @test nnz(r) == 0
        @test collect(r) == abs.(collect(v))
        @test collect(r) == fill(1.0, 5)
        # An entry that stays distinct from the new default is kept. Note that an
        # entry equal to the array's *own* default is not stored in the first
        # place — `setindex` drops it — so a broadcast is the only way to test
        # this rule.
        w = setindex(u, 5.0, 2)
        @test nnz(w) == 1
        @test nnz(w .+ 1) == 1
        @test Air.defaultvalue(w .+ 1) == 0.0
        @test nnz(setindex(u, -1.0, 2)) == 0
    end

    @testset "two persistent arrays" begin
        a = setindex(PVector{Float64}(1.0, (6,)), 5.0, 2)
        b = setindex(PVector{Float64}(2.0, (6,)), 7.0, 4)
        r = a .+ b
        @test r isa PVector{Float64}
        @test Air.defaultvalue(r) == 3.0
        @test nnz(r) == 2                   # the union of the two, and only that
        # at position 2 only `a` has an entry, and at 4 only `b` does
        @test r[1] == 3.0 && r[2] == 7.0 && r[4] == 8.0
        @test collect(r) == collect(a) .+ collect(b)
    end

    @testset "a dense operand makes a dense result" begin
        u = setindex(PVector{Float64}(1.0, (4,)), 9.0, 2)
        r = u .+ [0.0, 0.0, 0.0, 0.0]
        @test r isa PVector{Float64}
        @test collect(r) == collect(u) .+ zeros(4)
        @test Air.defaultvalue(r) === undef     # no default, so every entry is stored
        @test nnz(r) == 4
        # A PArray built with `undef` is dense in the same way. Every position
        # must be set before it is read: such an array has no default, so an
        # unset position raises rather than yielding a value.
        @test_throws ErrorException PArray{Float64,1}(undef, (4,))[1]
        d = PArray{Float64,1}(undef, (4,))
        for k in 1:4
            d = setindex(d, float(k), k)
        end
        rd = d .+ 1
        @test nnz(rd) == 4
        @test collect(rd) == collect(d) .+ 1
    end

    @testset "dimensions" begin
        m = setindex(PMatrix{Float64}(1.0, (2, 3)), 4.0, 2, 2)
        rm = m .* 2.0
        @test rm isa PMatrix{Float64}
        @test size(rm) == (2, 3)
        @test Air.defaultvalue(rm) == 2.0
        @test collect(rm) == collect(m) .* 2.0
        # operands of different shapes, where one stored entry of the result
        # would correspond to several positions
        col = PVector{Float64}(10.0, (2,))
        r2 = m .+ col
        @test r2 isa PMatrix{Float64}
        @test size(r2) == (2, 3)
        @test collect(r2) == collect(m) .+ collect(col)
    end

    @testset "element types follow Base" begin
        u = PVector{Int}(0, (3,))
        @test (u .+ 1.0) isa PVector{Float64}
        @test (u ./ 2) isa PVector{Float64}
        @test (u .== 0) isa PVector{Bool}
        @test collect(u .== 0) == (collect(u) .== 0)
        @test (u .+ u) isa PVector{Int}
        # the default's type follows too
        @test Air.defaultvalue(u .+ 1.0) == 1.0
    end

    @testset "broadcasting into a PArray is refused" begin
        # A `PArray` is immutable, so `.=` has nowhere to write. It must raise
        # rather than quietly materialise a mutable copy.
        u = PVector{Int}(0, (3,))
        @test_throws Exception (u .= 1)
    end

    @testset "a transient operand makes a transient result" begin
        t = transient(PVector{Float64}(1.0, (5,)))
        t[2] = 5.0
        r = t .+ 1
        @test r isa TVector{Float64}
        @test !(r isa PArray)
        @test Air.defaultvalue(r) == 2.0
        @test collect(r) == collect(t) .+ 1
        # the operands are left alone
        @test length(t) == 5 && t[2] == 5.0
        # the result owns its own tree: persisting it hands back a collection
        # with no node still marked as owned
        @test _owned_count(persistent!(r)) == 0
        # a persistent operand alongside a transient one is still a transient
        # result, since the result is what a batch update continues to update
        p = PVector{Float64}(3.0, (5,))
        r2 = t .+ p
        @test r2 isa TVector{Float64}
        @test collect(r2) == collect(t) .+ collect(p)
        # writing into the result does not reach the transient it came from
        r2[1] = -1.0
        @test t[1] == 1.0
        @test length(t) == 5
    end
end
