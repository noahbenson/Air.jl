################################################################################
# arrayops.jl
#
# Tests for the array operations that answer a persistent array with a
# persistent array.
#
# The property under test is that an operation over a persistent collection does
# not quietly stop being persistent, and that a sparse operand stays sparse: each
# operation carries the operand's default along. Values are compared against the
# same operation on the equivalent `Array`.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using LinearAlgebra: triu, rotl90

@testset "array operations keep the kind" begin
    u = PVector(collect(1:5))
    # a `PArray` is immutable, so sharing is a copy
    @test copy(u) === u
    @test copy(u) isa PVector{Int}
    @test map(x -> x + 1, u) isa PVector{Int}
    @test collect(map(x -> x + 1, u)) == collect(u) .+ 1

    # a sparse operand keeps its default through every one of these, so the
    # positions it does not store stay unstored
    s = setindex(setindex(PVector{Float64}(2.0, (6,)), 9.0, 2), 9.0, 5)
    @test nnz(s) == 2
    @test Air.defaultvalue(s) == 2.0
    m = map(x -> x + 1, s)
    @test m isa PVector{Float64}
    @test Air.defaultvalue(m) == 3.0
    @test nnz(m) == 2
    @test collect(m) == collect(s) .+ 1

    @testset "filter" begin
        f = filter(isodd, u)
        @test f isa PVector{Int}
        @test collect(f) == [1, 3, 5]
        @test collect(filter(iseven, u)) == filter(iseven, collect(u))
        # keeping an entry equal to the default does not store it
        keen = filter(x -> true, s)
        @test keen isa PVector{Float64}
        @test Air.defaultvalue(keen) == 2.0
        @test nnz(keen) == 2
        @test length(keen) == 6
        @test collect(keen) == collect(s)
        @test collect(filter(iseven, s)) == filter(iseven, collect(s))
        @test isempty(filter(x -> false, u))
    end

    @testset "reverse" begin
        r = reverse(u)
        @test r isa PVector{Int}
        @test collect(r) == reverse(collect(u))
        rs = reverse(s)
        @test rs isa PVector{Float64}
        @test Air.defaultvalue(rs) == 2.0
        @test nnz(rs) == 2
        @test length(rs) == 6
        @test collect(rs) == reverse(collect(s))
        @test reverse(reverse(s)) == s
    end

    @testset "indexing with a vector of positions" begin
        @test u[1:2] isa PVector{Int}
        @test collect(u[1:2]) == [1, 2]
        @test collect(u[[1, 3, 5]]) == [1, 3, 5]
        @test collect(u[5:-1:3]) == [5, 4, 3]
        @test length(u[2:4]) == 3
        @test collect(Int[]) == collect(u[Int[]])
        # a selected element equal to the default is not stored
        si = s[[1, 2, 6]]
        @test si isa PVector{Float64}
        @test Air.defaultvalue(si) == 2.0
        @test nnz(si) == 1
        @test collect(si) == collect(s)[[1, 2, 6]]
        # out-of-range positions are an error, as for an `Array`
        @test_throws BoundsError u[1:6]
        @test_throws BoundsError u[[0]]
    end

    @testset "indexing a matrix with vector positions" begin
        m = setindex(PMatrix(0.0, (3, 4)), 7.0, 2, 3)
        s = m[1:2, [1, 3]]
        @test s isa PArray{Float64,2}
        @test size(s) == (2, 2)
        @test Air.defaultvalue(s) == 0.0
        @test nnz(s) == 1                     # only the one entry it selects
        @test s[2, 2] == 7.0
        @test collect(s) == collect(m)[1:2, [1, 3]]
        # a scalar index drops its dimension, as it does for an `Array`
        r = m[1:3, 3]
        @test r isa PVector{Float64}
        @test length(r) == 3
        @test nnz(r) == 1
        @test collect(r) == collect(m)[1:3, 3]
        c = m[2, 1:4]
        @test c isa PVector{Float64}
        @test collect(c) == collect(m)[2, 1:4]
        # selecting positions that are not set stores nothing, and so does
        # selecting a position that holds the default
        @test nnz(m[1:2, [1, 2]]) == 0
        @test nnz(m[1:1, [1]]) == 0
        # agreement with Base, both for a sparsely-set matrix and a dense one
        d = PMatrix(0.0, (3, 3))
        for k in 1:3
            d = setindex(d, float(k), k, k)
        end
        @test collect(d[[3, 1], [2, 3]]) == collect(d)[[3, 1], [2, 3]]
        @test collect(d[2:3, 2:3]) == collect(d)[2:3, 2:3]
        @test collect(d[3, [1, 3]]) == collect(d)[3, [1, 3]]
        @test collect(d[[2], 1:2]) == collect(d)[[2], 1:2]
        # an out-of-range position is an error
        @test_throws BoundsError m[1:4, 1:4]
        @test_throws BoundsError m[1:2, [6]]
        # a transient still answers with a transient
        t = transient(m)
        @test t[1:2, [1, 3]] isa TArray{Float64,2}
        @test collect(t[1:2, [1, 3]]) == collect(m)[1:2, [1, 3]]
    end

    @testset "a transient answers with a transient" begin
        t = transient(u)
        @test map(x -> x + 1, t) isa TVector{Int}
        @test copy(t) isa TVector{Int}
        @test collect(map(x -> x + 1, t)) == collect(u) .+ 1
        @test length(t) == 5                 # the transient is unchanged
    end

    @testset "concatenation" begin
        a = PVector(collect(1:3))
        b = PVector(collect(4:6))
        c = vcat(a, b)
        @test c isa PVector{Int}
        @test length(c) == 6
        @test collect(c) == vcat(collect(a), collect(b))
        @test collect(vcat(a)) == collect(a)
        @test collect(hcat(PVector([1, 2]), PVector([3, 4]))) ==
              hcat([1, 2], [3, 4])
        # the first operand's default becomes the result's, so concatenating
        # sparse vectors stays sparse
        s = setindex(PVector{Float64}(2.0, (4,)), 9.0, 3)
        t = PVector{Float64}(2.0, (3,))
        sc = vcat(s, t)
        @test sc isa PVector{Float64}
        @test Air.defaultvalue(sc) == 2.0
        @test nnz(sc) == 1
        @test length(sc) == 7
        @test collect(sc) == vcat(collect(s), collect(t))
        # an operand with a different default is materialised: its default value
        # is an explicit value of the result
        u2 = PVector{Float64}(5.0, (3,))
        su = vcat(s, u2)
        @test Air.defaultvalue(su) == 2.0
        @test length(su) == 7
        @test nnz(su) == 4                   # s's one entry, plus u2's three
        @test collect(su) == vcat(collect(s), collect(u2))
        # a plain vector operand is materialised the same way
        @test collect(vcat(s, [1.0, 2.0])) == vcat(collect(s), [1.0, 2.0])
        # element types promote as they do for `Base`
        @test vcat(a, PVector([1.5])) isa PVector{Float64}
        @test collect(vcat(a, PVector([1.5]))) == [1, 2, 3, 1.5]
        # matrices: hcat stacks columns, vcat stacks rows
        m1 = setindex(PMatrix(0.0, (2, 2)), 1.0, 1, 1)
        m2 = setindex(PMatrix(0.0, (2, 3)), 2.0, 2, 2)
        h = hcat(m1, m2)
        @test h isa PArray{Float64,2}
        @test size(h) == (2, 5)
        @test Air.defaultvalue(h) == 0.0
        @test nnz(h) == 2
        @test collect(h) == hcat(collect(m1), collect(m2))
        # stacking rows maps positions through the operand's row count, so this
        # is the case where the shift is not constant
        v = vcat(m1, m1)
        @test v isa PArray{Float64,2}
        @test size(v) == (4, 2)
        @test Air.defaultvalue(v) == 0.0
        @test collect(v) == vcat(collect(m1), collect(m1))
        @test_throws DimensionMismatch hcat(m1, PMatrix(0.0, (3, 1)))
        @test_throws DimensionMismatch vcat(m1, PMatrix(0.0, (2, 3)))
        # a transient operand answers with a transient
        t1 = transient(a)
        @test vcat(t1, b) isa TVector{Int}
        @test collect(vcat(t1, b)) == collect(vcat(a, b))
        @test length(t1) == 3                # the transient is unchanged
    end

    @testset "arithmetic" begin
        # a sparse vector: default 1.0, with one stored entry
        a = setindex(PVector{Float64}(1.0, (4,)), 5.0, 2)
        b = setindex(PVector{Float64}(2.0, (4,)), 7.0, 4)
        @test collect(a) == [1.0, 5.0, 1.0, 1.0]
        s = a + b
        @test s isa PVector{Float64}
        @test Air.defaultvalue(s) == 3.0          # the defaults are added too
        @test s[2] == 7.0 && s[4] == 8.0          # 5 + 2, and 1 + 7
        @test nnz(s) == 2
        @test collect(s) == collect(a) + collect(b)
        for (r, expect) in ((a - b, collect(a) - collect(b)),)
            @test r isa PVector{Float64}
            @test collect(r) == expect
        end
        @test (-a) isa PVector{Float64}
        @test collect(-a) == -collect(a)
        @test Air.defaultvalue(-a) == -1.0
        # scaling and division by a number are elementwise
        for (r, expect) in (
            (2a, 2 .* collect(a)), (a * 2, collect(a) .* 2), (a / 2, collect(a) ./ 2)
        )
            @test r isa PVector{Float64}
            @test collect(r) == expect
        end
        @test Air.defaultvalue(2a) == 2.0
        # with a scalar, on either side
        for (r, expect) in (
            (a + 1, collect(a) .+ 1), (1 + a, 1 .+ collect(a)),
            (a - 1, collect(a) .- 1), (1 - a, 1 .- collect(a)),
        )
            @test r isa PVector{Float64}
            @test collect(r) == expect
        end
        # with a plain array, on either side: materialised, as in a broadcast
        for (r, expect) in (
            (a + [1.0, 1, 1, 1], collect(a) + [1, 1, 1, 1]),
            ([1.0, 1, 1, 1] + a, [1, 1, 1, 1] + collect(a)),
        )
            @test r isa PVector{Float64}
            @test collect(r) == expect
        end
        # matrices are elementwise the same
        m = setindex(PMatrix(0.0, (2, 2)), 1.0, 1, 1)
        @test (m + m) isa PMatrix{Float64}
        @test collect(m + m) == collect(m) + collect(m)
        @test (2m) isa PMatrix{Float64}
        # ... but multiplication by a matrix is not elementwise, so it stays
        # Base's, which answers with a mutable `Matrix`
        @test (m * m) isa Matrix{Float64}
        @test collect(m * m) == collect(m) * collect(m)
        # a transient operand makes a transient result
        t = transient(a)
        @test (t + t) isa TVector{Float64}
        @test (2t) isa TVector{Float64}
        @test collect(t + t) == collect(a) + collect(a)
        # the operands are unchanged
        @test collect(a) == [1.0, 5.0, 1.0, 1.0]
        @test length(a) == 4 && Air.defaultvalue(a) == 1.0
    end

    @testset "reshape" begin
        u = setindex(PVector{Float64}(1.0, (6,)), 9.0, 4)
        m = reshape(u, (2, 3))
        @test m isa PArray{Float64,2}
        @test size(m) == (2, 3)
        @test Air.defaultvalue(m) == 1.0
        @test nnz(m) == 1
        @test collect(m) == reshape(collect(u), 2, 3)
        @test m[2, 2] == 9.0                    # the linear order is preserved
        # reshape changes only the shape: it is the same tree
        @test Air.getfield(m, :_tree) === Air.getfield(u, :_tree)
        @test collect(reshape(m, (6,))) == collect(u)
        # dimensions by tuple and by vararg
        @test size(reshape(u, (6,))) == (6,)
        @test size(reshape(u, 2, 3)) == (2, 3)
        # `:` takes up the slack, as it does for an Array
        for dims in ((2, 3), (3, :), (:, 2), (1, 6), (6, 1), (2, 3, 1), (1, :, 2))
            @test size(reshape(u, dims)) == size(reshape(collect(u), dims))
            @test collect(reshape(u, dims)) == reshape(collect(u), dims)
        end
        # a shape that does not fit is an error, as for an Array
        for dims in ((4, 2), (:, :), (5, :), (0, :), (2, 2))
            @test_throws DimensionMismatch reshape(u, dims)
        end
        # The empty cases are asserted directly rather than against Base, because
        # Base is not consistent across the versions Air supports: Julia 1.10 and
        # 1.11 raise a `DivideError` for `reshape(Int[], (0, :))` and 1.13 gives
        # `(0, 0)`. Air follows the latter on every version — a zero among the
        # known dimensions pins the colon to zero, provided there is nothing to
        # hold — so its own answers are what is checked.
        e = PVector{Int}()
        @test size(reshape(e, (0,))) == (0,)
        @test size(reshape(e, (0, 1))) == (0, 1)
        @test size(reshape(e, (0, :))) == (0, 0)
        @test size(reshape(e, (:, 0))) == (0, 0)
        @test size(reshape(e, (1, :))) == (1, 0)
        # a transient keeps Base's mutable view, as an Array does
        t = transient(u)
        r = reshape(t, (2, 3))
        @test parent(r) === t
        r[1, 1] = -1.0
        @test t[1] == -1.0
    end

@testset "repetition and cat" begin
    v = PVector([1, 2, 3])
    u = collect(v)
    @test collect(repeat(v, 2)) == repeat(u, 2)
    @test repeat(v, 2) isa PVector{Int}
    @test collect(repeat(v; inner = 2)) == repeat(u; inner = 2)
    @test collect(repeat(v; outer = 3)) == repeat(u; outer = 3)
    @test collect(repeat(v; inner = 2, outer = 3)) == repeat(u; inner = 2, outer = 3)
    @test collect(repeat(v, 0)) == Int[]
    @test collect(repeat(v, 1)) == u
    @test_throws ArgumentError repeat(v; inner = -1)
    @test_throws ArgumentError repeat(v, -1)

    @test collect(cat(v, v; dims = 1)) == cat(u, u; dims = 1)
    @test cat(v, v; dims = 1) isa PVector{Int}
    @test collect(cat(v, v, v; dims = 1)) == cat(u, u, u; dims = 1)
    m = setindex(PMatrix(0.0, (2, 3)), 1.0, 1, 1)
    @test collect(cat(m, m; dims = 1)) == cat(collect(m), collect(m); dims = 1)
    @test collect(cat(m, m; dims = 2)) == cat(collect(m), collect(m); dims = 2)
    @test cat(m, m; dims = 2) isa PMatrix{Float64}
    @test_throws ArgumentError cat(m, m; dims = (1, 2))
    @test_throws ArgumentError cat(m, m; dims = 3)

    # a sparse vector keeps its default and its length through a repeat
    s = setindex(PVector{Float64}(0.0, (3,)), 5.0, 2)
    r = repeat(s, 2)
    @test length(r) == 6
    @test Air.defaultvalue(r) == 0.0
    @test nnz(r) == 2
    @test collect(r) == repeat(collect(s), 2)
end

@testset "reordering the values" begin
    v = PVector([3, 1, 4, 1, 5, 9, 2, 6])
    u = collect(v)

    @test collect(sort(v)) == sort(u)
    @test sort(v) isa PVector{Int}
    @test collect(sort(v; rev = true)) == sort(u; rev = true)
    @test collect(sort(v; by = x -> -x)) == sort(u; by = x -> -x)

    @test collect(unique(v)) == unique(u)
    @test unique(v) isa PVector{Int}

    @test collect(circshift(v, 2)) == circshift(u, 2)
    @test collect(circshift(v, -3)) == circshift(u, -3)
    @test collect(circshift(v, 0)) == u
    @test collect(circshift(v, 8)) == u              # a full turn

    @test collect(deleteat(v, 3)) == [u[1:2]; u[4:end]]
    @test collect(deleteat(v, 1)) == u[2:end]
    @test collect(deleteat(v, length(v))) == u[1:(end - 1)]
    @test collect(splice(v, 2, 3)) == [u[1]; u[5:end]]
    @test collect(splice(v, 1, 0)) == u              # nothing removed
    @test_throws BoundsError deleteat(v, 0)
    @test_throws BoundsError deleteat(v, 9)
    @test_throws BoundsError splice(v, 5, 5)
    @test collect(circshift(PVector{Int}(), 3)) == Int[]

    # A sparse vector keeps its default *and* its length: the positions it does
    # not store are positions, not absences, so sorting must move them with
    # everything else rather than collapse the result to the stored entries.
    s = setindex(setindex(PVector{Float64}(0.0, (6,)), 3.0, 2), 1.0, 5)
    @test length(sort(s)) == 6
    @test Air.defaultvalue(sort(s)) == 0.0
    @test nnz(sort(s)) == 2
    @test collect(sort(s)) == sort(collect(s))
    @test collect(unique(s)) == unique(collect(s))
    @test collect(circshift(s, 2)) == circshift(collect(s), 2)
    @test Air.defaultvalue(circshift(s, 2)) == 0.0
    @test length(circshift(s, 2)) == 6
    @test collect(deleteat(s, 2)) == [collect(s)[1:1]; collect(s)[3:end]]
    @test length(deleteat(s, 2)) == 5
end

@testset "rotation and the triangular part" begin
    # `rotl90` both rotates and transposes, so the result's shape is the
    # argument's reversed. `triu` keeps the shape and is checked at every `k`,
    # because the diagonal is where its two regimes meet.
    a = [1 2 3; 4 5 6]
    m = PMatrix(a)

    r = rotl90(m)
    @test r isa PMatrix{Int}
    @test size(r) == (3, 2)
    @test collect(r) == rotl90(a)
    @test Air.defaultvalue(r) === UndefInitializer()
    @test collect(rotl90(permutedims(m))) == rotl90(permutedims(a))

    t = triu(m)
    @test t isa PMatrix{Int}
    @test size(t) == (2, 3)
    @test collect(t) == triu(a)
    @test Air.defaultvalue(t) === UndefInitializer()
    for k in -3:3
        @test collect(triu(m, k)) == triu(a, k)
    end

    # A sparse matrix with a zero default. What lies below the diagonal already
    # reads as the default, so `triu` only has to drop the entries that are
    # there, and the operation costs the entries rather than the shape.
    s = setindex(setindex(PMatrix(0.0, (3, 3)), 1.0, 2, 1), 2.0, 1, 2)
    st = triu(s)
    @test st isa PMatrix{Float64}
    @test collect(st) == triu(collect(s))
    @test Air.defaultvalue(st) == 0.0
    @test nnz(st) == 1                     # the below-diagonal 1.0 is gone
    # rotation carries both the entries and the default without densifying
    sr = rotl90(s)
    @test collect(sr) == rotl90(collect(s))
    @test Air.defaultvalue(sr) == 0.0
    @test nnz(sr) == nnz(s)

    # A non-zero default is the other regime. `triu` *zeroes* the entries below
    # the diagonal, and zero is not what an unset position reads as, so each of
    # those positions has to be stored explicitly — this is the costly case, and
    # it is inherent to the operation rather than to the representation.
    n = setindex(PMatrix(7.0, (3, 3)), 1.0, 1, 1)
    nt = triu(n)
    @test collect(nt) == triu(collect(n))
    @test Air.defaultvalue(nt) == 7.0
    @test nt[2, 1] == 0.0
    @test nt[3, 3] == 7.0                  # still the default, so still unset
    @test nnz(nt) == 4                     # the 1.0 and three explicit 0.0s
    for k in -3:3
        @test collect(triu(n, k)) == triu(collect(n), k)
    end
end
end
