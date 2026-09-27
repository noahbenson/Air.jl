# Tests for the PArray type.
# Author: Noah C. Benson <n@nben.net>

@testset "PArray" begin
    numops = 100
    @testset "1D" begin
        a = Real[]
        p = Air.PVector{Real}()
        ops = [:push, :pushfirst, :pop, :popfirst, :set, :get]
        for ii in 1:numops
            q = rand(ops)
            if q == :push
                x = rand(Float64)
                push!(a, x)
                p = Air.push(p, x)
            elseif q == :pop && length(a) > 0
                x = pop!(a)
                @test x == p[end]
                p = Air.pop(p)
            elseif q == :pushfirst
                x = rand(Float64)
                pushfirst!(a, x)
                p = Air.pushfirst(p, x)
            elseif q == :popfirst && length(a) > 0
                x = popfirst!(a)
                @test x == p[1]
                p = Air.popfirst(p)
            elseif q == :set && length(a) > 0
                k = rand(1:length(a))
                v = rand(Float64)
                a[k] = v
                p = Air.setindex(p, v, k)
            elseif q == :get && length(a) > 0
                k = rand(1:length(a))
                @test a[k] == p[k]
            end
            @test a == p
            @test size(a) == size(p)
        end
    end
    @testset "2D" begin
        # These should be more fleshed out, but for now, we can do just a
        # few simple tests.
        a = Array(reshape(1:200, (20, 10)))
        p = PArray(a)
        @test p == a
        @test p[5, 8] == a[5, 8]
        @test p[:, 4] == a[:, 4]
        @test p[9, :] == a[9, :]
        p1 = Air.setindex(p, -5, 15, 8)
        @test p1 != a
        @test p1[:, 4] == a[:, 4]
        @test p1[9, :] == a[9, :]
        @test p1[:, 8] != a[:, 8]
        @test p1[15, :] != a[15, :]
    end
    @testset "3D" begin
        a = Array(reshape(1:1000, (20, 10, 5)))
        p = PArray(a)
        @test p == a
        @test p[5, 8, 2] == a[5, 8, 2]
        @test p[:, 4, 1] == a[:, 4, 1]
        @test p[9, :, 3] == a[9, :, 3]
        @test p[:, 2, 4] == a[:, 2, 4]
        p1 = Air.setindex(p, -5, 15, 8, 5)
        @test p1 != a
        @test p1[9, :, :] == a[9, :, :]
        @test p1[:, 3, :] == a[:, 3, :]
        @test p1[:, :, 4] == a[:, :, 4]
        @test p1[15, :, :] != a[15, :, :]
        @test p1[:, 8, :] != a[:, 8, :]
        @test p1[:, :, 5] != a[:, :, 5]
    end

@testset "generators" begin
    # `psparse` from positions and values, with a default of zero
    m = psparse([1, 2], [2, 3], [10.0, 20.0], 2, 3)
    @test m isa PArray{Float64,2}
    @test size(m) == (2, 3)
    @test m[1, 2] == 10.0 && m[2, 3] == 20.0 && m[1, 1] == 0.0
    @test nnz(m) == 2                       # the zeros are not stored
    # and from an array: the entries that are not its default
    @test collect(psparse([1, 0, 2])) == [1, 0, 2]
    # a plain `Array` has no default, so nothing is droppable and every entry is
    # kept; naming a default gives the `sparse(A)` behaviour
    @test nnz(psparse([1, 0, 2])) == 3
    @test nnz(psparse([1, 0, 2]; default = 0)) == 2
    @test Air.defaultvalue(psparse([1, 0, 2])) === undef
    @test Air.defaultvalue(psparse([1, 0, 2]; default = 0)) == 0
    @test nnz(psparse(setindex(PVector{Int}(0, (3,)), 5, 2))) == 1

    @test collect(pdiagm([1, 2, 3])) == [1 0 0; 0 2 0; 0 0 3]
    @test collect(pdiagm(1 => [1, 2])) == [0 1 0; 0 0 2; 0 0 0]
    @test nnz(pdiagm([1, 2, 3])) == 3

    # `permutedims` carries the default, so a sparse array stays sparse
    u = setindex(PMatrix(0.0, (2, 3)), 1.0, 2, 1)
    @test size(permutedims(u)) == (3, 2)
    @test collect(permutedims(u)) == permutedims(collect(u))
    @test Air.defaultvalue(permutedims(u)) == 0.0
    @test nnz(permutedims(u)) == 1
    @test collect(permutedims(u, (2, 1))) == permutedims(collect(u), (2, 1))
    @test_throws ArgumentError permutedims(u, (1, 1))

    b = SparseArrays.blockdiag(PMatrix(1.0, (2, 2)), PMatrix(2.0, (2, 2)))
    @test size(b) == (4, 4)
    # the first operand's default is the result's, so the off-diagonal positions
    # read as 1.0 rather than as zero
    @test collect(b) == [1 1 1 1; 1 1 1 1; 1 1 2 2; 1 1 2 2]
    @test Air.defaultvalue(b) == 1.0
    @test_throws DimensionMismatch SparseArrays.blockdiag(PVector([1, 2]), PMatrix(0.0, (2, 2)))

    # the random generators: a dense form, and a sparse one whose sparsity is the
    # probability each position is set
    Random.seed!(0x5eed)
    @test size(prand(2, 3)) == (2, 3)
    @test size(prandn(2, 3)) == (2, 3)
    @test size(prand(4, 4, 0.5)) == (4, 4)
    @test size(prandn(4, 4, 0.5)) == (4, 4)
    @test nnz(prand(4, 4, 0.0)) == 0
    @test nnz(prandn(4, 4, 0.0)) == 0
    @test nnz(prand(4, 4, 1.0)) == 16
    @test Air.defaultvalue(prand(4, 4, 0.5)) == 0.0
end
end
