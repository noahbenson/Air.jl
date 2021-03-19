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
        a = Array(reshape(1:200, (20,10)))
        p = PArray(a)
        @test p == a
        @test p[5,8] == a[5,8]
        @test p[:,4] == a[:,4]
        @test p[9,:] == a[9,:]
        p1 = Air.setindex(p, -5, 15, 8)
        @test p1 != a
        @test p1[:,4] == a[:,4]
        @test p1[9,:] == a[9,:]
        @test p1[:,8] != a[:,8]
        @test p1[15,:] != a[15,:]
    end
    @testset "3D" begin
        a = Array(reshape(1:1000, (20,10,5)))
        p = PArray(a)
        @test p == a
        @test p[5,8,2] == a[5,8,2]
        @test p[:,4,1] == a[:,4,1]
        @test p[9,:,3] == a[9,:,3]
        @test p[:,2,4] == a[:,2,4]
        p1 = Air.setindex(p, -5, 15, 8, 5)
        @test p1 != a
        @test p1[9,:,:] == a[9,:,:]
        @test p1[:,3,:] == a[:,3,:]
        @test p1[:,:,4] == a[:,:,4]
        @test p1[15,:,:] != a[15,:,:]
        @test p1[:,8,:] != a[:,8,:]
        @test p1[:,:,5] != a[:,:,5]
    end
end
