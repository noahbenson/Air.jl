################################################################################
# variables.jl
#
# Tests for task-local `Var` objects and the `@var` macro.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "Var" begin
    @testset "construction and access" begin
        v = Var{Symbol}(:initval)
        @test v isa Var{Symbol}
        @test v[] == :initval
        # A Var may only be assigned through `withvars`/`setvars`.
        @test_throws DomainError (v[] = :other)
        @test v[] == :initval
    end

    @testset "withvars" begin
        v = Var{Symbol}(:initval)
        @test withvars(() -> v[], v => :newval) == :newval
        # The original binding is restored afterwards.
        @test v[] == :initval
        # Pair arguments and a dict argument are both accepted.
        @test withvars(() -> v[], IdDict{Var,Any}(v => :dictval)) == :dictval
        @test withvars() do
            v[] == :initval
        end
        # Nested bindings shadow and then unwind.
        withvars(v => :outer) do
            @test v[] == :outer
            withvars(v => :inner) do
                @test v[] == :inner
            end
            @test v[] == :outer
        end
    end

    @testset "setvars" begin
        v = Var{Symbol}(:initval)
        @test setvars(() -> v[], IdDict{Var,Any}(v => :newval)) == :newval
        @test v[] == :initval
    end

    @testset "vars" begin
        v = Var{Symbol}(:initval)
        @test vars() isa AbstractDict
        withvars(v => :bound) do
            @test vars()[v] == :bound
        end
    end

    @testset "task locality" begin
        v = Var{Int}(0)
        withvars(v => 42) do
            @test v[] == 42
            # A spawned task does not inherit the current task's bindings.
            @test fetch(Threads.@spawn v[]) == 0
        end
        @test v[] == 0
    end

    @testset "wrapwithvars" begin
        v = Var{Int}(0)
        f = wrapwithvars((x -> v[] + x), v => 10)
        @test f(5) == 15
        @test v[] == 0
    end
end

@testset "@var" begin
    @eval module VarMacroTest
        using Air
        @var mu = :start_sym
        @var mu_typed = :start_sym::Any
    end
    @test VarMacroTest.mu isa Var{Symbol}
    @test VarMacroTest.mu[] == :start_sym
    @test VarMacroTest.mu_typed isa Var{Any}
    @test VarMacroTest.mu_typed[] == :start_sym
    # `@var` declares a constant binding.
    @test isconst(VarMacroTest, :mu)
    # The error paths use `ArgumentError` (they previously referenced the
    # non-existent `ArgumentException`).
    @test_throws LoadError @eval @var 3 = 1
end
