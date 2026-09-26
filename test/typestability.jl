################################################################################
# typestability.jl
#
# Type-stability baselines. Each check asserts that a call's return type can be
# inferred, which is what lets the compiler generate specialised (fast) code.
#
# There are no `@test_broken` calls left here: every call below is inferred. The
# structural non-concreteness that does remain — an unparameterised `Union`
# field, an abstractly-typed `Function` field — is pinned in the "deferred:
# representation changes" testset as `fieldtype` assertions rather than as
# broken inference tests, since fixing it means changing the layout of those
# types.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

@testset "type stability" begin
    @testset "collections" begin
        @test @inferred(PDict(:a => 1)) isa PDict{Symbol,Int}
        @test @inferred(PDict{Symbol,Int}()) isa PDict{Symbol,Int}
        @test @inferred(push(PDict{Symbol,Int}(), :a => 1)) isa PDict{Symbol,Int}
        @test @inferred(PDict(:a => 1)[:a]) == 1
        @test @inferred(delete(PDict(:a => 1), :a)) isa PDict{Symbol,Int}

        @test @inferred(PSet([:a, :b])) isa PSet{Symbol}
        @test @inferred(push(PSet{Symbol}(), :a)) isa PSet{Symbol}
        @test @inferred(:a in PSet([:a])) === true
        @test @inferred(delete(PSet([:a, :b]), :a)) isa PSet{Symbol}

        @test @inferred(PVector([1, 2, 3])) isa PVector{Int}
        @test @inferred(push(PVector{Int}(), 1)) isa PVector{Int}
        @test @inferred(PVector([1, 2, 3])[2]) == 2
        @test @inferred(setindex(PVector([1, 2, 3]), 9, 2)) isa PVector{Int}
        @test @inferred(pop(PVector([1, 2, 3]))) isa PVector{Int}
        @test @inferred(PArray(0.0, (2, 2))) isa PArray{Float64,2}

        @test @inferred(PWSet{Symbol}(:a => 1.0)) isa PWSet{Symbol,Float64}
        @test @inferred(first(PWSet{Symbol}(:a => 1.0))) === :a
        @test @inferred(push(PWSet{Symbol,Float64}(), :a => 1.0)) isa PWSet{Symbol,Float64}
        @test @inferred(PWDict()) isa PWDict

        @test @inferred(LazyDict{Symbol,Int}()) isa LazyDict{Symbol,Int}
        @test @inferred(Delay{Int}(() -> 1)[]) === 1
    end

    @testset "@p macro" begin
        # `@inferred` needs a call expression, so these go through small
        # functions rather than being written inline.
        p1() = @p{:a => 1}
        p2() = @p[1, 2]
        p3() = @p(:a, :b)
        @test @inferred(p1()) isa PDict{Symbol,Int}
        @test @inferred(p2()) isa PVector{Int}
        @test @inferred(p3()) isa PSet{Symbol}
    end

    @testset "task-local variables" begin
        @test @inferred(Var{Int}(0)[]) === 0
        v = Var{Int}(0)
        @test @inferred(withvars(() -> v[], v => 1)) === 1
    end

    @testset "transactions" begin
        # Reading a volatile outside a transaction is inferred. This only
        # became true once `_volatile_getindex` was given type assertions: the
        # data stored for a `Volatile{T}` lives in dictionaries typed
        # `IdDict{Volatile,VolatileData}`, so without them the read inferred as
        # `Any` and every read was boxed.
        @test @inferred(Volatile{Int}(0)[]) === 0

        # Reading inside a transaction is inferred, and so is the transaction's
        # own result. `tx` used to return `Union{Nothing,T}`: the result was
        # assigned to a `res` initialised to `nothing` before the retry loop, so
        # every transaction's value was a small union. Binding the result with
        # the `try` expression fixed it (the `catch` arm never falls through).
        readintx() = begin
            v = Volatile{Int}(0)
            tx() do
                v[]
            end
        end
        @test @inferred(readintx()) === 0
        # The nested case returns the value too, through the enclosing
        # transaction rather than a new one.
        readnested() = begin
            v = Volatile{Int}(0)
            tx() do
                tx() do
                    v[]
                end
            end
        end
        @test @inferred(readnested()) === 0
        # A transaction whose body returns a non-Int is equally concrete, so the
        # fix is not specific to the Integer case.
        readstr() = begin
            v = Volatile{String}("x")
            tx() do
                v[]
            end
        end
        @test @inferred(readstr()) == "x"
    end

    @testset "deferred: representation changes" begin
        # These assertions are not failures; they pin down the structural
        # reasons the calls above cannot be inferred, so that a later change is
        # visible here. Fixing each means changing the layout of the type
        # involved, which is deliberately out of scope for this pass.
        #
        # 1. `PTree{T}.cells` is a three-way `Union`, so every traversal has to
        #    re-assert which vector it is holding.
        @test !isconcretetype(fieldtype(Air.PTree{Int}, :cells))
        # 2. `VolatileData`'s filter/finalize fields are typed `Function`, so
        #    calling them is a dynamic dispatch. Parameterizing the type would
        #    fix that, but it would ripple into `Volatile` and `Transaction` for
        #    a saving measured at a small fraction of the ~130 ns a filtered
        #    write costs over a plain one; and with no filter installed — the
        #    default — neither field is ever called. Left deliberately. (Note
        #    `PHeap` does *not* have this problem: its comparison function is a
        #    type parameter.)
        @test fieldtype(Air.VolatileData{Int}, :filter) === Union{Nothing,Function}
        @test fieldtype(Air.VolatileData{Int}, :finalize) === Union{Nothing,Function}
        # 3. `ActorMsg` stores its function in an abstractly-typed field.
        @test fieldtype(Air.ActorMsg, :fn) === Function
        # 4. The transaction maps are keyed by `Volatile` and valued by
        #    unparameterised `VolatileData`, so lookups are `Any`-shaped.
        @test valtype(fieldtype(Air.Transaction, :reads)) === Air.VolatileData
    end
end
