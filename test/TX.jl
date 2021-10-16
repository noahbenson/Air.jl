# Tests for the Transaction system.
# Author: Noah C. Benson <n@nben.net>

@testset "TX" begin

    @testset "tx1" begin
        vals = Int[]
        vols = Vector{Volatile{Int}}[]
        log = Actor{typeof(stdout)}(stdout)
        worker() = begin
            global vals
            global vols
            while vols[end][1][] < 0
                # Pick a random layer and cell to calculate
                l = rand((2,3,4,5))
                layer = vols[l]
                k = rand(1:length(layer))
                # synconize for the rest.
                tx() do
                    #println("$l $k $(current_tx() === nothing) $(getfield(_current_tx, :values)))")
                    res = layer[k]
                    if res[] < 0
                        up = vols[l - 1]
                        a = up[k*2 - 1]
                        b = up[k*2]
                        if a[] >= 0 && b[] >= 0
                            res[] = a[] + b[]
                            send(log) do log
                                println("Layer $l, item $k complete $(objectid(current_task())).")
                                log
                            end
                        end
                    end
                end
            end
            send(log) do log
                println("Ending Woker $(objectid(current_task())).")
                log
            end
        end
        tx1_run(n::Int) = begin
            global vals
            global vols
            vals = rand(0:1000, 16)
            vols = Vector{Volatile{Int}}[
                [Volatile{Int}(u)   for u in vals],
                [Volatile{Int}(-1) for _ in 1:8],
                [Volatile{Int}(-1) for _ in 1:4],
                [Volatile{Int}(-1) for _ in 1:2],
                [Volatile{Int}(-1)]]
            # Spawn n threads.
            threads = [(Threads.@spawn worker()) for _ in 1:n]
            for th in threads
                Threads.wait(th)
            end
            println("  ---->  $(vols[end][1][]) == $(sum(vals))")
            return (vols[end][1][] == sum(vals))
        end
        @test tx1_run(1)
        @test tx1_run(2)
        @test tx1_run(3)
        @test tx1_run(4)
        @test tx1_run(5)
        @test tx1_run(6)
        @test tx1_run(7)
        @test tx1_run(8)
        @test tx1_run(9)
        @test tx1_run(10)
        @test tx1_run(11)
    end
    
end
