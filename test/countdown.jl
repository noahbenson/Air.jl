
# (1) The gold-standard solution:
function st_solve(k::Int, state::Symbol=:okay)
    r = Tuple{Symbol,Int}[]
    while k > 0
        push!(r, (state, k))
        if state == :okay
            k -= 3
            state = :ping
        elseif state == :ping
            k ÷= 2
            state = :pong
        elseif state == :pong
            k -= 4
            state = :okay
        end
    end
    return (r, k)
end

# (2) How to set up the data:
struct ProblemData
    value::Air.Volatile{Int}
    state::Air.Volatile{Symbol}
    actor::Air.Actor{Air.PVector{Tuple{Symbol,Int}}}
    function ProblemData(k::Int, state::Symbol=:okay)
        return new(Air.Volatile{Int}(k), Air.Volatile(state),
                   Air.Actor(PVector{Tuple{Symbol,Int}}()))
    end
end

# (3) Task functions:
function mt_task(pd::ProblemData, onstate, nextstate, op)
    go = true
    while go
        go = tx() do
            (pd.value[] <= 0) && return false
            yield()
            if pd.state[] == onstate
                let val = pd.value[]
                    Air.send(a -> Air.push(a, (onstate, val)), pd.actor)
                end
                pd.value[] = op(pd.value[])
                pd.state[] = nextstate
            end
            return true
        end
    end
end

# (4) How to solve the problem.
function mt_solve(k::Int, state::Symbol=:okay)
    pd = ProblemData(k, state)
    t1 = Threads.@spawn mt_task($pd, :pong, :okay, x -> x - 4)
    t2 = Threads.@spawn mt_task($pd, :ping, :pong, x -> x ÷ 2)
    t3 = Threads.@spawn mt_task($pd, :okay, :ping, x -> x - 3)
    wait(t1)
    wait(t2)
    wait(t3)
    sleep(0.1)
    return (pd.actor[], pd.value[])
end

function countdown_test(k::Int)
    (st_steps, st_val) = st_solve(k)
    (mt_steps, mt_val) = mt_solve(k)
    @test mt_val == st_val
    @test mt_steps == st_steps
end

@testset "CountDown" begin
    countdown_test(5)
    countdown_test(7)
    countdown_test(15)
    countdown_test(21)
    countdown_test(23)
    countdown_test(351)
    countdown_test(5016)
    countdown_test(18622)
    countdown_test(195122)
end
