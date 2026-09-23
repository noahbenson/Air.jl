using Air
using Test
using Random.Random
import Base.IdSet
import Base.delete!
import Base.isready

# Several testsets compare against randomized `Base` collections; seed the RNG
# so failures are reproducible.
Random.seed!(0x5eed)

@testset "Air.jl" begin
    include("api.jl")
    include("util.jl")
    include("regressions.jl")
    include("parray.jl")
    include("pset.jl")
    include("pdict.jl")
    include("lazydict.jl")
    include("pwset.jl")
    include("pwdict.jl")
    include("pheap.jl")
    include("variables.jl")
    include("TX.jl")
    include("countdown.jl")
    include("typestability.jl")
    include("aqua.jl")
end
