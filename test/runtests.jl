using Air
using Test
using Random.Random
# `Air` extends SparseArrays' functions (`nnz`, for one), so the tests that check
# that interop need the functions themselves in scope, not just Air.
using SparseArrays: nnz
import Base.IdSet
import Base.delete!
import Base.isready

# Several testsets compare against randomized `Base` collections; seed the RNG
# so failures are reproducible.
Random.seed!(0x5eed)

@testset "Air.jl" begin
    include("api.jl")
    include("ptree.jl")
    include("util.jl")
    include("regressions.jl")
    include("parray.jl")
    include("pset.jl")
    include("pdict.jl")
    include("lazydict.jl")
    include("pwset.jl")
    include("pwdict.jl")
    include("pheap.jl")
    include("transient.jl")
    # after `transient.jl`, which defines `_owned_count`
    include("broadcast.jl")
    include("variables.jl")
    include("TX.jl")
    include("countdown.jl")
    include("iteration.jl")
    include("typestability.jl")
    include("aqua.jl")
end
