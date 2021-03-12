using Air
using Test
using Random.Random
import Base.IdSet
import Base.delete!
import Base.isready

@testset "Air.jl" begin
    include("API.jl")
    include("PArray.jl")
    include("PSet.jl")
    include("PDict.jl")
    include("LazyDict.jl")
    include("PWSet.jl")
    include("PWDict.jl")
    include("TX.jl")
end
