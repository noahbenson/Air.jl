using Air
using Test
using Random.Random
import Base.IdSet
import Base.delete!
import Base.isready

@testset "Air.jl" begin
    include("api.jl")
    include("parray.jl")
    include("pset.jl")
    include("pdict.jl")
    include("lazydict.jl")
    include("pwset.jl")
    include("pwdict.jl")
    include("TX.jl")
end
