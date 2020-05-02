module Air

include("API.jl")
include("PTree.jl")
#include("PArray.jl")
#include("PMap.jl")
#include("PSet.jl")

# Exports; sorted by file:
#      PArray exports:
export PArray, PVector, TArray, TVector, PMatrix, TMatrix,
       push, pop, pushfirst, popfirst, setindex #, insert, delete,

end # module
