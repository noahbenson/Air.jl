module Air

include("Util.jl")
include("API.jl")
include("PTree.jl")
include("PList.jl")
include("PSet.jl")
include("PDict.jl")
include("LazyDict.jl")
include("PArray.jl")
include("PHeap.jl")
include("PWDict.jl")
include("PWSet.jl")
include("TX.jl")

# Exports; sorted by file:
#      PArray exports:
export PList, PSet, PDict, PArray, PVector, PMatrix, LazyDict,
    push, pop, pushfirst, popfirst, setindex, delete, #, insert
    Mutability, Immutable, Mutable, mutability, isimmtype,
    isequiv, equivhash, equalfn, hashfn, EquivSet, EquivDict,
    Var, Volatile, Actor, Source, tx, ReentrantRef,
    TransactionalRef, AbstractSourceKernel, getfilter,
    getfinalize, setfilter!, setfinalize!, send, geterror,
    receive, reset, PWDict, PWSet, getweight, setweight

end # module
