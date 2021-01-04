module Air

include("Util.jl")
include("API.jl")
include("PTree.jl")
include("PList.jl")
include("PLinearSet.jl")
include("PLinearDict.jl")
include("PSet.jl")
include("PDict.jl")
include("LazyDict.jl")
include("PArray.jl")
include("PHeap.jl")
include("PWDict.jl")
include("PWSet.jl")
include("TX.jl")

include("PTrie.jl")

# Exports; sorted by file:
#      PArray exports:
export PList, PSet, PDict, PArray, PVector, PMatrix, LazyDict,
    PEqualSet, PEqualDict, PIdSet, PIdDict,
    AbstractPDict, AbstractPArray, PLinearSet, PLinearDict, assoc,
    push, pop, pushfirst, popfirst, delete, #, insert, setindex
    Mutability, Immutable, Mutable, mutability, isimmtype,
    isequiv, equivhash, equalfn, hashfn, EquivSet, EquivDict,
    Delay, Var, @var, Volatile, Actor, Source, tx, @tx, ReentrantRef,
    TransactionalRef, AbstractSourceKernel, getfilter,
    getfinalize, setfilter!, setfinalize!, send, geterror,
    receive, reset, PWDict, PWSet, getweight, setweight,
    iscoll, iscolltype, issingle, issingletype, @memoize, @delay

end # module
