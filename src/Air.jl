################################################################################
# Air.jl
#
# The Air library.
# Functional collections and tools for Julia.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

module Air

include("util.jl")
include("api.jl")
include("ptree.jl")
include("pset.jl")
include("pdict.jl")
include("lazydict.jl")
include("parray.jl")
include("pheap.jl")
include("pwdict.jl")
include("pwset.jl")

include("TX.jl")

export Var, @var, Volatile, Actor, Source, tx, @tx, ReentrantRef,
    TransactionalRef, AbstractSourceKernel, getfilter, getfinalize,
    setfilter!, setfinalize!, send, geterror, receive, reset

end 
