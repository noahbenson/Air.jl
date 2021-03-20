#! /usr/bin/env julia

import Pkg

(pwd() == @__DIR__) || cd(@__DIR__)

Pkg.activate(pwd())
Pkg.instantiate()

using Air


jobs = Symbol[]
if length(ARGS) == 0
    push!(jobs, :test)
    push!(jobs, :docs)
else
    for a in ARGS
        a = lowercase(a)
        push!(jobs, Symbol(a))
    end
end
for j in jobs
    if j == :test
        Pkg.test("Air", coverage=true)
        Pkg.add("Documenter")
        import Documenter
        Documenter.doctest(Air)
    elseif j == :codecov
        Pkg.add("Coverage")
        using Coverage
        Codecov.submit(Codecov.process_folder())
    elseif j == :docs
        Pkg.add("Documenter")
        include("docs/make.jl")
    else
        error("Unknown make command: $j")
    end
end

exit(0)


