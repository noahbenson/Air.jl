#! /usr/bin/env julia

import Pkg

(pwd() == @__DIR__) || cd(@__DIR__)

push!(LOAD_PATH, pwd())

using Air

#Pkg.activate(@__DIR__)

air_package = Pkg.PackageSpec(path=pwd())
Pkg.develop(air_package)

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
        Pkg.test("Air")
    elseif j == :codecov
        Pkg.add("Coverage")
        using Coverage
        Codecov.submit(Codecov.process_folder())
        Pkg.rm("Coverage")
    elseif j == :docs
        Pkg.add("Documenter")
        Pkg.instantiate()
        include("docs/make.jl")
    else
        error("Unknown make command: $j")
    end
end

exit(0)

