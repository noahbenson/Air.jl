################################################################################
# runtests.jl
#
# Runs the Air benchmark suite and compares the results against a stored
# baseline.
#
#   julia --project=bench bench/runtests.jl
#   AIR_UPDATE_BASELINE=1 julia --project=bench bench/runtests.jl
#
# Benchmark timings are noisy on shared or virtualized machines, so the
# comparison uses generous tolerances and is intended to catch order-of-magnitude
# regressions rather than small fluctuations. The baseline is only written when
# AIR_UPDATE_BASELINE=1, so ordinary runs never modify the repository.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__)))

using BenchmarkTools

include(joinpath(@__DIR__, "benchmarks.jl"))

const BASELINE = joinpath(@__DIR__, "baseline.json")

# A baseline is machine-specific, so none is committed: create one locally with
# AIR_UPDATE_BASELINE=1 and comparisons will be made against it. The file is
# gitignored so that it is not accidentally shared between machines.
results = run(SUITE; verbose=true, seconds=1.0, samples=200)

if get(ENV, "AIR_UPDATE_BASELINE", "") == "1"
    BenchmarkTools.save(BASELINE, median(results))
    @info "Wrote benchmark baseline" BASELINE
elseif isfile(BASELINE)
    baseline = BenchmarkTools.load(BASELINE)[1]
    judge(median(results), baseline)
else
    @info "No baseline found; run with AIR_UPDATE_BASELINE=1 to create one" BASELINE
end

show(stdout, MIME"text/plain"(), results)
println()
