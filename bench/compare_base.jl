################################################################################
# compare_base.jl
#
# Compares each persistent collection against the mutable `Base` collection it
# stands in for. The README makes claims about this ("nearly as fast as native
# Dict objects"), and the only way to keep such a claim honest is to measure it.
#
#   julia --project=bench bench/compare_base.jl
#
# Persistent data structures are expected to be slower than mutable ones: every
# update copies the path from the root to the changed node, which is what buys
# the persistence. The useful question is how much slower, and whether the gap
# grows or shrinks as the library changes. Run this before and after a change to
# the collection internals.
#
# Two caveats when reading the ratios: Base's loops over an `Array` are
# SIMD-vectorised, so the vector rows compare a scalar loop against a vectorised
# one rather than comparing the index operations themselves; and Base's mutable
# collections update in place, so the construct/push rows compare an O(1)
# mutation against an O(log n) path copy, which is the price of persistence and
# not a defect.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using Air, BenchmarkTools, Random

Random.seed!(0x5eed)

const N = 1000
const KEYS = [Symbol("k", i) for i in 1:N]
const MID = KEYS[N ÷ 2]

function pdict_built()
    d = PDict{Symbol,Int}()
    for (i, k) in enumerate(KEYS)
        d = push(d, k => i)
    end
    return d
end
function dict_built()
    d = Dict{Symbol,Int}()
    for (i, k) in enumerate(KEYS)
        d[k] = i
    end
    return d
end

const PD = pdict_built()
const DD = dict_built()
const PS = PSet(KEYS)
const BS = Set(KEYS)
const PV = PVector(collect(1:N))
const BV = collect(1:N)

# Accumulating rather than discarding the reads keeps the compiler from
# eliminating the loop entirely (which would make the comparison meaningless).
sumidx(v) = (s = 0; for i in 1:N; s += v[i]; end; s)

function report(name, air, base, n)
    a, b = minimum(air), minimum(base)
    println(
        rpad(name, 22),
        "Air ", lpad(round(a.time / n, digits=1), 7), " ns, ",
        lpad(a.allocs, 6), " allocs   |   Base ",
        lpad(round(b.time / n, digits=1), 7), " ns, ",
        lpad(b.allocs, 6), " allocs   |  ",
        round(a.time / max(b.time, 1e-9), digits=1), "x",
    )
end

println(rpad("operation", 22), "Air (per op)                |   Base (per op)                | ratio")
println("-"^95)

report("dict lookup",
    @benchmark(for k in $KEYS; $PD[k]; end),
    @benchmark(for k in $KEYS; $DD[k]; end), N)
report("dict construct",
    @benchmark(pdict_built()),
    @benchmark(dict_built()), N)
report("set membership",
    @benchmark(for k in $KEYS; k in $PS; end),
    @benchmark(for k in $KEYS; k in $BS; end), N)
report("vector index",
    @benchmark(sumidx($PV)),
    @benchmark(sumidx($BV)), N)
report("vector push",
    @benchmark(push($PV, 0)),
    @benchmark(push!($BV, 0)), 1)
report("dict iterate",
    @benchmark(for x in $PD; end),
    @benchmark(for x in $DD; end), N)
