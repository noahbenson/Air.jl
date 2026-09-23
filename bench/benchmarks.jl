################################################################################
# benchmarks.jl
#
# The BenchmarkTools suite for Air. Run it with bench/runtests.jl.
#
# Each group covers the persistent collections and the transaction system. The
# operations chosen are the ones on the hot path for typical use: constructing a
# collection, updating it (which is where the persistent data structures do
# their path copying), and reading from it.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

using BenchmarkTools, Air, Random

# Fixed sizes so that runs are comparable across machines, and a fixed seed so
# that the data is reproducible.
const NDICT = 1000
const NVEC = 1000
const NHEAP = 1000
const NTX = 100

Random.seed!(0x5eed)

const SUITE = BenchmarkGroup()

# ==============================================================================
# Persistent dictionaries

function _pdict_pairs(n::Int)
    return [Symbol("k", i) => i for i in 1:n]
end
function _pdict(n::Int)
    d = PDict{Symbol,Int}()
    for p in _pdict_pairs(n)
        d = push(d, p)
    end
    return d
end

SUITE["pdict"] = BenchmarkGroup(["construct", "read", "update"])
SUITE["pdict"]["construct"] = @benchmarkable _pdict($NDICT)
SUITE["pdict"]["getindex"] = @benchmarkable $( _pdict(NDICT) )[:k1]
SUITE["pdict"]["push"] = @benchmarkable push($( _pdict(NDICT) ), :new => 0)
SUITE["pdict"]["delete"] = @benchmarkable delete($( _pdict(NDICT) ), :k1)
SUITE["pdict"]["iterate"] = @benchmarkable sum(values($( _pdict(NDICT) )))

# ==============================================================================
# Persistent sets

function _pset(n::Int)
    s = PSet{Symbol}()
    for i in 1:n
        s = push(s, Symbol("k", i))
    end
    return s
end

SUITE["pset"] = BenchmarkGroup(["construct", "read", "update"])
SUITE["pset"]["construct"] = @benchmarkable _pset($NDICT)
SUITE["pset"]["in"] = @benchmarkable (Symbol("k1") in $( _pset(NDICT) ))
SUITE["pset"]["push"] = @benchmarkable push($( _pset(NDICT) ), :new)
SUITE["pset"]["iterate"] = @benchmarkable foreach(identity, $( _pset(NDICT) ))

# ==============================================================================
# Persistent vectors/arrays

function _pvector(n::Int)
    v = PVector{Int}()
    for i in 1:n
        v = push(v, i)
    end
    return v
end

SUITE["pvector"] = BenchmarkGroup(["construct", "read", "update"])
SUITE["pvector"]["construct"] = @benchmarkable _pvector($NVEC)
SUITE["pvector"]["getindex"] = @benchmarkable $( _pvector(NVEC) )[NVEC ÷ 2]
SUITE["pvector"]["push"] = @benchmarkable push($( _pvector(NVEC) ), 0)
SUITE["pvector"]["setindex"] = @benchmarkable setindex($( _pvector(NVEC) ), 0, NVEC ÷ 2)
SUITE["pvector"]["iterate"] = @benchmarkable sum($( _pvector(NVEC) ))

# ==============================================================================
# Weighted collections

function _pwset(n::Int)
    s = PWSet{Symbol}()
    for i in 1:n
        s = push(s, Symbol("k", i) => float(i))
    end
    return s
end

SUITE["pwset"] = BenchmarkGroup(["construct", "read", "update"])
SUITE["pwset"]["construct"] = @benchmarkable _pwset($NHEAP)
SUITE["pwset"]["first"] = @benchmarkable first($( _pwset(NHEAP) ))
SUITE["pwset"]["push"] = @benchmarkable push($( _pwset(NHEAP) ), :new => 1.0)
SUITE["pwset"]["pop"] = @benchmarkable pop($( _pwset(NHEAP) ))
SUITE["pwset"]["rand"] = @benchmarkable rand($( _pwset(NHEAP) ))
SUITE["pwset"]["setweight"] = @benchmarkable setweight(
    $( _pwset(NHEAP) ), :k1, 1.0
)

# ==============================================================================
# Transactions

SUITE["tx"] = BenchmarkGroup(["commit", "read"])
SUITE["tx"]["commit"] = @benchmarkable begin
    local v = Volatile{Int}(0)
    tx() do
        for _ in 1:$NTX
            v[] = v[] + 1
        end
    end
    v[]
end
SUITE["tx"]["read"] = @benchmarkable begin
    local v = Volatile{Int}(0)
    tx() do
        v[]
    end
end
# Reading a volatile many times within one transaction isolates the per-read
# cost from the fixed cost of starting and committing a transaction.
SUITE["tx"]["readloop"] = @benchmarkable begin
    local v = Volatile{Int}(0)
    tx() do
        local acc = 0
        for _ in 1:$NTX
            acc += v[]
        end
        acc
    end
end

# ==============================================================================
# Task-local variables

SUITE["var"] = BenchmarkGroup(["read"])
SUITE["var"]["getindex"] = @benchmarkable begin
    local v = Var{Int}(0)
    withvars(() -> v[], v => 1)
end

# ==============================================================================
# Delays

SUITE["delay"] = BenchmarkGroup(["read"])
SUITE["delay"]["first"] = @benchmarkable begin
    local d = Delay{Int}(() -> 1)
    d[]
end
SUITE["delay"]["repeat"] = @benchmarkable $( Delay{Int}(1) )[]
