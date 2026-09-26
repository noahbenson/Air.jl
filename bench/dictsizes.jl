################################################################################
# dictsizes.jl
#
# Per-operation cost of lookup, insert and delete for Dict, PDict and TDict, as a
# function of dictionary size.
#
# Fairness notes, since these are easy to get wrong:
#
#  * A `PArray`-style persistent insert allocates one new version; a `TDict`
#    insert is only meaningful in a batch, so the TDict numbers are the cost of
#    the whole batch — `transient`, every update, and `persistent!` — divided by
#    the number of updates. That is the comparison a transient has to win.
#  * Every size starts from a dictionary of that size, so "insert" adds new keys
#    (n+1, n+2, ...) rather than overwriting, and "delete" removes existing ones.
#  * A size-10 dictionary cannot lose 100 keys, so the delete batch is
#    min(OPS, n) for every implementation.
#  * Mutating benchmarks use `setup` with a copy, so no sample sees the result of
#    the previous one.
#  * Every batch is checked once against its expected result before being timed,
#    so a benchmark cannot quietly measure a no-op.
#
# Run:  julia --project=bench bench/dictsizes.jl
#
# This is a tool rather than a test: it prints a table for a person to read, and
# it is not part of any gate. It reports both the mean and the best case, which
# matters — the insert and delete batches allocate 10-16 objects per operation,
# so one GC pause can inflate the mean by 2-3x, while the allocation-free rows
# agree to two decimal places. The best case is what repeats across runs.
#
# It is what found the dictionary read-allocation bug: `PDict` and `PWDict`
# lookups allocated a box per read from about a thousand entries up, because they
# had no `getindex` or `haskey` of their own and fell through to Base's, whose
# sentinel comparison stops splitting a union at that size. See `src/pdict.jl`.
################################################################################
using BenchmarkTools, Air, Random

Random.seed!(0x5eed)

const SIZES = (10, 100, 1_000, 10_000)
const OPS = 100           # updates per insert batch
const SAMPLES = 30

# ---- the collections ---------------------------------------------------------

function _pdict(n::Int)
    d = PDict{Int,Int}()
    for i in 1:n
        d = push(d, i => i)
    end
    return d
end

# ---- the measured operations -------------------------------------------------

function _lookup(c, n::Int)
    s = 0
    for i in 1:n
        s += c[i]
    end
    return s
end

function _lookup_transient(c, n::Int)
    t = transient(c)
    s = 0
    for i in 1:n
        s += t[i]
    end
    return s
end

function _insert_dict!(d, n::Int, ops::Int)
    for i in 1:ops
        d[n + i] = i
    end
    return d
end

function _insert_pdict(d, n::Int, ops::Int)
    for i in 1:ops
        d = push(d, (n + i) => i)
    end
    return d
end

function _insert_tdict(d, n::Int, ops::Int)
    t = transient(d)
    for i in 1:ops
        t[n + i] = i
    end
    return persistent!(t)
end

function _delete_dict!(d, n::Int, ops::Int)
    for i in 1:ops
        delete!(d, i)
    end
    return d
end

function _delete_pdict(d, n::Int, ops::Int)
    for i in 1:ops
        d = delete(d, i)
    end
    return d
end

function _delete_tdict(d, n::Int, ops::Int)
    t = transient(d)
    for i in 1:ops
        delete!(t, i)
    end
    return persistent!(t)
end

# ---- measuring ---------------------------------------------------------------

# The mean and the best-case time per operation, in nanoseconds, with the
# allocation *count* per operation. (`@allocated` reports bytes; `allocs` is a
# count.) `evals` is left to BenchmarkTools: forcing `evals=1` on an operation
# faster than the clock's own resolution makes the best case read as zero.
function _perop(trial, divisor::Int)
    return (mean(trial).time / divisor, minimum(trial).time / divisor,
            minimum(trial).allocs / divisor)
end

function _bench_read(f, c, n::Int)
    trial = @benchmark $f($c, $n) samples = SAMPLES
    return _perop(trial, n)
end

# The functional implementations do not touch their input, so they can be timed
# directly; the mutating one gets a fresh copy outside the measured region.
function _bench_fn(f, c, n::Int, ops::Int)
    trial = @benchmark $f($c, $n, $ops) samples = SAMPLES
    return _perop(trial, ops)
end

function _bench_mut!(f, c, n::Int, ops::Int)
    trial = @benchmark $f(c0, $n, $ops) setup = (c0 = copy($c)) samples = SAMPLES
    return _perop(trial, ops)
end

# ---- the analysis ------------------------------------------------------------

rows = Vector{Any}[]
println("size,operation,impl,mean_ns,best_ns,allocs")
for n in SIZES
    d = Dict{Int,Int}(i => i for i in 1:n)
    p = _pdict(n)
    ndel = min(OPS, n)

    # Sanity: each measured batch must do what its name says.
    @assert _lookup(d, n) == n * (n + 1) ÷ 2
    @assert _lookup(p, n) == n * (n + 1) ÷ 2
    @assert _lookup_transient(p, n) == n * (n + 1) ÷ 2
    @assert length(_insert_dict!(copy(d), n, OPS)) == n + OPS
    @assert length(_insert_pdict(p, n, OPS)) == n + OPS
    @assert length(_insert_tdict(p, n, OPS)) == n + OPS
    @assert length(_delete_dict!(copy(d), n, ndel)) == n - ndel
    @assert length(_delete_pdict(p, n, ndel)) == n - ndel
    @assert length(_delete_tdict(p, n, ndel)) == n - ndel

    for (op, f, g, h) in (
        ("lookup", _lookup, _lookup, _lookup_transient),
        ("insert", _insert_dict!, _insert_pdict, _insert_tdict),
        ("delete", _delete_dict!, _delete_pdict, _delete_tdict),
    )
        dop = op == "insert" ? OPS : ndel
        for (impl, fn) in (("Dict", f), ("PDict", g), ("TDict", h))
            coll = impl == "Dict" ? d : p
            m, b, a =
                op == "lookup" ? _bench_read(fn, coll, n) :
                (impl == "Dict" ? _bench_mut!(fn, coll, n, dop) :
                 _bench_fn(fn, coll, n, dop))
            push!(rows, [n, op, impl, m, b, a])
            println("$n,$op,$impl,$(round(m, digits=2)),$(round(b, digits=2)),$(round(a, digits=2))")
        end
    end
end

# ---- a compact table ---------------------------------------------------------

println()
print(rpad("size", 8))
for op in ("lookup", "insert", "delete"), impl in ("Dict", "PDict", "TDict")
    print(rpad("$op/$impl", 18))
end
println()
for n in SIZES
    print(rpad(n, 8))
    for op in ("lookup", "insert", "delete"), impl in ("Dict", "PDict", "TDict")
        r = rows[findfirst(x -> x[1] == n && x[2] == op && x[3] == impl, rows)]
        print(rpad(string(round(r[4], digits=1)) * " (" * string(round(r[6], digits=1)) * " a)", 18))
    end
    println()
end
