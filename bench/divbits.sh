#!/usr/bin/env bash
################################################################################
# divbits.sh
#
# Sweeps the trie's branching geometry — PTREE_NODE_SHIFT and PTREE_TWIG_SHIFT,
# i.e. how many bits of the key each internal node and each twig consume — and
# records each combination's correctness and timings.
#
# The constants are compile-time: they decide the width of the occupancy bitmap,
# the depth encoding in a node id, and the tree's depth. So each combination has
# to be compiled in its own process, which is what bench/divbits_payload.jl runs
# as. This script rewrites the two constants, runs the payload, and restores the
# source afterwards (including on failure).
#
# Usage:  bench/divbits.sh > divbits.txt
#         (JULIA=/path/to/julia to pick an interpreter)
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2021 Noah C. Benson

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
src="$root/src/ptree.jl"
julia=${JULIA:-julia}

backup=$(mktemp)
cp "$src" "$backup"
trap 'cp "$backup" "$src"; rm -f "$backup"' EXIT

set_constants() {
    python3 - "$src" "$1" "$2" <<'PY'
import re, sys
path, node_shift, twig_shift = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
s, n1 = re.subn(r"const PTREE_NODE_SHIFT\s*=\s*\d+",
                f"const PTREE_NODE_SHIFT    = {node_shift}", s)
s, n2 = re.subn(r"const PTREE_TWIG_SHIFT\s*=\s*\d+",
                f"const PTREE_TWIG_SHIFT    = {twig_shift}", s)
if n1 != 1 or n2 != 1:
    sys.exit(f"divbits.sh: expected one of each constant, found {n1} and {n2}")
open(path, "w").write(s)
PY
}

for node_shift in 4 5 6; do
    for twig_shift in 4 5 6; do
        set_constants "$node_shift" "$twig_shift"
        if ! "$julia" --project="$root/bench" --startup-file=no \
                "$root/bench/divbits_payload.jl" 2>/dev/null | grep '^SHIFTS'; then
            echo "SHIFTS node=$node_shift twig=$twig_shift failed=true"
        fi
    done
done
