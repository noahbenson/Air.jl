################################################################################
# PHeap.jl
# The Persistent heap (priority queue) type and related types such as THeap (the
# transient heap type).
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020 Noah C. Benson

import Random

# ==============================================================================
# #TODO List

# THeap
# random(pheap) to select weighted


# ==============================================================================
# PHeap definition.
"""
    PHeap{T,W,F}

A PHeap object is a priority queue to which objects can be pushed and given a
particular weight. PHeaps can then be iterated, popped, and the first element
may be obtained in a priority-weighted order.

PHeap objects may be initialized with a single argument of type <:Function, in
which case this function is used to compare the weights for ordering. By default
this is > such that high weights are retrieved first by the first() function and
removed first by the pop() function.
"""
struct PHeap{T,W<:Number,F<:Function,D<:AbstractPDict{T,Int}}
    _heap::PVector{Tuple{T,W,W}} # value, weight, total subtree weight
    _index::D
    _compare::F
end
PHeap{T,W,F,D}(f::F) where {T, W<:Number, F<:Function, D<:AbstractPDict} = PHeap{T,W,F,D}(
    PVector{Tuple{T,W,W}}(),
    D{T,Int}(),
    f)
PHeap{T,W,F,D}(f::F) where {T, W<:Number, F<:Function, D<:AbstractPDict{T,Int}} = PHeap{T,W,F,D}(
    PVector{Tuple{T,W,W}}(),
    D(),
    f)
PHeap{T,W,D}(f::F) where {T, W<:Number, F<:Function, D<:AbstractPDict} = PHeap{T,W,F,D}(
    PVector{Tuple{T,W,W}}(),
    D{T,Int}(),
    f)
PHeap{T,W,D}(f::F) where {T, W<:Number, F<:Function, D<:AbstractPDict{T,Int}} = PHeap{T,W,F,D}(
    PVector{Tuple{T,W,W}}(),
    D(),
    f)
PHeap{T,W,D}()   where {T, W<:Number, D<:AbstractPDict} = PHeap{T,W,D}(>)
PHeap{T,W}(f::F) where {T, W<:Number, F<:Function} = PHeap{T,W,PDict}(f)
PHeap{T,W}()     where {T, W<:Number} = PHeap{T,W,PDict}(>)
PHeap{T,D}(f::F) where {T, F<:Function, D<:AbstractPDict} = PHeap{T,Float64,D}(f)
PHeap{T,D}()     where {T, D<:AbstractPDict} = PHeap{T,Float64,D}(>)
PHeap{T}(f::F)   where {T, F<:Function} = PHeap{T,Float64,F,PDict}(f)
PHeap{T}()       where {T} = PHeap{T,Float64,PDict}(>)
PHeap(f::F)      where {F <: Function} = PHeap{Any,Float64,F,PDict}(f)
PHeap() = PHeap{Any,Float64,PDict}(>)
_pheap_swap(heap::PVector{Tuple{T,W,W}}, index::D, ii::Int, jj::Int, subii::Bool, subjj::Bool) where {T,W,D<:AbstractPDict{T,Int}} = begin
    (tii,wii,totii) = heap[ii]
    (tjj,wjj,totjj) = heap[jj]
    newtotii = totjj + wii - (subjj ? wjj : 0)
    newtotjj = totii + wjj - (subii ? wii : 0)
    heap = setindex(heap, (tjj, wjj, newtotjj), ii)
    heap = setindex(heap, (tii, wii, newtotii), jj)
    index = setindex(index, ii, tjj)
    index = setindex(index, jj, tii)
    return (heap, index)
end
_pheap_fix_tot(heap::PVector{Tuple{T,W,W}}, ii::Int, dw::W) where {T,W} = begin
    dw == 0 && return heap
    while ii > 0
        node = heap[ii]
        heap = setindex(heap, (node[1], node[2], node[3] + dw), ii)
        ii = div(ii, 2)
    end
    return heap
end
_pheap_fix_up(heap::PVector{Tuple{T,W,W}}, index::D, cmp::F, ii::Int, oldw::W) where {T,W,F<:Function,D<:AbstractPDict{T,Int}} = begin
    # This change affects the weight totals of nodes above us by this abount:
    newnode = heap[ii]
    neww = newnode[2]
    dw = neww - oldw
    zW = zero(W)
    # if ii is 1, then we don't need to worry about the parent nodes and
    # swapping on the way up in order to fix things.
    while ii > 1
        pp = div(ii, 2)
        parent = heap[pp]
        # Break when we don't need to percolate up anymore.
        cmp(neww, parent[2]) || break
        # We need to swap with the parent
        (heap, index) = _pheap_swap(heap, index, ii, pp, true, false)
        ii = pp
    end
    # At this point we may have gone all the way up the heap, but we might have
    # only traversed up some of it. In the latter case, if the change in the
    # weight is not 0, then we need to fix weight totals up the heap.
    heap = _pheap_fix_tot(heap, div(ii, 2), dw)
    return (heap, index)
end
_pheap_fix_down(heap::PVector{Tuple{T,W,W}}, index::D, cmp::F, ii::Int) where {T,W,F<:Function,D<:AbstractPDict{T,Int}} = begin
    # We can assume when this function is called that the totals are correct for
    # the entire heap at this point.
    n = length(heap)
    zW = zero(W)
    while true
        lch = ii*2
        # If we've reached the end, we can return the (updated?) heap and index.
        (lch > n) && break
        # Otherwise, we need to check the children
        rch = lch + 1
        node = heap[ii]
        lnode = heap[lch]
        if rch > n
            # Only check the left...
            if cmp(lnode[2], node[2])
                # Swap these two
                (heap, index) = _pheap_swap(heap, index, ii, lch, true, true)
            end
            # There can't be children beyond this point.
            break
        else
            rnode = heap[rch]
            if cmp(lnode[2], rnode[2])
                cmp(node[2], lnode[2]) && break
                # lnode < node and lnode < rnode
                (heap, index) = _pheap_swap(heap, index, ii, lch, true, true)
                ii = lch
            else
                cmp(node[2], rnode[2]) && beak
                # rnode < node and rnode < lnode
                (heap, index) = _pheap_swap(heap, index, ii, rch, true, true)
                ii = rch
            end
        end
    end
    return (heap, index)
end
_pheap_fixw(heap::PVector{Tuple{T,W,W}}, index::D, cmp::F, ii::Int, neww::W) where {T,W,F<:Function,D<:AbstractPDict{T,Int}} = begin
    node = heap[ii]
    oldw = node[2]
    dw = neww - oldw
    (dw == 0) && return (heap, index)
    # Fix this in the current node.
    heap = setindex(heap, (node[1], neww, node[3] + dw), ii)
    # If dw is positive then there is no chance that we will need to
    # swap with the children; on the flip side, if it is negative,
    # we still need to fix totals up the heap. So we can do that
    # now.
    (heap, index) = _pheap_fix_up(heap, index, cmp, ii, oldw)
    # Now, we only need to check the lower heap nodes if the weight
    # is now lower. If it is now lower, then ii will still point
    # to the same node after the _pheap_fix_up() call.
    if cmp(oldw, neww)
        (heap, index) = _pheap_fix_down(heap, index, cmp, ii)
    end
    # That is all.
    return (heap, index)
end
push(p::PHeap{T,W,F,D}, tw::Tuple{S,X}) where {T,W,F,D,S,X<:Number} = begin
    (t,w) = tw
    (w <= 0) && throw(ArgumentError("Weights must be positive"))
    cmp = p._compare
    heap = p._heap
    index = p._index
    # If t is already in the heap, we are just editing the weight.
    ii = get(index, t, 0)
    if ii == 0
        # Since t isn't already in the heap, we need to add it to the end
        # with a weight of 0 (which keeps the tree valid and all of the
        # weight totals in tact) then fix its weight.
        zW = zero(W)
        heap = push(heap, (t,zW,zW))
        ii = length(heap)
        index = setindex(index, ii, t)
    end
    (heap, index) = _pheap_fixw(heap, index, cmp, ii, w)
    (p._heap === heap && p._index === index) && return p
    return PHeap{T,W,F,D}(heap, index, cmp)
end
push(p::PHeap{T,W,F,D}, tw::Pair{S,X}) where {T,W,F,D,S,X<:Number} = begin
    return push(p, (tw[1],tw[2]))
end
"""
    setweight(h, x, w)

Sets the weight of the value x in the given persistent heap or persistent
weighted collection object h to be w and yields the new updated version of h.
If x is not already in the object h, then an error is thrown.
"""
setweight(p::PHeap{T,W,F,D}, t::S, w::X) where {T,W,F,D,S,X<:Number} = begin
    (w <= 0) && throw(ArgumentError("Weights must be positive"))
    ii = get(p._index, t, 0)
    (ii == 0) && error("cannot setweight on item that is not in the PHeap")
    (heap, index) = _pheap_fixw(p._heap, p._index, p._compare, ii, w)
    (heap === p._heap && index === p._index) && return p
    return PHeap{T,W,F,D}(heap, index, p._compare)
end
"""
    getweight(h, x)

Yields the weight for the given value x in the given persistent heap or
persistent weighted collection h. If h is a weighted dictionary, then x
should be the key.

If the key or object x is not found in the collection h, then 0 is returned.
"""
getweight(h::PHeap{T,W,F,D}, t::S) where {T,W,F,D,S} = begin
    ii = get(h._index, t, 0)
    (ii == 0) && return 0
    return h._heap[ii][2]
end
_pheap_delete(p::PHeap{T,W,F,D}, ii::Int) where {T,W,F,D} = begin
    n = length(p)
    cmp = p._compare
    (n == 1) && return PHeap{T,W,F,D}(pop(p._heap), empty(p._index), cmp)
    # Grab the very last node, then we can prep the removal by updating its
    # weight to be zero (this subtracts its weight from the tree's totals).
    (tn,wn,totn) = p._heap[n]
    (heap, index) = _pheap_fixw(p._heap, p._index, cmp, n, zero(W))
    # Now we replace the node in slot ii with the last node, but keep the
    # old weight so that the tree's totals are still valid.
    (tii, wii, totii) = heap[ii]
    heap = setindex(pop(heap), (tn, wii, totii), ii)
    index = setindex(delete(index, tii), ii, tn)
    # Then set the weight at that node to be that of its old value!
    (heap, index) = _pheap_fixw(heap, index, cmp, ii, wn)
    # Thwt will leave us in a valid state; return the new heap.
    (heap === p._heap && index === p._index) && return p
    return PHeap{T,W,F,D}(heap, index, cmp)
end
pop(p::PHeap{T,W,F,D}) where {T,W,F,D} = begin
    n = length(p)
    (n == 0) && throw(ArgumentError("PHeap must be non-empty"))
    return _pheap_delete(p, 1)
end
delete(p::PHeap{T,W,F,D}, t::S) where {T,W,F,D,S} = begin
    # Find this particular node we need to remove
    ii = get(p._index, t, 0)
    return ii == 0 ? p : _pheap_delete(p, ii)
end
Base.length(p::PHeap{T,W,F,D}) where {T,W,F,D} = length(p._index)
Base.iterate(p::PHeap{T,W,F,D}) where {T,W,F,D} = iterate(p, p)
Base.iterate(::PHeap{T,W,F,D}, p::PHeap{T,W,F,D}) where {T,W,F,D} = begin
    return (length(p) == 0 ? nothing : (first(p), pop(p)))
end
Base.first(p::PHeap{T,W,F,D}) where {T,W,F,D} = begin
    (length(p) == 0) && throw(ArgumentError("PHeap must be non-empty"))
    return p._heap[1][1]
end
Base.in(x::S, p::PHeap{T,W,F,D}) where {S,T,W,F,D} = begin
    return (get(p._index, x, 0) != 0)
end
Random.rand(p::PHeap{T,W,F,D}) where {T,W,F,D} = begin
    n = length(p)
    (n == 0) && throw(ArgumentError("PHeap must be non-empty"))
    x = Random.rand(Float64) * p._heap[1][3]
    ii = 1
    (t,w,tot) = p._heap[ii]
    while true
        (x <= w) && return t
        x -= w
        ii *= 2
        (t,w,tot) = p._heap[ii]
        (x <= tot) && continue
        x -= tot
        ii += 1
        (t,w,tot) = p._heap[ii]
    end
    error("invalid state reached")
end
