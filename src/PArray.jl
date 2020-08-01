################################################################################
# PArrray.jl
# The Persistent Array type and related types such as TArray (the transient
# array type).
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2019 Noah C. Benson

#struct PArray{T,N} <: AbstractArray{T,N}
#    _i0::Int
#    _size::NTuple{N,Int}
#    _tree::PTree{T}
#end



# ==============================================================================
# Private
#
# This section contains private constants and function definitions intendef for
# use only in this library or occasionally within the library, but not outside
# of the current module.
#
# # Functions
#
#   Function to dice up arrays.
"""
    _array_dice(a)
Yields an array with the same number of dimensions as a but in which each
element in a has been assigned to a super-array storing a chunck of the 
original array.

The optional arguments cpy (default: true) specifies whether the array
should be copied or used verbatim (false), and the optional output type
specifies that the result should consist of arrays of type S instead of
the eltype of a.
"""
_array_dice(a::AbstractArray{T,N}) where {T,N} = begin
    sz = size(a)
    osz = Tuple([Int(ceil(n / _PTREE_COUNT)) for n in sz])
    out = Array{Array{T,N}}(undef, osz...)
    idx = CartesianIndices(osz)
    for ii in idx
        jj = ((1 + _PTREE_COUNT*(k - 1)):min(mx,k*_PTREE_COUNT)
              for (k,mx) in zip(Tuple(ii), sz))
        out[ii] = getindex(a, jj...)
    end
    return out
end
_array_dice(a::AbstractArray{T,N}, S::Type) where {T,N} = begin
    sz = size(a)
    osz = Tuple([Int(ceil(n / _PTREE_COUNT)) for n in sz])
    out = Array{Array{S,N}}(undef, osz...)
    idx = CartesianIndices(osz)
    for ii in idx
        jj = ((1 + _PTREE_COUNT*(k - 1)):min(mx,k*_PTREE_COUNT)
              for (k,mx) in zip(Tuple(ii), sz))
        out[ii] = getindex(a, jj...)
    end
    return out
end


# ==============================================================================
# Structures / Types
#
# Definitions of structures and types that should come before the various
# methods and constructors that operate on them.
#
# #PArray ======================================================================
"""
    PArray{T, N}
The PArray class is a mimic of the Array class with the exception that it is
designed for immutable operations. A PArray object can act like and be accessed
like an Array object, but its values cannot be set. Instead, they support
efficient push and setindex operations (note the lack of the trailing !) that do
not modify the original object but instead yield a new object with the given
updates that is also persistent.

Warning: All of PArray's fields are to be considered private. In particular
until immutable arrays are supported, PArrays are implemented using mutable
arrays that are immutable only so long as the fields of all PArray objects
remain unchanged. Because persistent collections employ heavy data-sharing by
reference, changing one object's data may change the data for many objects.
"""
struct PArray{T, N} <: AbstractArray{T, N}
    # The size tuple stores the size of the array: the length of each dim.
    _size::NTuple{N, Int}
    # The index0 indicates if the first node(s) contain(s) fewer than _PTREE_COUNT
    # elements; this allows for efficient shifting and unshifting.
    # An index0 of x indicates that element 1 maps to element x in the data-tree
    _index0::NTuple{N, Int}
    # The tree that actually stores the values
    _tree::PTree{T}
    # The default value for elements not set in the tree.
    _default::T
end
# PArray is immutable, even though it employs a mutable field (Array)
mutability(::Type{PArray}) = Immutable
mutability(::Type{PArray{T,N}}) where {T,N} = Immutable
#
# PArrray and TArray aliases
const PVector = PArray{T,1} where {T}
const PMatrix = PArray{T,2} where {T}
#
# PArray and TArray are very similar, so I have a macro for when to duplicate
# behavior across functions for both:
macro _ptmethod(asym::Symbol, ex::Expr)
    h = ex.head
    fnsign = nothing
    if h === :function
        fnsign = ex.args[1]
        (fnsign.head === :where) && (fnsign = fnsign.args[1])
    elseif h === :(=) || isa(ex.args[1], Expr)
        ea = ex.args[1]
        if ea.head === :call
            fnsign = ex.args[1]
        elseif ea.head === :where && ea.args[1] isa Expr && ea.args[1].head == :call
            fnsign = ea.args[1].args[1]
        end
    end
    (fnsign === nothing) && throw(ErrorException("_ptmethods must be function definitions"))
    fnname = fnsign.args[1]
    # convert the array symbol to PArray or TArray
    q = quote
        # First, a version for PArray...
        let $asym = PArray
            global $fnname
            $ex
        end
        let $asym = TArray
            global $fnname
            $ex
        end
    end
    return esc(q)
end
#
# Example of the above macro; a function for building up an array from another
# array or from sequences of already-full lower-height arrays...
@_ptmethod IArray function _built_size(u::Array{IArray{T,N},N}) where {T,N}
    return u[1]._size .* (size(u) .- 1) .+ u[end]._size
end
@_ptmethod IArray function _build_iarray(tt::Array{IArray{T,N},N}) where {T,N}
    ph = Tuple(0 for n in tt[1]._size)
    h0 = tt[1]._height
    while true
        h0 += 1
        tsz = size(tt)
        maxsz = max(tsz...)
        (maxsz <= _PTREE_COUNT) && return IArray{T,N}(_built_size(tt), ph, h0, tt)
        div = _array_dice(tt)
        tt = IArray{T,N}[IArray{T,N}(_built_size(u), ph, h0, u) for u in div]
    end
end
@_ptmethod IArray function _build_iarray(tt::Type{IArray},
                                         a::AbstractArray{T,N}) where {T,N}
    sz = size(a)
    ph = Tuple(0 for n in sz)
    maxsz = max(sz...)
    # if the array is small enough, we can store it as one array right here
    (maxsz <= _PTREE_COUNT) && return IArray{T,N}(sz, ph, 0, Array{T,N}(a))
    # start by dividing up this array along each dimension
    div = _array_dice(a)
    # make that into an array of P/TArrays...
    tt = IArray{T,N}[IArray{T,N}(size(u), ph, 0, u) for u in div]
    # and build this up...
    return _build_iarray(tt)
end
@_ptmethod IArray function _build_iarray(tt::Type{IArray},
                                         a::AbstractArray{T,0}) where {T}
    return IArray{T,0}((), (), 0, convert(Array{T,0}, fill(a[], ())))
end
@_ptmethod IArray function _build_iarray(tt::Type{IArray},
                                         a::AbstractArray{T,N},
                                         S::Type) where {T,N}
    sz = size(a)
    ph = Tuple(0 for n in sz)
    maxsz = max(sz...)
    # if the array is small enough, we can store it as one array right here
    (maxsz <= _PTREE_COUNT) && return IArray{S,N}(sz, ph, 0, Array{S,N}(a))
    # start by dividing up this array along each dimension
    div = _array_dice(a)
    # make that into an array of P/TArrays...
    tt = IArray{S,N}[IArray{S,N}(size(u), ph, 0, u) for u in div]
    # and build this up...
    return _build_iarray(sz, ph, tt)
end

################################################################################
# Constructors
#
# PArray and TArray Constructors.
#
# #PArray #Constructors ========================================================
#
# Make empty PArrays
PArray() = PArray{Any,1}(NTuple{0,Int}(), NTuple{0,Int}(), 0, Any[])
PArray{T}() where {T} = PArray{T,1}(NTuple{0,Int}(), NTuple{0,Int}(), 0, T[])
PArray{T,N}() where {T,N} = begin
    sz = NTuple{N,Int}([0, [1 for x in 2:N]...])
    ph = NTuple{N,Int}([0 for x in 1:N])
    return PArray{T,N}(sz, ph, 0, Array{T,N}(undef, sz...))
end
                
#
# Make PArrays from other PArrays (just return them!)
PArray{T,N}(a::PArray{T,N}) where {T,N} = a
PArray(a::PArray{T,N}) where {T,N} = a
#
# Make PArrays from TArrays; note that this doesn't reset the TArray (the way
# that calling persist! does); rather it copies the TArray elements and leaves
# the TArray untouched.
PArray{T,N}(a::TArray{T,N}) where {T,N} = begin
    if a._root isa Array{T,N}
        return PArray{T,N}(a._size, a._phase, a._height, copy(a._root))
    else
        root = PArray{T,N}[PArray{T,N}(u) for u in a._root]
        return PArray{T,N}(a._size, a._phase, a._height, root)
    end
end
PArray(a::TArray{T,N}) where {T,N} = PArray{T,N}(a)
#
# Make PArrays from any kind of Array
PArray{T,N}(a::AbstractArray{T,N}) where {T,N} = _build_iarray(PArray, a)
PArray{S,N}(a::AbstractArray{T,N}) where {T,S,N} = _build_iarray(PArray, a, S)
PArray(a::AbstractArray{T,N}) where {T,N} = PArray{T,N}(a)
#
#
# #TArray #Constructors ========================================================
#
# Make empty TArrays
TArray() = TArray{Any,1}(NTuple{0,Int}(), NTuple{0,Int}(), 0, Any[])
TArray{T}() where {T} = TArray{T,1}(NTuple{0,Int}(), NTuple{0,Int}(), 0, T[])
#
# Make TArrays from other TArrays
TArray{T,N}(a::TArray{T,N}) where {T,N} = begin
    # we recursively duplicate the arrays that aren't parrays
    root = copy(a._root)
    if !isa(a._root, Array{T,N})
        for (ii,el) in enumerate(root)
            (el isa TArray{T,N}) && (root[ii] = TArray{T,N}(el))
        end
    end
    return TArray{T,N}(a._size, a._phase, a._height, root)
end
TArray(a::TArray{T,N}) where {T,N} = TArray{T,N}(a)
#
# Make TArrays from PArrays
TArray{T,N}(a::PArray{T,N}) where {T,N} = begin
    if a._root isa Array{T,N}
        root = copy(a._root)
    else
        root = convert(Array{Union{PArray{T,N}, TArray{T,N}}, N}, a._root)
    end
    return TArray{T,N}(a._size, a._phase, a._height, root)
end
TArray(a::PArray{T,N}) where {T,N} = TArray{T,N}(a)
#
# Make TArrays from any kind of Array
TArray{T,N}(a::AbstractArray{T,N}) where {T,N} = _build_iarray(TArray, a)
TArray{S,N}(a::AbstractArray{T,N}) where {T,S,N} = _build_iarray(TArray, a, S)
TArray(a::AbstractArray{T,N}) where {T,N} = TArray{T,N}(a)
#
#
# Additional Conversion Functions ==============================================
#
# Make a persistent array into a transient array.
transient(u::PArray{T,N}) where {T,N} = TArray{T,N}(u)
# Make a persistent array from a transient array (without editing the transient)
persistent(u::TArray{T,N}) where {T,N} = PArray{T,N}(u)
persistent(u::PArray{T,N}) where {T,N} = u

################################################################################
# Base interface methods
#
# Methods such as Base.length and Base.getindex that are implemented by Parray
# and TArray.
#
# Helper methods
_subindices(ph::NTuple{N,Int}, h::Int, ii::CartesianIndex{N}) where {T,N} = begin
    ii = Tuple(ii)
    jj = ii .+ ph .- 1
    sh = h * _PTREE_SHIFT
    mm = (1 << sh) - 1
    ii = (jj .>> sh) .+ 1
    jj = jj .- (ph .* (ii .== 1))
    jj = (jj .& mm) .+ 1
    return (CartesianIndex(ii), CartesianIndex(jj))
end
_locdelta(ph::NTuple{N,Int}, sz::NTuple{N,Int}, arr::Array{T,N}) where {T,N} = begin
    return ph #size(arr) .- (sz .- ph)
end
_locindices(ph::NTuple{N,Int}, sz::NTuple{N,Int}, arr::Array{T,N}, ii::CartesianIndex{N}) where {T,N} = begin
    return CartesianIndex(Tuple(ii) .+ _locdelta(ph, sz, arr))
end
#
# #getindex ====================================================================
# Flattening -------------------------------------------------------------------
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ::Colon) where {T,N} = begin
    if length(u._size) == 1
        return u
    else
        return reshape(u, :)
    end
end
# Linear indexing --------------------------------------------------------------
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ii::Int) where {T,N} = begin
    return u[CartesianIndices(u._size)[ii]]
end
@_ptmethod IArray Base.getindex(u::IArray{T,N}, iis::AbstractArray{Int,1}) where {T,N} = begin
    return u[CartesianIndices(u._size)[iis]]
end
# Logical indexing -------------------------------------------------------------
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ib::AbstractArray{Bool,1}) where {T,N} = begin
    return u[CartesianIndices(u._size)[ib]]
end
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ib::AbstractArray{Bool,N}) where {T,N} = begin
    return u[CartesianIndices(u._size)[ib]]
end
# Cartesian indexing -----------------------------------------------------------
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ii::CartesianIndex{N}) where {T,N} = begin
    if u._root isa Array{T,N}
        # we take from the back...
        jj = _locindices(u._phase, u._size, u._root, ii)
        return u._root[jj]
    else
        (ii, jj) = _subindices(u._phase, u._height, ii)
        return getindex(u._root[ii], jj)
    end
end
@_ptmethod IArray Base.getindex(u::IArray{T,N}, ii::Vararg{Int,N}) where {T,N} = u[CartesianIndex(ii)]
@_ptmethod IArray Base.getindex(u::IArray{T,N}, iis::AbstractArray{CartesianIndex{N},M}) where {T,N,M} = begin
    return IArray{T,M}(T[u[ii] for ii in iis])
end
# Multidimensional indexing ----------------------------------------------------
# Simple approximation of efficient multidimensional indexing is to use a view;
# since PArrays can't be overwritten anyway, and we want slices to share memory
# when possible, we can just use view for PArray.
# For TArray, the custom multidimensional index, in which an ordinary array is
# returned, would be fine, but for keeping things a TArray; so we do that.
Base.getindex(u::PArray{T,N}, I...) where {T,N} = view(u, I...)
Base.getindex(u::TArray{T,N}, I...) where {T,N} = TArray(view(u, I...))
#
# #size ========================================================================
Base.size(u::PArray) = u._size
Base.size(u::TArray) = u._size
#
# #length ======================================================================
Base.length(u::PArray) = prod(u._size)
Base.length(u::TArray) = prod(u._size)
#
# #setindex ====================================================================
"""
    setindex(parray, val, indices...) 
Yields a copy of the given persistent array with the value or values at the 
given indices changed to val. The setindex(parray, ...) function is the 
persistent/immutable version of the setindex! function.
"""
setindex(u::PArray{T,N}, el::S, I...) where {T,N,S} = begin
    t = transient(u)
    # we use transient arrays for complex update operations
    setindex!(t, el, I...)
    return persistent(t)
end
setindex(u::PArray{T,N}, el::S, ii::Int) where {T,N,S} = begin
    idx = CartesianIndices(u._size)
    return setindex(u, idx[ii])
end
setindex(u::PArray{T,N}, el::S, ii::Vararg{Int,N}) where {T,N,S} = begin
    return setindex(u, el, CartesianIndex(ii))
end
setindex(u::PArray{T,N}, el::S, ii::CartesianIndex{N}) where {T,N,S} = begin
    root = copy(u._root)
    if u._root isa Array{T,N}
        # we take from the back (if phase isn't 0)
        jj = _locindices(u._phase, u._size, u._root, ii)
        root[jj] = el
    else
        (ii, jj) = _subindices(u._phase, u._height, ii)
        uu = root[ii]
        root[ii] = setindex(uu, el, jj)
    end
    return PArray{T,N}(u._size, u._phase, u._height, root)
end
#
# #setindex! ===================================================================
Base.setindex!(u::PArray{T,N}, el::S, args...) where {T,S,N} = begin
    throw(ErrorException("object of type $(typeof(u)) is immutable"))
end
Base.setindex!(u::TArray{T,N}, el::S, ii::Int) where {T,N,S} = begin
    idx = CartesianIndices(u._size)
    return setindex!(u, el, idx[ii])
end
Base.setindex!(u::TArray{T,N}, el::S, ii::Vararg{Int,N}) where {T,N,S} = begin
    return setindex!(u, el, CartesianIndex(ii))
end
Base.setindex!(u::TArray{T,N}, el::S, ii::CartesianIndex{N}) where {T,N,S} = begin
    root = u._root
    if root isa Array{T,N}
        # we take from the back (if phase isn't 0)
        jj = _locindices(u._phase, u._size, u._root, ii)
        root[jj] = el
    else
        (ii, jj) = _subindices(u._phase, u._height, ii)
        uu = root[ii]
        if uu isa PArray{T,N}
            uu = TArray{T,N}(uu)
            root[jj] = uu
        end
        setindex!(uu, el, jj)
    end
    return el
end
#
# #push ========================================================================
push(u::PArray{T,1}, el::S) where {T,S} = begin
    if u._root isa Array{T,1}
        n0 = u._size[1]
        ph = u._phase[1]
        n = n0 + ph
        if n == _PTREE_COUNT
            return _build_iarray([u, PArray{T,1}((1,), (0,), 0, T[el])])
        else
            root = Array{T,1}(undef, n0+1)
            copyto!(root, 1, u._root, size(u._root)[1] - (n0 - ph) + 1, n0)
            root[end] = el
            return PArray{T,1}(u._size .+ 1, u._phase, 0, root)
        end
    else
        uu = u._root[end]
        vv = push(uu, el)
        # if we went off the end...
        if vv._height == u._height
            # we might be able to append...
            if length(u._root) < _PTREE_COUNT
                root = copy(u._root)
                push!(root, vv._root[2])
                return PArray{T,1}(u._size .+ 1, u._phase, u._height, root)
            else
                v = PArray{T,1}((1,), (0,), u._height, PArray{T,1}[vv._root[end]])
                uv = PArray{T,1}[u, v]
                return PArray{T,1}(u._size .+ 1, u._phase, u._height + 1, uv)
            end
        else # otherwise we just update root
            root = copy(u._root)
            root[end] = vv
            return PArray{T,1}(u._size .+ 1, u._phase, u._height, root)
        end
    end
end
#
# #push! ========================================================================
Base.push!(u::TArray{T,1}, el::S) where {T,S} = begin
    if u._root isa Array{T,1}
        n0 = u._size[1]
        ph = u._phase[1]
        n = n0 + ph
        if n == _PTREE_COUNT
            uu = TArray{T,1}(u._size, u._phase, 0, u._root)
            vv = TArray{T,1}((1,),    (0,),     0, T[el])
            u._size = u._size .+ 1
            u._height = 1
            u._root = PArray{T,1}[uu, vv]
        else
            push!(u._root, el)
        end
    else
        u._size = u._size .+ 1
        uu = u._root[end]
        if uu isa PArray{T,1}
            uu = TArray{T,1}(uu)
            u._root[end] = uu
        end
        push!(uu, el)
        u._size = u._size .+ 1
        # if we went off the end...
        if uu._height == u._height
            u._root[end] = uu._root[1]
            vv = uu._root[2]
            # we might be able to append...
            if length(u._root) < _PTREE_COUNT
                push!(u._root, vv)
            else
                vv = TArray{T,1}((1,), (0,), u._height, TArray{T,1}[vv])
                uu._root = u._root
                uu._size = u._size .- 1
                u._root = TArray{T,1}[uu, vv]
                u._height = u._height + 1
            end
        end
    end
    return u
end
#
# #pop =========================================================================
pop(u::PVector{T}) where {T} = begin
    if u._root isa Array{T,1}
        return PVector{T}(u._size .- 1, u._phase, 0, u._root[1:end-1])
    else
        uu = u._root[end]
        vv = pop(uu)
        if vv._size[1] == 0
            if length(u._root) == 2
                return u._root[1]
            else
                return PVector{T}(u._size .- 1, u._phase, u._height, u._root[1:end-1])
            end
        else
            root = copy(u._root)
            root[end] = vv
            return PVector{T}(u._size .- 1, u._phase, u._height, root)
        end
    end
end
#
# #pop! ========================================================================
Base.pop!(u::TVector{T}) where {T} = begin
    if u._root isa Vector{T}
        x = pop!(u._root)
        u._size = u._size .- 1
        return x
    else
        uu = u._root[end]
        if uu isa PVector{T}
            uu = TVector{T}(uu)
            u._root[end] = uu
        end
        x = pop!(uu)
        u._size = u._size .- 1
        # we might have to pop this sub-array off (if it's now empty)
        if uu._size[1] == 0
            pop!(u._root)
            # we might even have to simplify u
            if length(u._root) == 1
                u._height -= 1
                u._root = u._root[1]._root
            end
        end
        return x
    end
end
#
# #pushfirst ===================================================================
pushfirst(u::PVector{T}, el::S) where {T,S} = begin
    if u._root isa Vector{T}
        n0 = u._size[1]
        ph = u._phase[1]
        if ph > 0
            root = copy(u._root)
            pushfirst!(root, el)
            return PArray{T,1}((n0 + 1,), (ph - 1,), 0, root)
        elseif n0 == _PTREE_COUNT
            u0 = PArray{T,1}((1,), (0,), 0, T[el])
            root = PArray{T,1}[u0, u]
            return PArray{T,1}((n0 + 1,), (_PTREE_COUNT - 1,), 1, root)
        else
            root = copy(u._root)
            pushfirst!(root, el)
            return PArray{T,1}((n0 + 1,), (0,), 0, root)
        end
    else
        if u._phase[1] > 0
            root = copy(u._root)
            root[1] = pushfirst(root[1], el)
            return PArray{T,1}(u._size .+ 1, u._phase .- 1, u._height, root)
        elseif length(u._root) < _PTREE_COUNT
            st = _height_to_stride(u._height)
            root = copy(u._root)
            pushfirst!(root, PArray{T,1}((1,), (0,), 1, [el]))
            return PArray{T,1}(u._size .+ 1, (st - 1,), u._height, root)
        else
            st = _height_to_stride(u._height + 1)
            v = PArray{T,1}((1,), (0,), 0, T[el])
            return PArray{T,1}(u._size .+ 1, (st - 1,), u._height + 1, PArray{T,1}[v, u])
        end
    end
end
#
# #pushfirst! ==================================================================
Base.pushfirst!(u::TArray{T,1}, el::S) where {T,S} = begin
    if u._root isa Array{T,1}
        n0 = u._size[1]
        ph = u._phase[1]
        n = n0 + ph
        if n == _PTREE_COUNT
            uu = TArray{T,1}(u._size, u._phase, 0, u._root)
            vv = TArray{T,1}((1,),    (0,),     0, T[el])
            u._size = u._size .+ 1
            u._height = 1
            u._root = PArray{T,1}[uu, vv]
        else
            push!(u._root, el)
        end
    else
        u._size = u._size .+ 1
        uu = u._root[end]
        if uu isa PArray{T,1}
            uu = TArray{T,1}(uu)
            u._root[end] = uu
        end
        push!(uu, el)
        u._size = u._size .+ 1
        # if we went off the end...
        if uu._height == u._height
            u._root[end] = uu._root[1]
            vv = uu._root[2]
            # we might be able to append...
            if length(u._root) < _PTREE_COUNT
                push!(u._root, vv)
            else
                vv = TArray{T,1}((1,), (0,), u._height, TArray{T,1}[vv])
                uu._root = u._root
                uu._size = u._size .- 1
                u._root = TArray{T,1}[uu, vv]
                u._height = u._height + 1
            end
        end
    end
    return u
end
#
# #popfirst ====================================================================
popfirst(u::PVector{T}) where {T} = begin
    if u._root isa Vector{T}
        (u._size[1] == 0) && throw(ArgumentError("array must be non-empty"))
        return PVector{T}(u._size .- 1, u._phase .+ 1, 0, u._root)
    else
        uu = u._root[1]
        vv = popfirst(uu)
        if vv._size[1] == 0
            if length(u._root) == 2
                return u._root[2]
            else
                return PVector{T}(u._size .- 1, (0,), u._height, u._root[2:end])
            end
        else
            root = copy(u._root)
            root[1] = vv
            return PVector{T}(u._size .- 1, u._phase .+ 1, u._height, root)
        end
    end
end
#
# #popfirst! ===================================================================
Base.popfirst!(u::TVector{T}) where {T} = begin
    if u._root isa Vector{T}
        x = popfirst!(u._root)
        u._size = u._size .- 1
        return x
    else
        uu = u._root[1]
        if uu isa PVector{T}
            uu = TVector{T}(uu)
            u._root[1] = uu
        end
        x = popfirst!(uu)
        if uu._size[1] == 0
            if length(u._root) == 2
                vv = u._root[2]
                u._root = vv._root
                u._height = vv._height
                u._phase = (0,)
                u._size = vv._size
            else
                popfirst!(u._root)
                u._phase = (0,)
                u._size = u._size .- 1
            end
        else
            u._phase = u._phase .+ 1
            u._size = u._size .- 1
        end
        return x
    end
end
#
# #isequiv =====================================================================
@_ptmethod IArray function isequiv(u::IArray{T,N}, v::IArray{S,N}) where {T,S,N}
    (u._size == v._size) || return false
    (u === v) && return true
    for (a,b) in zip(u, v)
        isequiv(a, b) || return false
    end
    return true
end
#
# #equivhash ===================================================================
@_ptmethod IArray function equivhash(u::IArray{T,N}) where {T,N}
    h = hash(u._size)
    (IArray === TArray) && (h += 1)
    for x in u
        h += 31 * equivhash(x)
    end
    return h
end
