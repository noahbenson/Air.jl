
import Base.setindex
import Base.delete!

################################################################################
# Abstract Types
"""
    AbstractPDict{K,V}

AbstractPDict is a subtype of AbstractDict that is extended only by persistent
dictionary types such as PDict and LazyDict.
"""
abstract type AbstractPDict{K,V} <: AbstractDict{K,V} end
"""
    AbstractPSet{T}

AbstractPSet is an abstract type extended by persistent set types such as PSet,
PIdSet, and PEqualSet, as well as the weighted persistent set types.
"""
abstract type AbstractPSet{T} <: AbstractSet{T} end
"""
    AbstractPArray{T,N}

The Abstract persistent array type represents any Array-like type that is
persistent. Reified types include PArray{T,N} and LazyArray{T,N}.
"""
abstract type AbstractPArray{T,N} <: AbstractArray{T,N} end


################################################################################
# Mutability
"""
    Mutability
The Mutability type is an abstract trait type with children Mutable, Immutable,
and MutabilityUnknown.
"""
abstract type Mutability end
struct Mutable <: Mutability end
struct Immutable <: Mutability end
struct UnknownMutability <: Mutability end
const _MUT_TYPE = Mutable()
const _IMM_TYPE = Immutable()
const _UNK_TYPE = UnknownMutability()
"""
    mutability(obj)
Yields an object of type Mutable, Immutable, or UnknownMutability depending on 
whether the given type object is mutable, immutable, or unknown.
"""
mutability(T::Type) = begin
    if typeof(T) === Union
        let mA = mutability(T.a)
            return (mA === mutability(T.b) ? mA : _UNK_TYPE)
        end
    else
        return _UNK_TYPE
    end
end
# Implementations for mutability below
mutability(T::DataType) = begin
    if isbitstype(T)
        return _IMM_TYPE
    elseif !T.isconcretetype
        return _UNK_TYPE
    elseif T.mutable
        return _MUT_TYPE
    else
        return _IMM_TYPE
    end
end
mutability(::Core.TypeofBottom) = _UNK_TYPE
mutability(T::UnionAll) = mutability(T{T.var.ub})
# these two come out wrong in the above code, according to the documentation:
mutability(::Type{String}) = _IMM_TYPE
mutability(::Type{Symbol}) = _IMM_TYPE
# These also seem to come out wrong:
mutability(::Type{Set}) = _MUT_TYPE
mutability(::Type{Set{T}}) where {T} = _MUT_TYPE
mutability(::Type{Base.IdSet}) = _MUT_TYPE
mutability(::Type{Base.IdSet{T}}) where {T} = _MUT_TYPE
mutability(::Type{Dict}) = _MUT_TYPE
mutability(::Type{Dict{K,V}}) where {K,V} = _MUT_TYPE
mutability(::Type{IdDict}) = _MUT_TYPE
mutability(::Type{IdDict{K,V}}) where {K,V} = _MUT_TYPE


################################################################################
# #ismut, #isimm, #isimmtype, #ismuttype functions for testing immutabilities
_ismuttype(::Mutable) = true
_ismuttype(::Immutable) = false
_ismuttype(::UnknownMutability) = false
"""
    ismuttype(T)
Yields true if the given type T is known to be mutable and false otherwise.
"""
ismuttype(T::Type) = _ismuttype(mutability(T))
_isimmtype(::Mutable) = false
_isimmtype(::Immutable) = true
_isimmtype(::UnknownMutability) = false
"""
    isimmtype(T)
Yields true if the given type T is known to be immutable and false otherwise.
"""
isimmtype(T::Type) = _isimmtype(mutability(T))
"""
    ismut(x)
Yields true if x is mutable and false if x is immutable or has unknown mutability.
"""
ismut(::T) where {T} = _ismuttype(mutability(T))
"""
    isimm(x)
Yields true if x is immutable and false if x is mutable or has unknown mutability.
"""
isimm(::T) where {T} = _isimmtype(mutability(T))

################################################################################
# #mutable_cotype #immutable_cotype #abstract_cotype

_mutable_cotype(::Mutable, x) = error(
    "Mutable type $(typeof(x)) cannot have a mutable cotype")
_mutable_cotype(::Immutable, x) = error(
    "Immutable type $(typeof(x)) has not overloaded mutable_cotype()")
_mutable_cotype(::UnknownMutability, x) = error(
    "Type $(typeof(x)) with unknown mutability does not have a mutable cotype")
"""
    mutable_cotype(T)
Yields the mutable equivalent of the immutable or abstract type T. If T is
already a mutable type or there is no mutable cotype of T then an error is
raised.
"""
mutable_cotype(T::Type) = _mutable_cotype(mutability(T))

_immutable_cotype(::Immutable, x) = error(
    "Immutable type $(typeof(x)) cannot have an immutable cotype")
_immutable_cotype(::Mutable, x) = error(
    "Mutable type $(typeof(x)) has not overloaded immutable_cotype()")
_immmutable_cotype(::UnknownMutability, x) = error(
    "Type $(typeof(x)) with unknown mutability does not have an immutable cotype")
"""
    immutable_cotype(T)
Yields the immutable equivalent of the mutable or abstract type T. If T is
already an immutable type or there is no immutable cotype of T then an error is
raised.
"""
immutable_cotype(T::Type) = _immutable_cotype(mutability(T), x)

_abstract_cotype(::Mutable, x) = error(
    "Mutable type $(typeof(x)) has not overloaded abstract_cotype()")
_abstract_cotype(::Immutable, x) = error(
    "Immutable type $(typeof(x)) has not overloaded abstract_cotype()")
_abstract_cotype(::UnknownMutability, x) = error(
    "Type $(typeof(x)) with unknown mutability cannot have an abstract cotype")
"""
    abstract_cotype(T)
Yields the abstract equivalent of the immutable or mutable type T. If T is
already an abstract type or there is no mutable cotype of T then an error is
raised.
"""
abstract_cotype(T::Type) = _abstract_cotype(mutability(T), x)

################################################################################
# #aimtypes
"""
    @aimtypes A I M
Declares that the given types are related types; specifically:
 * A is an abstract type,
 * I <: A is an immutable type,
 * M <: A is a mutable type,
 * freeze(x::M) yields an object of type I, and
 * thaw(x::I) yields an object of type M.
"""
macro aimtypes(A::Symbol, I::Symbol, M::Symbol)
    return quote
        $(esc(:abstract_cotype))(x::Type{$I}) = $A
        $(esc(:abstract_cotype))(x::Type{$M}) = $A
        $(esc(:immutable_cotype))(x::Type{$A}) = $I
        $(esc(:immutable_cotype))(x::Type{$M}) = $I
        $(esc(:mutable_cotype))(x::Type{$A}) = $M
        $(esc(:mutable_cotype))(x::Type{$I}) = $M
        $(esc(:mutability))(::Type{$I}) = Immutable()
        $(esc(:mutability))(::Type{$M}) = Mmutable()
        $(esc(:copy))(x::$I) = x

    end
end

################################################################################
# #freeze, #thaw, #deepfreeze, #deepthaw

_freeze(::Immutable, x) = x
_freeze(::UnknownMutability, x) = error(
    "type $(typeof(x)) has not overloaded the freeze() method")
_freeze(::Mutable, x) = error(
    "type $(typeof(x)) has not overloaded the freeze() method")
"""
    freeze(x)
Yields x if x is immutable or a shallow immutable copy of x if x is not
immutable. If x is of a non-immutable type that has not overloaded the freeze()
method, then an error is raised.
"""
freeze(x::T) where {T} = _freeze(mutability(T), x)

"""
    deepfreeze(x)
Yields a deep copy of x in which x and all of its sub-values have been converted
from mutable into immutable objects. Immutable objects are left as-is unless
their children are changed.

Like the deepcopy function, deepfreeze() can be customized; however, unlike
deepcopy() which employs a deepcopy_internal() co-routine, deepfreeze() should
be overloaded directly for the particular type.
"""
deepfreeze(@nospecialize(x)) = begin
    T = typeof(x)::DataType
    isbitstype(T) && return x
    # There are a few possibilities here:
    # (1) the object is mutable, in which case we'll need to freeze() it
    # (2) the object is immutable, in which case it's fine except maybe its
    #     children;
    # (3) the object may have children whose identities change, in which case
    #     we need to make a new object even if it's immutable
    # The easiest case to catch is the mutable case, in which case we freeze it
    # then deepfreeze the result (to ensure that children are also frozen).
    isimmtype(T) || return deepfreeze(freeze(x))
    # For an immutable object, we start by freezing the children/fields:
    n = nfields(x)
    changed = false
    field_defined = Vector{Bool}()
    field_changed = Vector{Bool}()
    new_field_val = Vector{Any}()
    for i in 1:n
        if isdefined(x, i)
            old_obj = getfield(x, i)
            new_obj = deepfreeze(old_obj)
            push!(new_field_val, new_obj)
            push!(field_changed, old_obj !== new_obj)
            push!(field_defined, true)
            old_obj === new_obj || (changed = true)
        else
            push!(new_field_val, nothing)
            push!(field_changed, false)
            push!(field_defined, false)
        end
    end
    # We now have the updated fields stored as these three arrays; if the
    # field was not defined in the old object then it should not be set in any
    # new object.
    if !changed
        # In this case, none of the children have changed and the object is
        # known to be immutable; we do not even have to call freeze on it!
        return x
    else
        # The object is immutable, but we need to freeze its children:
        y = ccall(:jl_new_struct_uninit, Any, (Any,), T)
        ftypes = fieldtypes(T)
        for i in 1:n
            isdef = field_defined[i]
            isnew = field_changed[i]
            new_val = new_field_val[i]
            if isdef
                # if it's a new object and the types don't match, we should
                # throw an error; otherwise we can set it
                if isnew && !(typeof(new_val) <: ftypes[i])
                    error("field $i of type $T was not frozen into" *
                          " a sub-type of $(ftypes[i])")
                else
                    ccall(:jl_set_nth_field,
                          Cvoid, (Any, Csize_t, Any), y, i-1,
                          new_val)
                end
            end
        end
        return y
    end
end
const TERMINAL_TYPES = Union{Symbol,String,Core.MethodInstance,Method,GlobalRef,
                             DataType,Union,UnionAll,Task,Regex}
deepfreeze(x::TERMINAL_TYPES) = x
deepfreeze(x::Tuple) = ntuple(i->deepfreeze(x[i]), length(x))
deepfreeze(x::Module) = error("deepfreeze of Modules not supported")
deepfreeze(x::Nothing) = nothing

_thaw(::Mutable, x) = x
_thaw(::UnknownMutability, x) = error(
    "type $(typeof(x)) has not overloaded the freeze() method")
_thaw(::Immutable, x) = error(
    "type $(typeof(x)) has not overloaded the freeze() method")
"""
    thaw(x)
Yields a mutable version of the object x if x is an immutable object; otherwise
if x is already mutable yields x unchanged. To get a mutable copy of an object
regardless of whether it is immutable, use thawcopy(). If x cannot be converted
to an immutable object, an error is raised.
"""
thaw(x::T) where {T} = _thaw(mutability(T), x)

"""
    deepthaw(x)
Yields a mutable version of the object x after also deepthawing any relevant
children of x. Note that deepthaw(x) should never affect x or any of its
children whether x is mutable or not (i.e., when the children of x change, a
copy of x with the new children is returned even if x is already mutable).

Unlike thaw(), freeze(), and deepfreeze(), deepthaw() does not raise errors when
it encounters immutable objects that it cannot thaw. This is because it is
presumed that immutable terminals (like String and Int64) should be left as they
are.

If you wish to obtain a mutable deepcopy of an object regardless of whether the
original object is mutable or immutable, use the deepthawcopy() function.
"""
deepthaw(@nospecialize(x)) = begin
    T = typeof(x)::DataType
    isbitstype(T) && return x
    # There are a few possibilities here:
    # (1) the object is immutable, in which case we'll need to thaw() it
    # (2) the object is mutable, in which case it's fine except maybe its
    #     children;
    # (3) the object may have children whose identities change when we deepthaw
    #     them, in which case we need to copy the object then update
    # The easiest case to catch is the immutable case, in which case we thaw it
    # then deepthaw the result (to ensure that children are also deepthawed).
    if !ismuttype(T)
        try
            x = thaw(x)
        catch e
            return x
        end
    end
    # We start by deepthaw'ing the children/fields:
    n = nfields(x)
    changed = false
    field_defined = Vector{Bool}()
    field_changed = Vector{Bool}()
    new_field_val = Vector{Any}()
    for i in 1:n
        if isdefined(x, i)
            old_obj = getfield(x, i)
            new_obj = deepthaw(old_obj)
            push!(new_field_val, new_obj)
            push!(field_changed, old_obj !== new_obj)
            push!(field_defined, true)
            old_obj === new_obj || (changed = true)
        else
            push!(new_field_val, nothing)
            push!(field_changed, false)
            push!(field_defined, false)
        end
    end
    # We now have the updated fields stored as these three arrays; if the
    # field was not defined in the old object then it should not be set in any
    # new object.
    if !changed
        # In this case, none of the children have changed and the object is
        # known to be mutable already
        return x
    else
        # The object is mutable, but we need to update its children:
        y = ccall(:jl_new_struct_uninit, Any, (Any,), T)
        ftypes = fieldtypes(T)
        for i in 1:n
            isdef = field_defined[i]
            isnew = field_changed[i]
            new_val = new_field_val[i]
            if isdef
                # if it's a new object and the types don't match, we should
                # throw an error; otherwise we can set it
                if isnew && !(typeof(new_val) <: ftypes[i])
                    error("field $i of type $T was not thawed into" *
                          " a sub-type of $(ftypes[i])")
                else
                    ccall(:jl_set_nth_field,
                          Cvoid, (Any, Csize_t, Any), y, i-1,
                          new_val)
                end
            end
        end
        return y
    end
end
deepthaw(x::TERMINAL_TYPES) = x
deepthaw(x::Tuple) = ntuple(i->deepthaw(x[i]), length(x))
deepthaw(x::Module) = x
deepthaw(x::Nothing) = nothing

_thawcopy(::Mutable, x) = copy(x)
_thawcopy(::Immutable, x) = thaw(x)
_thawcopy(::UnknownMutability, x) = copy(x)
"""
    thawcopy(x)
Yields thaw(x) if x is immutable and copy(x) if x is mutable.
"""
thawcopy(x::T) where {T} = _thawcopy(mutability(T), x)

"""
    deepthawcopy(x)
Yields deepthaw(deepcopy(x)).
"""
deepthawcopy(x) = deepthaw(deepcopy(x))

# isequiv should always default to === unnless both objects are immutable
# or the mutabiilities are different (in which case the are not equiv).
_isequiv(::UnknownMutability, ::UnknownMutability, t, s) = (t === s)
_isequiv(::Mutable,           ::Mutable,           t, s) = (t === s)
_isequiv(::UnknownMutability, ::Mutable,           t, s) = false
_isequiv(::UnknownMutability, ::Immutable,         t, s) = false
_isequiv(::Mutable,           ::UnknownMutability, t, s) = false
_isequiv(::Mutable,           ::Immutable,         t, s) = false
_isequiv(::Immutable,         ::UnknownMutability, t, s) = false
_isequiv(::Immutable,         ::Mutable,           t, s) = false
_isequiv(::Immutable,         ::Immutable,         t::T, s::S) where {T,S} = false
_isequiv(::Immutable,         ::Immutable,         t::T, s::T) where {T} = begin
    t === s && return true
    isbitstype(T) && return isequal(t, s)
    # call isequiv on each of the fields
    let n = nfields(t)
        n == 0 && return isequal(t, s)
        for k in 1:n
            (isdefined(t,k) != isdefined(s,k)) && return false
            !isequiv(getfield(t,k), getfield(s,k)) && return false
        end
        # All fields are equivalent, so return true
        return true
    end
end
function _isequiv(
    ::Immutable, ::Immutable, t::T, s::S
) where {T<:Number, S<:Number}
    return isequal(t, s)
end
function _isequiv(
    ::Immutable, ::Immutable, t::T, s::S
) where {T<:AbstractString, S<:AbstractString}
    return isequal(t, s)
end
function _isequiv(
    ::Immutable, ::Immutable, t::T, s::S
) where {
    KT,VT, T<:AbstractDict{KT,VT},
    KS,VS, S<:AbstractDict{KS,VS}
}
    # Dictionaries can only be equivalent if they use the same equal function.
    (equalfn(T) === equalfn(S)) || return false
    # Dictionaries are equivalent when immutable and when they share equivalent
    # key/value pairs.
    (t === s) && return true
    (length(t) == length(s)) || return false
    for pair in t
        in(pair, s, isequiv) || return false
    end
    return true
end
function _isequiv(
    ::Immutable, ::Immutable, t::T, s::S
) where {
    TT, T<:AbstractSet{TT},
    SS, S<:AbstractSet{SS}
}
    # Sets can only be equivalent if they use the same equal function.
    (equalfn(T) === equalfn(S)) || return false
    # Dictionaries are equivalent when immutable and when they share equivalent
    # key/value pairs.
    (t === s) && return true
    (length(t) == length(s)) || return false
    for el in t
        in(el, s, isequiv) || return false
    end
    return true
end
function _isequiv(
    ::Immutable, ::Immutable, t::T, s::S
) where {
    TT, NT, T<:AbstractArray{TT,NT},
    SS, NS, S<:AbstractArray{SS,NS}
}
    # Arrays are equivalent when immutable and when they share equivalent
    # key/value pairs.
    (t === s) && return true
    (size(t) == size(s)) || return false
    for (xt,xs) in zip(t,s)
        isequiv(xt, xs) || return false
    end
    return true
end

"""
    isequiv(x, y)
Yields true if x and y are equivalent and false otherwise.
Equivalents is similar to equality, but it distinguishes
between mutable itens that are equal (on the grounds that
they may not be equal in the future). The rules for 
deciding equivalence are:

1. If both x and y are immutable, then isequiv(x,y) yields
   isequal(x, y)
2. Otherwise, yields (x === y)

Note that immutable types should overload isequiv() in order
to ensure that they are considered equivalent only when
their subparts (like keys and values) are considered
equivalent rather than equal. By default, for a structure
type, all 
"""
isequiv(t::T, s::S) where {T, S} = _isequiv(mutability(T), mutability(S), t, s)
# Some builtin equivs:
isequiv(t::Pair{KT,VT}, s::Pair{KS,VS}) where {KT,VT,KS,VS} =
    isequiv(t.first, s.first) && isequiv(t.second, s.second)
isequiv(t::Tuple, s::Tuple) = begin
    (length(t) == length(s)) || return false
    for (xt,xs) in zip(t, s)
        isequiv(xt,xs) || return false
    end
    return true
end
isequiv(t::Symbol, s::Symbol) = isequal(t, s)
isequiv(t::String, s::String) = isequal(t, s)

_equivhash(::UnknownMutability, x) = objectid(x)
_equivhash(::Mutable, x) = objectid(x)
_equivhash(::Immutable, x::T) where {T} = begin
    isbitstype(T) && return hash(x)
    # call isequiv on each of the fields
    let n = nfields(x), h = objectid(T) + 41
        n == 0 && return objectid(x)
        for k in 1:n
            isdefined(x, k) || continue
            h += (~k) * equivhash(getfield(x, k))
        end
        return h
    end
end
_equivhash(::Immutable, x::T) where {T <: Number} = hash(x)
_equivhash(::Immutable, x::T) where {T <: AbstractString} = hash(x)
# Special versions for dictionaries, sets, and arrays...
_equivhash(::Immutable, u::DD) where {K,V,DD<:AbstractDict{K,V}} = begin
    h = UInt(0x0)
    for (k,v) in u
        h += equivhash(k) ⊻ equivhash(v)
    end
    return h
end
_equivhash(::Immutable, u::SS) where {T,SS<:AbstractSet{T}} = begin
    h = UInt(0x0)
    for x in u
        h += equivhash(x)
    end
    return h
end
_equivhash(::Immutable, u::AA) where {T,D,AA<:AbstractArray{T,D}} = begin
    h = UInt(sum(size(u) .* collect(1:D)))
    for (k,v) in enumerate(u)
        h += UInt(k) * equivhash(v)
    end
    return h
end
"""
    equivhash(x)
Yields a hash value appropriate for the isequiv() equality
function. Generally, the hasheq of an immutable object is
its hash; the equivhash of any other object is its objectid.
"""
equivhash(x::T) where {T} = _equivhash(mutability(T), x)
equivhash(x::Symbol) = hash(x)
equivhash(x::String) = hash(x)
# Tuples get equivhash'ed also, but slightly differently than arrays
equivhash(u::Tuple) = begin
    h = UInt(length(u))
    for x in u
        h = (0x1 + h) * equivhash(x)
    end
    return h
end
equivhash(u::Pair{K,V}) where {K,V} =
    (0x1f + equivhash(u.first)) * equivhash(u.second)

"""
    equalfn(x)
If x is an object (such as a PSet or Dict) that has an opinion about equality,
equalfn(x) returns the function that it uses.
"""
equalfn(x::T) where {T} = equalfn(T)
equalfn(x::Type{T}) where {T} = throw(ArgumentError("no equalfn for type $T"))
equalfn(x::Type{Dict}) = isequal
equalfn(x::Type{Set}) = isequal
equalfn(x::Type{IdDict}) = (===)
equalfn(x::Type{Base.IdSet}) = (===)
equalfn(x::Type{Dict{K,V}}) where {K,V} = isequal
equalfn(x::Type{Set{T}}) where {T} = isequal
equalfn(x::Type{IdDict{K,V}}) where {K,V} = (===)
equalfn(x::Type{Base.IdSet{T}}) where {T} = (===)
"""
    hashfn(x)
If x is an object (such as a PSet or Dict) that has an opinion about how it
hashes objects, hashfn(x) returns the function that it uses. It is
sufficient in almost all circumstances to define `equalfn(T)`; the `hashfn`
should always match the `equalfn` regardless.
"""
hashfn(x::T) where {T} = hashfn(equalfn(T))
hashfn(::Type{T}) where {T} = hashfn(equalfn(T))

# equalfn and hashfn can also be used on their respective functions
equalfn(::typeof(objectid))  = (===)
equalfn(::typeof(equivhash)) = isequiv
equalfn(::typeof(hash))      = isequal
hashfn(::typeof(===))     = objectid
hashfn(::typeof(isequiv)) = equivhash
hashfn(::typeof(isequal)) = hashcode

#
# #setindex ====================================================================
# The setindex for tuples from StaticArrays is better than this version:
#"""
#    setindex(tup, val, index) 
#Yields a copy of the given n-tuple with the value at the given indiex changed to
#val.
#"""
#setindex(u::NTuple{N,T}, el::S, I...) where {T,N,S} = begin
#    # make a temporary array...
#    ar = T[u...]
#    ar[I...] = el
#    return NTuple{N,T}(ar)
#end
"""
    setindex(arr, val, index) 
Yields a copy of the given array with the value at the given indiex changed to
val.

Note that setindex() *always* returns a copy of the array arr, so unless arr is
a persistent array (PArray), then this operation is not very efficient.
"""
setindex(u::Array{T,N}, x::S, I...) where {T,N,S} =
    Base.setindex!(copy(u), x, I...)
"""
    setindex(d, val, key)
Yields a copy of the given dictionary d with the given key associated with the
given value val.
"""
setindex(d::Dict{K,V}, v::U, k::J) where {K,V,U,J} =
    Base.setindex!(copy(d), v, k)
setindex(d::IdDict{K,V}, v::U, k::J) where {K,V,U,J} =
    Base.setindex!(copy(d), v, k)
# non-typed or redirected calls:
setindex(u::SubArray, args...) = Base.setindex!(copy(u), args...)
setindex(u::Base.ReshapedArray, args...) = Base.setindex!(copy(u), args...)
setindex(u::Base.ReinterpretArray, args...) =
    Base.setindex!(copy(u), args...)
setindex(u::BitArray, args...) = Base.setindex!(copy(u), args...)
#
# #push ========================================================================
"""
    push(tup, val)
Yields a copy of the given n-tuple with the given val appended.
"""
push(u::NTuple{N,T}, val::S) where {T,N,S} = begin
    return NTuple{N+1,T}(T[u..., val])
end
"""
    push(arr, val)
Yields a copy of the given (1D) array with the given val appended.

Note that push() *always* returns a copy of the array arr, so unless arr is a
persistent array (PArray), then this operation is not very efficient.
"""
push(u::Vector{T}, val::S) where {T,N,S} = T[u..., val]
"""
    push(d, k => v)
Yields a copy of the given dictionary d with the given key/value pair added.
This is equivalent to setindex(d, v, k).
"""
push(d::DT, kv::Pair{K,V}) where {K,V,DT<:AbstractDict} = begin
    return setindex(d, kv[2], kv[1])
end
"""
    push(s, u)
Yields a copy of the given set or bit-set object s with the given object u
appended to the set.
"""
push(s::Set{T}, u::S) where {T,S} = push!(copy(s), u)
push(s::Base.IdSet{T}, u::S) where {T,S} = push!(copy(s), u)
# non-typed
#push(A, a) = push!(copy(A), a)
push(A, a, b...) = push(push(A, a), b...)
#
# #pop =========================================================================
"""
    pop(tup)
Yields a tuple (last, most) where last is the last element of the given tuple
tup and most is a duplicate tuple of all but the last element of tup.
"""
pop(u::NTuple{N,T}) where {N,T} = (u[end], u[1:end-1])
pop(u::Tuple{}) = throw(ArgumentError("n-tuple must be non-empty"))
"""
    pop(arr)
Yields a tuple (last, most) where last is the final element of the given (1D)
array arr and most is a copy of the given array without its last element.

Note that this function *always* makes a copy of the array arr, so unless arr is
a persistent array (PArray) this operation is not very efficient.
"""
pop(u::AbstractArray{T,1}) where {T} = (u[end], copy(u[1:end-1]))
"""
    pop(d)
Yields a tuple (k=>v, rem) where k=>v is an element that has been removed from
the copy of dictionary d, rem (i.e., remainder of d).
"""
pop(d::Dict{K,V}) where {K,V} = begin
    rem = copy(d)
    kv = pop!(rem)
    return (kv, rem)
end
pop(d::IdDict{K,V}) where {K,V} = begin
    rem = copy(d)
    kv = pop!(rem)
    return (kv, rem)
end
"""
    pop(d, k)
Yields a tuple (v, rem) where v is the value, mapped to k in d that has been
removed from the copy of dictionary d, rem (i.e., remainder of d).
"""
pop(d::Dict{K,V}, k::J) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k)
    return (v, rem)
end
pop(d::IdDict{K,V}, k::J) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k)
    return (v, rem)
end
"""
    pop(d, k, dv)
Yields a tuple (v, rem) where v is the value, mapped to k in d that has been
removed from the copy of dictionary d, rem (i.e., remainder of d). If the key
k is not found in d, then rem will be a copy of d and v will be equal to dv.
"""
pop(d::Dict{K,V}, k::J, dv) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k, dv)
    return (v, rem)
end
pop(d::IdDict{K,V}, k::J, dv) where {K,V,J} = begin
    rem = copy(d)
    v = pop!(rem, k, dv)
    return (v, rem)
end
"""
    pop(s)
Yields a tuple (u, rem) where u is an element from the set s and rem is the
remainder of the set. Note that this always makes a copy of s.
"""
pop(s::Set{T}) where {T} = begin
    rem = copy(s)
    u = pop!(rem)
    return (u, rem)
end
pop(s::Base.IdSet{T}) where {T} = begin
    rem = copy(s)
    u = pop!(rem)
    return (u, rem)
end
pop(s::BitSet) where {T} = begin
    rem = copy(s)
    u = pop!(rem)
    return (u, rem)
end
"""
    pop(s, u)
Yields a tuple (u, rem) where u is an element in s and rem is a copy of s
with element u removed.
"""
pop(s::Set{T}, u::S) where {T,S} = begin
    rem = copy(s)
    u = pop!(rem, u)
    return (u, rem)
end
pop(s::Base.IdSet{T}, u::S) where {T,S} = begin
    rem = copy(s)
    u = pop!(rem, u)
    return (u, rem)
end
pop(s::BitSet, u) = begin
    rem = copy(s)
    u = pop!(rem, u)
    return (u, rem)
end
#
# #insert ======================================================================
"""
    insert(arr, idx, val)
Yields a copy of the given vector arr with the given value inserted at the given
index. Equivalent to insert!(copy(arr), idx, va).
"""
insert(a::Vector{T}, k::K, val::S) where {T,K<:Integer,S} = begin
    n = length(a) + 1
    (k > n) && throw(
        ArgumentError("insert: $k is out of range for Vector of length $(n-1)"))
    out = Vector{T}(undef, n)
    (k > 1) && copyto!(out, 1, a, 1, k - 1)
    @inbounds out[k] = val
    (k < n) && copyto!(out, k + 1, a, k, n - k)
    return out
end
insert(a::BitArray{T}, idx::II, val::S) where {T,S,II<:Integer} =
    insert!(copy(a), idx, val)
insert(a::Tuple, idx::II, val::T) where {II<:Integer,T} =
    (a[1:idx]..., val, a[idx:end]...)
insert(a::Tuple{}, idx::II, val::T) where {II<:Integer,T} = begin
    if idx == 1
        return (val,)
    else
        throw(ArgumentError("invalid index $idx for tuple of size 0"))
    end
end
# For tuples, we want to generate versions of this for NTuples up to size 64;
# beyond that we can use a generic function.
macro _tuple_insert_gencode(N::Int)
    # We'll want to refer to the tuple elements:
    els = [:(tup[$k]) for k in 1:N]
    # Build up the if-elseif-else expression, starting with the else:
    ifexpr = :(throw(ArgumentError(
        "insert: $k if out of range for a Tuple of length $(length(tup))")))
    ifexpr = Expr(:elseif, :(k == $(N+1)), :(push(tup, el)), ifexpr) 
    for k in N:-1:1
        ifexpr = Expr(k == 1 ? :if : :elseif,
                      :(k == $k),
                      :(($(els[1:k-1]...), el, $(els[k:end]...))),
                      ifexpr)
    end
    return quote
        insert(tup::NTuple{$N,T}, k::K, el::S) where {T,K<:Integer,S} = $ifexpr
    end |> esc
end
# Now generate functions for up to 64:
(@_tuple_insert_gencode  1)
(@_tuple_insert_gencode  2)
(@_tuple_insert_gencode  3)
(@_tuple_insert_gencode  4)
(@_tuple_insert_gencode  5)
(@_tuple_insert_gencode  6)
(@_tuple_insert_gencode  7)
(@_tuple_insert_gencode  8)
(@_tuple_insert_gencode  9)
(@_tuple_insert_gencode 10)
(@_tuple_insert_gencode 11)
(@_tuple_insert_gencode 12)
(@_tuple_insert_gencode 13)
(@_tuple_insert_gencode 14)
(@_tuple_insert_gencode 15)
(@_tuple_insert_gencode 16)
(@_tuple_insert_gencode 17)
(@_tuple_insert_gencode 18)
(@_tuple_insert_gencode 19)
(@_tuple_insert_gencode 20)
(@_tuple_insert_gencode 21)
(@_tuple_insert_gencode 22)
(@_tuple_insert_gencode 23)
(@_tuple_insert_gencode 24)
(@_tuple_insert_gencode 25)
(@_tuple_insert_gencode 26)
(@_tuple_insert_gencode 27)
(@_tuple_insert_gencode 28)
(@_tuple_insert_gencode 29)
(@_tuple_insert_gencode 30)
(@_tuple_insert_gencode 31)
(@_tuple_insert_gencode 32)
(@_tuple_insert_gencode 33)
(@_tuple_insert_gencode 34)
(@_tuple_insert_gencode 35)
(@_tuple_insert_gencode 36)
(@_tuple_insert_gencode 37)
(@_tuple_insert_gencode 38)
(@_tuple_insert_gencode 39)
(@_tuple_insert_gencode 40)
(@_tuple_insert_gencode 41)
(@_tuple_insert_gencode 42)
(@_tuple_insert_gencode 43)
(@_tuple_insert_gencode 44)
(@_tuple_insert_gencode 45)
(@_tuple_insert_gencode 46)
(@_tuple_insert_gencode 47)
(@_tuple_insert_gencode 48)
(@_tuple_insert_gencode 49)
(@_tuple_insert_gencode 50)
(@_tuple_insert_gencode 51)
(@_tuple_insert_gencode 52)
(@_tuple_insert_gencode 53)
(@_tuple_insert_gencode 54)
(@_tuple_insert_gencode 55)
(@_tuple_insert_gencode 56)
(@_tuple_insert_gencode 57)
(@_tuple_insert_gencode 58)
(@_tuple_insert_gencode 59)
(@_tuple_insert_gencode 60)
(@_tuple_insert_gencode 61)
(@_tuple_insert_gencode 62)
(@_tuple_insert_gencode 63)
(@_tuple_insert_gencode 64)
#
# #delete ======================================================================
"""
    delete(d, k)
Yields a copy of the dictionary d with the given key k deleted.
"""
delete(d::Dict{K,V}, k::J) where {K,V,J} = delete!(copy(d), k)
"""
    delete(v, k)
Yields a copy of the vector v with index k deleted.
"""
delete(u::Vector{T}, k::K) where {T,K<:Integer} = begin
    n = length(u)
    (n == 0 || k < 0 || k > n) && throw(
        ArgumentError("delete: $k is out of range for Vector of length $n"))
    out = Vector{T}(undef, n - 1)
    (k > 1) && copyto!(out, 1, u, 1, k - 1)
    (k < n) && copyto!(out, k, u, k + 1, n - k)
    return out
end
delete(a::Tuple, idx::K) where {K<:Integer,T} = (a[1:idx-1]..., a[idx+1:end]...)
delete(a::Tuple{}, idx::II, val::T) where {II<:Integer,T} =
    throw(ArgumentError("delete: $idx is out of range for a Tuple{}"))
# For tuples, we want to generate versions of this for NTuples up to size 64;
# beyond that we can use a generic function.
macro _tuple_delete_gencode(N::Int)
    # We'll want to refer to the tuple elements:
    els = [:(tup[$k]) for k in 1:N]
    # Build up the if-elseif-else expression, starting with the else:
    ifexpr = :(throw(ArgumentError(
        "delete: $k if out of range for a Tuple of length $(length(tup))")))
    for k in N:-1:1
        ifexpr = Expr(k == 1 ? :if : :elseif,
                      :(k == $k),
                      :(($(els[1:k-1]...), $(els[k+1:end]...))),
                      ifexpr)
    end
    return quote
        delete(tup::NTuple{$N,T}, k::K) where {T,K<:Integer} = $ifexpr
    end |> esc
end
# Now generate functions for up to 64:
(@_tuple_delete_gencode  1)
(@_tuple_delete_gencode  2)
(@_tuple_delete_gencode  3)
(@_tuple_delete_gencode  4)
(@_tuple_delete_gencode  5)
(@_tuple_delete_gencode  6)
(@_tuple_delete_gencode  7)
(@_tuple_delete_gencode  8)
(@_tuple_delete_gencode  9)
(@_tuple_delete_gencode 10)
(@_tuple_delete_gencode 11)
(@_tuple_delete_gencode 12)
(@_tuple_delete_gencode 13)
(@_tuple_delete_gencode 14)
(@_tuple_delete_gencode 15)
(@_tuple_delete_gencode 16)
(@_tuple_delete_gencode 17)
(@_tuple_delete_gencode 18)
(@_tuple_delete_gencode 19)
(@_tuple_delete_gencode 20)
(@_tuple_delete_gencode 21)
(@_tuple_delete_gencode 22)
(@_tuple_delete_gencode 23)
(@_tuple_delete_gencode 24)
(@_tuple_delete_gencode 25)
(@_tuple_delete_gencode 26)
(@_tuple_delete_gencode 27)
(@_tuple_delete_gencode 28)
(@_tuple_delete_gencode 29)
(@_tuple_delete_gencode 30)
(@_tuple_delete_gencode 31)
(@_tuple_delete_gencode 32)
(@_tuple_delete_gencode 33)
(@_tuple_delete_gencode 34)
(@_tuple_delete_gencode 35)
(@_tuple_delete_gencode 36)
(@_tuple_delete_gencode 37)
(@_tuple_delete_gencode 38)
(@_tuple_delete_gencode 39)
(@_tuple_delete_gencode 40)
(@_tuple_delete_gencode 41)
(@_tuple_delete_gencode 42)
(@_tuple_delete_gencode 43)
(@_tuple_delete_gencode 44)
(@_tuple_delete_gencode 45)
(@_tuple_delete_gencode 46)
(@_tuple_delete_gencode 47)
(@_tuple_delete_gencode 48)
(@_tuple_delete_gencode 49)
(@_tuple_delete_gencode 50)
(@_tuple_delete_gencode 51)
(@_tuple_delete_gencode 52)
(@_tuple_delete_gencode 53)
(@_tuple_delete_gencode 54)
(@_tuple_delete_gencode 55)
(@_tuple_delete_gencode 56)
(@_tuple_delete_gencode 57)
(@_tuple_delete_gencode 58)
(@_tuple_delete_gencode 59)
(@_tuple_delete_gencode 60)
(@_tuple_delete_gencode 61)
(@_tuple_delete_gencode 62)
(@_tuple_delete_gencode 63)
(@_tuple_delete_gencode 64)
#
# #getpair =====================================================================
let _private_symbol = gensym("private_getpair_c61d3fee7deae4d3d2237de9044247b2")
"""
    getpair(d, k)

If the key k is found in the dictionary d, yields the pair (k => d[k]); 
otherwise yields missing.
"""
global getpair(d::AbstractDict, k) = begin
    v = get(d, k, _private_symbol)
    return (v === _private_symbol) ? missing : (k => v)
end
end
 

################################################################################
# EquivRef and equiv-type collections.
"""
    EquivRef{T}
An EquivRef wraps an object and mimics its equivhash() and equiv() functions on
the hash() and isequal() functions.
"""
struct EquivRef{T} <: Ref{T}
    object::T
    function EquivRef{T}(u::T) where {T}
        return new(u)
    end
    function EquivRef{T}(u::EquivRef) where {T}
        return new(u.object)
    end
    function EquivRef{T}(u) where {T}
        return new(u)
    end
end
EquivRef(u::EquivRef{T}) where {T} = EquivRef{T}(u.object)
EquivRef(u::T) where {T} = EquivRef{T}(u)
Base.hash(r::EquivRef) = equivhash(r.object)
Base.isequal(r::EquivRef, s::EquivRef) = isequiv(r.object, s.object)
Base.:(==)(r::EquivRef, s::EquivRef) = isequiv(r.object, s.object)
Base.getindex(r::EquivRef{T}) where {T} = r.object
"""
    EquivSet{T}
The EquivSet type is like the Set and IdSet types, except that it works on the
isequiv() and equivhash() functions instead of the isequal() and hash()
functions or the === and objectid() functions.
"""
struct EquivSet{T} <: AbstractSet{T}
    set::Set{EquivRef{T}}
end
mutability(::Type{EquivSet}) = Mutable()
mutability(::Type{EquivSet{T}}) where {T} = Mutable()
EquivSet{T}() where {T} = EquivSet{T}(Set{EquivRef{T}}())
EquivSet() = EquivSet{Any}()
EquivSet{T}(arr::AbstractArray{S,1}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    return EquivSet{T}(Set{EquivRef{T}}(ers))
end
EquivSet{T}(arr::AbstractSet{S}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    return EquivSet{T}(Set{EquivRef{T}}(ers))
end
EquivSet(arr::AbstractArray{S,1}) where {T,S<:T} = begin
    let T = typejoin([typeof(x) for x in arr]...), ers = [EquivRef{T}(x) for x in arr]
        return EquivSet{T}(Set{EquivRef{T}}(ers))
    end
end
EquivSet(arr::AbstractSet{S}) where {T,S<:T} = let ers = [EquivRef{T}(x) for x in arr]
    let T = typejoin([typeof(x) for x in arr]...), ers = [EquivRef{T}(x) for x in arr]
        return EquivSet{T}(Set{EquivRef{T}}(ers))
    end
end
Base.length(s::EquivSet) = length(s.set)
Base.in(x::S, s::EquivSet{T}) where {T, S<:T} = in(EquivRef{T}(x), s.set)
Base.in(x::S, s::EquivSet{T}, f) where {T, S<:T} = in(EquivRef{T}(x), s.set, f)
Base.iterate(s::EquivSet{T}) where {T} = let u = iterate(s.set)
    (u === nothing) && return nothing
    return (u[1].object, u[2])
end
Base.iterate(s::EquivSet{T}, it) where {T} = let u = iterate(s.set, it)
    (u === nothing) && return nothing
    return (u[1].object, u[2])
end
Base.push!(s::EquivSet{T}, x::S) where {T, S<:T} = begin
    push!(s.set, EquivRef{T}(x))
    return s
end
Base.delete!(s::EquivSet{T}, x::S) where {T, S<:T} = begin
    delete!(s.set, EquivRef{T}(x))
    return s
end
equalfn(::Type{EquivSet}) = isequiv
equalfn(::Type{EquivSet{T}}) where {T} = isequiv
hashfn(::Type{EquivSet}) = equivhash
hashfn(::Type{EquivSet{T}}) where {T} = equivhash
# A utility function for turning a list of pairs/tuples into a list of pairs.
_to_pairs(kvs) = begin
    if length(kvs) == 0
        K = Any
        V = Any
        if Base.IteratorEltype(kvs) isa Base.HasEltype
            ET = Base.eltype(itr)
            if isa(ET, DataType)
                if ET <: Pair
                    K = ET.parameters[1]
                    V = ET.parameters[2]
                elseif ET <: Tuple && length(ET.parameters) == 2
                    K = ET.parameters[1]
                    V = ET.parameters[2]
                end
            end
        end
        return Pair{K,V}[]
    else
        ks = []
        vs = []
        for kv in kvs
            if kv isa Pair || (kv isa Tuple && length(kv) == 2)
                push!(ks, kv[1])
                push!(vs, kv[2])
            else
                msg = "EquivDict: arg must be iterator of tuples or pairs"
                throw(ArgumentError(msg))
            end
        end
        K = typejoin(map(typeof, ks)...)
        V = typejoin(map(typeof, vs)...)
        return Pair{K,V}[Pair{K,V}(k,v) for (k,v) in zip(ks,vs)]
    end
end
"""
    EquivDict{K,V}
The EquivDict type is like the Dict and IdDict types, except that it works on
the isequiv() and equivhash() functions instead of the isequal() and hash()
functions or the === and objectid() functions.
"""
struct EquivDict{K,V} <: AbstractDict{K,V}
    dict::Dict{EquivRef{K},V}
    function EquivDict{K,V}() where {K,V}
        return new(Dict{EquivRef{K},V}())
    end
    function EquivDict{K,V}(::Tuple{}) where {K,V}
        return new(Dict{EquivRef{K},V}())
    end
    function EquivDict{K,V}(d::EquivDict) where {K,V}
        return new(Dict{EquivRef{K},V}(d.dict))
    end
    function EquivDict{K,V}(d::AbstractDict{EquivRef{K},V}) where {K,V}
        return new(Dict{EquivRef{K},V}(d))
    end
    function EquivDict{K,V}(d::AbstractDict) where {K,V}
        return new(Dict{EquivRef{K},V}(d))
    end
    function EquivDict{K,V}(ps::Pair...) where {K,V}
        return reduce(push!, ps, init=new(Dict{EquivRef{K},V}()))
    end
    function EquivDict{K,V}(itr) where {K,V}
        return reduce(push!, itr, init=new(Dict{EquivRef{K},V}()))
    end
end
const EquivPair{K,V} = Pair{EquivRef{K},V} where {K,V}
const Equiv2Tuple{K,V} = Tuple{EquivRef{K},V} where {K,V}
# EquivDict is not immutable despite being a struct.
mutability(::Type{EquivDict}) = Mutable()
mutability(::Type{EquivDict{K,V}}) where {K,V} = Mutable()
# Constructors.
EquivDict() = EquivDict{Any,Any}()
EquivDict(kv::Tuple{}) = EquivDict{Any,Any}()
EquivDict(p::Pair{K,V}) where {K,V} = EquivDict{K,V}(p)
EquivDict(d::AbstractDict{EquivRef{K},V}) where {K,V} = EquivDict{K,V}(d)
EquivDict(d::AbstractDict{K,V}) where {K,V} = EquivDict{K,V}(d)
EquivDict(ps::Pair...) = begin
    ps = _to_pairs(ps)
    T = Base.eltype(ps)
    K = T.parameters[1]
    V = T.parameters[2]
    return EquivDict{K,V}(ps...)
end
EquivDict(itr) = begin
    ps = _to_pairs(itr)
    T = Base.eltype(ps)
    K = T.parameters[1]
    V = T.parameters[2]
    return EquivDict{K,V}(ps...)
end
# Base methods.
Base.length(s::EquivDict) = length(s.dict)
Base.get(x::EquivDict{K,V}, k::KK, df) where {K,V,KK} = get(x.dict, EquivRef{K}(k), df)
Base.in(x::Pair{KK,VV}, s::EquivDict{K,V}) where {K,V,KK,VV} = begin
    kv = Pair{EquivRef{K},VV}(EquivRef{K}(x.first), x.second)
    return in(kv, s.dict)
end
Base.in(x::Pair{KK,VV}, s::EquivDict{K,V}, f) where {K,V,KK,VV} = begin
    kv = Pair{EquivRef{K},VV}(EquivRef{K}(x.first), x.second)
    return in(kv, s.dict, f)
end
Base.iterate(s::EquivDict{K,V}) where {K,V} = let u = iterate(s.dict)
    (u === nothing) && return nothing
    p = u[1]
    return (Pair{K,V}(p[1].object, p[2]), u[2])
end
Base.iterate(s::EquivDict{K,V}, it) where {K,V} = let u = iterate(s.dict, it)
    (u === nothing) && return nothing
    p = u[1]
    return (Pair{K,V}(p[1].object, p[2]), u[2])
end
Base.push!(s::EquivDict{K,V}, kv::Pair{KK,VV}) where {K,V,KK,VV} = begin
    push!(s.dict, EquivPair{K,V}(EquivRef{K}(kv[1]), kv[2]))
    return s
end
Base.delete!(s::EquivDict{K,V}, k::KK) where {K,V,KK} = begin
    delete!(s.dict, EquivRef{K}(k))
    return s
end
equalfn(::Type{EquivDict}) = isequiv
equalfn(::Type{EquivDict{K,V}}) where {K,V} = isequiv
hashfn(::Type{EquivDict}) = equivhash
hashfn(::Type{EquivDict{K,V}}) where {K,V} = equivhash

# Given the new category of equivalence, we actually need to rewrite the isequal
# function for dictionaries and sets--this is because at least some versions of
# the Julia base use the assumption that there is either Dict or IdDict, and
# that is just not correct any more.
function _isequiv(l::DL, r::DR, MUTL, MUTR) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}}
    return (l === r)
end
function _isequiv(l::DL, r::DR, ::Immutable, ::Immutable) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}}
    # Identical dictionaries are always equal.
    (l === r) && return true
    # Equality functions must match.
    (equalfn(DL) === equalfn(DR)) || return false
    # Dictionaries must be the same length to be equal.
    (length(l) == length(r)) || return false
    # All pairs from one must be in the other.
    for pair in l
        in(pair, r, isequiv) || return false
    end
    return true
end
function isequiv(l::DL, r::DR) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}}
    return _isequiv(l, r, mutability(DL), mutability(DR))
end
function _isequiv(l::SL, r::SR, MUTL, MUTR) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}}
    return (l === r)
end
function _isequiv(l::SL, r::SR, ::Immutable, ::Immutable) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}}
    # Sets must have the same equality base to be equal.
    (equalfn(SL) === equalfn(SR)) || return false
    # Identical sets are always equal.
    (l === r) && return true
    # Sets must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all elements must be found in both.
    for el in l
        in(el, r) || return false
    end
    return true
end
function isequiv(l::SL, r::SR) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}}
    return _isequiv(l, r, mutability(SL), mutability(SR))
end
function Base.isequal(
    l::DL, r::DR
) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}
}
    # Identical dictionaries are always equal.
    (l === r) && return true
    # Equality functions must match.
    (equalfn(DL) === equalfn(DR)) || return false
    # Dictionaries must be the same length to be equal.
    (length(l) == length(r)) || return false
    # All pairs from one must be in the other.
    for pair in l
        in(pair, r, isequal) || return false
    end
    return true
end
function Base.isequal(
    l::SL, r::SR
) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}
}
    # Sets must have the same equality base to be equal.
    (equalfn(SL) === equalfn(SR)) || return false
    # Identical sets are always equal.
    (l === r) && return true
    # Sets must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all elements must be found in both.
    for el in l
        in(el, r) || return false
    end
    return true
end
function Base.:(==)(
    l::DL, r::DR
) where {
    KL,VL, DL <: AbstractDict{KL,VL},
    KR,VR, DR <: AbstractDict{KR,VR}
}
    # Dictionaries must have the same equality base to be equal.
    (equalfn(DL) === equalfn(DR)) || return false
    # Identical dictionaries are always equal.
    (l === r) && return true
    # Dictionaries must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all pairs must be found in both.
    anymissing = false
    for pair in l
        isin = in(pair, r)
        if ismissing(isin)
            anymissing = true
        else
            isin || return false
        end
    end
    return anymissing ? missing : true
end
function Base.:(==)(
    l::SL, r::SR
) where {
    TL, SL <: AbstractSet{TL},
    TR, SR <: AbstractSet{TR}
}
    # Sets must have the same equality base to be equal.
    (equalfn(SL) === equalfn(SR)) || return false
    # Identical sets are always equal.
    (l === r) && return true
    # Sets must be the same length to be equal.
    (length(l) == length(r)) || return false
    # And all elements must be found in both.
    for el in l
        in(el, r) || return false
    end
    return true
end
