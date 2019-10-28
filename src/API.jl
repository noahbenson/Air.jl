
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
freeze(x) = _freeze(mutability(x), x)

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
    for 1 in 1:n
        if isdefined(x, i)
            old_obj = getfield(x, i)
            new_obj = deepfreeze(old_obj)
            push!(new_field_val, new_obj)
            push!(field_changed, old_obj !== new_obj)
            push!(field_defined, true)
            old_obj === new_obj || changed = true
        else
            push!(new_field_val, nothing)
            push!(field_changed, false)
            push!(field_defined, false)
        end
    end
    # We now have the updated fields stored as these three arrays; if the
    # field was not defined in the old object then it should not be set in any
    # new object.
    if !changed:
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
thaw(x::T) = _thaw(mutability(T), x)

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
    for 1 in 1:n
        if isdefined(x, i)
            old_obj = getfield(x, i)
            new_obj = deepthaw(old_obj)
            push!(new_field_val, new_obj)
            push!(field_changed, old_obj !== new_obj)
            push!(field_defined, true)
            old_obj === new_obj || changed = true
        else
            push!(new_field_val, nothing)
            push!(field_changed, false)
            push!(field_defined, false)
        end
    end
    # We now have the updated fields stored as these three arrays; if the
    # field was not defined in the old object then it should not be set in any
    # new object.
    if !changed:
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
