################################################################################
# Meta Data
# This file is currently experimental and under development.

################################################################################
# SimpleDict
# Simple dictionaries are dictionaries whose keys are Symbols. Simple dicts can
# be used in argument processing if their keys are valid tokens.
const AbstractSimpleDict{T} = AbstractDict{Symbol, T} where {T}
const SimpleDict{T} = Dict{Symbol, T} where {T}
const PSimpleDict{T} = PDict{Symbol, T} where {T}
SimpleDict() = SimpleDict{Any}()
PSimpleDict() = PSimpleDict{Any}()
# How we parse arguments in simpleflat.
_simplebuild!(d::SimpleDict{T}, k::Symbol, v::S, args...) where {T,S} = begin
    d[k] = v
    return args
end
_simplebuild!(d::SimpleDict{T}, p::Pair{Symbol,S}, args...) where {T,S} = begin
    d[p[1]] = p[2]
    return args
end
_simplebuild!(d::SimpleDict{T}, p::Tuple{Symbol,S}) where {T,S} = begin
    d[p[1]] = p[2]
    return args
end
_simplebuild!(d::SimpleDict{T}, p::AbstractSimpleDict{S}, args...) where {T,S} = begin
    for (k,v) in p
        d[k] = v
    end
    return args
end
_simplebuild!(d::SimpleDict{T}, a0, args...) where {T} =
    error("invalid simpleflat argument type: $(typeof(a0))")
"""
    simpleflat(args...)
    simpleflat(T, args...)

Yields a simple dictionary that results from flattening the arguments. This
flattening occurs as follows:
 1. First, the arguments are parsed. The arguments may include simple
    dictionaries or paired key-value arguments (e.g., the following call is
    valid: `simpleflat(sym1, val1, dict, sym2, val2)`. These arguments are
    collapsed right-to-left. Additionally, Pair{Symbol,_} or Tuple{Symbol,_}
    objects may be passed in.
 2. Next, the keyword values are merged into the result.

The return type is always SimpleDict{Any} unless the type argument T is
supplied as the first argument, in which case, it is of type SimpleDict{T}.
If the supplied type T is Union{}, then the type used is the result of
typejoining across all the values.
"""
simpleflat(::Type{T0}, args...; kws...) where {T0} = begin
    tjoin = false
    T = T0 === Union{} ? Any : T0
    # Create the buffer dictionary into which we gather these.
    if length(args) > 0 && isa(args[1], AbstractSimpleDict)
        d = SimpleDict{T}(args[1])
        args = args[2:end]
    else
        d = SimpleDict{T}()
    end
    # Parse each argument using the above _simplebuild!() methods.
    aa = args
    while !isempty(aa)
        aa = _simplebuild!(d, aa...)
    end
    # Merge the keywords in last.
    for (k,v) in kws
        d[k] = v
    end
    # Check if we still need to join types.
    if T0 === Union{}
        T = typejoin([typeof(v) for v in values(d)])
        (T === Any) || return SimpleDict{T}(d)
    end
    return d
end
simpleflat(args...; kw...) = simpleflat(Union{}, args...; kw...)
simpleflat() = SimpleDict{Any}()


################################################################################
# MetaData
# Now that we've defined PDict, we can also fill out the meta-data concept for
# Air. Meta-data is an optional feature of a type that can be tracked and used
# Details for handling meta-data and declaring types that carry meta-data.
const AbstractMetaData = AbstractSimpleDict
const empty_metadata = PDict{Symbol,Any}()
# A trait for tracking if an object has meta-data or not.
"""
    MetaDataTrait

The MetaDataTrait allows types to indicate at compile time that their objects
support a meta-data dictionary. Meta-data dictionaries are Dict or PDict objects
that map symbols to anything.

A type T that supports meta-data will overload the following methods:
 * metatrait(T) yields a trait-object of type WithMetaData <: MetaDataTrait
 * metadata(::WithMetaData, t::T) yields the meta-data dictionary of t
 * Optionally, one of:
   * setmetadata!(::WithMetaData, t::T, dict), which changes the meta-data
     dictionary of to the dictionary dict, or
   * withmetadata(::WithMetaData, t::T, dict), which yields a duplicate of t
     with the new meta-data dictionary dict.
"""
abstract type MetaDataTrait end
struct WithMetaData <: MetaDataTrait end
struct WoutMetaData <: MetaDataTrait end
# #metatrait
"""
    metatrait(T)

Yields the meta-data trait (WithMetaData() or WoutMetaData(), both of abstract
type MetaDataTrait) of the type T. Types that yield WithMetaData() track a
dictionary whose type T <: AbstractMetaData.

Any type T that declares itself a meta-data type by overloading metatrait(T) to
return WithMetaData() must also implement the metadata() function, which yields
the meta-data dictionary of the object passed to it.
"""
metatrait(::Type{T}) where {T} = WoutMetaData() # Default: No metadata tracked.
# #metadata
"""
    metadata(t)

Yields the meta-data dictionary of the given object t. If t does not implement
meta-data (see ?MetaDataTrait) then an empty persistent dictionary is returned.
"""
metadata(t::T) where {T} = metadata(metatrait(T), t)
metadata(::WoutMetaData, t::T) where {T} = empty_metadata
# #meta
"""
    meta(object, key, default)
    meta(object, key)

Yields the value associated with the given meta-data key in the given object.
If the key is not found in the object's meta-data (or if the object does not
track meta-data) yields default. If default is not provided, then an error is
raised when it would be returned.
"""
meta(t::T, k::Symbol, dflt) where {T} = get(metadata(metatrait(T), t), k, dflt)
meta(t::T, k::Symbol) where {T} = begin
    md = metadata(metatrait(T), t)
    v = get(md, k, nothing)
    (v === nothing) && !in(k, md) && throw(KeyError(k))
    return v
end
# #setmetadata!
"""
    setmetadata!(object, dict)

Sets the meta-data dictionary for the giveb object to be equal to the given
dictionary dict. If the object is not an object that tracks its meta-data, then
an exception is thrown.

On success, the new value of the meta-data dictionary is returned.
"""
setmetadata!(t::T, d::D) where {T, D <: AbstractMetaData} =
    setmetadata!(metatrait(T), t, d)
setmetadata!(::WoutMetaData, t::T, d::D) where {T, D <: AbstractMetaData} =
    error("type $(T) does not support meta-data")
setmetadata!(::WithMetaData, t::T, d::D) where {T, D <: AbstractMetaData} =
    error("type $(T) does not support setting meta-data of type $(D)")
# #withmetadata
"""
    withmetadata(object, dict)

Yields a duplicate of the given object except that its meta-data dictionary
will have been updated to be equal to dict.
"""
withmetadata(t::T, d::D) where {T, D <: AbstractMetaData} =
    withmetadata(metatrait(T), t, d)
withmetadata(::WoutMetaData, t::T, d::D) where {T, D <: AbstractMetaData} =
    error("type $(T) does not support meta-data")
withmetadata(::WithMetaData, t::T, d::D) where {T, D <: AbstractMetaData} =
    error("type $(T) does not support clone-updating meta-data of type $(D)")
# #setmeta!
"""
    withmeta(obj, k, v)
    withmeta(obj, k1, v1, k2, v2, ...)
    withmeta(obj, <key>=v)
    withmeta(obj, <key1>=v1, <key2>=v2, ...)
    withmeta(obj, k1a, v1a, k2a, v2a, ..., <key1b>=v1b, <key2b>=v2b, ...)
    withmeta(T, obj, args...)

Yields a duplicate of the given object in which the given key or keys have been
updated to be associated with new values. In order for this function to work,
the object obj must support the `withmetadata()` method and it must have a
dictionary type that wortks with Air's `setindex()` and `delete()` methods,
otherwise behavior is undefined (likely an exception will be raised).

Note that matched `k, v` argument pairs may be replaced with or interspersed
with arguments such as `k => v` and `(k,v)` so long as the associated keys
are symbols.
"""
withmeta(t::T, args...; kw...) where {T} =
    withmeta(metatrait(T), t, args...; kw...)
withmeta(::WoutMetaData, t::T, args...; kw...) where {T} =
    error("type $(T) does not support meta-data")
_buildmeta(md0::AbstractSimpleDict{V}, args...; kw...) where {V} =
    simpleflat(V, md0, args...; kw...)
withmeta(Q::WithMetaData, t::T, args...; kw...) where {T} = begin
    # Get the meta-data that we start with
    md0 = metadata(Q, t)
    # Process the arguments
    md = _buildmeta(md0, args...; kw...)
    # Make sure the type is correct
    MD0 = typeof(md0)
    (MD0 === typeof(md)) || (md = MD0(md))
    # And return the updated version.
    return withmetadata(Q, t, md)
end
"""
    woutmeta(obj, k)
    woutmeta(obj, k1, k2...)

Yields a duplicate of the given object in which the given meta-data keys have
been removed. If obj does not support meta-data or does not support
clone-updating the meta-data, an exception will be thrown.
"""
woutmeta(t::T, args...) where {T} = woutmeta(metatrait(T), t, args...)
woutmeta(::WoutMetaData, t::T, args...) where {T} =
    error("type $(T) does not support meta-data")
_trimmeta(md0::AbstractSimpleDict{V}, args...) where {V} = begin
    for arg in args
        md0 = delete(md0, arg)
    end
    return md0
end
woutmeta(Q::WithMetaData, t::T, args...) where {T} = begin
    # Get the meta-data that we start with
    md0 = metadata(Q, t)
    # Process the arguments
    md = _trimmeta(md0, args...)
    # Make sure the type is correct
    MD0 = typeof(md0)
    (MD0 === typeof(md)) || (md = MD0(md))
    # And return the updated version.
    return withmetadata(Q, t, md)
end
