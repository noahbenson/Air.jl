################################################################################
# form.jl
#
# `Form`: an all-purpose persistent metadata object, in the spirit of JSON.
#
# A `Form` is a hybrid of a sequence and a map. It holds a `seq` of values, keyed
# by position, and a `map` of values keyed by symbol — so it reads like a function
# argument list, with positional arguments first and keyword arguments after.
# `Form("widget", meta=Form(name="abc"))` is a form whose first positional value
# is the string "widget" and whose `:meta` key holds another form.
#
# Values are drawn from a closed set, `FormLeaf`, so that every form can be
# written as JSON and read back. A value that is not already a leaf is converted
# when the conversion is easy and unambiguous — a `Char` becomes a `String`,
# another `Integer` becomes an `Int64`, a narrower real becomes a `Float64` — and
# otherwise construction raises.
# That closedness is the point: a form is metadata that can always be serialized,
# so it cannot hold an arbitrary object.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

"""
    FormLeaf

The values a [`Form`](@ref) can hold directly: `String`, `Symbol`, `Bool`,
`Int64`, `Float64`, `ComplexF64` and `Nothing`. A `Form` holds either one of
these or another `Form`.

A `String` and a `Symbol` are both JSON strings, and the first character inside
the quotes tells them apart: `'` for a string and `:` for a symbol.
"""
const FormLeaf = Union{String,Symbol,Bool,Int64,Float64,ComplexF64,Nothing}

"""
    Form(seq...; keys...)

A persistent metadata object holding both a sequence and a map: positional values
in `seq` and symbol-keyed values in `map`, so it reads like a function argument
list. Values are drawn from [`FormLeaf`](@ref) or are other forms.

Values are converted to a leaf where the conversion is unambiguous — a `Char`
becomes a `String`, another `Integer` becomes an `Int64`, a narrower real becomes
a `Float64` — and construction raises otherwise, since a form must always be
serializable.

A `Vector`, `PVector` or transient argument becomes a form with only a `seq`, and
a `Dict`, `PDict` or transient argument becomes a form with only a `map`, so
nested data can be written as ordinary Julia collections.

"""
struct Form
    # The union is written out rather than named: `Form` is not yet defined at
    # this point, so an alias for `Union{Form,FormLeaf}` cannot be either. The
    # alias is defined just below, for use everywhere else.
    seq::PVector{Union{Form,FormLeaf}}
    map::PDict{Symbol,Union{Form,FormLeaf}}
    # This inner constructor exists to *suppress* the one the compiler would
    # otherwise generate. That default takes any two arguments and converts them
    # to the field types, so it would capture `Form(a, b)` — which must mean two
    # positional values — and fail with a conversion error instead. Defining any
    # inner constructor prevents the default from being generated at all, and
    # this one is the fast path for a caller who already has the collections.
    Form(seq::PVector{Union{Form,FormLeaf}}, map::PDict{Symbol,Union{Form,FormLeaf}}) =
        new(seq, map)
end

"""
    FormValue

What a [`Form`](@ref) holds: either a [`FormLeaf`](@ref) or another `Form`.
"""
const FormValue = Union{Form,FormLeaf}

# ---- converting a value into the closed set ----------------------------------

# The error message names the type, because the usual cause is a value the caller
# assumed was acceptable — a `Set` where a leaf is wanted, or a dictionary with
# non-symbol keys.
_formerror(x, why) =
    throw(ArgumentError("a Form cannot hold $(typeof(x)): $why; " *
                        "see FormLeaf for what it can hold"))

# `Bool` before `Integer`, since `Bool <: Integer` and `true` must stay `true`.
_formvalue(x::Form) = x
_formvalue(x::FormLeaf) = x
_formvalue(x::Bool) = x
# The conversion is guarded so that an integer too large for an `Int64` reports
# itself the way every other unacceptable value does, rather than as a bare
# `InexactError` that says nothing about forms. A `try` costs nothing when it does
# not throw.
function _formvalue(x::Integer)
    return try
        Int64(x)
    catch
        _formerror(x, "it does not fit in an Int64, which is the only integer a Form holds")
    end
end
# A narrower real widens, which is the conversion that cannot lose anything.
_formvalue(x::Real) = Float64(x)
_formvalue(x::Complex) = ComplexF64(x)
_formvalue(x::Char) = string(x)
_formvalue(x::AbstractString) = String(x)
function _formvalue(x::AbstractVector)
    return Form(PVector(FormValue[_formvalue(v) for v in x]), PDict{Symbol,FormValue}())
end
function _formvalue(x::AbstractDict)
    pairs = Pair{Symbol,FormValue}[_formkey(k) => _formvalue(v) for (k, v) in x]
    return Form(PVector(FormValue[]), PDict{Symbol,FormValue}(pairs...))
end
_formvalue(x) = _formerror(x, "it is not a value type this can represent")

# A dictionary key must become a symbol for the map. Only the conversions that
# cannot lose information are made.
_formkey(k::Symbol) = k
_formkey(k::AbstractString) = Symbol(k)
_formkey(k::Char) = Symbol(string(k))
_formkey(k) = _formerror(k, "a map key must be a Symbol, a String or a Char")

# ---- construction ------------------------------------------------------------

function Form(args...; kwargs...)
    seq = PVector(FormValue[_formvalue(a) for a in args])
    map = PDict{Symbol,FormValue}(Pair{Symbol,FormValue}[k => _formvalue(v) for (k, v) in kwargs]...)
    return Form(seq, map)
end

# ---- reading -----------------------------------------------------------------

Base.length(f::Form) = length(f.seq) + length(f.map)

Base.getindex(f::Form, i::Integer) = f.seq[i]
Base.getindex(f::Form, k::Symbol) = f.map[k]

Base.get(f::Form, i::Integer, default) = get(f.seq, i, default)
Base.get(f::Form, k::Symbol, default) = get(f.map, k, default)

Base.haskey(f::Form, i::Integer) = i in eachindex(f.seq)
Base.haskey(f::Form, k::Symbol) = haskey(f.map, k)

# ---- updating ----------------------------------------------------------------

setindex(f::Form, v, i::Integer) = Form(setindex(f.seq, _formvalue(v), i), f.map)
setindex(f::Form, v, k::Symbol) = Form(f.seq, push(f.map, k => _formvalue(v)))

push(f::Form, v) = Form(push(f.seq, _formvalue(v)), f.map)

function pop(f::Form)
    isempty(f.seq) && throw(ArgumentError("a Form with an empty seq has nothing to pop"))
    return Form(pop(f.seq), f.map)
end

# `deleteat`, not `delete`: the latter removes a *value* from a vector, and here
# the argument is a position.
delete(f::Form, i::Integer) = Form(deleteat(f.seq, i), f.map)
delete(f::Form, k::Symbol) = Form(f.seq, delete(f.map, k))

# ---- equality ----------------------------------------------------------------

Base.isequal(a::Form, b::Form) = isequal(a.seq, b.seq) && isequal(a.map, b.map)
Base.:(==)(a::Form, b::Form) = isequal(a, b)
Base.hash(f::Form, h::UInt) = hash(f.map, hash(f.seq, hash(:Form, h)))

# ---- printing ----------------------------------------------------------------

# A nested form prints without the `Form` prefix, so that a form reads like a
# function call with its arguments: `Form("widget", meta=(name="abc"))`.
function _showvalue(io::IO, f::Form)
    print(io, "(")
    _showcontents(io, f)
    return print(io, ")")
end
_showvalue(io::IO, x) = show(io, x)

function _showcontents(io::IO, f::Form)
    sep = ""
    for x in f.seq
        print(io, sep)
        sep = ", "
        _showvalue(io, x)
    end
    # The keys are sorted so that what is printed is a property of the form and
    # not of the map's iteration order, which is an implementation detail of the
    # HAMT. Note that this means keyword values do not print in the order they
    # were given: the map is a `PDict`, so it does not record one.
    for k in sort!(collect(keys(f.map)))
        print(io, sep, k, "=")
        sep = ", "
        _showvalue(io, f.map[k])
    end
    return nothing
end

function Base.show(io::IO, f::Form)
    print(io, "Form(")
    _showcontents(io, f)
    return print(io, ")")
end

export Form, FormLeaf, FormValue

# ==============================================================================
# JSON

# A form is always written as the two-element list `[seq, map]`, with the symbol
# keys of the map as JSON strings. That is what makes the format self-describing:
# a list is a form and nothing else, so `from_JSON` never has to guess.
#
# Two encodings are worth stating because they are not forced by JSON itself.
#
#  * A `String` and a `Symbol` are both JSON strings, so the first character
#    inside the quotes says which: `'` for a string and `:` for a symbol.
#  * A `Float64` is written as a plain number but always with a decimal point or
#    an exponent, since a bare integer would read back as an `Int64`.
#  * A complex number has no JSON equivalent, so it is written as the object
#    `{"re": …, "im": …}` — always, even when the imaginary part is zero, since
#    `ComplexF64(1, 0)` and `Float64(1)` are different values.

const _JSON_ESCAPE = Dict{UInt8,String}(
    UInt8('"') => "\\\"", UInt8('\\') => "\\\\", UInt8('\b') => "\\b",
    UInt8('\f') => "\\f", UInt8('\n') => "\\n", UInt8('\r') => "\\r",
    UInt8('\t') => "\\t",
)

# `marker` is the first character inside the quotes, and is what says whether the
# value is a string or a symbol: JSON has only one string type, and a form holds
# both.
function _jsonstring(io::IO, s::AbstractString, marker::Char)
    print(io, '"', marker)
    for c in s
        u = UInt32(c)
        if u < 0x20
            haskey(_JSON_ESCAPE, UInt8(u)) ? print(io, _JSON_ESCAPE[UInt8(u)]) :
                print(io, "\\u", string(u, base = 16, pad = 4))
        elseif c == '"' || c == '\\'
            print(io, '\\', c)
        else
            print(io, c)
        end
    end
    return print(io, '"')
end

_jsonvalue(io::IO, x::Nothing) = print(io, "null")
_jsonvalue(io::IO, x::Bool) = print(io, x ? "true" : "false")
_jsonvalue(io::IO, x::Int64) = print(io, x)
_jsonvalue(io::IO, x::AbstractString) = _jsonstring(io, x, '\'')
_jsonvalue(io::IO, x::Symbol) = _jsonstring(io, string(x), ':')
# `print`, not `show`: `show(1.5)` writes `1.5` but `show(1.5f0)` would write
# `1.5f0`, which is not JSON. A float is always written with a point or an
# exponent, so it never reads back as an integer.
_jsonvalue(io::IO, x::Float64) = print(io, x)
# Always the object, even when the imaginary part is zero: `ComplexF64(1, 0)` and
# `Float64(1)` are different values, and a plain number would read back as the
# latter.
function _jsonvalue(io::IO, x::ComplexF64)
    print(io, "{\"re\": ", real(x), ", \"im\": ", imag(x), "}")
    return nothing
end

function _jsonvalue(io::IO, f::Form)
    # Two brackets: the form's own, then the seq array's.
    print(io, "[[")
    sep = ""
    for x in f.seq
        print(io, sep)
        sep = ", "
        _jsonvalue(io, x)
    end
    print(io, "], {")
    sep = ""
    for k in sort!(collect(keys(f.map)))
        print(io, sep)
        sep = ", "
        _jsonvalue(io, k)
        print(io, ": ")
        _jsonvalue(io, f.map[k])
    end
    return print(io, "}]")
end

"""
    to_JSON(f::Form) -> String

The form as a JSON string. A form is always written as the two-element list
`[seq, map]`, so the format is self-describing: `from_JSON` reads a list as a
form and nothing else.

A `String` and a `Symbol` are both written as JSON strings, with the first
character inside the quotes saying which: `'` for a string and `:` for a symbol.
A `Float64` is a plain number, always with a decimal point or an exponent so that
it reads back as a float rather than as an integer. A complex number has no JSON
equivalent and is written as `{"re": …, "im": …}`, always, even when its
imaginary part is zero.

"""
function to_JSON(f::Form)
    io = IOBuffer()
    _jsonvalue(io, f)
    return String(take!(io))
end

# ---- parsing -----------------------------------------------------------------

# A hand-written recursive-descent reader for the JSON subset, rather than a
# dependency: the values a form can hold are a closed set, so the reader has
# nothing to do that a general one would.
mutable struct _JSONReader
    s::String
    i::Int
end
_JSONReader(s::AbstractString) = _JSONReader(String(s), 1)

_jsonfail(r::_JSONReader, msg) =
    throw(ArgumentError("invalid JSON at byte $(r.i): $msg"))
_atend(r::_JSONReader) = r.i > ncodeunits(r.s)

function _skipws!(r::_JSONReader)
    while !_atend(r) && (r.s[r.i] == ' ' || r.s[r.i] == '\t' ||
                         r.s[r.i] == '\n' || r.s[r.i] == '\r')
        r.i = nextind(r.s, r.i)
    end
    return nothing
end

function _expect!(r::_JSONReader, c::Char)
    _skipws!(r)
    (_atend(r) || r.s[r.i] != c) && _jsonfail(r, "expected '$c'")
    r.i = nextind(r.s, r.i)
    return nothing
end

function _literal!(r::_JSONReader, word::String, value)
    _skipws!(r)
    startswith(SubString(r.s, r.i), word) || _jsonfail(r, "expected $word")
    for _ in 1:length(word)
        r.i = nextind(r.s, r.i)
    end
    return value
end

function _parsestring!(r::_JSONReader)
    _expect!(r, '"')
    io = IOBuffer()
    while true
        _atend(r) && _jsonfail(r, "unterminated string")
        c = r.s[r.i]
        if c == '"'
            r.i = nextind(r.s, r.i)
            return String(take!(io))
        elseif c != '\\'
            (UInt32(c) < 0x20) && _jsonfail(r, "unescaped control character")
            print(io, c)
            r.i = nextind(r.s, r.i)
        else
            r.i = nextind(r.s, r.i)
            _atend(r) && _jsonfail(r, "unterminated escape")
            e = r.s[r.i]
            r.i = nextind(r.s, r.i)
            if e == 'u'
                u = _hex4!(r)
                # A surrogate pair is one character, encoded as two escapes.
                if 0xd800 <= u <= 0xdbff
                    (startswith(SubString(r.s, r.i), "\\u")) ||
                        _jsonfail(r, "unpaired surrogate")
                    r.i = nextind(r.s, r.i)   # the backslash
                    r.i = nextind(r.s, r.i)   # the u
                    v = _hex4!(r)
                    (0xdc00 <= v <= 0xdfff) || _jsonfail(r, "unpaired surrogate")
                    u = 0x10000 + ((u - 0xd800) << 10) + (v - 0xdc00)
                end
                print(io, Char(u))
            elseif e == 'n'
                print(io, '\n')
            elseif e == 't'
                print(io, '\t')
            elseif e == 'r'
                print(io, '\r')
            elseif e == 'b'
                print(io, '\b')
            elseif e == 'f'
                print(io, '\f')
            elseif e == '"' || e == '\\' || e == '/'
                print(io, e)
            else
                _jsonfail(r, "unknown escape \\$e")
            end
        end
    end
end

function _hex4!(r::_JSONReader)
    (r.i + 3 <= ncodeunits(r.s)) || _jsonfail(r, "truncated \\u escape")
    hex = SubString(r.s, r.i, r.i + 3)
    all(c -> c in '0':'9' || c in 'a':'f' || c in 'A':'F', hex) ||
        _jsonfail(r, "bad \\u escape")
    r.i += 4
    return parse(UInt32, hex; base = 16)
end

# An integer literal reads back as an `Int64`; anything with a point or an
# exponent reads back as the `Float64` the writer would have produced for one.
# That is what makes the encoding above round-trip.
function _parsenumber!(r::_JSONReader)
    _skipws!(r)
    start = r.i
    while !_atend(r) && (isdigit(r.s[r.i]) || r.s[r.i] in ('-', '+', '.', 'e', 'E'))
        r.i = nextind(r.s, r.i)
    end
    (r.i == start) && _jsonfail(r, "expected a number")
    text = SubString(r.s, start, prevind(r.s, r.i))
    if !occursin('.', text) && !occursin('e', text) && !occursin('E', text)
        v = tryparse(Int64, text)
        (v === nothing) && _jsonfail(r, "integer out of range: $text")
        return v
    end
    v = tryparse(Float64, text)
    (v === nothing) && _jsonfail(r, "bad number: $text")
    return v
end

# A JSON string carries its own kind in its first character: `'` for a string and
# `:` for a symbol. Anything else is a mistake rather than a value, since the
# writer never produces one.
function _parseleafstring!(r::_JSONReader)
    text = _parsestring!(r)
    isempty(text) &&
        _jsonfail(r, "a JSON string must begin with ' (a string) or : (a symbol)")
    marker = text[1]
    rest = SubString(text, nextind(text, 1))
    (marker == '\'') && return String(rest)
    (marker == ':') && return Symbol(rest)
    _jsonfail(r, "a JSON string must begin with ' (a string) or : (a symbol), " *
                 "not '$marker'")
end

# A map key is always a symbol, so it is written with the `:` marker. The marker
# is required rather than assumed, since a key without one is a mistake the writer
# never makes — but it is checked here rather than while reading the object, so
# that the `{"re", "im"}` object a complex number is written as can carry plain
# keys without them meaning symbols.
function _parsekey(r::_JSONReader, k::AbstractString)
    isempty(k) && _jsonfail(r, "a map key cannot be empty")
    (k[1] == ':') || _jsonfail(r, "a map key must be a symbol, written \":name\"")
    return Symbol(SubString(k, nextind(k, 1)))
end

function _parsearray!(r::_JSONReader)
    _expect!(r, '[')
    out = FormValue[]
    _skipws!(r)
    if !_atend(r) && r.s[r.i] == ']'
        r.i = nextind(r.s, r.i)
        return out
    end
    while true
        push!(out, _parsevalue!(r))
        _skipws!(r)
        _atend(r) && _jsonfail(r, "unterminated array")
        c = r.s[r.i]
        r.i = nextind(r.s, r.i)
        (c == ']') && return out
        (c == ',') || _jsonfail(r, "expected ',' or ']'")
    end
end

function _parseobject!(r::_JSONReader)
    _expect!(r, '{')
    out = Pair{String,FormValue}[]
    _skipws!(r)
    if !_atend(r) && r.s[r.i] == '}'
        r.i = nextind(r.s, r.i)
        return out
    end
    while true
        _skipws!(r)
        k = _parsestring!(r)
        _expect!(r, ':')
        push!(out, k => _parsevalue!(r))
        _skipws!(r)
        _atend(r) && _jsonfail(r, "unterminated object")
        c = r.s[r.i]
        r.i = nextind(r.s, r.i)
        (c == '}') && return out
        (c == ',') || _jsonfail(r, "expected ',' or '}'")
    end
end

# A list is a form and nothing else, so it is read as exactly the two elements
# the encoding always writes: a sequence array and a map object. The two are read
# directly rather than through `_parsevalue!`, since routing the seq back through
# it would try to read the seq itself as a form.
function _parseform!(r::_JSONReader)
    _expect!(r, '[')
    seq = _parsearray!(r)
    _expect!(r, ',')
    pairs = _parseobject!(r)
    _expect!(r, ']')
    kv = Pair{Symbol,FormValue}[_parsekey(r, first(p)) => last(p) for p in pairs]
    return Form(PVector(FormValue[seq...]), PDict{Symbol,FormValue}(kv...))
end

function _parsevalue!(r::_JSONReader)
    _skipws!(r)
    _atend(r) && _jsonfail(r, "unexpected end of input")
    c = r.s[r.i]
    (c == 'n') && return _literal!(r, "null", nothing)
    (c == 't') && return _literal!(r, "true", true)
    (c == 'f') && return _literal!(r, "false", false)
    (c == '"') && return _parseleafstring!(r)
    (c == '[') && return _parseform!(r)
    if c == '{'
        pairs = _parseobject!(r)
        # `{"re": …, "im": …}` is how a complex number is written; any other
        # object is a form with only a map.
        if length(pairs) == 2 && Set(first.(pairs)) == Set(["re", "im"])
            re = pairs[findfirst(p -> first(p) == "re", pairs)].second
            im = pairs[findfirst(p -> first(p) == "im", pairs)].second
            (re isa Float64 && im isa Float64) ||
                _jsonfail(r, "re and im must be numbers")
            return ComplexF64(re, im)
        end
        kv = Pair{Symbol,FormValue}[_parsekey(r, first(p)) => last(p) for p in pairs]
        return Form(PVector(FormValue[]), PDict{Symbol,FormValue}(kv...))
    end
    return _parsenumber!(r)
end

"""
    from_JSON(s::AbstractString) -> Form

Reads a [`Form`](@ref) from the JSON that [`to_JSON`](@ref) writes. A form is
always a two-element list `[seq, map]`, and that is the only shape accepted as
one; a bare object is accepted as a form with only a map, since that is a
convenient way to write one by hand.

"""
function from_JSON(s::AbstractString)
    r = _JSONReader(s)
    v = _parsevalue!(r)
    _skipws!(r)
    _atend(r) || _jsonfail(r, "trailing data after the form")
    v isa Form || _jsonfail(r, "the JSON must describe a Form, not a $(typeof(v))")
    return v
end

export to_JSON, from_JSON
