```@meta
DocTestSetup = quote
    using Air
end
```

# Persistent Forms

A [`Form`](@ref) is a persistent metadata object: a hybrid of a sequence and a
map, so that one container can carry both positional and keyword values. It reads
like a function argument list, and prints like one.

```julia
julia> f = Form("widget", 3, meta=Form(name="abc"), hidden=true);

julia> f[1]                      # by position
"widget"

julia> f[:meta]                  # or by symbol
Form(name="abc")

julia> f
Form("widget", 3, hidden=true, meta=(name="abc"))
```

Note that a nested form prints without the `Form` prefix — `meta=(name="abc")`
rather than `meta=Form(name="abc")` — so that the whole reads as one call with its
arguments, however deep it goes.

## A closed set of values

A form holds a [`FormLeaf`](@ref) — a `String`, `Symbol`, `Bool`, `Int64`,
`Float64`, `ComplexF64` or `Nothing` — or another form, and nothing else. That closedness is the point: a
form is metadata that can always be written out, so it cannot come to hold an
object that has no serialized form.

Values that are not already leaves are converted when the conversion is
unambiguous, and refused otherwise:

```julia
julia> Form(:sym)[1]             # a Symbol is a leaf of its own
:sym

julia> Form(Int32(7))[1]         # a narrower Integer -> Int64
7

julia> Form(1.5f0)[1]            # a narrower real -> Float64
1.5

julia> Form(Set([1]))
ERROR: ArgumentError: a Form cannot hold Set{Int64}: it is not a value type this can represent; see FormLeaf for what it can hold
```

## Building nested data from collections

A `Vector` argument becomes a form with only a `seq`, and a `Dict` argument one
with only a `map`, so nested data can be written as ordinary Julia collections:

```julia
julia> Form([1, 2, 3])
Form((1, 2, 3))

julia> Form(Dict(:a => 1))
Form((a=1))
```

Two `Vector` arguments are two positional values, not a sequence and a map.

## JSON

A form is always written as the two-element list `[seq, map]`, with the symbol
keys of the map as JSON strings. That is what makes the format self-describing: a
list is a form and nothing else, so [`from_JSON`](@ref) never has to guess.

```julia
julia> to_JSON(Form("widget", meta=Form(name="abc")))
"[[\"'widget\"], {\":meta\": [[], {\":name\": \"'abc\"}]}]"

julia> from_JSON("[[\"'widget\"], {}]")
Form("widget")
```

Three encodings are worth knowing because JSON does not force them:

- A `String` and a `Symbol` are both JSON strings, so the first character inside
  the quotes says which: **`'` for a string and `:` for a symbol**. A map key is
  always a symbol, so it carries the `:` marker too.
- A `Float64` is a plain number, but **always with a decimal point or an
  exponent**, since a bare integer would read back as an `Int64`.
- A complex number has no JSON equivalent and is written as the object
  `{"re": …, "im": …}` — always, even when its imaginary part is zero, since
  `ComplexF64(1, 0)` and `Float64(1)` are different values.

Together these make `from_JSON(to_JSON(f)) == f` hold for every value the type
admits.

`from_JSON` accepts a bare JSON object as a form with only a map, since that is a
convenient way to write one by hand, but it refuses a bare array: a form is
always `[seq, map]`.

## Reference

The full docstring for each is in the API reference page, along with the rest of
Air's public surface: [`Form`](@ref), [`FormLeaf`](@ref), [`FormValue`](@ref),
[`to_JSON`](@ref) and [`from_JSON`](@ref).
