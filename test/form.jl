################################################################################
# form.jl
#
# Tests for `Form`, the sequence-and-map metadata object.
#
# Three properties are worth more than the rest and are tested first, because
# everything else follows from them:
#
#  * a form holds only the closed value set, so it can always be serialized —
#    which means construction must convert what it can and *refuse* the rest
#    rather than storing something that would fail later;
#  * a nested form prints without the `Form` prefix, so a form reads like a
#    function call with its arguments;
#  * `from_JSON(to_JSON(f)) == f` for every value the type admits, which is what
#    the JSON encoding exists for and what pins the two encodings that JSON does
#    not force — a real number written with a decimal point so that it reads back
#    as the complex it is stored as, and a complex written as `{"re","im"}`.
#
# @author Noah C. Benson
#
# MIT License
# Copyright (c) 2020-2026 Noah C. Benson

@testset "Form" begin
    @testset "construction and the closed value set" begin
        @test Form() isa Form
        @test Form() == Form()
        @test length(Form()) == 0

        f = Form("widget", 3, meta=Form(name="abc"), hidden=true)
        @test f isa Form
        @test length(f) == 4
        @test f[1] == "widget"
        @test f[2] == 3
        @test f[:meta] == Form(name="abc")
        @test f[:hidden] === true

        # A form holds only `FormLeaf` and other forms, so the field types are the
        # guarantee: an element is either a leaf or a form, never anything else.
        @test eltype(f.seq) == FormValue
        @test eltype(f.map) == Pair{Symbol,FormValue}
        for x in f.seq
            @test x isa Form || x isa FormLeaf
        end

        # The conversions that cannot lose information are made.
        @test Form(:sym)[1] === :sym           # a Symbol is a leaf of its own
        @test Form('c')[1] === "c"             # Char -> String
        @test Form(SubString("sub"))[1] === "sub"
        @test Form(Int32(7))[1] === Int64(7)   # narrower Integer -> Int64
        @test Form(UInt8(3))[1] === Int64(3)
        @test Form(true)[1] === true           # Bool stays Bool, not Int64
        @test Form(false)[1] === false
        @test Form(nothing)[1] === nothing
        @test Form(1.5)[1] === 1.5             # a real is a Float64
        @test Form(1.5f0)[1] === 1.5           # and a narrower one widens
        @test Form(2//4)[1] === 0.5
        @test Form(ComplexF64(1, 2))[1] === ComplexF64(1, 2)
        @test Form(ComplexF32(1, 2))[1] === ComplexF64(1, 2)   # widened
        @test Form(1.0 + 2im)[1] === ComplexF64(1, 2)
        # A real is a Float64 and *not* a complex with a zero imaginary part:
        # they are different leaves, and JSON keeps them apart.
        @test Form(1.0)[1] === 1.0
        @test Form(1.0)[1] isa Float64
        @test Form(1.0 + 0im)[1] isa ComplexF64

        # A form is never mutated: every update returns a new one.
        g = Form("a", x=1)
        @test push(g, "b") == Form("a", "b", x=1)
        @test g == Form("a", x=1)
    end

    @testset "what a form will not hold" begin
        # The error names the type, because the usual cause is a value the caller
        # assumed was acceptable.
        @test_throws ArgumentError Form((a=1,))              # a NamedTuple
        @test_throws ArgumentError Form(Dict(:a => 1), Set([1]))  # a Set
        @test_throws ArgumentError Form("x", y=Set([1]))
        @test_throws ArgumentError Form(big(2)^100)          # will not fit an Int64
        @test_throws ArgumentError Form(:a => 1)             # a bare Pair
        # ... but the message should say which value was refused.
        err = try
            Form(Set([1]))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Set", sprint(showerror, err))
        @test occursin("FormLeaf", sprint(showerror, err))
    end

    @testset "collections become nested forms" begin
        # A vector becomes a form with only a seq, and a dictionary one with only
        # a map, so nested data can be written as ordinary Julia collections.
        @test Form([1, 2, 3]) == Form(Form(1, 2, 3))
        @test Form(Dict(:a => 1, :b => 2)) == Form(Form(a=1, b=2))
        @test Form("x", [1, 2]) == Form("x", Form(1, 2))
        # `String` keys and `Int32` values convert, since both are unambiguous.
        @test Form(Dict("k" => Int32(1))) == Form(Form(k=1))
        # A `Vector{Any}` is fine as long as its elements convert.
        @test Form(Any[1, 2]) == Form(Form(1, 2))
        # A vector with a value that cannot convert is refused.
        @test_throws ArgumentError Form(Any[Set([1])])
        # The Air collections work the same way.
        @test Form(PVector([1, 2])) == Form(Form(1, 2))
        @test Form(PDict(:a => 1)) == Form(Form(a=1))
        # Two vectors are two positional values, not a seq and a map.
        @test Form([1], [2]) == Form(Form(1), Form(2))
    end

    @testset "reading" begin
        f = Form("a", "b", x=1, y=2)
        @test f[1] == "a"
        @test f[2] == "b"
        @test f[:x] == 1
        @test f[:y] == 2
        @test length(f) == 4
        @test haskey(f, 1) && haskey(f, 2) && !haskey(f, 3)
        @test haskey(f, :x) && !haskey(f, :z)
        @test get(f, 1, nothing) == "a"
        @test get(f, 9, "dflt") == "dflt"
        @test get(f, :x, nothing) == 1
        @test get(f, :z, "dflt") == "dflt"
        # Out of range is a BoundsError and a missing key a KeyError, as for the
        # collections they come from.
        @test_throws BoundsError f[9]
        @test_throws KeyError f[:z]
    end

    @testset "updating" begin
        f = Form("a", "b", x=1)

        @test push(f, "c") == Form("a", "b", "c", x=1)
        @test push(f, 5)[3] === 5
        @test pop(f) == Form("a", x=1)
        @test_throws ArgumentError pop(Form(x=1))   # nothing in the seq to pop

        @test setindex(f, "z", 1) == Form("z", "b", x=1)
        @test setindex(f, 9, :x) == Form("a", "b", x=9)
        @test setindex(f, 9, :new) == Form("a", "b", x=1, new=9)

        @test delete(f, 1) == Form("b", x=1)
        @test delete(f, :x) == Form("a", "b")
        @test delete(f, :absent) == f

        # The originals are untouched.
        @test f == Form("a", "b", x=1)
    end

    @testset "equality" begin
        @test Form("a", x=1) == Form("a", x=1)
        @test Form("a", x=1) != Form("a", x=2)
        @test Form("a") != Form("b")
        @test Form(1, 2) != Form(2, 1)          # order matters in the seq
        @test isequal(Form("a"), Form("a"))
        @test hash(Form("a", x=1)) == hash(Form("a", x=1))
        @test hash(Form("a", x=1)) != hash(Form("a", x=2))
        @test Form(1.5) == Form(1.5f0)                 # both are Float64
        @test Form(:a) != Form("a")                    # a symbol is not a string
        # A float and a complex with a zero imaginary part are *equal*, because
        # that is Julia's own numeric equality and a form compares its values.
        # They are still stored, and written to JSON, differently.
        @test Form(1.0) == Form(1.0 + 0im)
        @test to_JSON(Form(1.0)) != to_JSON(Form(1.0 + 0im))
    end

    @testset "printing" begin
        # Sequential values first, keyword values after, like a function call.
        @test string(Form()) == "Form()"
        @test string(Form("a")) == "Form(\"a\")"
        @test string(Form("a", 3, x=1)) == "Form(\"a\", 3, x=1)"
        # A nested form prints without the `Form` prefix.
        @test string(Form("widget", meta=Form(name="abc"))) ==
              "Form(\"widget\", meta=(name=\"abc\"))"
        @test string(Form(Form(1, 2), Form(a=1))) == "Form((1, 2), (a=1))"
        # A string and a symbol print differently, as Julia prints them.
        @test string(Form("s", :y)) == "Form(\"s\", :y)"
        @test string(Form(1.5)) == "Form(1.5)"
        @test string(Form(ComplexF64(1, 2))) == "Form(1.0 + 2.0im)"
        @test string(Form(nothing)) == "Form(nothing)"
        @test string(Form(true)) == "Form(true)"
        # The map's keys are sorted, so what is printed is a property of the form
        # rather than of the HAMT's iteration order.
        @test string(Form(b=1, a=2)) == "Form(a=2, b=1)"
    end

    @testset "JSON round-trips" begin
        cases = Any[
            Form(),
            Form("a"),
            Form(1, 2, 3),
            Form(x=1, y=2),
            Form("s", true, false, nothing, 42, -7, 1.5),
            Form(ComplexF64(1, 2)),
            Form(1.0 + 0im),
            Form(:sym),
            Form("str", :sym, k=:v),
            Form([1, 2], Dict(:a => 1)),
            Form("deep", inner=Form("deeper", inner=Form("deepest"))),
            Form(quote_test="he said \"hi\"\\ and\nnewline\ttab"),
            Form(unicode="héllo → 世界 😀"),
            Form(empty=""),
            Form(control="\u0001\u001f"),
        ]
        for c in cases
            j = to_JSON(c)
            @test from_JSON(j) == c
            @test to_JSON(from_JSON(j)) == j   # and the text is stable
        end
    end

    @testset "the JSON encoding is the documented one" begin
        # A form is always the two-element list [seq, map].
        @test to_JSON(Form()) == "[[], {}]"
        # A string carries a leading `'` and a symbol a leading `:`, since JSON
        # has only one string type and a form holds both.
        @test to_JSON(Form("a")) == "[[\"'a\"], {}]"
        @test to_JSON(Form(:a)) == "[[\":a\"], {}]"
        @test from_JSON("[[\"'a\"], {}]")[1] === "a"
        @test from_JSON("[[\":a\"], {}]")[1] === :a
        # Map keys are symbols, so they carry the `:` marker too.
        @test to_JSON(Form(x=1)) == "[[], {\":x\": 1}]"
        @test to_JSON(Form("a", x=1)) == "[[\"'a\"], {\":x\": 1}]"
        @test to_JSON(Form(name="abc")) == "[[], {\":name\": \"'abc\"}]"
        @test to_JSON(Form("widget", meta=Form(name="abc"))) ==
              "[[\"'widget\"], {\":meta\": [[], {\":name\": \"'abc\"}]}]"
        # A float is written with a point, so it reads back as a float rather
        # than as an integer; an integer has none.
        @test to_JSON(Form(1.5)) == "[[1.5], {}]"
        @test to_JSON(Form(1.0)) == "[[1.0], {}]"
        @test to_JSON(Form(2)) == "[[2], {}]"
        @test from_JSON("[[1.0], {}]")[1] === 1.0
        @test from_JSON("[[1], {}]")[1] === Int64(1)
        # A complex number has no JSON equivalent, so it is an object — always,
        # even with a zero imaginary part, since that is a different value from a
        # float.
        @test to_JSON(Form(ComplexF64(1, 2))) == "[[{\"re\": 1.0, \"im\": 2.0}], {}]"
        @test to_JSON(Form(1.0 + 0im)) == "[[{\"re\": 1.0, \"im\": 0.0}], {}]"
        @test from_JSON("[[{\"re\": 1.0, \"im\": 2.0}], {}]")[1] === ComplexF64(1, 2)
        @test from_JSON(to_JSON(Form(1.0 + 0im)))[1] === ComplexF64(1, 0)
    end

    @testset "from_JSON refuses what it cannot read" begin
        # A bare array is not a form: a form is always [seq, map].
        @test_throws ArgumentError from_JSON("[1,2,3]")
        @test_throws ArgumentError from_JSON("[[1],2]")     # the map must be an object
        @test_throws ArgumentError from_JSON("[[1]]")
        @test_throws ArgumentError from_JSON("[")
        @test_throws ArgumentError from_JSON("")
        @test_throws ArgumentError from_JSON("nul")
        @test_throws ArgumentError from_JSON("[[],{}]x")    # trailing data
        @test_throws ArgumentError from_JSON("{\"a\":1,}")  # trailing comma
        @test_throws ArgumentError from_JSON("[[\"'\\q\"],{}]")  # unknown escape
        @test_throws ArgumentError from_JSON("[[\"'unterminated],{}]")
        # A string with no marker is a mistake, not a value: the writer never
        # produces one.
        @test_throws ArgumentError from_JSON("[[\"bare\"],{}]")
        @test_throws ArgumentError from_JSON("[[\"'a\"],{\"k\":1}]")   # key with no marker
        # An empty symbol is legal, if odd, so `":"` is a key rather than a
        # mistake.
        empty_key = from_JSON("[[\"'a\"],{\":\":1}]")
        @test empty_key[1] == "a"
        @test empty_key[Symbol("")] == 1
        @test_throws ArgumentError from_JSON("[[99999999999999999999],{}]")  # not an Int64
        # A bare object is accepted as a form with only a map, since that is a
        # convenient way to write one by hand — but its keys are symbols, so they
        # carry the marker like any other.
        @test from_JSON("{\":a\": 1}") == Form(a=1)
        # A top-level value that is not a form is refused.
        @test_throws ArgumentError from_JSON("5")
        @test_throws ArgumentError from_JSON("\"s\"")
    end

    @testset "escapes and unicode" begin
        @test from_JSON("[[\"'a\\nb\"],{}]")[1] == "a\nb"
        @test from_JSON("[[\"'a\\tb\"],{}]")[1] == "a\tb"
        @test from_JSON("[[\"'a\\\"b\"],{}]")[1] == "a\"b"
        @test from_JSON("[[\"'a\\\\b\"],{}]")[1] == "a\\b"
        @test from_JSON("[[\"'a\\/b\"],{}]")[1] == "a/b"
        @test from_JSON("[[\"'\\u00e9\"],{}]")[1] == "é"
        @test from_JSON("[[\"'\\ud83d\\ude00\"],{}]")[1] == "😀"   # a surrogate pair
        # What we write, we can read: round-trip the awkward ones.
        for s in ("quote\"here", "back\\slash", "new\nline", "tab\there", "é", "😀", "\u0007")
            @test from_JSON(to_JSON(Form(v=s)))[:v] == s
        end
    end
end
