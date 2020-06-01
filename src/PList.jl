################################################################################
# PList.jl
# Persistent lists / sequences for Julia.
# These are not really intended to be used in common code, but they are used
# by PDict and PSet to track elements with identical hashes.

"""
PLink{T}
    PLink objects make up the persistent links in persistnet linked lists.
"""
struct PLink{T}
    first::T
    rest::Union{Nothing,PLink{T}}
end

struct PList{T}
    _n::Int
    _data::Union{Nothing,PLink{T}}
end
PList{T}() where {T} = PList{T}(0, nothing)
PList() = PList{Any}()
PList{T}(els::Vararg{T,N}) where {T,N} = begin
    p = PLink{T}(els[end], nothing)
    for k = 1:length(els)-1
        p = PLink{T}(els[end-k], p)
    end
    return PList{T}(N, p)
end
PList(els...) = begin
    T = typejoin(map(typeof, els)...)
    return PList{T}(els...)
end
PList{T}(u::PList{S}) where {T,S} = PList{T}(u...)
PList(u::PList{T}) where {T} = u

Base.length(l::PList) = l._n
Base.in(k, l::PList{T}, eqfn::Function) where {T} = begin
    (l._n == 0) && return false
    q = l._data::PLink{T}
    while true
        eqfn(k, q.first) && return true
        #(q.rest === nothing) && break
        if q.rest === nothing
            break
        else
            q = q.rest
        end
    end
    return false
end
Base.in(k, l::PList{T}) where {T} = in(k, l, isequiv)
Base.first(l::PList{T}) where {T} = begin
    (l._data === nothing) && throw(BoundsError([l], [l]))
    return l._data.first
end
pushfirst(l::PList{T}, el::S) where {T,S} = begin
    return PList{T}(l._n + 1, PLink{T}(el, l._data))
end
popfirst(l::PList{T}) where {T} = begin
    (l._data === nothing) && throw(ArgumentError("list must be non-empty"))
    return PList{T}(l._n - 1, l._data.rest)
end
_delete(l::PLink{T}, s::S, eqfn::Function) where {T, S} = begin
    if eqfn(s, l.first)
        (l.rest === nothing) && return (1, nothing)
        (n, ll) = _delete(l.rest, s, eqfn)
        return (n+1, ll)
    elseif l.rest === nothing
        return (0, l)
    else
        (n, q) = _delete(l.rest, s, eqfn)
        return (n, q === l.rest ? l : PLink{T}(l.first, q))
    end
end
delete(l::PList{T}, s::S, eqfn::Function) where {T, S} = begin
    (l._data === nothing) && return l
    (n, q) = _delete(l._data, s, eqfn)
    return q === l._data ? l : PList{T}(l._n - n, q)
end
delete(l::PList{T}, s::S) where {T, S} = delete(l, s, isequiv)
Base.show(io::IO, l::PList{T}) where {T} = begin
    print(io, "$(typeof(l))(")
    k = l._data
    i = 0
    while k !== nothing
        (i > 0) && print(io, ", ")
        print(io, "$(k.first)")
        k = k.rest
        i += 1
        if i == 10
            print(io, "... <$(l._n - 10) more>")
            break
        end
    end
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", l::PList{T}) where {T} = begin
    println(io, "$(l._n)-element $(typeof(l)):")
    k = l._data
    i = 0
    while k !== nothing
        println(io, " $(k.first)")
        k = k.rest
        i += 1
        if i == 10
            println(io, " ... <$(l._n - 10) more>")
            break
        end
    end
end
Base.iterate(l::PList{T}) where {T} = let d = l._data
    d === nothing ? nothing : (d.first, d.rest)
end
Base.iterate(l::PList{T}, u::PLink{T}) where {T} = (u.first, u.rest)
Base.iterate(l::PList{T}, ::Nothing) where {T} = nothing
Base.convert(::Type{PList{T}}, u::PList{S}) where {T,S} = PList{T}(u...)

            
