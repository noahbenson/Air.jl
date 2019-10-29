using Documenter, Air

makedocs(;
    modules=[Air],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/noahbenson/Air.jl/blob/{commit}{path}#L{line}",
    sitename="Air.jl",
    authors="Noah C. Benson",
    assets=String[],
)

deploydocs(;
    repo="github.com/noahbenson/Air.jl",
)
