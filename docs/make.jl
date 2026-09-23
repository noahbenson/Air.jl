using Documenter, Air

makedocs(;
    modules=[Air],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "Persistent Arrays" => "parray.md",
        "Persistent Dictionaries" => "pdict.md",
        "Persistent Sets" => "pset.md",
        "Persistent Weighted Dictionaries" => "pwdict.md",
        "Persistent Weighted Sets" => "pwset.md",
        "Persistent Lazy Dictionaries" => "lazydict.md",
        "Persistent Heaps" => "pheap.md",
        "Transactions (STM)" => "stm.md",
        "Task-Local Variables" => "var.md",
        "Utilities" => "util.md",
        "API Reference" => "API.md",
    ],
    repo="https://github.com/noahbenson/Air.jl/blob/{commit}{path}#L{line}",
    sitename="Air.jl",
    authors="Noah C. Benson",
    # Several internal docstrings cross-reference each other by unqualified
    # name, so strict cross-reference checking is left off for now.
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/noahbenson/Air.jl",
    deploy_config=Documenter.GitHubActions(),
    push_preview=false,
)
