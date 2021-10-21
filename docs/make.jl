using Documenter, Air

makedocs(;
         modules=[Air],
         format=Documenter.HTML(),
         pages=[
             "Home" => "index.md",
             "Persistent Arrays" => "parray.md",
             "Persistent Dictionaries" => "pdict.md",
             "Persistent Sets" => "pset.md",
             "Persistent Weighted Dictionaries" => "pwdict.md",
             "Persistent Weighted Sets" => "pwset.md",
             "Persistent Lazy Dictionaries" => "lazydict.md"],
         repo="https://github.com/noahbenson/Air.jl/blob/{commit}{path}#L{line}",
         sitename="Air.jl",
         authors="Noah C. Benson",
         strict=true,
         #assets=String[],
)

deploydocs(;
           repo="github.com/noahbenson/Air.jl",
           deploy_config=Documenter.GitHubActions(),
           push_preview=false
)
