using Documenter, Air

makedocs(;
    modules=[Air],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        # The default branch, named explicitly. Documenter otherwise asks git for
        # the remote HEAD and falls back to a `master` that this repository has
        # not had since the rename, which would leave every "edit this page" link
        # pointing at a branch that does not exist.
        edit_link="main",
        # `API.md` is one reference page holding every public docstring, so it is
        # large by construction: the default warning threshold of 100 KiB and
        # error threshold of 200 KiB are sized for prose pages, and this page is
        # most of the way to the error already. Raising them is a decision about
        # this page, not a way of ignoring a problem.
        size_threshold_warn=250 * 1024,
        size_threshold=500 * 1024,
    ),
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
    # A `Remotes.GitHub` rather than a URL string. Documenter derives the navbar
    # link, the "edit this page" link and the commit of a source link from it, and
    # with a bare string it warns that it cannot — and then falls back to a
    # `master` default that this repository does not have.
    repo=Documenter.Remotes.GitHub("noahbenson", "Air.jl"),
    sitename="Air.jl",
    authors="Noah C. Benson",
    # `:exports` requires that every exported symbol's docstring appear somewhere
    # in these pages. `:all` — which would also demand that every *internal*
    # docstring appear — is deliberately not used: those are written for people
    # working on Air, not for people reading the manual.
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/noahbenson/Air.jl",
    deploy_config=Documenter.GitHubActions(),
    push_preview=false,
    # The repository's default branch is `main`. `deploydocs` otherwise infers
    # the dev branch by asking git for the remote HEAD and silently falls back to
    # "master" when it cannot, which would leave the dev docs never updated.
    devbranch="main",
)
