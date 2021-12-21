using FemAdjoint
using Documenter

DocMeta.setdocmeta!(FemAdjoint, :DocTestSetup, :(using FemAdjoint); recursive=true)

makedocs(;
    modules=[FemAdjoint],
    authors="Edouard <edouard.oudet@imag.fr> and contributors",
    repo="https://github.com/edljk/FemAdjoint.jl/blob/{commit}{path}#{line}",
    sitename="FemAdjoint.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://edljk.github.io/FemAdjoint.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/edljk/FemAdjoint.jl",
    devbranch="main",
)
