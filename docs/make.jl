using Documenter
using FluorescenceExclusion

DocMeta.setdocmeta!(FluorescenceExclusion, :DocTestSetup, :(using FluorescenceExclusion); recursive=true)
makedocs(
    modules=[FluorescenceExclusion], 
    sitename="FluorescenceExclusion.jl",
    authors="Tamas Nagy",
    pages = [
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/tlnagy/FluorescenceExclusion.jl.git",
)