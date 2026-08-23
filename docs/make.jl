using SnellEnvelope
using Documenter

DocMeta.setdocmeta!(SnellEnvelope, :DocTestSetup, :(using SnellEnvelope); recursive=true)

makedocs(;
    modules=[SnellEnvelope],
    authors="Sebastian Schlenkrich",
    sitename="SnellEnvelope.jl",
    format=Documenter.HTML(;
        canonical="https://sschlenkrich.github.io/SnellEnvelope.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sschlenkrich/SnellEnvelope.jl",
    devbranch="main",
)
