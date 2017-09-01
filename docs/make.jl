using Documenter, NaturalSelection

makedocs(
    modules = [NaturalSelection.dNdS],
    format = :html,
    sitename = "NaturalSelection.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "dNdS" => "man/dNdS.md",
        ],
        "Contributing" => "contributing.md"
    ],
    authors = "Ben J. Ward, The BioJulia organisation, and contributors."
)

deploydocs(
    repo = "github.com/BioJulia/NaturalSelection.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
