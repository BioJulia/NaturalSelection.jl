using Documenter, NaturalSelection

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/NaturalSelection.jl.git",
    julia = "0.6",
    osname = "linux",
)
