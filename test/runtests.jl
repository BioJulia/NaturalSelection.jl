module TestNaturalSelection

using Base.Test

using BioSequences, NaturalSelection

include("dNdS.jl")
include("tajima.jl")
include("mkt.jl")

end
