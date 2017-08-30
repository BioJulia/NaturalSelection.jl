# NaturalSelection.jl
# ===================
#
# A julia package that provides methods for detecting the presence, strength
# and effects of natural selection, in biological data.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

__precompile__()

module NaturalSelection

export
    dNdS

include("dnds/dnds.jl")
end
