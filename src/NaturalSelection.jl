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
    dNdS_NG86,
    S_N_NG86,
    DS_DN_NG86,
    tajimad

import BioSequences:
    BioSequences,
    BioSequence,
    Kmer,
    NucAlphs,
    GeneticCode,
    ispurine

import GeneticVariation:
    NL79,
    Segregating,
    avg_mut

include("dnds/dnds.jl")
include("tajima.jl")
include("mkt.jl")
end
