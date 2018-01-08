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
    tajimad,
    mkt,
    make_edge_reference

import BioSequences:
    BioSequences,
    DNA,
    BioSequence,
    DNASequence,
    NucleicAcid,
    Kmer,
    NucAlphs,
    GeneticCode,
    ispurine,
    @kmer_str,
    ncbi_trans_table,
    standard_genetic_code

import GeneticVariation:
    NL79,
    Segregating,
    avg_mut

using PairwiseListMatrices

const DEFAULT_TRANS = BioSequences.ncbi_trans_table[1]
const Codon{T} = BioSequences.Kmer{T, 3}

@inline bitindex(x::Kmer{T,K}, i::Integer) where {T,K} = 2 * (K - i)

include("codons/codon_set.jl")

include("dnds/dnds.jl")
include("tajima.jl")
end
