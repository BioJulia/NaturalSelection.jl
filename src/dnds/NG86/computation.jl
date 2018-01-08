# NG86/computation.jl
# ===================
#
# Core computations for the NG86 method of computing dN and dS statistics.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

"""
    S_N_NG86(codon::C, code::GeneticCode) where {C <: CDN}

Enumerate the number of synonymous (S) and non-synonymous (N) sites in a codon,
using the method used by [Nei and Gojobori (1986)](https://www.ncbi.nlm.nih.gov/pubmed/3444411).

Returns a tuple where S is the first element and N is the second (S, N).

Each site in a codon may be both partially synonymous and non-synonymous.
"""
function S_N_NG86(codon::C, code::GeneticCode) where {C <: Codon}
    cdn_bits = UInt64(codon)
    aa = code[codon]
    S = N = 0
    for (pos, msk) in enumerate(CDN_POS_MASKS)
        bidx = bitindex(codon, pos)
        for base in 0:3
            # Create the neighbor codon.
            neighbor = C((cdn_bits & msk) | (base << bidx))
            if codon == neighbor # Codon created is not a neighbor: should happen 3 times.
                continue
            end
            # See if the protein changes between codon and neighbor, and update
            # N and S counts accordingly.
            neighbor_aa = code[neighbor]
            if neighbor_aa == BioSequences.AA_Term
                N += 1
            elseif neighbor_aa == aa
                S += 1
            else
                N += 1
            end
        end
    end
    normalization = (N + S) / 3
    return (S / normalization), (N / normalization)
end

"""
    DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN

Compute the number of synonymous (DS) and non-synonymous (DN) mutations between
two codons, using the all paths method used by the [Nei and Gojobori (1986)](https://www.ncbi.nlm.nih.gov/pubmed/3444411).
"""
function DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: Codon
    return DS_DN_enumerator(AllPaths, x, y, code)
end

@inline function d_(p::Float64)
    if p < 3 / 4
        return - 3 / 4 * log(1 - 4.0 / 3 * p)
    else
        return -1.0
    end
end

function dNdS_NG86_kernel(x, y,
    S::Float64, N::Float64,
    DS::Float64, DN::Float64,
    snlookup::S_N_NG86_LOOKUP, dsdnlookup::DS_DN_NG86_LOOKUP)

    # Iterate over every pair of codons.
    @inbounds for (i, j) in zip(x, y)
        si, ni = snlookup[i]
        sj, nj = snlookup[j]
        S += (si + sj)
        N += (ni + nj)
        DSi, DNi = dsdnlookup[i, j]
        DS += DSi
        DN += DNi
    end

    S = S / 2.0
    N = N / 2.0
    pN = DN / N
    pS = DS / S
    dN = d_(pN)
    dS = d_(pS)
    return dN, dS
end

function _dNdS_NG86(x, y, addone::Bool, snlookup::S_N_NG86_LOOKUP, dsdnlookup::DS_DN_NG86_LOOKUP, ::Type{C}, ::Type{C}) where C <: Codon
    # Expected no. of syn and nonsyn sites.
    S = N = 0.0
    # Observed no. of syn and nonsyn mutations.
    DS = ifelse(addone, 1.0, 0.0)
    DN = 0.0
    return dNdS_NG86_kernel(x, y, S, N, DS, DN, snlookup, dsdnlookup)
end

function _dNdS_NG86(x, y, addone::Bool, snlookup::S_N_NG86_LOOKUP, dsdnlookup::DS_DN_NG86_LOOKUP, ::Type{N}, ::Type{N}) where N <: NucleicAcid
    xcdns, ycdns = aligned_codons(x, y)
    _dNdS_NG86(xcdns, ycdns, addone, snlookup, dsdnlookup, eltype(xcdns), eltype(ycdns))
end
