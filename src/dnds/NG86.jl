# NG86.jl
# =======
#
# dNdS computation using the NG86 method.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

"""
    dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode)

Compute dN and dS, using the [Nei and Goborjei 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

This function requires two iterables `x` and `y`, which yield `DNACodon` or
`RNACodon` type variables. These two types are defined in the BioSequences
package.
"""
function dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS, addone::Bool = false)
    return _dNdS_NG86(x, y, k, code, addone, eltype(x), eltype(y))
end

"""
    dNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}

Compute dN and dS, using the [Nei and Goborjei 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

This method adds conveinience when working with DNA or RNA sequences, by taking
two sequences, and creating two vectors of aligned codons from them. These two
iterables are then passed into the generic NG86 method.
"""
function dNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS, addone::Bool = false) where {A <: NucAlphs}
    xcdns, ycdns = aligned_codons(x, y)
    return dNdS_NG86(xcdns, ycdns, k, code, addone)
end

function pairwise_dNdS_NG86(x, opt...)
    n = length(x)
    @assert n >= 2 "At least two sequences are required."
    results = Matrix{Tuple{Float64, Float64}}(n, n)
    for i in 1:n
        results[i,i] = 0.0, 0.0
        for j in (i + 1):n
            results[i,j] = results[j,i] = dNdS_NG86(x[i], x[j], opt...)
        end
    end
    return results
end

function _dNdS_NG86(x, y, k::Float64, code::GeneticCode, addone::Bool, xtype::Type{C}, ytype::Type{C}) where C <: Codon
    # Expected no. of syn and nonsyn sites.
    S = N = 0.0
    # Observed no. of syn and nonsyn mutations.
    DS = ifelse(addone, 1.0, 0.0)
    DN = 0.0
    # Iterate over every pair of codons.
    @inbounds for (i, j) in zip(x, y)
        si, ni = S_N_NG86(i, k, code)
        sj, nj = S_N_NG86(j, k, code)
        S += (si + sj)
        N += (ni + nj)
        DSi, DNi = DS_DN_NG86(i, j, code)
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

"""
    S_N_NG86(codon::C, k::Float64, code::GeneticCode) where {C <: CDN}

Enumerate the number of synonymous (S) and non-synonymous (N) sites in a codon,
using the method used by [Nei and Goborjei (1986)](https://www.ncbi.nlm.nih.gov/pubmed/3444411).

Returns a tuple where S is the first element and N is the second (S, N).

Each site in a codon may be both partially synonymous and non-synonymous.
"""
function S_N_NG86(codon::C, k::Float64, code::GeneticCode) where {C <: Codon}
    cdn_bits = UInt64(codon)
    aa = code[codon]
    S = N = 0.0
    for (pos, msk) in enumerate(CDN_POS_MASKS)
        bidx = bitindex(codon, pos)
        @inbounds for base in 0:3
            # Create the neighbor codon.
            neighbor = C((cdn_bits & msk) | (base << bidx))
            if codon == neighbor # Codon created is not a neighbor: should happen 3 times.
                continue
            end
            # See if the mutation is transition or transversion.
            cdn_purine = ispurine(codon[pos])
            neighbor_purine = ispurine(neighbor[pos])
            istransition = (cdn_purine && neighbor_purine) || (!cdn_purine && !neighbor_purine)
            # See if the protein changes between codon and neighbor, and update
            # N and S counts accordingly.
            inc = ifelse(istransition, 1.0, k)
            neighbor_aa = code[neighbor]
            if neighbor_aa == BioSequences.AA_Term
                N += inc
            elseif neighbor_aa == aa
                S += inc
            else
                N += inc
            end
        end
    end
    normalization = (N + S) / 3
    return (S / normalization), (N / normalization)
end

"""
    DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN

Compute the number of synonymous (DS) and non-synonymous (DN) mutations between
two codons, using the all paths method used by the [Nei and Goborjei (1986)](https://www.ncbi.nlm.nih.gov/pubmed/3444411).
"""
function DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: Codon
    return DS_DN_enumerator(AllPaths, x, y, code)
end

@inline function d_(p::Float64)
    return - 3 / 4 * log(1 - 4.0 / 3 * p)
end
