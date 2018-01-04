
const S_N_NG86_LOOKUP = CodonLookupTable{1, Tuple{Float64, Float64}}
const DS_DN_NG86_LOOKUP = CodonLookupTable{2, Tuple{Float64, Float64}}

include("computation.jl")
include("lookups.jl")

"""
    dNdS_NG86(x, y, addone::Bool = true, code::GeneticCode = BioSequences.standard_genetic_code)

Compute dN and dS, using the [Nei and Gojobori 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

If the genetic code is one of those defined in `ncbi_trans_table` then
the correct D, N, DS, and DN lookup tables for NG86 will be used, and if it is
a unique user defined genetic code then D, N, DS, and DN lookup tables will be
generated.

This function requires two iterables `x` and `y`.
If these iterables yield `Codon{DNA}` or `Codon{RNA}` type variables. Then it is
assumed that `x` and `y` are iterables that yield a sequence of aligned codons.
If the iterables produce `DNA` or `RNA` type variables, then it is assumed `x`
and `y` iterables that conform to the behaviour of DNA or RNA sequences as
defined in the BioSequences package. In this case, a new `x` and `y` that do
have an element type of `Codon{DNA}` or `Codon{RNA}`.
"""
function dNdS_NG86(x, y; addone::Bool = true, code::GeneticCode = BioSequences.standard_genetic_code)
    snlookup = get(S_N_NG86_LOOKUPS, make_S_N_NG86_table(code))
    dsdnlookup = get(DS_DN_NG86_LOOKUPS, make_DS_DN_NG86_table(code))
    return dNdS_NG86(x, y, addone, snlookup, dsdnlookup)
end

function dNdS_NG86(x, y, addone::Bool = true, snlookup::SN_NG86_LOOKUP, dsdnlookup::DSDN_NG86_LOOKUP)
    return _dNdS_NG86(x, y, addone, snlookup, dsdnlookup, eltype(x), eltype(y))
end

function pairwise_dNdS_NG86(x, addone::Bool = true, snlookup::SN_NG86_LOOKUP, dsdnlookup::DSDN_NG86_LOOKUP)
    n = length(x)
    @assert n >= 2 "At least two sequences are required."
    results = Matrix{Tuple{Float64, Float64}}(n, n)
    for i in 1:n
        results[i,i] = 0.0, 0.0
        for j in (i + 1):n
            results[i,j] = results[j,i] = dNdS_NG86(x[i], x[j], addone, snlookup, dsdnlookup)
        end
    end
    return results
end

function pairwise_dNdS_NG86(x; addone::Bool = true, code::GeneticCode = BioSequences.standard_genetic_code)
    snlookup = get(S_N_NG86_LOOKUPS, make_S_N_NG86_table(code))
    dsdnlookup = get(DS_DN_NG86_LOOKUPS, make_DS_DN_NG86_table(code))
    return pairwise_dNdS_NG86(x, addone, snlookup, dsdnlookup)
end
