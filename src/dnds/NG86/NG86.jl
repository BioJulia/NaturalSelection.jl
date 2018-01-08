
const S_N_NG86_LOOKUP = SingleCodonLookup{Tuple{Float64, Float64}}
const DS_DN_NG86_LOOKUP = PairwiseCodonLookup{Tuple{Float64, Float64}}

include("computation.jl")
include("lookups.jl")

"""
    dNdS_NG86(x, y, addone::Bool = true, code::Int = 1)

Compute dN and dS, using the [Nei and Gojobori 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

The genetic code that is used, is defined according to the numbering of
`ncbi_trans_table`. Code 1 is the standard genetic code.

This function requires two iterables `x` and `y`.
If these iterables yield `Codon{DNA}` or `Codon{RNA}` type variables. Then it is
assumed that `x` and `y` are iterables that yield a sequence of aligned codons.
If the iterables produce `DNA` or `RNA` type variables, then it is assumed `x`
and `y` iterables that conform to the behaviour of DNA or RNA sequences as
defined in the BioSequences package. In this case, a new `x` and `y` that do
have an element type of `Codon{DNA}` or `Codon{RNA}`.

NG86 is a counting method of computing dN/dS and is typically safer to use on
sequence data where codon usage, (esp. at 3rd position), is uniform,
the sequences are not very divergent, and transition/transversion rates, are similar.
"""
function dNdS_NG86(x, y; addone::Bool = false, code::Int = 1)
    snlookup = S_N_NG86_LOOKUPS[code]
    dsdnlookup = DS_DN_NG86_LOOKUPS[code]
    return _dNdS_NG86(x, y, addone, snlookup, dsdnlookup, eltype(x), eltype(y))
end

function pairwise_dNdS_NG86(x, addone::Bool = false, code::Int = 1)
    n = length(x)
    @assert n >= 2 "At least two sequences are required."
    snlookup = S_N_NG86_LOOKUPS[code]
    dsdnlookup = DS_DN_NG86_LOOKUPS[code]
    results = Matrix{Tuple{Float64, Float64}}(n, n)
    for i in 1:n
        results[i,i] = 0.0, 0.0
        for j in (i + 1):n
            results[i,j] = results[j,i] = _dNdS_NG86(x[i], x[j], addone, snlookup, dsdnlookup, eltype(x[i]), eltype(x[j]))
        end
    end
    return results
end
