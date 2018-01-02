
include("computation.jl")
include("lookups.jl")

@inline function d_(p::Float64)
    return - 3 / 4 * log(1 - 4.0 / 3 * p)
end

"""
    dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode)

Compute dN and dS, using the [Nei and Gojobori 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

This function requires two iterables `x` and `y`, which yield `Codon{DNA}` or
`Codon{RNA}` type variables. These two types are defined in the BioSequences
package.
"""
function dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS, addone::Bool = false)
    return _dNdS_NG86(x, y, k, code, addone, eltype(x), eltype(y))
end

"""
    dNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}

Compute dN and dS, using the [Nei and Gojobori 1986](https://www.ncbi.nlm.nih.gov/pubmed/3444411) method.

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
