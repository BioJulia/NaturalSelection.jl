# mkt.jl
# ======
#
# An implementation of the McDonald Kreitman test.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

include("lookups.jl")

function codondiff(sequences::Vector{Vector{Codon{DNA}}}, position::Integer)
    s = CodonSet{DNA}()
    @inbounds for seq in sequences
        s |= seq[position]
    end
    return s
end

function codondiff(sequences::Vector{Vector{Codon{DNA}}})
    out = Vector{CodonSet{DNA}}(length(sequences[1]))
    @inbounds for site in eachindex(sequences[1])
        out[site] = codondiff(sequences, site)
    end
    return out
end



function multi_short_path(codonset::CodonSet{T},
                          codoncache::Vector{Codon{T}},
                          edgecache::Vector{Tuple{Codon{T}, Codon{T}, Int}},
                          rankref::PairwiseListMatrix) where T <: NucleicAcid

    # Fill the codon cache with codons from the codon set.
    ncodons = length(codonset)
    copy!(codoncache, codonset)
    nedges = div(ncodons * (ncodons - 1), 2)

    e = 1
    for i in 1:ncodons
        @inbounds for j in (i + 1):ncodons
            x = codoncache[i]
            y = codoncache[j]
            edgecache[e] = (lookup(rankref, x, y), x, y)
            e += 1
        end
    end

    println(codoncache)
    println(edgecache)

end
