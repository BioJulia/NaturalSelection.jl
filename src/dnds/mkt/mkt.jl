# mkt.jl
# ======
#
# An implementation of the McDonald Kreitman test.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

include("refs_and_graphs/lookups.jl")
include("refs_and_graphs/codon_graphs.jl")
include("mst.jl")

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

function most_parsimonious_evolution(cs::CodonSet{T}, cg::CodonGraph{Codon{T}}, msts::MSTState{Codon{T}}) where T <: NucleicAcid
    reset!(cg, cs)
    return mst(cg, msts)
end
