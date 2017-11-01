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

function simplest_mutation_path(cs::CodonSet{T}, cg::CodonGraph{Codon{T}}, msts::MSTState{Codon{T}}) where T <: NucleicAcid
    reset!(cg, cs)
    return mst(cg, msts)
end

function div_count(csa::CodonSet{T}, csb::CodonSet{T}, ref::CodonGraphReference{Codon{T}}) where T <: NucleicAcid
    DS = DN = 0
    R = 10
    if (csa & csb) == CodonSet{T}()
        for i in csa
            for j in csb
                edge = lookup_edge(ref, i, j)
                DS_i = edge[3]
                DN_i = edge[4]
                R_i = rankof(DS_i, DN_i)
                lowrank = R_i < R
                DS = ifelse(lowrank, DS_i, DS)
                DN = ifelse(lowrank, DN_i, DN)
                R = ifelse(lowrank, R_i, R)
            end
        end
    end
    return DS, DN
end

function site_count(csa::CodonSet{T}, csb::CodonSet{T}, cg::CodonGraph{Codon{T}}, msts::MSTSState{Codon{T}}) where T <: NucleicAcid
    a = simplest_mutation_path(cs, cg, msts)

end

function mkt(sa::Vector{Vector{Codon{DNA}}}, sb::Vector{Vector{Codon{DNA}}})


end
