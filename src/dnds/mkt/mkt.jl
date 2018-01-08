# mkt.jl
# ======
#
# An implementation of the McDonald Kreitman test.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

include("refs_and_graphs/codon_graphs.jl")
include("mst.jl")

function CodonSet{T}(sequences::Vector{Vector{Codon{T}}},
                     position::Integer) where T <: NucleicAcid
    s = CodonSet{T}()
    @inbounds for seq in sequences
        s |= seq[position]
    end
    return s
end

function simplest_mutation_path(cs::CodonSet{T},
                                cg::CodonGraph{Codon{T}},
                                msts::MSTState{Codon{T}}) where T <: NucleicAcid
    reset!(cg, cs)
    return mst(cg, msts)
end

function poly_count(cs::CodonSet{T},
                    cg::CodonGraph{Codon{T}},
                    msts::MSTState{Codon{T}}) where T <: NucleicAcid

    return simplest_mutation_path(cs, cg, msts)
end

function div_count(csa::CodonSet{T},
                   csb::CodonSet{T},
                   ref::CodonGraphReference{Codon{T}}) where T <: NucleicAcid
    DS = DN = 0
    R = 10
    if (csa & csb) == CodonSet{T}()
        for i in csa
            for j in csb
                edge = lookup_edge(ref, i, j)
                @inbounds DS_i = edge[3]
                @inbounds DN_i = edge[4]
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

function mkt_PSPN(x::CodonSet{T},
                  y::CodonSet{T},
                  cg::CodonGraph{Codon{T}},
                  msts::MSTState{Codon{T}}) where T <: NucleicAcid

    a = poly_count(x, cg, msts)
    b = poly_count(y, cg, msts)
    return a[1] + b[1], a[2] + b[2]
end

mkt_α(PS, PN, DS, DN) = 1 - (DS * PN) / (DN * PS)


"""
    mkt(x, y, ref)

Compute McDonald Kreitman statistics for two sets of codons.
"""
function mkt(x::Vector{Vector{Codon{T}}},
             y::Vector{Vector{Codon{T}}},
             ref::CodonGraphReference{Codon{T}} = DEFAULT_CODON_GRAPH_REFERENCE) where T <: NucleicAcid

    n = min(minimum(length(xi) for xi in x), minimum(length(yi) for yi in y))
    graph = CodonGraph{Codon{T}}(ref)
    msts = MSTState{Codon{T}}()
    PS = PN = DS = DN = 0
    for i in 1:n
        xset = CodonSet{T}(x, i)
        yset = CodonSet{T}(y, i)

        PS_i, PN_i = mkt_PSPN(xset, yset, graph, msts)
        DS_i, DN_i = div_count(xset, yset, ref)

        PS += PS_i
        PN += PN_i
        DS +- DS_i
        DN += DN_i
    end
    α = mkt_α(PS, PN, DS, DN)
    return PS, PN, DS, DN, α
end
