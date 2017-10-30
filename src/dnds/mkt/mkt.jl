# mkt.jl
# ======
#
# An implementation of the McDonald Kreitman test.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

include("refs_and_graphs/lookups.jl")
include("refs_and_graphs/codon_graphs.jl")

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



function multi_short_path(cs::CodonSet{T}, cg::CodonGraph{Codon{T}}, msts::MSTState{Codon{T}}) where T <: NucleicAcid
    reset!(cg, cs)
    mst(vertices, edges)
    mst(cg, msts)
end

function find(v, p)
    if p[v] != v
        p[v] = find(v, p[v])
    end
    return p[v]
end

function union(v1, v2, p, r)
    r1 = find(v1, p)
    r2 = find(v2, p)
    if r1 != r2
        if r[r1] > r[r2]
            p[r2] = r1
        else
            p[r1] = r2
            if r[r1] == r[r2]
                r[r2] += 1
            end
        end
    end
end

function mst(V, E)
    parents = copy(V)
    ranks = zeros(V)
    sort!(V)
    for (weight, v1, v2) in E
        if find(v1, parents) != find(v2, parents)
            union(v1, v2, parents, ranks)
            # add edge to tree
        end
    end
    return tree
end

function mst(codongraph::CodonGraph, state::MSTState)
    reset!(state)
    for (v1, v2, weight) in codongraph
        if find(v1, parents) != find(v2, parents)
            union(v1, v2, parents, ranks)
            # add edge to tree
        end
    end
    return tree
end
