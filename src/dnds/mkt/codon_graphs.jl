
const DSDN_RANK_LOOKUP = [0 2 5 9;
                          1 4 8 0;
                          3 7 0 0;
                          6 0 0 0;]

@inline function rankof(DS::Integer, DN::Integer)
    @inbounds return DSDN_RANK_LOOKUP[DS + 1, DN + 1]
end

struct CodonGraphReference{C <: Codon}
    edges::Vector{Tuple{C, C, Int}}
    edge_permutation::Vector{Int}
    edge_order::Vector{Int}
end

function CodonGraphReference{C}(code::GeneticCode) where C <: Codon
    edgetable = make_edge_reference(code)
    edgetable_perm = sortperm(edgetable, by = x -> x[3])
    edgetable_order = zeros(edgetable_perm)
    for i in eachindex(edgetable_perm)
        perm = edgetable_perm[i]
        edgetable_order[perm] = i
    end
    return CodonGraphReference{C}(edgetable, edgetable_perm, edgetable_order)
end

function make_edge_reference(code::GeneticCode{C}) where C <: Codon
    edges = Vector{Tuple{C, C, Int}}(2016)
    k = 1
    for i in 0x00:0x3F
        for j in (i + 1):0x3F
            x = C(UInt64(i))
            y = C(UInt64(j))
            DS, DN = DS_DN_enumerator(ShortestPath, x, y, code)
            edges[k] = (x, y, rankof(Integer(DS), Integer(DN)))
            k += 1
        end
    end
    return edges
end

@inline function ij_to_k(i::Integer, j::Integer)
    nelements = 64
    x = nelements - i
    return div(nelements * (nelements - 1) - (x * (x - 1)), 2) - nelements + j
end

@inline function codons_to_k(i::T, j::T) where T <: Codon
    return ij_to_k(UInt64(i) + 1, UInt64(j) + 1)
end

@inline function lookup_edge(table::CodonGraphReference{T}, i::T, j::T) where T <: Codon
    @inbounds return table.edges[codons_to_k(i, j)]
end

@inline function lookup_perm(table::CodonGraphReference{T}, i::T, j::T) where T <: Codon
    @inbounds return table.edge_permutation[codons_to_k(i, j)]
end

@inline function lookup_order(table::CodonGraphReference{T}, i::T, j::T) where T <: Codon
    @inbounds return table.edge_order[codons_to_k(i, j)]
end

const DEFAULT_CODON_GRAPH_REFERENCE = CodonGraphReference{Codon{DNA}}(DEFAULT_TRANS)

mutable struct CodonGraph{C<:Codon}
    ref::CodonGraphReference{C}
    vertices::Vector{C}
    edges::BitVector
    nv::Int
end

function CodonGraph{C}(ref::CodonGraphReference{C}) where C <: Codon
    return CodonGraph{C}(ref, Vector{C}(64), flases(2016))
end

function CodonGraph{C}(code::GeneticCode) where C <: Codon
    reference = CodonGraphReference{C}(code)
    return CodonGraph{C}(reference)
end

function reset!(cg::CodonGraph, cs::CodonSet)
    # Fill the vertices buffer with codons from the codon set.
    cg.nv = length(cs)
    copy!(cg.vertices, cs)
    # Reset and turn on the bits indicating which edges we have.
    # Order of bits is not the same as cg.ref.edges, but is dictated by
    # cg.ref.edge_order. This removes the need to sort a lot, removing a key
    # bottleneck in Kruskal's algorithm.
    fill!(cg.edges, false)
    for i in 1:cg.nv
        x = cg.vertices[i]
        @inbounds for j in (i + 1):(cg.nv)
            y = cg.vertices[j]
            e = lookup_order(cg.ref, x, y)
            cg.edges[e] = true
        end
    end
end

@inline Base.start(x::CodonGraph) = findnext(x, 1)
@inline function Base.next(x::CodonGraph, state::Int)
    @inbounds e = x.ref.edges[x.ref.edge_permutation[state]]
    s = findnext(x, state + 1)
    return e, s
end
@inline Base.done(x::CodonGraph, state::Int) = state == 0
Base.eltype(x::CodonGraph{C}) where C <: Codon = Tuple{C, C, Int}
Base.iteratorsize(::Type{CodonGraph{C}}) where C <: Codon = HasLength()
Base.iteratoreltype(::Type{CodonGraph{C}}) where C <: Codon = HasEltype()
Base.length(x::CodonGraph) = count(x.edges)

struct MSTState{C <: Codon}
    parents::Vector{C}
    ranks::Vector{Int}
end

function MSTState{C}() where C <: Codon
    MSTState{C}(collect(1:64), zeros(Int, 64))
end

@inline function reset!(x::MSTState{C}) where C <: Codon
    fill!(x.ranks, 0)
    copy!(x.parents, 1:64)
end
