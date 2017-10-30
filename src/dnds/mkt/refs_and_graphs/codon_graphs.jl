
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
    @inbounds for i in 1:cg.nv
        x = cg.vertices[i]
        for j in (i + 1):(cg.nv)
            y = cg.vertices[j]
            e = lookup_order(cg.ref, x, y)
            cg.edges[e] = true
        end
    end
end

@inline Base.start(x::CodonGraph) = findnext(x, 1)
@inline function Base.next(x::CodonGraph, state::Int)
    e = x.ref.edges[x.ref.edge_permutation[state]]
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
    MSTState{C}(collect((C(i) for i in 0x00:0x3F)), zeros(Int, 64))
end

@inline function reset!(x::MSTState{C}) where C <: Codon
    fill!(x.ranks, 0)
    copy!(x.parents, (C(i) for i in 0x00:0x3F))
end
