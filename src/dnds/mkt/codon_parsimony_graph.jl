struct PGEdge{C<:Codon}
    rank::Int
    I::C
    J::C
end

Base.isless(x::PGEdge{C}, y::PGEdge{C}) where c <: Codon = isless(x.rank, y.rank)

mutable struct ParsimonyGraph{C<:Codon}
    nvertices::Int
    nedges::Int
    vertices::Vector{C}
    lookup::PairwiseListMatrix{Int,false,Vector{Int}}
end
