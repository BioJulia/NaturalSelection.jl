
struct CodonLookupTable{T}
    table::T
end

const SingleCodonLookup{T} = CodonLookupTable{Vector{T}}
const PairwiseCodonLookup{T} = CodonLookupTable{PairwiseListMatrix{T, true, Vector{T}}}

function SingleCodonLookup{T}(f::Function) where T
    v = SingleCodonLookup{T}(Vector{T}(64))
    cdn = kmer"AAA"
    @inbounds while cdn < kmer"TTT"
        v[cdn] = f(cdn)
        cdn += 1
    end
    v[cdn] = f(cdn)
    return v
end

function PairwiseCodonLookup{T}(f::Function) where T
    v = PairwiseCodonLookup{T}(PairwiseListMatrix(Vector{T}(2080), true))
    ci = kmer"AAA"
    while ci < kmer"TTT"
        cj = ci
        @inbounds while cj < kmer"TTT"
            v[ci, cj] = f(ci, cj)
            cj += 1
        end
        v[ci, cj] = f(ci, cj)
        ci += 1
    end
    return v
end

@inline function Base.getindex(tbl::SingleCodonLookup, cdn::C) where C <: Codon
    return tbl.table[UInt64(cdn) + 1]
end
@inline function Base.setindex!(tbl::SingleCodonLookup{T}, val::T, cdn::C) where {C <: Codon, T}
    return tbl.table[UInt64(cdn) + 1] = val
end

@inline function Base.getindex(tbl::PairwiseCodonLookup, i::C, j::C) where C <: Codon
    return tbl.table[UInt64(i) + 1, UInt64(j) + 1]
end
@inline function Base.setindex!(tbl::PairwiseCodonLookup{T}, val::T, i::C, j::C) where {C <: Codon, T}
    return tbl.table[UInt64(i) + 1, UInt64(j) + 1] = val
end









## To Sort out!

function make_edge_reference(::Type{C}, code::GeneticCode) where C <: Codon
    edges = Vector{Tuple{C, C, Int, Int}}(2016)
    k = 1
    for i in 0x00:0x3F
        for j in (i + 1):0x3F
            x = C(UInt64(i))
            y = C(UInt64(j))
            DS, DN = DS_DN_enumerator(ShortestPath, x, y, code)
            edges[k] = (x, y, Integer(DS), Integer(DN))
            k += 1
        end
    end
    return edges
end

const DSDN_RANK_LOOKUP = [0 2 5 9;
                          1 4 8 0;
                          3 7 0 0;
                          6 0 0 0;]

@inline function rankof(DS::Integer, DN::Integer)
    @inbounds return DSDN_RANK_LOOKUP[DS + 1, DN + 1]
end

struct CodonGraphReference{C <: Codon}
    edges::Vector{Tuple{C, C, Int, Int}}
    permutation::Vector{Int}
    order::Vector{Int}
end

function CodonGraphReference{C}(code::GeneticCode) where C <: Codon
    edgetable = make_edge_reference(C, code)
    edgetable_perm = sortperm(edgetable, by = x -> rankof(x[3], x[4]))
    edgetable_order = zeros(edgetable_perm)
    for i in eachindex(edgetable_perm)
        perm = edgetable_perm[i]
        edgetable_order[perm] = i
    end
    return CodonGraphReference{C}(edgetable, edgetable_perm, edgetable_order)
end

@inline function ij_to_k(i::Integer, j::Integer)
    i, j = ifelse(i > j, (j, i), (i, j))
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

@inline function lookup_edge(table::CodonGraphReference, i::Int)
    @inbounds return table.edges[table.permutation[i]]
end

@inline function lookup_perm(table::CodonGraphReference{T}, i::T, j::T) where T <: Codon
    @inbounds return table.permutation[codons_to_k(i, j)]
end

@inline function lookup_order(table::CodonGraphReference{T}, i::T, j::T) where T <: Codon
    @inbounds return table.order[codons_to_k(i, j)]
end

const DEFAULT_CODON_GRAPH_REFERENCE = CodonGraphReference{Codon{DNA}}(DEFAULT_TRANS)
