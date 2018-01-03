
# This is an internal type to abstract codon lookup tables for various dN/dS
# related operations. It has many similarities with the PairwiseListMatrix type
# from the PairwiseListMatrices.jl package, but it does not support the storage
# of diagonal values and many normal matrix operations, as it is simply used for
# lookup of precomputed values, and not for any other purpose.

struct CodonLookupTable{N, T}
    table::Vector{T}
end

function CodonLookupTable{1, T}() where T
    return CodonLookupTable{1, T}(Vector{T}(64))
end

function CodonLookupTable{2, T}() where T
    return CodonLookupTable{2, T}(Vector{T}(2016))
end

include("indexing.jl")

function setuplookup!(tbl::CodonLookupTable{1, T}, f::Function) where T
    cdn = kmer"AAA"
    @inbounds while cdn < kmer"TTT"
        tbl[cdn] = f(cdn)
        cdn += 1
    end
    tbl[cdn] = f(cdn)
end

function setuplookup!(tbl::CodonLookupTable{2, T}, f::Function) where T
    ci = kmer"AAA"
    while ci < kmer"TTT"
        cj = ci + 1
        @inbounds while cj < kmer"TTT"
            tbl[ci, cj] = f(ci, cj)
            cj += 1
        end
        tbl[ci, cj] = f(ci, cj)
        ci += 1
    end
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
