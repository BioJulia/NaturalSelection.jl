
const DSDN_RANK_LOOKUP = [0 2 5 9;
                          1 4 8 0;
                          3 7 0 0;
                          6 0 0 0;]

@inline function rankof(DS::Integer, DN::Integer)
    @inbounds return DSDN_RANK_LOOKUP[DS + 1, DN + 1]
end

function make_rank_lookup(code::GeneticCode)
    table = Vector{Int}(2016)
    k = 1
    for i in 0x00:0x3F
        for j in (i + 1):0x3F
            x = BioSequences.DNACodon(UInt64(i))
            y = BioSequences.DNACodon(UInt64(j))
            DS, DN = DS_DN_enumerator(ShortestPath, x, y, code)
            table[k] = rankof(Integer(DS), Integer(DN))
            k += 1
        end
    end
    return table
end

const DEFAULT_RANK_TABLE = make_rank_lookup(DEFAULT_TRANS)

@inline function ij_to_k(i, j)
    nelements = 64
    x = nelements - i
    return div(nelements * (nelements-1) - (x * (x - 1)), 2) - nelements + j
end

@inline function codons_to_k(i::T, j::T) where T <: Codon
    return ij_to_k(UInt64(i) + 1, UInt64(j) + 1)
end

@inline function lookup(table, i::T, j::T) where T <: Codon
    @inbounds return table[codons_to_k(i, j)]
end
