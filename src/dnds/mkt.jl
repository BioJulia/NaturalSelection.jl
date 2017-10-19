# mkt.jl
# ======
#
# An implementation of the McDonald Kreitman test.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function codondiff(sequences::Vector{Vector{DNACodon}}, position::Integer)
    s = CodonSet{DNA}()
    @inbounds for seq in sequences
        s |= seq[position]
    end
    return s
end

function codondiff(sequences::Vector{Vector{DNACodon}})
    out = Vector{CodonSet{DNA}}(length(sequences[1]))
    @inbounds for site in eachindex(sequences[1])
        out[site] = codondiff(sequences, site)
    end
    return out
end

function multi_short_path(codonset)
    for i in 1:length(codonset)

    end
end

function short_path(cdnx, cdny)

end
