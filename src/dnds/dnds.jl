# dnds.jl
# =======
#
# dNdS statistics computation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

const CDN_POS_MASKS = (0xFFFFFFFFFFFFFFCF, 0xFFFFFFFFFFFFFFF3, 0xFFFFFFFFFFFFFFFC)
const SITE_PERMUTATIONS = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

function aligned_codons(x::BioSequence{T}, y::BioSequence{T}, start::Int = 1) where T <: NucAlphs
    xcdns = Vector{Kmer{eltype(T), 3}}()
    ycdns = Vector{Kmer{eltype(T), 3}}()
    pos = start
    while pos + 2 ≤ min(endof(x), endof(y))
        cdnx, okx = BioSequences.extract_kmer_impl(x, pos, 3)
        cdny, oky = BioSequences.extract_kmer_impl(y, pos, 3)
        if okx && oky
            push!(xcdns, convert(Kmer{eltype(T), 3}, cdnx))
            push!(ycdns, convert(Kmer{eltype(T), 3}, cdny))
        end
        pos += 3
    end
    return xcdns, ycdns
end

function aligned_codons(v::Vector{DNASequence}, start::Int = 1)
    nv = endof(v)
    cdns = Vector{Vector{Codon{DNA}}}(nv)
    oks = Vector{Vector{Bool}}(nv)
    @inbounds for i in eachindex(v)
        cdns[i] = Vector{Codon{DNA}}()
        oks[i] = Vector{Bool}()

        pos = start
        while pos + 2 ≤ minimum(endof.(v))
            cdn, ok = BioSequences.extract_kmer_impl(v[i], pos, 3)
            push!(cdns[i], convert(Codon{DNA}, cdn))
            push!(oks[i], ok)
            pos += 3
        end
    end
    finalbool = falses(oks[1])
    @inbounds for ok in oks
        for i in eachindex(ok)
            finalbool[i] |= ok[i]
        end
    end
    return Vector{Codon{DNA}}[cdn[finalbool] for cdn in cdns]
end

include("enumerator/enumerator.jl")
include("enumerator/allpaths.jl")
include("enumerator/shortestpath.jl")

include("lookups/lookups.jl")

include("NG86/NG86.jl")

include("mkt/mkt.jl")
