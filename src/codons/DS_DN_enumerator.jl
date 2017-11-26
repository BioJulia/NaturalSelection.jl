# DS_DN.jl
# ========
#
# Parameterized procedure for identifying and calculating the number of
# synonymous and nonsynonymous mutations between codons, accounting for the
# different possible pathways of substitution.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

abstract type DS_DN_Counter end

const ONEPATH = Val{1}()
const TWOPATHS = Val{2}()
const THREEPATHS = Val{3}()

function splice_into(x::C, y::C, pos::Integer) where C <: Codon
    mask = UInt64(3) << bitindex(x, pos)
    return C((UInt64(x) & ~mask) | (UInt64(y) & mask))
end

function find_differences(x::Codon, y::Codon)
    diffs = 0x00
    @inbounds for pos in 1:3
        diffs = (x[pos] != y[pos]) | (diffs << 1)
    end
    return diffs, count_ones(diffs)
end

@inline function classify_mutation(x::Codon, y::Codon, code::GeneticCode, weight::N = N(1)) where N <: Real
    @inbounds aresame = code[x] == code[y]
    DS = ifelse(aresame, weight, N(0))
    DN = ifelse(aresame, N(0), weight)
    return DS, DN
end

function DS_DN_enumerator(::Type{T}, x::Codon, y::Codon, code::GeneticCode = DEFAULT_TRANS) where T<:DS_DN_Counter
    if x == y # Early escape, codons are the same, no syn or nonsyn mutations.
        return 0.0, 0.0
    else
        diff_positions, n_diffs = find_differences(x, y) # Which positions are different.
        if n_diffs == 1
            DS, DN = classify_mutation(x, y, code, weighting(T, ONEPATH))
            # One site in the two codons is different. It is obvious and simple
            # then to count whether it is a synonymous or nonsynonymous mutation.
        elseif n_diffs == 2
            DS, DN = DS_DN_init(T)
            # For two changes, the number of synonymous and non-synonymous
            # differences per codon, sum to 2, there are two pathways,
            # each possible pathway having two steps.
            # For example, comparing CTA and GTT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            # 2: CTA (L) -> CTT (L) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            @inbounds for pos in 1:3
                if ((diff_positions >> (3 - pos)) & 0x01) == 0x01
                    temp_cdn = splice_into(x, y, pos)

                    DS_a, DN_a = classify_mutation(x, temp_cdn, code, weighting(T, TWOPATHS))
                    DS_b, DN_b = classify_mutation(temp_cdn, y, code, weighting(T, TWOPATHS))
                    DS_i = DS_a + DS_b
                    DN_i = DN_a + DN_b

                    DS, DN = accumulate(T, DS, DN, DS_i, DN_i)

                end
            end
        elseif n_diffs == 3
            DS, DN = DS_DN_init(T)
            # For two changes, there are 6 pathways, each with three steps.
            # For example, comparing CTA and GAT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 2: CTA (L) -> GTA (V) -> GTT (V) -> GAT (D) : 2 nonsynonymous and 1 synonymous change.
            # 3: CTA (L) -> CAA (Q) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 4: CTA (L) -> CAA (Q) -> CAT (H) -> GAT (D) : 3 nonsynonymous changes.
            # 5: CTA (L) -> CTT (L) -> GTT (V) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.
            # 6: CTA (L) -> CTT (L) -> CAT (H) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.
            @inbounds for path in SITE_PERMUTATIONS
                tmp_cdn_a = splice_into(x, y, path[1])
                tmp_cdn_b = splice_into(tmp_cdn_a, y, path[2])

                DS_a, DN_a = classify_mutation(x, tmp_cdn_a, code, weighting(T, THREEPATHS))
                DS_b, DN_b = classify_mutation(tmp_cdn_a, tmp_cdn_b, code, weighting(T, THREEPATHS))
                DS_c, DN_c = classify_mutation(tmp_cdn_b, y, code, weighting(T, THREEPATHS))
                DS_i = DS_a + DS_b + DS_c
                DN_i = DN_a + DN_b + DN_c

                DS, DN = accumulate(T, DS, DN, DS_i, DN_i)

            end
        end
        return DS, DN
    end
end
