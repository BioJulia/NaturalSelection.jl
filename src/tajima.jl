
function harmonic(n::Integer)
    H = 0.0
    for i in n:-1:1
        H += 1 / i
    end
    return H
end

function harmonic(n::Integer, m::Integer)
    H = 0.0
    for i in n:-1:1
        H += 1 / (i ^ m)
    end
    return H
end

td_a1(n::Integer) = harmonic(n - 1)

td_a2(n::Integer) = harmonic(n - 1, 2)

td_b1(n::Integer) = (n + 1) / (3 * (n - 1))

function td_b2(n::Integer)
    return (2 * ((n ^ 2) + n + 3)) / ((9 * n) * (n - 1))
end

td_c1(a1::AbstractFloat, b1::AbstractFloat) = b1 - 1 / a1

function td_c2(n::Integer, a1::AbstractFloat, a2::AbstractFloat, b2::AbstractFloat)
    return b2 - ((n + 2) / (a1 * n)) + a2 / (a1 ^ 2)
end

td_e1(a1::AbstractFloat, c1::AbstractFloat) = c1 / a1

td_e2(a1::AbstractFloat, a2::AbstractFloat, c2::AbstractFloat) = c2 / (a1 ^ 2 + a2)

function tajimad(π::AbstractFloat, S::Integer, a1::AbstractFloat, e1::AbstractFloat, e2::AbstractFloat)
    return (π - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))
end

"""
    tajimad(π::AbstractFloat, S::Integer, n::Integer)

Compute Tajima's D from:

* `π`: The average number of SNPs found in (n choose 2) pairwise comparisons of
       a sample of sequences.

* `S`: The number of segregating sites in a sample of sequences.

* `n`: The number of sequences in your sample.

*Example*

```julia
tajimad(3.88888, 16, 10)
```

"""
function tajimad(π::AbstractFloat, S::Integer, n::Integer)
    a1 = td_a1(n)
    a2 = td_a2(n)
    b1 = td_b1(n)
    b2 = td_b2(n)
    c1 = td_c1(a1, b1)
    c2 = td_c2(n, a1, a2, b2)
    e1 = td_e1(a1, c1)
    e2 = td_e2(a1, a2, c2)
    return tajimad(π, S, a1, e1, e2)
end

"""
    tajimad(seqs)

Compute Tajima's D from a collection of BioSequences{DNAAlphabet{n}} (n = 2 or 4).

This will estimate the `π`, `S`, and `n` parameters from the sequences and use
those parameters to estimate Tajima's D.

*Example*

```julia

sample = [dna"ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA",
          dna"AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
          dna"AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA",
          dna"AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA",
          dna"AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA",
          dna"AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA",
          dna"AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
          dna"AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA",
          dna"AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA",
          dna"AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA"]

tajimad(sample)
```

"""
function tajimad(seqs)
    π = avg_mut(seqs)
    S = count(Segregating, seqs)
    n = length(seqs)
    return tajimad(π, S[1], n)
end
