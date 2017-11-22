@testset "McDonald Kreitman" begin

    @testset "mkt" begin
        testseqs = [[dna"AAAGGGCCC", dna"AGAGGGCGC"],
                    [dna"AACGGGCCG", dna"ACAGGGCCC"]]

        @test mkt(aligned_codons(NaturalSelection.testseqs[1]),
                  aligned_codons(NaturalSelection.testseqs[2])) == (2, 3, 0, 1, 1.0)
    end

end
