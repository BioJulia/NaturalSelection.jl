@testset "McDonald Kreitman" begin

    @testset "mkt" begin
        testseqs = [[dna"AAAGGGCCC", dna"AGAGGGCGC"],
                    [dna"AACGGGCCG", dna"ACAGGGCCC"]]

        @test mkt(NaturalSelection.aligned_codons(testseqs[1]),
                  NaturalSelection.aligned_codons(testseqs[2])) == (2, 3, 0, 1, 1.0)
    end

end
