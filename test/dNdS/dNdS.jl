@testset "dNdS" begin
    @testset "NG86" begin
        codonsA = [kmer"ATG",
                   kmer"AAA",
                   kmer"CCC",
                   kmer"GGG",
                   kmer"TTT",
                   kmer"TAA",
                   kmer"GGG"]

        codonsB = [kmer"ATG",
                   kmer"AAA",
                   kmer"CGC",
                   kmer"GGC",
                   kmer"TAC",
                   kmer"TAA",
                   kmer"GGG"]

        include("NG86/computation.jl")

        @testset "dN/dS" begin
            @test dNdS_NG86(codonsA, codonsB)[1] ≈ 0.125 atol=0.001
            @test dNdS_NG86(codonsA, codonsB)[2] ≈ 0.974 atol=0.001
        end

    end
end
