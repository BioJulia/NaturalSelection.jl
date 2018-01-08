@testset "dNdS" begin
    @testset "NG86" begin
        codonsA = [kmer"ATG",
                   kmer"AAA",
                   kmer"CCC",
                   kmer"GGG",
                   kmer"TTT",
                   kmer"GGG",
                   kmer"TAA"]

        codonsB = [kmer"ATG",
                   kmer"AAA",
                   kmer"CGC",
                   kmer"GGC",
                   kmer"TAC",
                   kmer"GGG",
                   kmer"TAA"]

        @testset "computation" begin
            @testset "expected" begin
                #all_codons = convert(Vector{Kmer{DNA,3}}, collect(UInt64(0):UInt64(63)))
                n_ans = [3.0, 2.666, 2.0, 2.0, 2.666, 2.0, 3.000]
                s_ans = [0.0, 0.333, 1.0, 1.0, 0.333, 1.0, 0.000]
                for i in 1:endof(codonsA)
                    cdn = codonsA[i]
                    @test S_N_NG86(cdn, ncbi_trans_table[1])[1] ≈ s_ans[i] atol=0.001
                    @test S_N_NG86(cdn, ncbi_trans_table[1])[2] ≈ n_ans[i] atol=0.001
                end
            end

            @testset "observed" begin
                function testobserved(a, b, ans)
                    @test DS_DN_NG86(a, b, ncbi_trans_table[1])[1] ≈ ans[1] atol=0.001
                    @test DS_DN_NG86(a, b, ncbi_trans_table[1])[2] ≈ ans[2] atol=0.001
                end
                answers = [(0.0, 0.0),
                           (0.0, 0.0),
                           (0.0, 1.0),
                           (1.0, 0.0),
                           (1.0, 1.0),
                           (0.0, 0.0),
                           (0.0, 0.0)]
                for i in 1:length(answers)
                    testobserved(codonsA[i], codonsB[i], answers[i])
                end
                testobserved(kmer"TTT", kmer"GAC", (1.0, 2.0))
            end
        end

        @testset "lookups" begin
            firstcdn = UInt64(0)
            lastcdn = UInt64(63)

            for i in firstcdn:lastcdn
                cdn = DNACodon(i)
                @test S_N_NG86(cdn, ncbi_trans_table[1])[1] == NaturalSelection.make_S_N_NG86_table(BioSequences.standard_genetic_code)[cdn][1]
                @test S_N_NG86(cdn, ncbi_trans_table[1])[2] == NaturalSelection.make_S_N_NG86_table(BioSequences.standard_genetic_code)[cdn][2]
            end

            for i in firstcdn:lastcdn
                for j in (i + 1):lastcdn
                    a = DNACodon(i)
                    b = DNACodon(j)
                    @test DS_DN_NG86(a, b, ncbi_trans_table[1])[1] == NaturalSelection.make_DS_DN_NG86_table(BioSequences.standard_genetic_code)[a, b][1]
                    @test DS_DN_NG86(a, b, ncbi_trans_table[1])[2] == NaturalSelection.make_DS_DN_NG86_table(BioSequences.standard_genetic_code)[a, b][2]
                end
            end
        end

        @testset "dN/dS" begin
            @test dNdS_NG86(codonsA[1:end-1], codonsB[1:end-1])[1] ≈ 0.154 atol=0.001
            @test dNdS_NG86(dna"ATGAAACCCGGGTTTGGG", dna"ATGAAACGCGGCTACGGG")[1] ≈ 0.154 atol=0.001
            @test dNdS_NG86(codonsA[1:end-1], codonsB[1:end-1])[2] ≈ 0.974 atol=0.001
            @test dNdS_NG86(dna"ATGAAACCCGGGTTTGGG", dna"ATGAAACGCGGCTACGGG")[2] ≈ 0.974 atol=0.001
        end

    end
end
