@testset "Tajima's D" begin
    n = (10, 77, 72)
    a1 = (2.828968, 4.914514, 4.846921)
    a2 = (1.539768, 1.631862, 1.630948)
    b1 = (0.407407, 0.342105, 0.342723)
    b2 = (0.279012, 0.228184, 0.228612)
    c1 = (0.053922, 0.138626, 0.136406)
    c2 = (0.047227, 0.086985, 0.085989)
    e1 = (0.019061, 0.028208, 0.028143)
    e2 = (0.004949, 0.003374, 0.003423)
    π = (3.888889, 8.438483, 15.339984)
    S = (16, 103, 88)
    D = (-1.446172, -2.021749, -0.525801)

    for i in eachindex(n)
        @test NaturalSelection.td_a1(n[i]) ≈ a1[i] atol=10e-5
        @test NaturalSelection.td_a2(n[i]) ≈ a2[i] atol=10e-5
        @test NaturalSelection.td_b1(n[i]) ≈ b1[i] atol=10e-5
        @test NaturalSelection.td_b2(n[i]) ≈ b2[i] atol=10e-5
        @test NaturalSelection.td_c1(a1[i], b1[i]) ≈ c1[i] atol=10e-5
        @test NaturalSelection.td_c2(n[i], a1[i], a2[i], b2[i]) ≈ c2[i] atol=10e-5
        @test NaturalSelection.td_e1(a1[i], c1[i]) ≈ e1[i] atol=10e-5
        @test NaturalSelection.td_e2(a1[i], a2[i], c2[i]) ≈ e2[i] atol=10e-5
        @test tajimad(π[i], S[i], n[i]) ≈ D[i] atol=10e-5
    end

    testpop = [dna"ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA",
               dna"AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
               dna"AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA",
               dna"AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA",
               dna"AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA",
               dna"AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA",
               dna"AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
               dna"AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA",
               dna"AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA",
               dna"AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA"]

    @test tajimad(testpop) ≈ -1.446172 atol=10e-5
end
