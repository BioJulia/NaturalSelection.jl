# lookup_tables.jl
# ================
#
# dNdS computation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

struct DNDSLookup
    tbl::PairwiseListMatrix{Tuple{Float64,Float64},false,Vector{Tuple{Float64,Float64}}}
    DNDSLookup() = new(PairwiseListMatrix(Tuple{Float64,Float64}, 64, false, (0.0,0.0)))
end








function dnds_lookup_generator()
    table = PairwiseListMatrix(Tuple{Float64,Float64}, 64, false, (0.0,0.0))

    for i in 0x00:0x3F

    end

end
