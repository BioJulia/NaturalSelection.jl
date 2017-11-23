# allpaths.jl
# ===========
#
# A DS_DN_Counter for use with the DS_DN_enumerator in DS_DN_enumerator.jl
# Counts DS and DN according to the all-paths method devised for use with
# Nei and Gojoborei's method of estimating dN and dS ratios.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md

struct AllPaths <: DS_DN_Counter end
@inline DS_DN_init(::Type{AllPaths}) = 0.0, 0.0
@inline weighting(::Type{AllPaths}, ::Val{1}) = 1.0
@inline weighting(::Type{AllPaths}, ::Val{2}) = 0.5
@inline weighting(::Type{AllPaths}, ::Val{3}) = 0.5 / 3
@inline function accumulate(::Type{AllPaths}, DS::Float64, DN::Float64, DS_i::Float64, DN_i::Float64)
    return DS + DS_i, DN + DN_i
end
