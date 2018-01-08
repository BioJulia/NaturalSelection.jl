# shortestpath.jl
# ===============
#
# A DS_DN_Counter for use with the DS_DN_enumerator in DS_DN_enumerator.jl
# Only keeps and returns the DS and DN values for 
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md


struct ShortestPath <: DS_DN_Counter end
@inline DS_DN_init(::Type{ShortestPath}) = 0.0, 3.0
@inline weighting(::Type{ShortestPath}, ::Val{T}) where T = 1.0
@inline function accumulate(::Type{ShortestPath}, DS::Float64, DN::Float64, DS_i::Float64, DN_i::Float64)
    switchout = DN_i < DN
    DN = ifelse(switchout, DN_i, DN)
    DS = ifelse(switchout, DS_i, DS)
    return DS, DN
end
