# codon_set.jl
# ============
#
# A tiny immutable set type for DNA and RNA codons.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/NaturalSelection.jl/blob/master/LICENSE.md


# Type definition
# ---------------

"""
CodonSet is a primitive type that represents a set of all possible 64 codons.

It is useful to compactly represent such a set in a single 64 bit word, when
one is working with variations between biological sequences. 
"""
primitive type CodonSet{T<:NucleicAcid} 64 end


# Conversion and construction
# ---------------------------

CodonSet{T}() where T <: NucleicAcid = reinterpret(CodonSet{T}, UInt64(0))


# Base operators
# --------------

@inline function Base.:|(x::CodonSet{T}, y::UInt64) where T<:NucleicAcid
    return reinterpret(CodonSet{T}, reinterpret(UInt64, x) | y)
end

@inline function Base.:|(x::CodonSet{T}, y::Kmer{T,3}) where T<:NucleicAcid
    return x | UInt64(1) << UInt64(y)
end

@inline function Base.:&(x::CodonSet{T}, y::CodonSet{T}) where T <: NucleicAcid
    return reinterpret(CodonSet{T}, reinterpret(UInt64, x) & reinterpret(UInt64, y))
end


# Iteration interface
# -------------------

@inline function getnext(x::CodonSet, start::Integer)
    return trailing_zeros(reinterpret(UInt64, x) >> start) + start
end

Base.iteratorsize(::Type{CodonSet}) = HasLength()
Base.iteratoreltype(::Type{CodonSet}) = HasEltype()
Base.eltype(::Type{CodonSet{T}}) where T<:NucleicAcid = Kmer{T,3}
Base.length(x::CodonSet) = count_ones(reinterpret(UInt64, x))
Base.start(x::CodonSet) = getnext(x, UInt64(0))

@inline function Base.next(x::CodonSet{T}, state::UInt64) where T<:NucleicAcid
    return Kmer{T,3}(state), getnext(x, state + UInt64(1))
end

@inline function Base.done(x::CodonSet, state::UInt64)
    return state > UInt64(63)
end
