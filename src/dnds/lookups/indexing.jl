
# Indicies conversions
# ====================

@inline function ij_to_k(i::Integer, j::Integer)
    i, j = ifelse(i > j, (j, i), (i, j))
    nelements = 64
    x = nelements - i
    return div(nelements * (nelements - 1) - (x * (x - 1)), 2) - nelements + j
end

@inline function codons_to_k(i::T, j::T) where T <: Codon
    return ij_to_k(UInt64(i) + 1, UInt64(j) + 1)
end


# Single index operations
# =======================

# Bounds checking
# ---------------

# Basic boundschecking for when indexing into a 1 dimensional lookup table with
# a single index.
@inline function Base.checkbounds(table::CodonLookupTable{1, T}, i::Integer) where T
    if i > 64 || i < 1
        throw(BoundsError(table, i))
    end
end

# Basic boundschecking for when indexing into a 2 dimensional lookup table with
# a single index.
@inline function Base.checkbounds(table::CodonLookupTable{2, T}, i::Integer) where T
    if i > 2016 || i < 1
        throw(BoundsError(table, i))
    end
end

# getindex overloads
# ------------------

# Indexing into both 1 and 2 dimensional lookup tables with a single integer index.
@inline function Base.getindex(table::CodonLookupTable{N, T}, i::Integer) where {N, T}
    @boundscheck checkbounds(table, i)
    @inbounds return x.table[i]
end

# Indexing into a 1-dimensional lookup table with a single codon.
@inline function Base.getindex(table::CodonLookupTable{1, T}, c::C) where {C <: Codon, T}
    # Do an unsafe indexing assuming inbounds, because you shouldn't ever have a
    # codon that is convertable to an out of bounds index.
    @inbounds return table[UInt64(c) + 1]
end

# setindex overloads
# ------------------

# Setting an element in both 1 and 2 dimensional lookup tables with a single index.
@inline function Base.setindex!(table::CodonLookupTable{N, T}, i::Integer, x::T) where {N, T}
    @boundscheck checkbounds(table, i)
    @inbounds table.table[i] = x
end

@inline function Base.setindex!(table::CodonLookupTable{1, T}, i::C, x::T) where {C <: Codon, T}
    # Do an unsafe indexing assuming inbounds, because you shouldn't ever have a
    # codon that is convertable to an out of bounds index.
    @inbounds table[UInt64(i) + 1] = x
end


# Two index operations
# ====================

# getindex overloads
# ------------------

@inline function Base.getindex(table::CodonLookupTable{2, T}, i::Integer, j::Integer) where T
    k = ij_to_k(i, j)
    @boundscheck checkbounds(table, k)
    @inbounds return table.table[k]
end

@inline function Base.getindex(table::CodonLookupTable{2, T}, i::C, j::C) where {C <: Codon, T}
    # Do an unsafe indexing assuming inbounds, because you shouldn't ever have
    # codons that are convertable to an out of bounds index.
    @inbounds return table[UInt64(i) + 1, UInt64(j) + 1]
end

# setindex overloads
# ------------------

@inline function Base.setindex!(table::CodonLookupTable{2, T}, i::Integer, j::Integer, x::T) where T
    k = ij_to_k(i, j)
    @boundscheck checkbounds(table, k)
    @inbounds table.table[k] = x
end

@inline function Base.setindex!(table::CodonLookupTable{2, T}, i::C, j::C, x::T) where {C <: Codon,T}
    # Do an unsafe indexing assuming inbounds, because you shouldn't ever have
    # codons that are convertable to an out of bounds index.
    @inbounds table[UInt64(i) + 1, UInt64(j) + 1] = x
end
