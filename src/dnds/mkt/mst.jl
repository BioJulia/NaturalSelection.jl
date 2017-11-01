struct MSTState{C <: Codon}
    parents::Vector{C}
    ranks::Vector{Int}
end

function MSTState{C}() where C <: Codon
    MSTState{C}(collect((C(i) for i in UInt64(0):UInt64(63))), zeros(Int, 64))
end

@inline function reset!(x::MSTState{C}) where C <: Codon
    fill!(x.ranks, 0)
    copy!(x.parents, (C(i) for i in UInt64(0):UInt64(63)))
end

cdn2i(c::Codon) = UInt64(c) + 1

@inline function parentof(vertex::C, state::MSTState{C}) where C <: Codon
    @inbounds return state.parents[cdn2i(vertex)]
end

@inline function setparentof(vertex::C, parent::C, state::MSTState{C}) where C <: Codon
    @inbounds state.parents[cdn2i(vertex)] = parent
end

@inline function rankof(vertex::C, state::MSTState{C}) where C <: Codon
    @inbounds return state.rank[cdn2i(vertex)]
end

@inline function incrank(vertex::C, state::MSTState{C}) where C <: Codon
    state.rank[cdn2i(vertex)] += 1
end

function find(vertex::C, state::MSTState{C}) where C <: Codon
    p = parentof(vertex, state)
    if p != vertex
        setparentof(vertex, find(p, state), state)
    end
    return parentof(vertex, state)
end

function union(v1::C, v2::C, state::MSTState{C}) where C <: Codon
    r1 = find(v1, p)
    r2 = find(v2, p)
    if r1 != r2
        if r[r1] > r[r2]
            p[r2] = r1
        else
            p[r1] = r2
            if r[r1] == r[r2]
                r[r2] += 1
            end
        end
    end
end

function union(vertexa::C, vertexb::C, state::MSTState{C}) where C <: Codon
    roota = find(vertexa, state)
    rootb = find(vertexb, state)
    if roota != rootb
        ranka = rankof(roota, state)
        rankb = rankof(rootb, state)
        if ranka > rankb
            setparentof(rootb, roota, state)
        else
            setparentof(roota, rootb, state)
            if ranka == rankb
                incrank(rootb, state)
            end
        end
    end
end

function mst(codongraph::CodonGraph, state::MSTState)
    DS = DN = 0
    reset!(state)
    for (v1, v2, DS_i, DN_i) in codongraph
        if find(v1, state) != find(v2, state)
            println("Edge: ", v1, " --> ", v2, " added to tree.")
            println("Edge DS cost: ", DS_i, ". Edge DN cost: ", DN_i)
            union(v1, v2, parents, ranks)
            DS += DS_i
            DN += DN_i
        end
    end
    return DS, DN
end
