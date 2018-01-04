
function make_S_N_NG86_table(code::GeneticCode)
    f = i -> S_N_NG86(i, code)
    return S_N_NG86_LOOKUP(f)
end

function make_DS_DN_NG86_table(code::GeneticCode)
    f = (i, j) -> DS_DN_NG86(i, j, code)
    return DS_DN_NG86_LOOKUP(f)
end

isinteractive() && info("Compiling lookup tables for NG86...")
const S_N_NG86_LOOKUPS = (function ()
    lookups = Dict{Int, S_N_NG86_LOOKUP}()
    for (i, code) in ncbi_trans_table.tables
        lookups[i] = make_S_N_NG86_table(code)
    end
    return lookups
end)()

const DS_DN_NG86_LOOKUPS = (function ()
    lookups = Dict{Int, DS_DN_NG86_LOOKUP}()
    for (i, code) in ncbi_trans_table.tables
        lookups[i] = make_DS_DN_NG86_table(code)
    end
    return lookups
end)()
