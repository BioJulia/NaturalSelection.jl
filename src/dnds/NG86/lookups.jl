
function make_S_N_NG86_table(code::GeneticCode)
    table = CodonLookupTable{1, Tuple{Float64, Float64}}()
    f = i -> S_N_NG86(i, code)
    setuplookup!(table, f)
    return table
end

function make_DS_DN_NG86_table(code::GeneticCode)
    table = CodonLookupTable{2, Tuple{Float64, Float64}}()
    f = (i, j) -> DS_DN_NG86(i, j, code)
    setuplookup!(table, f)
    return table
end

isinteractive() && info("Compiling lookup tables for NG86...")
const S_N_NG86_LOOKUPS = (function ()
    tables = Dict{GeneticCode, S_N_NG86_LOOKUP}()
    for table in values(ncbi_trans_table.tables)
        tables[table] = make_S_N_NG86_table(table)
    end
    return tables
end)()

const DS_DN_NG86_LOOKUPS = (function ()
    tables = Dict{GeneticCode, DS_DN_NG86_LOOKUP}()
    for table in values(ncbi_trans_table.tables)
        tables[table] = make_DS_DN_NG86_table(table)
    end
    return tables
end)()
