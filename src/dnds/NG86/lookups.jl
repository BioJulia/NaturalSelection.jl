
const S_N_NG86_LOOKUP = CodonLookupTable{1, Tuple{Float64, Float64}}
const DS_DN_NG86_LOOKUP = CodonLookupTable{2, Tuple{Float64, Float64}}

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

info("Compiling default lookup tables for NG86 dN dS...")
const S_N_NG86_LOOKUPS = (function ()
    tables = Dict{GeneticCode, S_N_NG86_LOOKUP}()
    for table in ncbi_trans_table.tables
        tables[table] = make_S_N_NG86_table(table)
    end
    return tables
end)()

const DS_DN_NG86_LOOKUPS = (function ()
    tables = Dict{GeneticCode, DS_DN_NG86_LOOKUP}()
    for table in ncbi_trans_table.tables
        tables[table] = make_DS_DN_NG86_table(table)
    end
    return tables
end)()
