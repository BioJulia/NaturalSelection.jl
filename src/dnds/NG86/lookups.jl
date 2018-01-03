
function make_S_N_NG86_table(code::GeneticCode, k::Float64 = 1.0)
    info("Setting up S_N_NG86")
    table = CodonLookupTable{1, Tuple{Float64, Float64}}()
    f = i -> S_N_NG86(i, k, code)
    setuplookup!(table, f)
    return table
end

function make_DS_DN_NG86_table(code::GeneticCode)
    info("Setting up DS_DN_NG86")
    table = CodonLookupTable{2, Tuple{Float64, Float64}}()
    f = (i, j) -> DS_DN_NG86(i, j, code)
    setuplookup!(table, f)
    return table
end

info("Compiling default lookup tables for NG86 dN dS...")
const DEFAULT_S_N_NG86_LOOKUP = make_S_N_NG86_table(DEFAULT_TRANS, 1.0)
const DEFAULT_DS_DN_NG86_LOOKUP = make_DS_DN_NG86_table(DEFAULT_TRANS)
