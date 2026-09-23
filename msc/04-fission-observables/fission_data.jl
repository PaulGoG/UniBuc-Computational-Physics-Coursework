# Loaders for the ²³⁵U(n_th,f) tables in `data/`. Every file is whitespace
# separated with at most one header line, and is read once into a typed
# container; nothing downstream filters a table by boolean mask.

"One row of the yield matrix: heavy-fragment mass and charge, TKE [MeV], yield and its error [%]."
struct YieldRow
    A_H::Int
    Z_H::Int
    TKE::Int
    Y::Float64
    σY::Float64
end

"Fields of each non-empty line of `path` after the first `skip` lines."
function table_fields(path; skip = 0)
    rows = Vector{Vector{SubString{String}}}()
    for (i, line) in enumerate(eachline(path))
        i <= skip && continue
        fields = split(line)
        isempty(fields) || push!(rows, fields)
    end
    return rows
end

"""
    load_yields(path) -> Vector{YieldRow}

The yield matrix Y(A_H, Z_H, TKE) in % per fission, heavy fragment only.
"""
function load_yields(path)
    return [YieldRow(parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                parse(Float64, f[4]), parse(Float64, f[5]))
            for f in table_fields(path; skip = 1)]
end

"""
    load_masses(path) -> Dict{Tuple{Int,Int},Tuple{Float64,Float64}}

Atomic mass excess and its uncertainty in keV, keyed on `(Z, A)`.
"""
function load_masses(path)
    return Dict((parse(Int, f[1]), parse(Int, f[2])) =>
                    (parse(Float64, f[4]), parse(Float64, f[5]))
    for f in table_fields(path))
end

"Möller–Nix ground-state quadrupole deformation β₂ keyed on `(Z, A)`."
function load_beta2(path)
    return Dict((parse(Int, f[1]), parse(Int, f[2])) => parse(Float64, f[3])
    for f in table_fields(path; skip = 1))
end

"""
    load_gilbert_cameron(path) -> Dict{Int,@NamedTuple{S_N::Float64, S_Z::Float64}}

Gilbert–Cameron shell corrections in MeV keyed on nucleon number: `S_N` when the
key counts neutrons, `S_Z` when it counts protons.
"""
function load_gilbert_cameron(path)
    return Dict(parse(Int, f[1]) => (S_N = parse(Float64, f[2]), S_Z = parse(Float64, f[3]))
    for f in table_fields(path; skip = 1))
end

"Three-column measurement file with one header line: abscissa, value, uncertainty."
function load_measurement(path)
    fields = table_fields(path; skip = 1)
    return (x = [parse(Float64, f[1]) for f in fields],
        y = [parse(Float64, f[2]) for f in fields],
        σ = [parse(Float64, f[3]) for f in fields])
end
