# Loading and indexing of the nuclear mass tables shared by this directory.
#
# The original scripts accessed these tables with expressions of the form
#
#     df.D[(df.A .== A) .& (df.Z .== Z)][1]
#
# a full boolean scan of a ~3000-row table, allocating two temporary vectors,
# executed several times per nuclide inside a double loop over all nuclides.
# Radionuclizi_3.jl called a function containing six such scans six times per
# (A, Z) pair, which is of order 2×10⁹ element comparisons — and because the
# table was a non-`const` global, every one of those accesses was also
# type-unstable. A dictionary keyed on (Z, A) replaces all of it.

using CSV, DataFrames

"One entry of a mass evaluation. Mass excess and its uncertainty in keV."
struct Nuclide
    Z::Int
    A::Int
    symbol::String
    Δ::Float64
    σΔ::Float64
end

"""
    load_masses(path)

Read an AME-format table — `Z A symbol Δ σΔ`, whitespace separated, no header —
and return a `Dict{Tuple{Int,Int},Nuclide}` keyed on (Z, A).
"""
function load_masses(path)
    df = CSV.read(path, DataFrame; header = ["Z", "A", "symbol", "Δ", "σΔ"],
                  delim = ' ', ignorerepeated = true, silencewarnings = true)
    table = Dict{Tuple{Int,Int},Nuclide}()
    for row in eachrow(df)
        (ismissing(row.Z) || ismissing(row.A) || ismissing(row.Δ)) && continue
        table[(row.Z, row.A)] = Nuclide(row.Z, row.A, string(row.symbol),
                                        row.Δ, ismissing(row.σΔ) ? 0.0 : row.σΔ)
    end
    return table
end

"""
    load_moller(path)

Möller–Nix FRDM microscopic (shell-plus-pairing) correction, `Z A δW` in MeV.
The file is the only one in the archive with CRLF line endings and a trailing
space on every row; both are absorbed here rather than left as a latent parser
dependency.
"""
function load_moller(path)
    table = Dict{Tuple{Int,Int},Float64}()
    for line in eachline(path)
        parts = split(strip(line))
        length(parts) < 3 && continue
        table[(parse(Int, parts[1]), parse(Int, parts[2]))] = parse(Float64, parts[3])
    end
    return table
end

"""
Mass excess of ¹H — the hydrogen **atom**, which is the correct reference for
atomic mass excesses. The original code named this `Dᵖ` and called it the proton;
the value was right, the name was not.
"""
Δ_H1(t) = t[(1, 1)].Δ
"Mass excess of the free neutron."
Δ_n(t) = t[(0, 1)].Δ

"""
    binding_energy(t, Z, A)

Mean binding energy per nucleon in keV, or `nothing` if the nuclide is absent:

    B/A = [Z·Δ(¹H) + N·Δ(n) − Δ(A,Z)] / A
"""
function binding_energy(t, Z, A)
    haskey(t, (Z, A)) || return nothing
    return (Z * Δ_H1(t) + (A - Z) * Δ_n(t) - t[(Z, A)].Δ) / A
end

"""
    separation_energy(t, Z, A, Z_x, A_x)

Energy to remove the particle (Z_x, A_x) from the nuclide (Z, A), in keV, or
`nothing` if either the parent, the ejectile or the residual is absent. The original
version returned `[true, S]` on success and `[false, 0]` on failure — two
different element types from one function, then indexed as `[1]` for a boolean
and `[2]` for an energy.
"""
function separation_energy(t, Z, A, Z_x, A_x)
    (haskey(t, (Z, A)) && haskey(t, (Z_x, A_x)) &&
     haskey(t, (Z - Z_x, A - A_x))) || return nothing
    return t[(Z - Z_x, A - A_x)].Δ + t[(Z_x, A_x)].Δ - t[(Z, A)].Δ
end
