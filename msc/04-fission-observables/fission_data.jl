# Loaders for the ²³⁵U(n_th,f) datasets used in this directory.
#
# All the originals accessed these tables by boolean masking inside nested
# loops. `Fisiune_2.jl:TKE_A` alone made roughly 8000 full passes over the
# 20 705-row yield table. Everything is indexed once here.

using CSV, DataFrames

"Straede triple-differential yield Y(A_H, Z_H, TKE) in % per fission."
function load_yields(path)
    CSV.read(path, DataFrame; header = ["A_H", "Z_H", "TKE", "Y", "σY"],
        skipto = 2, delim = ' ', ignorerepeated = true, silencewarnings = true,)
end

"AME mass excesses keyed on (Z, A), in keV."
function load_masses(path)
    df = CSV.read(path, DataFrame; header = ["Z", "A", "symbol", "Δ", "σΔ"],
        delim = ' ', ignorerepeated = true, silencewarnings = true,)
    Dict((r.Z, r.A) => Float64(r.Δ)
    for r in eachrow(df)
    if !ismissing(r.Z) && !ismissing(r.A) && !ismissing(r.Δ))
end

"Möller–Nix FRLDM ground-state β₂ keyed on (Z, A)."
function load_beta2(path)
    t = Dict{Tuple{Int,Int},Float64}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        p = split(strip(line))
        length(p) < 3 && continue
        t[(parse(Int, p[1]), parse(Int, p[2]))] = parse(Float64, p[3])
    end
    return t
end

"Gilbert–Cameron shell corrections S(N), S(Z) keyed on nucleon number."
function load_gilbert_cameron(path)
    S = Dict{Int,Tuple{Float64,Float64}}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        p = split(strip(line))
        length(p) < 3 && continue
        S[parse(Int, p[1])] = (parse(Float64, p[2]), parse(Float64, p[3]))
    end
    return S
end

"Two-column measurement file with a trailing header comment: x, y, σy."
function load_measurement(path)
    x = Float64[]
    y = Float64[]
    σ = Float64[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        p = split(strip(line))
        length(p) < 3 && continue
        push!(x, parse(Float64, p[1]))
        push!(y, parse(Float64, p[2]))
        push!(σ, parse(Float64, p[3]))
    end
    return (x = x, y = y, σ = σ)
end

"Mass excess lookup returning `nothing` rather than throwing."
Δ(masses, Z, A) = get(masses, (Z, A), nothing)
