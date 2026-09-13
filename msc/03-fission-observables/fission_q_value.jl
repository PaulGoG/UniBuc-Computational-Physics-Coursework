# Q-value of ²³⁵U(n_th,f) from the AME mass excesses,
#
#   Q(A_H, Z_H) = Δ(²³⁶U) − Δ(A_H, Z_H) − Δ(A_L, Z_L)
#
# averaged over the isobaric charge distribution about the most probable charge
# Z_p(A) = Z_UCD + ΔZ with ΔZ = −0.5 and an rms width of 0.6.
#
# The mass range, the polarisation and the width are all set by the assignment
# (Tudora, 10 Oct 2022): "A de la 76 la 160 (adica AH de la 118 la 160)", three
# charges per A at the nearest integer either side of Z_p(A) = Z_UCD(A) + ΔZ(A),
# polarisation 0.5 with + for the light fragment and − for the heavy, and a
# Gaussian isobaric charge distribution of rms 0.6. ΔZ = −0.5 and rms = 0.6 are
# the yield-weighted averages over the studied actinides, from the Wahl Z_p
# model as tabulated in Wagemans, The Nuclear Fission Process, Fig. 85.
#
# Ported from Fisiune_1.jl. The physics was right. The implementation had four
# boolean re-filters of the whole array per inner iteration, parsed its CSV
# inside the calculation rather than taking it as an argument — unlike every
# sibling file — and used `round(Z_p)` as a loop bound, producing a float range
# whose values were then pushed into an `Int[]` field. It also omitted the
# `q > 0` guard that the identical code in Fisiune_2, _3 and _4 carries.

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_data.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Mass number of the compound system ²³⁶U."
const A₀ = 236
"Proton number of the compound system."
const Z₀ = 92
"Charge polarisation of the isobaric charge distribution."
const ΔZ = -0.5
"Rms width of the isobaric charge distribution."
const σ_Z = 0.6
"""
Heavy-fragment mass range, as set by the assignment: A from 76 to 160, that is
A_H from 118 to 160. ΔZ = −0.5 and rms = 0.6 are set there too, both as the
yield-weighted averages over the studied actinides.
"""
const A_H_RANGE = 118:160

"Unchanged-charge-density most probable charge, with polarisation."
Z_p(A_H) = Z₀ * A_H / A₀ + ΔZ

"Gaussian weight of charge Z about the most probable value."
charge_weight(Z, A_H) = exp(-(Z - Z_p(A_H))^2 / (2σ_Z^2))

function main()
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    Δ₀ = Δ(masses, Z₀, A₀)
    @printf("Δ(²³⁶U) = %.1f keV, from %d tabulated nuclides\n", Δ₀, length(masses))

    A_H_range = A_H_RANGE
    Q = Float64[]; A_kept = Int[]
    for A_H in A_H_range
        A_L = A₀ - A_H
        num = 0.0; den = 0.0
        # three charge splits about the most probable charge, as in Fisiune_1.jl
        for Z_H in (round(Int, Z_p(A_H)) - 1):(round(Int, Z_p(A_H)) + 1)
            Z_L = Z₀ - Z_H
            δH = Δ(masses, Z_H, A_H); δL = Δ(masses, Z_L, A_L)
            (δH === nothing || δL === nothing) && continue
            q = (Δ₀ - δH - δL) / 1000        # MeV
            q > 0 || continue
            w = charge_weight(Z_H, A_H)
            num += w * q; den += w
        end
        den > 0 || continue
        push!(A_kept, A_H); push!(Q, num / den)
    end

    @printf("Q(A_H) over %d mass splits: mean %.2f MeV, range %.2f to %.2f\n",
            length(Q), mean(Q), minimum(Q), maximum(Q))
    @printf("maximum at A_H = %d\n", A_kept[argmax(Q)])
    sym = findfirst(==(118), A_kept)
    sym !== nothing && @printf("symmetric split A_H = 118: Q = %.2f MeV\n", Q[sym])
    @printf("the textbook figure for ²³⁵U(n,f) is about 200 MeV\n")

    fig = Figure(size = (800, 470))
    ax = Axis(fig[2, 1], xlabel = L"Heavy fragment mass $A_H$",
        ylabel = L"$Q$ [MeV]")
    l = lines!(ax, A_kept, Q, color = PALETTE.blue, linewidth = 1.8)
    h = hlines!(ax, [mean(Q)], color = PALETTE.red, linestyle = :dash, linewidth = 1.2)
    text!(ax, 0.04, 0.06; text = @sprintf("mean %.1f MeV", mean(Q)),
        space = :relative, align = (:left, :bottom), color = PALETTE.red, fontsize = 15)

    Legend(fig[1, 1], [l, h],
        [L"$Q(A_H)$, charge-averaged", "Mean over mass splits"],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 22)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("wrote ", savefigure(fig, FIGURES, "fission_q_value"))
end

main()
