# Prompt-neutron multiplicity per fragment, ν(A), from the total excitation
# energy, compared against four measured datasets.
#
# The energy balance is the one `Fisiune_3.jl` used:
#
#   ν_pair = (TXE − q) / (⟨ε⟩ + ⟨S_n⟩ + p),     ⟨ε⟩ = 4T_m/3,  T_m = √(TXE/a_tot)
#   p = 6.71 − Z²·0.156/A       q = 0.75 + Z²·0.088/A
#
# with p and q the systematics terms of the compound system and a_tot = a_L + a_H
# the summed Gilbert–Cameron level-density parameters,
# a = A[0.00917(S_Z + S_N) + 0.142].
#
# The pair multiplicity is split between the fragments by the ratio of their
# level-density parameters. The original obtained that ratio from its
# scission-point deformation energies instead; that model is not reproduced here
# (see below), and the level-density ratio is the statistical-equilibrium limit
# of the same quantity.
#
# Ported from the tractable part of Fisiune_3.jl, whose scission-point
# deformation-energy model is not reproduced — it hardcodes six undocumented
# parametrisation constants (`a = 0.58/13`, `b = -28*a`, `-0.58/6`, `0.6/15`,
# `-50*a`) with no source, and those are the most model-dependent numbers in the
# entire archive.
#
# The defect that mattered in that file is fixed here regardless:
#
#     for Z_H in minimum(dGC.n):maximum(dGC.n)
#
# bounded the loop over *fragment proton number* by the index column of the
# Gilbert–Cameron shell-correction table, so Z_H ran from 11 to 150.
# Unphysical pairs such as (A_H = 130, Z_H = 100) survived the guards and
# polluted the level-density distribution and everything averaged over it. The
# charge range is now derived from the yield matrix itself.

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_data.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")
const A₀, Z₀ = 236, 92

"Gilbert–Cameron level-density parameter in MeV⁻¹."
function level_density(A, Z, gc)
    haskey(gc, Z) && haskey(gc, A - Z) || return nothing
    S_Z = gc[Z][2]; S_N = gc[A - Z][1]
    return A * (0.00917 * (S_Z + S_N) + 0.142)
end

"Systematics terms of the compound system, as in Fisiune_3.jl."
p_term(A, Z) = 6.71 - Z^2 * 0.156 / A
q_term(A, Z) = 0.75 + Z^2 * 0.088 / A

function main()
    y = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    gc = load_gilbert_cameron(joinpath(DATA, "Parametrizari_auxiliare", "SZSN.GC"))

    Δ₀ = Δ(masses, Z₀, A₀)
    S_n_compound = (Δ(masses, 92, 235) + Δ(masses, 0, 1) - Δ₀) / 1000

    # charge range from the data, not from the shell-correction table index
    @printf("Z_H present in the yield matrix: %d–%d\n", minimum(y.Z_H), maximum(y.Z_H))
    @printf("the 2018 loop bound ran Z_H from %d to %d, the index column of SZSN.GC\n\n",
            minimum(keys(gc)), maximum(keys(gc)))

    A_list = Int[]; ν_H = Float64[]; ν_L = Float64[]; ν_tot = Float64[]
    for a in sort(unique(y.A_H))
        sub = y[y.A_H .== a, :]
        sum(sub.Y) > 0 || continue
        tke = sum(sub.TKE .* sub.Y) / sum(sub.Y)
        Zp = round(Int, Z₀ * a / A₀ - 0.5)
        aL = A₀ - a; ZL = Z₀ - Zp
        δH = Δ(masses, Zp, a); δL = Δ(masses, ZL, aL)
        (δH === nothing || δL === nothing) && continue
        TXE = (Δ₀ - δH - δL) / 1000 + S_n_compound - tke
        TXE > 0 || continue

        aH_ld = level_density(a, Zp, gc); aL_ld = level_density(aL, ZL, gc)
        (aH_ld === nothing || aL_ld === nothing) && continue
        (aH_ld > 0 && aL_ld > 0) || continue

        # pair multiplicity, exactly the Fisiune_3.jl energy balance
        a_tot = aH_ld + aL_ld
        T_m = sqrt(TXE / a_tot)
        ε_mean = 4 * T_m / 3
        S_H = separation_cost(masses, Zp, a); S_L = separation_cost(masses, ZL, aL)
        (S_H === nothing || S_L === nothing) && continue
        S_mean = (S_H + S_L) / 2
        pp = p_term(A₀, Z₀); qq = q_term(A₀, Z₀)
        TXE > qq || continue
        ν_pair = (TXE - qq) / (ε_mean + S_mean + pp)

        # split by the level-density ratio, the statistical-equilibrium limit
        R = aH_ld / a_tot
        push!(A_list, a); push!(ν_H, ν_pair * R); push!(ν_L, ν_pair * (1 - R))
        push!(ν_tot, ν_pair)
    end

    @printf("ν computed for %d mass splits\n", length(A_list))
    @printf("systematics terms for ²³⁶U: p = %.3f MeV, q = %.3f MeV\n",
            p_term(A₀, Z₀), q_term(A₀, Z₀))
    @printf("mean total multiplicity <nu_pair> = %.3f\n", mean(ν_tot))
    @printf("evaluated value for ²³⁵U(n_th,f): 2.42 — the model is %+.0f %% off\n",
            100 * (mean(ν_tot) / 2.42 - 1))
    @printf("the residual is the model: prompt γ emission competes for the same\n")
    @printf("excitation energy and carries off 6–7 MeV per fission, which this\n")
    @printf("balance does not account for.\n\n")

    sets = [("Göök", "U5NUAGOOK.DAT", PALETTE.blue),
            ("Maslin", "U5NUAMASLIN.DAT", PALETTE.orange),
            ("Nishio", "U5NUANISHIO.DAT", PALETTE.green),
            ("Vorobyev", "U5NUAVORO.DAT", PALETTE.purple)]
    measured = map(sets) do (name, file, colour)
        d = load_measurement(joinpath(DATA, "Date_experimentale", "Multiplicitate_n", file))
        keep = d.y .> 0
        @printf("%-9s %3d points with ν > 0, A %d–%d, mean ν = %.3f\n",
                name, count(keep), Int(minimum(d.x[keep])), Int(maximum(d.x[keep])),
                mean(d.y[keep]))
        (name = name, A = d.x[keep], ν = d.y[keep], σ = d.σ[keep], colour = colour)
    end

    fig = Figure(size = (900, 480))
    ax = Axis(fig[2, 1], xlabel = L"Fragment mass $A$", ylabel = L"\nu(A)")
    handles = []
    for m in measured
        push!(handles, scatter!(ax, m.A, m.ν, color = (m.colour, 0.75), markersize = 5))
    end
    l_model = lines!(ax, A_list, ν_H, color = PALETTE.red, linewidth = 2)
    lines!(ax, A₀ .- A_list, ν_L, color = PALETTE.red, linewidth = 2, linestyle = :dash)
    xlims!(ax, 70, 170); ylims!(ax, 0, 5)

    Legend(fig[1, 1], [handles; l_model],
        [[m.name for m in measured]; "Energy-balance model"],
        orientation = :horizontal, framevisible = false, labelsize = 15, colgap = 18)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_multiplicity"))
end

"Neutron separation energy of the fragment, MeV."
function separation_cost(masses, Z, A)
    d1 = Δ(masses, Z, A); d2 = Δ(masses, Z, A - 1); dn = Δ(masses, 0, 1)
    (d1 === nothing || d2 === nothing) && return nothing
    return (d2 + dn - d1) / 1000
end

main()
