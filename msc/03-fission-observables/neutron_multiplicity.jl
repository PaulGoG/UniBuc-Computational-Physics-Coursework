# Prompt-neutron multiplicity per fragment, ν(A), from the total excitation
# energy, compared against four measured datasets.
#
# The simplest energy-balance estimate divides the excitation energy available
# to a fragment by the cost of evaporating one neutron:
#
#   ν(A) ≈ E*(A) / (⟨S_n⟩ + ⟨ε⟩),        ⟨ε⟩ = 4T/3,  T = √(E*/a)
#
# with the level-density parameter from the Gilbert–Cameron systematics
# a = A[0.00917(S_Z + S_N) + 0.142]. The excitation energy is partitioned
# between the fragments in the ratio of their level-density parameters, which is
# the limit of statistical equilibrium at scission.
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

        E_H = TXE * aH_ld / (aH_ld + aL_ld)
        E_L = TXE * aL_ld / (aH_ld + aL_ld)
        function multiplicity(E, ald, Z, A)
            S = separation_cost(masses, Z, A)
            S === nothing && return nothing
            T = sqrt(E / ald)
            return E / (S + 4T/3)
        end
        nH = multiplicity(E_H, aH_ld, Zp, a)
        nL = multiplicity(E_L, aL_ld, ZL, aL)
        (nH === nothing || nL === nothing) && continue
        push!(A_list, a); push!(ν_H, nH); push!(ν_L, nL); push!(ν_tot, nH + nL)
    end

    @printf("ν computed for %d mass splits\n", length(A_list))
    @printf("mean total multiplicity <nu_pair> = %.3f\n", mean(ν_tot))
    @printf("evaluated value for ²³⁵U(n_th,f): 2.42 — the model is %.0f %% high\n",
            100 * (mean(ν_tot) / 2.42 - 1))
    @printf("that is the model, not a coding error: dividing the whole excitation\n")
    @printf("energy by the neutron cost ignores the competition with prompt γ\n")
    @printf("emission, which carries off roughly 6–7 MeV of TXE per fission, and\n")
    @printf("<ε> = 4T/3 understates the mean emitted neutron energy.\n\n")

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
