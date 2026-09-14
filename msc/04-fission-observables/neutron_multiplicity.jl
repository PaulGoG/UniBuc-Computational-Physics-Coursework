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
# scission-point deformation energies instead; both routes are computed here and
# compared, and the level-density ratio is the statistical-equilibrium limit of
# the same quantity.
#
# Ported from Fisiune_3.jl.
#
# **A retraction.** An earlier version of this header called the six constants of
# that file's scission-point deformation model — `a = 0.58/13`, `b = -28*a`,
# `-0.58/6`, `0.6/15`, `-50*a` — "undocumented, with no source", and declined to
# reproduce the model on that ground. They were a table. Tudora's notes on the
# partition of TXE give β at scission as a piecewise-linear function of fragment
# charge through five breakpoints,
#
#     Z    28    41    44    50    65
#     β     0   0.58  0.58    0    0.6
#
# and linear interpolation between them reproduces all six exactly: 0.58/13 and
# −28·(0.58/13) on 28→41, −0.58/6 on 44→50, and 0.6/15 with −50·(0.6/15) on
# 50→65. The parameterisation is implemented below as `scission_deformation` and
# plotted against the Möller–Nix ground-state β₂, as the assignment asks; the
# liquid-drop deformation energy from the same notes gives the alternative
# excitation-energy split, which is reported alongside the level-density one.
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

"Charge polarisation of the isobaric distribution, as in `fission_q_value.jl`."
const ΔZ_POL = -0.5
"Rms width of the isobaric charge distribution."
const σ_Z = 0.6

"""
    charge_averaged_q(masses, Δ₀, A_H)

Q value for the mass split `A_H`, averaged over the three charges nearest
Z_p(A) with a Gaussian weight of rms 0.6 — the treatment the assignment
specifies and that `fission_q_value.jl` and `fragment_yields.jl` both use. Using
the single most probable charge instead raises ⟨Q⟩, and so TXE, by about
0.5 MeV.
"""
function charge_averaged_q(masses, Δ₀, A_H)
    Zp = Z₀ * A_H / A₀ + ΔZ_POL
    centre = round(Int, Zp)
    num = 0.0; den = 0.0
    for z in (centre - 1, centre, centre + 1)
        δH = Δ(masses, z, A_H); δL = Δ(masses, Z₀ - z, A₀ - A_H)
        (δH === nothing || δL === nothing) && continue
        w = exp(-(z - Zp)^2 / (2σ_Z^2))
        num += w * (Δ₀ - δH - δL) / 1000
        den += w
    end
    return den > 0 ? num / den : nothing
end

"Gilbert–Cameron level-density parameter in MeV⁻¹."
function level_density(A, Z, gc)
    haskey(gc, Z) && haskey(gc, A - Z) || return nothing
    S_Z = gc[Z][2]; S_N = gc[A - Z][1]
    return A * (0.00917 * (S_Z + S_N) + 0.142)
end

"Systematics terms of the compound system. Tudora gives both."
p_term(A, Z) = 6.71 - Z^2 * 0.156 / A
q_term(A, Z) = 0.75 + Z^2 * 0.088 / A

"""
Breakpoints of the scission deformation parameterisation: fragment charge and
the corresponding β. Linearly interpolated between, constant outside.
"""
const SCISSION_BETA = [(28, 0.0), (41, 0.58), (44, 0.58), (50, 0.0), (65, 0.6)]

"""
    scission_deformation(Z)

Quadrupole deformation β of a fission fragment of charge `Z` at scission, by
linear interpolation of [`SCISSION_BETA`](@ref).

The two zeros sit at the closed shells: Z = 28 and Z = 50 are spherical at
scission, and deformation peaks in between. This is the parameterisation the
original file encoded as slopes and intercepts.
"""
function scission_deformation(Z)
    Z <= first(SCISSION_BETA)[1] && return first(SCISSION_BETA)[2]
    Z >= last(SCISSION_BETA)[1] && return last(SCISSION_BETA)[2]
    for i in 1:(length(SCISSION_BETA) - 1)
        (z1, b1) = SCISSION_BETA[i]; (z2, b2) = SCISSION_BETA[i+1]
        z1 <= Z <= z2 && return b1 + (b2 - b1) * (Z - z1) / (z2 - z1)
    end
    return 0.0
end

"""
    ldm_energy(A, Z, β)

Liquid-drop energy in MeV of a nucleus deformed to β, as the sum of the
volume–surface and electrostatic terms,

    η  = (A − 2Z)/A,      x = 1 − 1.7826η²,      α₂ = 5β²/4π
    X_vs = −x[15.4941A − 17.9439A^{2/3}(1 + 0.4α₂)]
    X_e  = Z²[0.7053(1 − 0.2α₂)/A^{1/3} − 1.1529/A]

The 0.4 and 0.2 are 2/5 and 1/5, the leading quadrupole corrections to the
surface and Coulomb terms, so the α₂ of these expressions is the square of the
usual expansion coefficient.
"""
function ldm_energy(A, Z, β)
    η = (A - 2Z) / A
    x = 1 - 1.7826 * η^2
    α₂ = 5 * β^2 / (4π)
    X_vs = -x * (15.4941 * A - 17.9439 * A^(2/3) * (1 + 0.4 * α₂))
    X_e = Z^2 * (0.7053 * (1 - 0.2 * α₂) / A^(1/3) - 1.1529 / A)
    return X_vs + X_e
end

"""
    extra_deformation_energy(A, Z, β_gs)

Energy in MeV stored in deforming a fragment from its ground-state β to its
scission β, `E_LDM(β_sciss) − E_LDM(β_gs)`.
"""
extra_deformation_energy(A, Z, β_gs) =
    ldm_energy(A, Z, scission_deformation(Z)) - ldm_energy(A, Z, β_gs)

function main()
    y = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    gc = load_gilbert_cameron(joinpath(DATA, "Parametrizari_auxiliare", "SZSN.GC"))
    β_gs = load_beta2(joinpath(DATA, "Parametrizari_auxiliare", "B2MOLLER.ANA"))

    Δ₀ = Δ(masses, Z₀, A₀)
    S_n_compound = (Δ(masses, 92, 235) + Δ(masses, 0, 1) - Δ₀) / 1000

    # charge range from the data, not from the shell-correction table index
    @printf("Z_H present in the yield matrix: %d–%d\n", minimum(y.Z_H), maximum(y.Z_H))
    @printf("the original loop bound ran Z_H from %d to %d, the index column of SZSN.GC\n\n",
            minimum(keys(gc)), maximum(keys(gc)))

    A_list = Int[]; ν_H = Float64[]; ν_L = Float64[]; ν_tot = Float64[]
    yield_A = Float64[]                    # mass yield, the weight for every average
    R_ld = Float64[]; R_def = Float64[]
    for a in sort(unique(y.A_H))
        sub = y[y.A_H .== a, :]
        sum(sub.Y) > 0 || continue
        tke = sum(sub.TKE .* sub.Y) / sum(sub.Y)
        Zp = round(Int, Z₀ * a / A₀ - 0.5)
        aL = A₀ - a; ZL = Z₀ - Zp
        q = charge_averaged_q(masses, Δ₀, a)
        q === nothing && continue
        TXE = q + S_n_compound - tke
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
        push!(ν_tot, ν_pair); push!(yield_A, sum(sub.Y))

        # the deformation route, for comparison: the same statistical split of
        # what is left after each fragment pays for its own extra deformation
        bH = get(β_gs, (Zp, a), nothing); bL = get(β_gs, (ZL, aL), nothing)
        if bH !== nothing && bL !== nothing
            ΔH = extra_deformation_energy(a, Zp, bH)
            ΔL = extra_deformation_energy(aL, ZL, bL)
            E_rest = TXE - ΔH - ΔL
            if E_rest > 0
                r = aL_ld / aH_ld
                push!(R_ld, R)
                push!(R_def, (ΔH + E_rest / (1 + r)) / TXE)
            end
        end
    end

    @printf("ν computed for %d mass splits\n", length(A_list))
    @printf("systematics terms for ²³⁶U: p = %.3f MeV, q = %.3f MeV\n",
            p_term(A₀, Z₀), q_term(A₀, Z₀))
    ν_weighted = sum(ν_tot .* yield_A) / sum(yield_A)
    @printf("mean total multiplicity <nu_pair> = %.3f, yield-weighted\n", ν_weighted)
    @printf("  unweighted over mass splits      = %.3f\n", mean(ν_tot))
    @printf("  the two differ because the symmetric splits carry TXE some 3.5 MeV\n")
    @printf("  above the yield-weighted mean and a yield around 10^-4 of the peak,\n")
    @printf("  so counting them equally inflates the average.\n")
    @printf("evaluated value for ²³⁵U(n_th,f): 2.42 — the model is %+.0f %% off\n",
            100 * (ν_weighted / 2.42 - 1))
    @printf("the residual is the model: prompt γ emission competes for the same\n")
    @printf("excitation energy and carries off 6–7 MeV per fission, which this\n")
    @printf("balance does not account for.\n\n")

    @printf("excitation-energy share of the heavy fragment, R = E*_H/TXE:\n")
    @printf("  level-density ratio alone      <R> = %.3f\n", mean(R_ld))
    @printf("  with the scission deformation  <R> = %.3f   (%d splits)\n",
            mean(R_def), length(R_def))
    @printf("They differ by %.3f, or %.0f %% — paying for the extra deformation\n",
            abs(mean(R_def) - mean(R_ld)),
            100 * abs(mean(R_def) - mean(R_ld)) / mean(R_ld))
    @printf("first moves the heavy fragment from just above half the excitation\n")
    @printf("energy to just below it. That is not negligible. ν(A) below is still\n")
    @printf("computed from the level-density ratio alone, because the deformation\n")
    @printf("route carries a parameterised β(Z) and a liquid-drop energy on top of\n")
    @printf("it, and there is nothing here to validate the extra model against.\n")
    @printf("The comparison is reported rather than chosen between.\n\n")

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

    fig = Figure(size = (1250, 500))
    ax = Axis(fig[2, 1], xlabel = L"Fragment mass $A$", ylabel = L"\nu(A)")
    handles = []
    for m in measured
        push!(handles, scatter!(ax, m.A, m.ν, color = (m.colour, 0.75),
              markersize = MARKERSIZE.cloud + 2))
    end
    l_model = lines!(ax, A_list, ν_H, color = PALETTE.red, linewidth = 2)
    # The light-fragment branch is the same model in a different line style and
    # had no legend entry, so the dashed orange curve was unexplained.
    l_model_L = lines!(ax, A₀ .- A_list, ν_L, color = PALETTE.red, linewidth = 2,
        linestyle = :dash)
    xlims!(ax, 70, 170); ylims!(ax, 0, 3.9)
    text!(ax, 0.98, 0.97;
        text = rich("⟨", it("ν"), subscript("pair"), @sprintf("⟩ = %.3f", ν_weighted),
                    " against 2.42 evaluated"),
        space = :relative, align = (:right, :top), fontsize = 15, color = PALETTE.red)

    Legend(fig[1, 1], [handles; [l_model, l_model_L]],
        [[m.name for m in measured];
         ["Energy balance, heavy fragment", "Energy balance, light fragment"]],
        orientation = :horizontal, framevisible = false, labelsize = 15,
        nbanks = 2, colgap = 18)

    # The assignment asks for β at scission and β in the ground state against Z.
    # Plotting them is also the check on the parameterisation: the two zeros
    # must fall on the closed shells at Z = 28 and Z = 50.
    ax2 = Axis(fig[2, 2], xlabel = L"Fragment charge $Z$",
        ylabel = L"Quadrupole deformation $\beta_2$")
    Zs = sort(unique(vcat(y.Z_H, Z₀ .- y.Z_H)))
    Zgrid = range(minimum(Zs), maximum(Zs), length = 400)
    gs_pts = [(Z, b) for ((Z, A), b) in β_gs if Z in Zs]
    l_gs = scatter!(ax2, first.(gs_pts), last.(gs_pts),
        color = (PALETTE.sky, 0.35), markersize = MARKERSIZE.cloud)
    l_sc = lines!(ax2, Zgrid, scission_deformation.(Zgrid),
        color = PALETTE.red, linewidth = 2)
    vlines!(ax2, [28, 50], color = PALETTE.black, linestyle = :dot, linewidth = 1.2)
    text!(ax2, 28, -0.22; text = rich(it("Z"), " = 28"),
        align = (:center, :bottom), fontsize = 13)
    text!(ax2, 50, -0.22; text = rich(it("Z"), " = 50"),
        align = (:center, :bottom), fontsize = 13)
    ylims!(ax2, -0.28, 0.88)
    axislegend(ax2,
        [MarkerElement(color = PALETTE.sky, marker = :circle,
                       markersize = MARKERSIZE.key),
         LineElement(color = PALETTE.red, linewidth = 2)],
        ["Ground state, Möller–Nix", "At scission"];
        position = :rt, framevisible = false, labelsize = 14, padding = 4)

    rowsize!(fig.layout, 2, Relative(0.86))
    colgap!(fig.layout, 1, 30)
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_multiplicity"))
end

"Neutron separation energy of the fragment, MeV."
function separation_cost(masses, Z, A)
    d1 = Δ(masses, Z, A); d2 = Δ(masses, Z, A - 1); dn = Δ(masses, 0, 1)
    (d1 === nothing || d2 === nothing) && return nothing
    return (d2 + dn - d1) / 1000
end

main()
