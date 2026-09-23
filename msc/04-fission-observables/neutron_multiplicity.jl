# Prompt-neutron multiplicity per fragment, ν(A), from the energy balance of
# Fisiune_3.jl on the `legacy` branch (Julia-Workflow-FFUB/Fisiune_M_2/),
# against four measured ν(A) sets.
#
#   ν_pair = (TXE − q) / (⟨ε⟩ + ⟨Sₙ⟩ + p),   ⟨ε⟩ = 4T_m/3,   T_m = √(TXE/(a_L + a_H))
#   p = 6.71 − 0.156 Z₀²/A₀,   q = 0.75 + 0.088 Z₀²/A₀      [MeV, compound system]
#
# with a = A[0.00917 (S_Z + S_N) + 0.142] the Gilbert–Cameron level-density
# parameter (Can. J. Phys. 43, 1446 (1965), doi:10.1139/p65-139) and ⟨Sₙ⟩ the
# mean neutron separation energy of the two fragments. The pair multiplicity is
# split between the fragments by the share R of the excitation energy taken by
# the heavy one. Fisiune_3.jl obtains R from the scission-point deformation:
# each fragment first pays the liquid-drop energy of deforming from its
# ground-state β₂ to a scission β(Z), and the remainder is shared in the ratio
# of the level-density parameters,
#
#   E*_H = ΔE_H + (TXE − ΔE_H − ΔE_L) / (1 + a_L/a_H),      R = E*_H / TXE.
#
# The statistical limit R = a_H/(a_L + a_H) is computed alongside and compared
# with the same data. Every quantity is evaluated per (A_H, Z_H) over the
# charges the mass table lists, averaged over Z_H with the Gaussian isobaric
# weights of rms 0.6 about Z_p = Z_UCD − 0.5, and the totals are weighted by the
# mass yield. The β(Z) breakpoints, the liquid-drop constants and the p, q
# terms are those of the fission course's TXE-partition application, as
# Fisiune_3.jl carried them; the README gives their sources.
#
# Fisiune_3.jl bounded its charge loop by the index column of the
# Gilbert–Cameron table, 11 to 150; charges far from Z_p enter with a Gaussian
# weight below 10⁻⁵ beyond |Z − Z_p| = 3, so the bound cost time, not results.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_core.jl"))

"""
Breakpoints (Z, β) of the scission deformation, the piecewise-linear
parametrisation of the course application that Fisiune_3.jl implemented; linear
between them and undefined outside 28 ≤ Z ≤ 65, its four segments reproduced.
"""
const SCISSION_BETA = [(28, 0.0), (41, 0.58), (44, 0.58), (50, 0.0), (65, 0.6)]

"""
Liquid-drop energy constants [MeV], the Myers–Swiatecki set of 1967 (Ark. Fys.
36, 343) as Fisiune_3.jl carried them:
E(β) = −χ [c_v A − c_s A^{2/3} (1 + 0.4α₂)] + Z² [c_c (1 − 0.2α₂)/A^{1/3} − c_x/A],
χ = 1 − κ I², α₂ = 5β²/4π.
"""
const LDM = (c_v = 15.4941, c_s = 17.9439, κ = 1.7826, c_c = 0.7053, c_x = 1.1529)

"Systematics terms of the energy balance for the compound system (Z₀, A₀) [MeV]."
const P_TERM = 6.71 - 0.156 * Z₀^2 / A₀
const Q_TERM = 0.75 + 0.088 * Z₀^2 / A₀

"Tolerance of the identity ν_H + ν_L = ν_pair per mass."
const SPLIT_TOLERANCE = 1e-12

"Gilbert–Cameron level-density parameter [MeV⁻¹], or `nothing` off the table."
function level_density(gc, Z, A)
    (haskey(gc, Z) && haskey(gc, A - Z)) || return nothing
    return A * (0.00917 * (gc[Z].S_Z + gc[A - Z].S_N) + 0.142)
end

"Scission deformation β(Z) by linear interpolation of `SCISSION_BETA`; `nothing` outside its range."
function scission_deformation(Z)
    (first(SCISSION_BETA)[1] <= Z <= last(SCISSION_BETA)[1]) || return nothing
    for i in 1:(length(SCISSION_BETA) - 1)
        (z₁, β₁), (z₂, β₂) = SCISSION_BETA[i], SCISSION_BETA[i + 1]
        z₁ <= Z <= z₂ && return β₁ + (β₂ - β₁) * (Z - z₁) / (z₂ - z₁)
    end
    return nothing
end

"Liquid-drop energy of the nucleus (Z, A) deformed to β [MeV]."
function ldm_energy(β, Z, A)
    I = (A - 2Z) / A
    χ = 1 - LDM.κ * I^2
    α₂ = 5β^2 / (4π)
    return -χ * (LDM.c_v * A - LDM.c_s * A^(2 / 3) * (1 + 0.4α₂)) +
           Z^2 * (LDM.c_c * (1 - 0.2α₂) / A^(1 / 3) - LDM.c_x / A)
end

"|E(β_scission) − E(β_ground)| of the fragment (Z, A), as Fisiune_3.jl takes it; `nothing` if either β is absent."
function deformation_cost(β_gs, Z, A)
    β_sc = scission_deformation(Z)
    β₀ = get(β_gs, (Z, A), nothing)
    (β_sc === nothing || β₀ === nothing) && return nothing
    return abs(ldm_energy(β_sc, Z, A) - ldm_energy(β₀, Z, A))
end

"Weighted-mean accumulator: Σwx and Σw."
mutable struct Accumulator
    num::Float64
    den::Float64
end
Accumulator() = Accumulator(0.0, 0.0)
add!(a::Accumulator, w, x) = (a.num += w * x; a.den += w; a)
value(a::Accumulator) = a.den > 0 ? a.num / a.den : nothing

"""
    energy_balance(A_H, masses, gc, β_gs, TKE, S_n) -> NamedTuple or nothing

The per-mass averages over Z_H of ν_pair, of the two partitions of it and of
their shares R, for the heavy-fragment mass `A_H` at total kinetic energy
`TKE`. The ν_pair average is taken over the charges where the partition is
defined, so that ν_H + ν_L = ν_pair holds exactly.
"""
function energy_balance(A_H, masses, gc, β_gs, TKE, S_n)
    A_L = A₀ - A_H
    acc = (ν = Accumulator(), ν_all = Accumulator(), TXE = Accumulator(),
        R_def = Accumulator(), R_ld = Accumulator(),
        ν_H_def = Accumulator(), ν_L_def = Accumulator(),
        ν_H_ld = Accumulator(), ν_L_ld = Accumulator(),)
    for Z_H in tabulated_charges(masses, A_H)
        Z_L = Z₀ - Z_H
        q = q_value(masses, Z_H, A_H)
        (q === nothing || q[1] <= 0) && continue
        TXE = q[1] + S_n - TKE
        TXE > 0 || continue
        w = charge_weight(Z_H, A_H)
        add!(acc.TXE, w, TXE)
        a_H, a_L = level_density(gc, Z_H, A_H), level_density(gc, Z_L, A_L)
        (a_H === nothing || a_L === nothing || a_H <= 0 || a_L <= 0) && continue
        s_H, s_L = separation_energy(masses, Z_H, A_H), separation_energy(masses, Z_L, A_L)
        (s_H === nothing || s_L === nothing) && continue
        S̄ = (s_H[1] + s_L[1]) / 2
        TXE >= S̄ || continue
        T_m = sqrt(TXE / (a_H + a_L))
        ν_pair = (TXE - Q_TERM) / (4T_m / 3 + S̄ + P_TERM)
        add!(acc.ν_all, w, ν_pair)
        R_ld = a_H / (a_H + a_L)

        ΔE_H, ΔE_L = deformation_cost(β_gs, Z_H, A_H), deformation_cost(β_gs, Z_L, A_L)
        (ΔE_H === nothing || ΔE_L === nothing) && continue
        ε = TXE - ΔE_H - ΔE_L
        ε > 0 || continue
        R_def = (ΔE_H + ε / (1 + a_L / a_H)) / TXE
        0 < R_def < 1 || error("R = $R_def outside (0, 1) at A_H = $A_H, Z_H = $Z_H")
        add!(acc.ν, w, ν_pair)
        add!(acc.R_def, w, R_def)
        add!(acc.R_ld, w, R_ld)
        add!(acc.ν_H_def, w, R_def * ν_pair)
        add!(acc.ν_L_def, w, (1 - R_def) * ν_pair)
        add!(acc.ν_H_ld, w, R_ld * ν_pair)
        add!(acc.ν_L_ld, w, (1 - R_ld) * ν_pair)
    end
    value(acc.ν) === nothing && return nothing
    return map(value, acc)
end

"χ²/N and rms of a model ν(A) against a measured set, over the masses both cover."
function agreement(model::Dict{Int,Float64}, d)
    r = [(d.y[i] - model[Int(d.x[i])]) / d.σ[i]
         for i in eachindex(d.x) if d.y[i] > 0 && d.σ[i] > 0 && haskey(model, Int(d.x[i]))]
    dev = [d.y[i] - model[Int(d.x[i])]
           for i in eachindex(d.x) if d.y[i] > 0 && haskey(model, Int(d.x[i]))]
    return (χ²N = mean(abs2, r), rms = sqrt(mean(abs2, dev)), n = length(dev))
end

function main()
    rows = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    gc = load_gilbert_cameron(joinpath(DATA, "Parametrizari_auxiliare", "SZSN.GC"))
    β_gs = load_beta2(joinpath(DATA, "Parametrizari_auxiliare", "B2MOLLER.ANA"))
    bins = tke_by_mass(yield_cells(rows))
    S_n = separation_energy(masses, Z₀, A₀)[1]

    scission_deformation(28) == 0 && scission_deformation(50) == 0 ||
        error("the scission deformation does not vanish at the closed shells")
    @printf("compound system ²³⁶U: p = %.3f MeV, q = %.3f MeV, Sₙ = %.3f MeV\n",
        P_TERM, Q_TERM, S_n)

    A_H = Int[]
    balance = NamedTuple[]
    for a in A_H_RANGE
        haskey(bins, a) || continue
        b = energy_balance(a, masses, gc, β_gs, bins[a].TKE, S_n)
        b === nothing && continue
        push!(A_H, a)
        push!(balance, b)
    end
    Y = [bins[a].Y for a in A_H]
    for b in balance
        abs(b.ν_H_def + b.ν_L_def - b.ν) < SPLIT_TOLERANCE &&
        abs(b.ν_H_ld + b.ν_L_ld - b.ν) < SPLIT_TOLERANCE ||
            error("ν_H + ν_L differs from ν_pair")
    end
    yw(f) = sum(f.(balance) .* Y) / sum(Y)
    @printf("energy balance on %d heavy-fragment masses, A_H %d–%d\n", length(A_H), A_H[1],
        A_H[end])
    @printf("⟨ν_pair⟩ = %.3f yield-weighted (%.3f unweighted over masses); ⟨TXE⟩ = %.2f MeV\n",
        yw(b -> b.ν), mean(b.ν for b in balance), yw(b -> b.TXE))
    @printf("heavy-fragment share of the excitation energy, yield-weighted: R = %.3f with the scission deformation, %.3f from the level densities alone\n",
        yw(b -> b.R_def), yw(b -> b.R_ld))

    # ν(A) over both fragments, for each partition
    curve(h, l) = merge(Dict(zip(A_H, h)), Dict(zip(A₀ .- A_H, l)))
    ν_def = curve([b.ν_H_def for b in balance], [b.ν_L_def for b in balance])
    ν_ld = curve([b.ν_H_ld for b in balance], [b.ν_L_ld for b in balance])

    sets = [("Göök", "U5NUAGOOK.DAT", PALETTE.blue, :circle),
        ("Maslin", "U5NUAMASLIN.DAT", PALETTE.orange, :rect),
        ("Nishio", "U5NUANISHIO.DAT", PALETTE.green, :utriangle),
        ("Vorobyev", "U5NUAVORO.DAT", PALETTE.purple, :diamond),]
    Y_of = Dict(a => bins[a].Y for a in keys(bins))
    println("\nmeasured ν(A): yield-weighted mean per fragment, ×2 for the pair; model against data")
    measured = map(sets) do (name, file, colour, marker)
        d = load_measurement(joinpath(DATA, "Date_experimentale", "Multiplicitate_n", file))
        keep = [d.y[i] > 0 && (haskey(Y_of, Int(d.x[i])) || haskey(Y_of, A₀ - Int(d.x[i])))
                for i in eachindex(d.x)]
        wts = [get(Y_of, Int(a), get(Y_of, A₀ - Int(a), 0.0)) for a in d.x[keep]]
        per_fragment = sum(wts .* d.y[keep]) / sum(wts)
        g_def, g_ld = agreement(ν_def, d), agreement(ν_ld, d)
        @printf("  %-9s %3d points  ν̄ = %.3f (×2 = %.3f)   deformation partition χ²/N = %6.1f, rms %.3f;   level densities χ²/N = %6.1f, rms %.3f\n",
            name, count(d.y .> 0), per_fragment, 2per_fragment, g_def.χ²N, g_def.rms,
            g_ld.χ²N, g_ld.rms)
        (; name, A = d.x[d.y .> 0], ν = d.y[d.y .> 0], σ = d.σ[d.y .> 0], colour, marker,
            rms_def = g_def.rms, rms_ld = g_ld.rms,)
    end

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1600, 720))
    ax = Axis(fig[2, 1], xlabel = L"Fragment mass $A$", ylabel = L"\nu(A)", xticks = 80:20:160)
    handles = Any[]
    for m in measured
        errorbars_unstroked!(
            ax, m.A, m.ν, m.σ, color = (m.colour, 0.45), linewidth = GUIDE_WIDTH,
            whiskerwidth = 0,)
        push!(handles,
            scatter!(ax, m.A, m.ν, color = m.colour, marker = m.marker,
                markersize = MARKERSIZE.dense,),)
    end
    A_curve = sort!(collect(keys(ν_def)))
    l_def = lines!(ax, A_curve, [ν_def[a] for a in A_curve], color = PALETTE.red)
    l_ld = lines!(
        ax, A_curve, [ν_ld[a] for a in A_curve], color = PALETTE.red, linestyle = :dash,)
    xlims!(ax, 70, 170)
    ylims!(ax, 0, 4.2)
    text!(ax, 0.03, 0.97;
        text = rich(
            rich("⟨", it("ν"), subscript("pair"), @sprintf("⟩ = %.3f", yw(b -> b.ν)),
                color = PALETTE.red,),
            "\n",
            @sprintf("rms against Göök: %.2f (deformation), %.2f (level densities)",
                measured[1].rms_def, measured[1].rms_ld)),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)

    # β₂ of the fragment nuclides in the ground state against the scission line
    ax2 = Axis(fig[2, 2], xlabel = L"Fragment charge $Z$",
        ylabel = L"Quadrupole deformation $\beta_2$",
        xticks = 30:10:60,)
    Zs = Set(vcat([r.Z_H for r in rows], [Z₀ - r.Z_H for r in rows]))
    gs = [(Z, b) for ((Z, A), b) in β_gs if Z in Zs]
    l_gs = scatter!(ax2, first.(gs), last.(gs), color = (PALETTE.sky, 0.4),
        markersize = MARKERSIZE.cloud,
        strokewidth = 0,)
    Z_grid = range(28, 65, length = 300)
    l_sc = lines!(ax2, Z_grid, scission_deformation.(Z_grid), color = PALETTE.red)
    vlines!(ax2, [28, 50], color = PALETTE.black, linestyle = :dot, linewidth = GUIDE_WIDTH)
    for Zm in (28, 50)
        text!(
            ax2, Zm + 0.6, -0.24; text = rich(it("Z"), " = $Zm"), align = (:left, :bottom),
            fontsize = ANNOTATION_SIZE,)
    end
    xlims!(ax2, 26, 67)
    ylims!(ax2, -0.28, 0.9)

    Legend(fig[1, 1:2],
        [handles, [l_def, l_ld],
            [
                MarkerElement(color = PALETTE.sky, marker = :circle, markersize = MARKERSIZE.key), l_sc,],],
        [[m.name for m in measured],
            ["Scission-deformation partition", "Level-density partition"],
            ["Ground state, FRDM", "At scission"],],
        [rich("Measured ", it("ν"), "(", it("A"), ")"),
            "Energy balance", rich(it("β"), subscript("2")),];
        titleposition = :left, nbanks = 2, groupgap = 30,)
    colsize!(fig.layout, 2, Relative(0.34))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_multiplicity"))
end

main()
