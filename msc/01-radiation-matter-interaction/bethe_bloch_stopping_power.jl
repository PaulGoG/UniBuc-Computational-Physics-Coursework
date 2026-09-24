# Bethe-Bloch mass electronic stopping power of a 5 MeV α particle in silicon.
#
#   -dE/dx = 4π N_A r_e² m_ec² ρ (Z/A) (z²/β²) [ ln(2m_ec²γ²β²W_max/I²)/2 - β² ]
#
# in the Leo form, with W_max ≈ 2m_ec²β²γ² and the mean excitation potential
# from the Sternheimer-Barkas parametrisation I = Z(9.76 + 58.8 Z^{-1.19}) eV,
# valid for Z ≥ 13.
#
# Ported from Calcul_Bethe_Bloch.jl on the `legacy` branch (Julia-Workflow-FFUB/IRM_M_1/).
# The formula was transcribed correctly. What
# it lacked was any statement of what it was computing: the target material is
# named nowhere in the file, A was written as 28 rather than the natural silicon
# molar mass 28.085, the electron rest energy appeared three times as a bare
# 0.511, and the log argument typed the same sub-expression twice instead of
# squaring it, obscuring that the second factor is W_max.
#
# The physics caveat the original did not make: at β = 0.052, which is where a
# 5 MeV α sits, this is near the low-energy validity limit of Bethe-Bloch. The
# shell correction C/Z and the Barkas effective-charge term are not negligible
# there, and neither is included. ASTAR (Berger, Coursey, Zucker and Chang,
# NIST Standard Reference Database 124, doi:10.18434/T4NC7P) gives
# 617.4 MeV cm² g⁻¹ at 5 MeV and a CSDA range of 5.651 × 10⁻³ g cm⁻²; both are
# asserted below within 5 %.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Avogadro constant [mol⁻¹]."
const N_A = 6.02214076e23
"Classical electron radius [m]."
const R_E = 2.8179403262e-15
"Electron rest energy [MeV]."
const MEC² = 0.51099895

"Silicon target."
const SILICON = (Z = 14, A = 28.085, ρ = 2.329e6)   # g m⁻³
"Alpha projectile."
const ALPHA = (z = 2, Mc² = 4.0015065 * 931.494)     # MeV

"ASTAR electronic stopping power of a 5 MeV α in silicon [MeV cm² g⁻¹], doi:10.18434/T4NC7P."
const ASTAR_STOPPING_5MEV = 617.4
"ASTAR CSDA ranges of a 0.5 MeV and a 5 MeV α in silicon [g cm⁻²]."
const ASTAR_CSDA = (at_0p5MeV = 5.286e-4, at_5MeV = 5.651e-3)
"Relative tolerance of the Bethe–Bloch results against ASTAR."
const ASTAR_TOLERANCE = 0.05

"Mean excitation potential in MeV, Sternheimer-Barkas, valid for Z ≥ 13."
excitation_potential(Z) = Z * (9.76 + 58.8 * Z^(-1.19)) * 1e-6

"""
    stopping_power(E_kin, projectile, target)

Electronic stopping power −dE/dx in MeV cm⁻¹.
"""
function stopping_power(E_kin, projectile, target)
    γ = E_kin / projectile.Mc² + 1
    β² = 1 - 1 / γ^2
    I = excitation_potential(target.Z)
    W_max = 2 * MEC² * β² * γ^2                     # maximum energy transfer
    prefactor = 2π * N_A * R_E^2 * MEC² * target.ρ * target.Z / target.A
    bracket = log(2 * MEC² * γ^2 * β² * W_max / I^2) - 2β²
    return prefactor * projectile.z^2 / β² * bracket / 100   # MeV m⁻¹ -> MeV cm⁻¹
end

function main()
    E = 5.0
    S = stopping_power(E, ALPHA, SILICON)
    γ = E / ALPHA.Mc² + 1
    β = sqrt(1 - 1/γ^2)
    @printf("target: silicon, Z = %d, A = %.3f, ρ = %.3f g/cm³\n",
        SILICON.Z, SILICON.A, SILICON.ρ / 1e6)
    @printf("mean excitation potential I = %.1f eV\n",
        excitation_potential(SILICON.Z) * 1e6)
    @printf("5 MeV α: γ = %.6f, β = %.5f\n", γ, β)
    @printf("dE/dx = %.1f MeV/cm = %.1f MeV cm²/g\n", S, S / (SILICON.ρ / 1e6))
    S_mass = S / (SILICON.ρ / 1e6)
    @printf("ASTAR at 5 MeV: %.1f MeV cm²/g, ratio Bethe–Bloch / ASTAR = %.3f\n",
        ASTAR_STOPPING_5MEV, S_mass / ASTAR_STOPPING_5MEV)
    abs(S_mass / ASTAR_STOPPING_5MEV - 1) < ASTAR_TOLERANCE ||
        error("Bethe–Bloch stopping power $S_mass MeV cm²/g is off ASTAR by more than $(100ASTAR_TOLERANCE) %")

    # independent check: integrate 1/(dE/dx) to a CSDA range
    E_low = 0.5
    range_cm, _ = quadgk(e -> 1 / stopping_power(e, ALPHA, SILICON), E_low, E)
    astar_range_cm = (ASTAR_CSDA.at_5MeV - ASTAR_CSDA.at_0p5MeV) / (SILICON.ρ / 1e6)
    @printf("\nCSDA range from %.1f to %.1f MeV = %.2f µm; ASTAR %.2f µm, ratio %.3f\n",
        E_low, E, range_cm * 1e4, astar_range_cm * 1e4, range_cm / astar_range_cm)
    abs(range_cm / astar_range_cm - 1) < ASTAR_TOLERANCE ||
        error("CSDA range $(range_cm * 1e4) µm is off ASTAR by more than $(100ASTAR_TOLERANCE) %")

    energies = 10 .^ range(log10(0.5), log10(200), length = 300)
    S_curve = stopping_power.(energies, Ref(ALPHA), Ref(SILICON))

    fig = Figure(size = (900, 600))
    ax = Axis(fig[1, 1],
        xlabel = L"Kinetic energy $E$ [MeV]",
        ylabel = L"$-\mathrm{d}E/\mathrm{d}x$ [MeV cm$^{-1}$]",
        xscale = log10, yscale = log10,
        xticks = logticks(-1, 2),
        # explicit: over a range narrower than a decade and a half Makie labels
        # this axis 10^{3.3}, 10^{3.0}, ..., which no stopping power is read in
        yticks = ([200, 500, 1000, 2000], [L"200", L"500", L"1000", L"2000"]),)
    lines!(ax, energies, S_curve, color = PALETTE.blue)
    vlines!(ax, [E], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    S_astar = ASTAR_STOPPING_5MEV * SILICON.ρ / 1e6
    scatter!(ax, [E], [S], color = PALETTE.blue, markersize = MARKERSIZE.emphasis)
    scatter!(ax, [E], [S_astar], color = PALETTE.black, marker = :diamond,
        markersize = MARKERSIZE.emphasis,)
    text!(ax, E * 1.25, S * 1.12;
        text = rich("Bethe–Bloch ", @sprintf("%.0f MeV cm", S), superscript("−1")),
        color = PALETTE.blue, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)
    # Left of the marker: the curve descends through the space to its right.
    text!(ax, E / 1.12, S_astar / 1.12;
        text = rich("ASTAR ", @sprintf("%.0f MeV cm", S_astar), superscript("−1")),
        color = PALETTE.black, align = (:right, :top), fontsize = ANNOTATION_SIZE,)
    text!(ax, 0.97, 0.95;
        text = "Shell and Barkas corrections omitted",
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,)
    println("wrote ", savefigure(fig, FIGURES, "bethe_bloch_stopping_power"))
end

main()
