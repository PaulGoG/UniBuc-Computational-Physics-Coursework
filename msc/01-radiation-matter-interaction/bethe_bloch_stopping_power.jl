# Bethe-Bloch mass electronic stopping power of a 5 MeV α particle in silicon.
#
#   -dE/dx = 4π N_A r_e² m_ec² ρ (Z/A) (z²/β²) [ ln(2m_ec²γ²β²W_max/I²)/2 - β² ]
#
# in the Leo form, with W_max ≈ 2m_ec²β²γ² and the mean excitation potential
# from the Sternheimer-Barkas parametrisation I = Z(9.76 + 58.8 Z^{-1.19}) eV,
# valid for Z ≥ 13.
#
# Ported from Calcul_Bethe_Bloch.jl. The formula was transcribed correctly. What
# it lacked was any statement of what it was computing: the target material is
# named nowhere in the file, A was written as 28 rather than the natural silicon
# molar mass 28.085, the electron rest energy appeared three times as a bare
# 0.511, and the log argument typed the same sub-expression twice instead of
# squaring it, obscuring that the second factor is W_max.
#
# The physics caveat the original did not make: at β = 0.052, which is where a
# 5 MeV α sits, this is near the low-energy validity limit of Bethe-Bloch. The
# shell correction C/Z and the Barkas effective-charge term are not negligible
# there, and neither is included. The result should therefore sit above the
# tabulated value, and it does.

using Printf, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

const N_A = 6.02214076e23        # mol⁻¹
const R_E = 2.8179403262e-15     # classical electron radius, m
const MEC² = 0.51099895          # MeV

"Silicon target."
const SILICON = (Z = 14, A = 28.085, ρ = 2.329e6)   # g m⁻³
"Alpha projectile."
const ALPHA = (z = 2, Mc² = 4.0015065 * 931.494)     # MeV

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
    @printf("mean excitation potential I = %.1f eV\n", excitation_potential(SILICON.Z) * 1e6)
    @printf("5 MeV α: γ = %.6f, β = %.5f\n", γ, β)
    @printf("dE/dx = %.1f MeV/cm = %.1f MeV cm²/g\n", S, S / (SILICON.ρ / 1e6))

    # independent check: integrate 1/(dE/dx) to a CSDA range
    E_low = 0.5
    range_cm, _ = quadgk(e -> 1 / stopping_power(e, ALPHA, SILICON), E_low, E)
    @printf("\nCSDA range from %.1f to %.1f MeV = %.2f µm\n", E_low, E, range_cm * 1e4)
    @printf("a 5 MeV α in silicon has a full range of about 25 µm, so the missing\n")
    @printf("piece is the sub-%.1f MeV portion. Below roughly that energy the Bethe\n", E_low)
    @printf("logarithm turns over and the formula stops describing the physics --\n")
    @printf("the shell and Barkas corrections it omits are no longer small.\n")

    energies = 10 .^ range(log10(0.5), log10(200), length = 300)
    S_curve = stopping_power.(energies, Ref(ALPHA), Ref(SILICON))

    fig = Figure(size = (760, 480))
    ax = Axis(fig[2, 1],
        xlabel = L"Kinetic energy $E$ [MeV]",
        ylabel = L"$-\mathrm{d}E/\mathrm{d}x$ [MeV cm$^{-1}$]",
        xscale = log10, yscale = log10,
        xticks = logticks(-1, 2),
        # explicit: over a range narrower than a decade and a half Makie labels
        # this axis 10^{3.3}, 10^{3.0}, ..., which no stopping power is read in
        yticks = ([200, 500, 1000, 2000], [L"200", L"500", L"1000", L"2000"]))
    l = lines!(ax, energies, S_curve, color = PALETTE.blue, linewidth = 1.8)
    vlines!(ax, [E], color = PALETTE.red, linestyle = :dash, linewidth = 1.2)
    scatter!(ax, [E], [S], color = PALETTE.red, markersize = MARKERSIZE.emphasis)
    # The units here follow the axis rather than switching to MeV/cm, and the
    # legend entry that repeated this line is gone.
    text!(ax, E * 1.3, S;
        text = rich(@sprintf("%.0f MeV cm", S), superscript("−1"), " at 5 MeV"),
        color = PALETTE.red, align = (:left, :center), fontsize = 16)
    text!(ax, 0.97, 0.93;
        text = "Shell and Barkas corrections omitted;\nBethe–Bloch breaks down below ≈0.5 MeV",
        space = :relative, align = (:right, :top), fontsize = 15)

    Legend(fig[1, 1], [l], [L"$\alpha$ in silicon"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("wrote ", savefigure(fig, FIGURES, "bethe_bloch_stopping_power"))
end

main()
