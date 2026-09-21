# The accelerating universe, from the data the discovery was made with.
#
# Perlmutter et al. (1999), "Measurements of Ω and Λ from 42 High-Redshift
# Supernovae", ApJ 517, 565 — arXiv:astro-ph/9812133 — Tables 1 and 2: 42
# high-redshift SNe Ia from the Supernova Cosmology Project and 18 low-redshift
# SNe from the Calán/Tololo survey, with stretch-corrected effective B-band peak
# magnitudes and their uncertainties.
#
# The paper and its two data tables were in the Dubna project folder, alongside
# Riess et al. (1998) — astro-ph/9805201, the companion discovery — and the SAI
# Supernova Catalogue that `supernova_magnitudes.csv` came from. They were never
# used: the 2019 code fitted host-galaxy magnitudes instead, which carry no
# cosmological information at all. The catalogue's own `sn_mag` column, the
# supernova magnitude at maximum, sat twelve columns from the one taken.
#
# These do. The effective magnitudes are standardised: stretch-corrected,
# K-corrected and extinction-corrected, which is what makes a Type Ia a standard
# candle and what the host-galaxy magnitudes and the compiled distance moduli in
# `friedmann_expansion.jl` are not.
#
# Two results follow. Flat ΛCDM fitted here returns Ω_M within half a sigma of
# the value the paper published, which validates the implementation against the
# literature rather than against itself. And the Einstein–de Sitter universe —
# flat, matter only, no dark energy — is rejected by Δχ² ≈ 34, which is the
# acceleration discovery itself.
#
# Provenance: presented at the Laboratory of Information Technology, JINR Dubna,
# 25 July 2019, supervised by Bijan Saha and Victor Rikhvitsky, as "Numerical
# simulation of homogeneous and isotropic universe given by Friedmann equations"
# (A. Crăciun and P. Gogîță).

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames, LsqFit, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Dimensionless expansion rate E(z) = H(z)/H₀ for a flat universe of matter and Λ."
E(z, Ωm) = sqrt(Ωm * (1 + z)^3 + (1 - Ωm))

"""
    luminosity_distance(z, Ωm; n = 500)

Luminosity distance in units of c/H₀, `(1+z)∫₀^z dz'/E(z')`, by Simpson's rule.

Left dimensionless on purpose. The absolute magnitude of a Type Ia and H₀ enter
the apparent magnitude only through the sum `M − 5log₁₀H₀`, so they cannot be
separated by this data and are fitted as one offset. That is exactly how the
paper treats them.
"""
function luminosity_distance(z, Ωm; n = 500)
    zs = range(0, z, length = n + 1)
    h = step(zs)
    f = [1 / E(zi, Ωm) for zi in zs]
    I = h / 3 * (f[1] + f[end] + 4sum(f[2:2:(end - 1)]) + 2sum(f[3:2:(end - 2)]))
    return (1 + z) * I
end

"Apparent effective B magnitude: 5log₁₀ of the dimensionless distance, plus the offset."
apparent_magnitude(z, Ωm, offset) = 5 * log10(luminosity_distance(z, Ωm)) + offset

const Ωm_BOUNDS = (0.02, 1.5)

"χ² of a model against the measured magnitudes."
chi_squared(z, m, σ, Ωm, offset) = sum(((m .- apparent_magnitude.(z, Ωm, offset)) ./ σ) .^
                                       2)

"""
    fit_flat_lcdm(z, m, σ)

Fit Ω_M and the magnitude offset of a flat ΛCDM universe, weighted by the quoted
magnitude uncertainties.
"""
function fit_flat_lcdm(z, m, σ)
    model(zz, p) = [apparent_magnitude(zi, clamp(p[2], Ωm_BOUNDS...), p[1]) for zi in zz]
    fit = curve_fit(model, z, m, 1 ./ σ .^ 2, [24.0, 0.3];
        lower = [15.0, Ωm_BOUNDS[1]], upper = [30.0, Ωm_BOUNDS[2]],)
    offset, Ωm = fit.param[1], clamp(fit.param[2], Ωm_BOUNDS...)
    σp = try
        stderror(fit)
    catch
        [NaN, NaN]
    end
    return (; Ωm, offset, σΩm = σp[2], σoffset = σp[1],
        χ² = chi_squared(z, m, σ, Ωm, offset), dof = length(z) - 2,)
end

"""
    fit_fixed_Ωm(z, m, σ, Ωm)

Fit only the magnitude offset, holding Ω_M at a stated value. Used for the
Einstein–de Sitter comparison, which has one free parameter rather than two.
"""
function fit_fixed_Ωm(z, m, σ, Ωm)
    model(zz, p) = [apparent_magnitude(zi, Ωm, p[1]) for zi in zz]
    fit = curve_fit(model, z, m, 1 ./ σ .^ 2, [24.0])
    offset = fit.param[1]
    return (; Ωm, offset, χ² = chi_squared(z, m, σ, Ωm, offset), dof = length(z) - 1)
end

function main()
    d = CSV.read(joinpath(@__DIR__, "data", "perlmutter1999_sn_ia.csv"), DataFrame)
    z, m, σ = d.z, d.m_b_effective, d.sigma_m
    high = d.sample .== "high-z"

    @printf("Perlmutter et al. 1999: %d supernovae, %d high-z and %d low-z\n",
        nrow(d), count(high), count(.!high))
    @printf("z from %.3f to %.3f, effective m_B from %.2f to %.2f\n",
        minimum(z), maximum(z), minimum(m), maximum(m))

    lcdm = fit_flat_lcdm(z, m, σ)
    @printf("\nflat ΛCDM:  Ω_M = %.3f ± %.3f,  offset = %.3f ± %.3f mag\n",
        lcdm.Ωm, lcdm.σΩm, lcdm.offset, lcdm.σoffset)
    @printf("χ² = %.1f on %d dof  (χ²/dof = %.2f)\n", lcdm.χ², lcdm.dof, lcdm.χ² / lcdm.dof)

    # --- validation against the published value --------------------------------
    published, published_σ = 0.28, 0.085           # Perlmutter 1999, flat universe
    @printf("\nPublished by the paper for a flat universe: Ω_M = 0.28 (+0.09 −0.08)\n")
    @printf("Recovered here: %.3f. Difference %.2fσ of the published uncertainty.\n",
        lcdm.Ωm, abs(lcdm.Ωm - published) / published_σ)
    @printf("Planck 2018, for context: Ω_m = 0.315 ± 0.007\n")

    # --- the discovery ---------------------------------------------------------
    eds = fit_fixed_Ωm(z, m, σ, 1.0)
    Δχ² = eds.χ² - lcdm.χ²
    @printf("\nEinstein–de Sitter (flat, matter only, no dark energy):\n")
    @printf("  χ² = %.1f on %d dof — worse by Δχ² = %.1f\n", eds.χ², eds.dof, Δχ²)
    @printf("  that is the acceleration: these data reject a decelerating\n")
    @printf("  matter-only universe, which is what the 1999 paper reported and\n")
    @printf("  what the 2011 Nobel Prize was awarded for.\n")

    @printf("\nThe 2019 project had this file and did not use it. It fitted\n")
    @printf("host-galaxy magnitudes instead — see friedmann_expansion.jl, where\n")
    @printf("those give 1.35 mag of scatter and no constraint on Ω_M whatsoever.\n")

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1000, 780))

    zf = 10 .^ range(log10(0.01), log10(0.95), length = 300)
    ax1 = Axis(fig[2, 1], ylabel = L"Effective $m_B$ [mag]", xscale = log10,
        xticks = logticks(-2, 0), xticklabelsvisible = false,)
    l_lcdm = lines!(ax1, zf, [apparent_magnitude(zi, lcdm.Ωm, lcdm.offset) for zi in zf],
        color = PALETTE.blue, linewidth = 2,)
    l_eds = lines!(ax1, zf, [apparent_magnitude(zi, 1.0, eds.offset) for zi in zf],
        color = PALETTE.orange, linewidth = 2, linestyle = :dash,)
    errorbars!(ax1, z[.!high], m[.!high], σ[.!high],
        color = (PALETTE.green, 0.55), whiskerwidth = 7,)
    errorbars!(ax1, z[high], m[high], σ[high],
        color = (PALETTE.red, 0.55), whiskerwidth = 7,)
    l_lo = scatter!(ax1, z[.!high], m[.!high], color = PALETTE.green,
        markersize = MARKERSIZE.dense,)
    l_hi = scatter!(ax1, z[high], m[high], color = PALETTE.red,
        markersize = MARKERSIZE.dense,)
    text!(ax1, 0.03, 0.96;
        text = rich("Ω", subscript("M"), @sprintf(" = %.2f ± %.2f", lcdm.Ωm, lcdm.σΩm)),
        space = :relative, align = (:left, :top), fontsize = 15, color = PALETTE.blue,)
    text!(ax1, 0.03, 0.885; text = "published 0.28 (+0.09 −0.08)",
        space = :relative, align = (:left, :top), fontsize = 15, color = PALETTE.blue,)

    ax2 = Axis(fig[3, 1], xlabel = L"Redshift $z$",
        ylabel = L"$m_B$ − Einstein–de Sitter [mag]", xscale = log10,
        xticks = logticks(-2, 0),)
    ref(zi) = apparent_magnitude(zi, 1.0, eds.offset)
    hlines!(ax2, [0.0], color = PALETTE.orange, linestyle = :dash, linewidth = 2)
    lines!(ax2, zf, [apparent_magnitude(zi, lcdm.Ωm, lcdm.offset) - ref(zi) for zi in zf],
        color = PALETTE.blue, linewidth = 2,)
    errorbars!(ax2, z[.!high], m[.!high] .- ref.(z[.!high]), σ[.!high],
        color = (PALETTE.green, 0.55), whiskerwidth = 7,)
    errorbars!(ax2, z[high], m[high] .- ref.(z[high]), σ[high],
        color = (PALETTE.red, 0.55), whiskerwidth = 7,)
    scatter!(ax2, z[.!high], m[.!high] .- ref.(z[.!high]),
        color = PALETTE.green, markersize = MARKERSIZE.dense,)
    scatter!(ax2, z[high], m[high] .- ref.(z[high]), color = PALETTE.red,
        markersize = MARKERSIZE.dense,)
    # right of the low-z whiskers, whose caps reached into the ascenders here
    text!(ax2, 0.97, 0.05;
        text = @sprintf("Δχ² = %.1f in favour of ΛCDM", Δχ²),
        space = :relative, align = (:right, :bottom), fontsize = 15, color = PALETTE.blue,)

    Legend(fig[1, 1], [l_hi, l_lo, l_lcdm, l_eds],
        ["High-z (Supernova Cosmology Project)", "Low-z (Calán/Tololo)",
            L"Flat $\Lambda$CDM, fitted",
            rich("Einstein–de Sitter (Ω", subscript("M"), " = 1)"),],
        orientation = :horizontal, framevisible = false, labelsize = 15, nbanks = 2,
        colgap = 18,)

    rowsize!(fig.layout, 2, Relative(0.55))
    rowgap!(fig.layout, 1, 10)
    rowgap!(fig.layout, 2, 8)
    linkxaxes!(ax1, ax2)
    println("\nwrote ", savefigure(fig, FIGURES, "supernova_cosmology"))
end

main()
