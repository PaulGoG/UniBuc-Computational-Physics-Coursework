# Flat ΛCDM fitted to the sixty Type Ia supernovae of Perlmutter et al. (1999),
# "Measurements of Ω and Λ from 42 High-Redshift Supernovae", ApJ 517, 565,
# doi:10.1086/307221 — Tables 1 and 2: 42 high-redshift supernovae of the
# Supernova Cosmology Project and 18 low-redshift ones of the Calán/Tololo
# survey, with stretch-, K- and extinction-corrected effective B-band peak
# magnitudes whose quoted uncertainties include 0.17 mag of intrinsic dispersion.
#
# The model is m = 5 log₁₀[(1+z) ∫₀^z dz′/E(z′)] + offset with
# E² = Ω_M(1+z)³ + 1 − Ω_M. The absolute magnitude and H₀ enter only through the
# offset, so two parameters are fitted, weighted by the tabulated uncertainties.
#
# Two samples are fitted: all sixty supernovae (the paper's Fit A) and the 54 of
# its primary Fit C, which drops two stretch outliers, two residual outliers and
# two probably reddened supernovae. The paper's flat-universe results are
# Ω_M = 0.29 (+0.09 −0.08) and 0.28 (+0.09 −0.08). Its fit is not this one: it
# has four parameters (Ω_M, Ω_Λ, the magnitude zero point and the slope α of the
# width–luminosity relation), propagates the redshift uncertainty, and reads
# Ω_M off the flat line of the (Ω_M, Ω_Λ) plane. The comparison is therefore a
# consistency check of the implementation, not a reproduction.
#
# The Einstein–de Sitter universe (flat, Ω_M = 1) is fitted with the offset
# alone. Only flat models are tested here, so what Δχ² rejects is flat
# matter-only expansion; open decelerating models need the curvature term,
# which this module does not carry.
#
# The 2019 analysis (`Friedmann_Magnitude.jl` on the `legacy` branch, joint work
# with Alexandru Crăciun at JINR Dubna) fitted host-galaxy magnitudes, which are
# not standard candles; `friedmann_expansion.jl` shows what those constrain.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "cosmology_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Supernovae excluded from the paper's primary Fit C (notes column of its Tables 1 and 2)."
const FIT_C_EXCLUDED = ("1992bo", "1992br", "1994H", "1997O", "1996cg", "1996cn")

"Flat-universe Ω_M published for Fits A and C (Table 3 of the paper), statistical errors."
const PUBLISHED = (fit_a = (Ωm = 0.29, up = 0.09, down = 0.08, χ² = 98, dof = 56),
    fit_c = (Ωm = 0.28, up = 0.09, down = 0.08, χ² = 56, dof = 50),)

function report(name, f, published)
    @printf("%-22s Ω_M = %.3f ± %.3f   offset = %.3f ± %.3f mag   χ²/dof = %.1f/%d = %.2f\n",
        name, f.Ωm, f.σΩm, f.offset, f.σoffset, f.χ², f.dof, f.χ² / f.dof)
    @printf("%-22s Ω_M = %.2f (+%.2f −%.2f)                        χ²/dof = %d/%d (four parameters)\n",
        "  published", published.Ωm, published.up, published.down, published.χ²,
        published.dof)
end

function main()
    d = CSV.read(joinpath(@__DIR__, "data", "perlmutter1999_sn_ia.csv"), DataFrame)
    z, m, σ = d.z, d.m_b_effective, d.sigma_m
    high = d.sample .== "high-z"
    kept = [!(name in FIT_C_EXCLUDED) for name in d.sn_name]
    count(.!kept) == length(FIT_C_EXCLUDED) ||
        error("expected all $(length(FIT_C_EXCLUDED)) Fit C exclusions in the table")

    @printf("%d supernovae: %d high-z, %d low-z, z = %.3f–%.3f\n",
        nrow(d), count(high), count(.!high), minimum(z), maximum(z))

    fit_a = fit_flat_lcdm(z, m, σ)
    fit_c = fit_flat_lcdm(z[kept], m[kept], σ[kept])
    report("Fit A sample (60)", fit_a, PUBLISHED.fit_a)
    report("Fit C sample (54)", fit_c, PUBLISHED.fit_c)
    @printf("Ω_M error scaled by √(χ²/dof): ±%.3f (Fit A), ±%.3f (Fit C)\n",
        fit_a.σΩm * sqrt(fit_a.χ² / fit_a.dof), fit_c.σΩm * sqrt(fit_c.χ² / fit_c.dof))
    @printf("Planck 2018: Ω_m = %.3f ± %.3f\n", PLANCK2018.Ωm, PLANCK2018.σΩm)

    eds_a = fit_flat_lcdm(z, m, σ; Ωm = 1.0)
    eds_c = fit_flat_lcdm(z[kept], m[kept], σ[kept]; Ωm = 1.0)
    Δχ²_a, Δχ²_c = eds_a.χ² - fit_a.χ², eds_c.χ² - fit_c.χ²
    @printf("Einstein–de Sitter: χ²/dof = %.1f/%d (Fit A), %.1f/%d (Fit C)\n",
        eds_a.χ², eds_a.dof, eds_c.χ², eds_c.dof)
    @printf("Δχ² against flat ΛCDM, one parameter: %.1f (Fit A), %.1f (Fit C)\n",
        Δχ²_a, Δχ²_c)

    # The published values are not the target of anything above; what is checked
    # is that the same data and model class land within the published interval.
    for (f, p) in ((fit_a, PUBLISHED.fit_a), (fit_c, PUBLISHED.fit_c))
        p.Ωm - p.down <= f.Ωm <= p.Ωm + p.up ||
            error("fitted Ω_M = $(f.Ωm) outside the published interval of $(p.Ωm)")
    end

    # --- figure: the primary (Fit C) sample, its six exclusions drawn hollow ----
    fig = Figure(size = (1000, 880))
    zf = 10 .^ range(-2, 0, length = 300)
    reference(zi) = magnitude(zi, 1.0, eds_c.offset)
    samples = ((high, PALETTE.orange, :circle), (.!high, PALETTE.purple, :diamond))
    decades = ([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0],
        [L"0.01", L"0.02", L"0.05", L"0.1", L"0.2", L"0.5", L"1"],)

    function draw_samples!(ax, y)
        return map(samples) do (sel, colour, marker)
            unstroked_errorbars!(ax, z[sel], y[sel], σ[sel], color = colour)
            scatter!(ax, z[sel .& .!kept], y[sel .& .!kept], color = :white,
                marker = marker, strokecolor = colour, markersize = MARKERSIZE.dense,)
            scatter!(ax, z[sel .& kept], y[sel .& kept], color = colour, marker = marker,
                strokecolor = darker(colour), markersize = MARKERSIZE.dense,)
        end
    end

    ax1 = Axis(fig[2, 1], ylabel = L"Effective $m_B$ [mag]", xscale = log10,
        xticks = decades, xticklabelsvisible = false,)
    l_lcdm = lines!(ax1, zf, magnitude.(zf, fit_c.Ωm, fit_c.offset), color = PALETTE.blue)
    l_eds = lines!(ax1, zf, reference.(zf), color = PALETTE.red, linestyle = :dash)
    l_data = draw_samples!(ax1, m)
    text!(ax1, 0.03, 0.95;
        text = rich(
            "Ω", subscript("M"), @sprintf(" = %.2f ± %.2f,   ", fit_c.Ωm, fit_c.σΩm),
            it("χ"), superscript("2"), @sprintf("/dof = %.1f/%d", fit_c.χ², fit_c.dof)),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,)

    ax2 = Axis(fig[3, 1], xlabel = L"Redshift $z$",
        ylabel = L"$m_B - m_B^{\mathrm{EdS}}$ [mag]", xscale = log10, xticks = decades,)
    hlines!(ax2, [0.0], color = PALETTE.red, linestyle = :dash)
    lines!(ax2, zf, magnitude.(zf, fit_c.Ωm, fit_c.offset) .- reference.(zf),
        color = PALETTE.blue,)
    draw_samples!(ax2, m .- reference.(z))
    text!(ax2, 0.03, 0.93;
        text = rich("Δ", it("χ"), superscript("2"), @sprintf(" = %.1f", Δχ²_c)),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,)

    hollow = MarkerElement(color = :white, marker = :circle, strokecolor = PALETTE.black,
        strokewidth = 1.5, markersize = MARKERSIZE.dense,)
    Legend(fig[1, 1], [[l_data..., hollow], [l_lcdm, l_eds]],
        [["SCP, high z", "Calán/Tololo, low z", "Not in Fit C"],
            ["ΛCDM fit", "Einstein–de Sitter"],],
        ["Supernovae", "Flat models"]; nbanks = 2, titleposition = :left,)

    linkxaxes!(ax1, ax2)
    xlims!(ax2, 0.01, 1.0)
    rowsize!(fig.layout, 3, Auto(0.5))
    println("wrote ", savefigure(fig, FIGURES, "supernova_cosmology"))
end

main()
