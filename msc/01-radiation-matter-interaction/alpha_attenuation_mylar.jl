# Energy loss of α particles through stacked Mylar absorbers. The residual energy
# ε is measured after each absorber setting x, and the constant-stopping-power
# (thin-absorber) model
#
#   ε(x) = ε₀ − S x
#
# is fitted by weighted least squares to all five measurements, each with its
# own uncertainty and independent of the others.
#
# Ported from Atenuare_alpha.jl on the `legacy` branch
# (`Julia-Workflow-FFUB/IRM_M_1/`). The original differenced the data first,
# Δε = ε(0) − ε(x), fitted the four differences unweighted with a free intercept,
# drew error bars that did not enter the fit, and never printed the
# coefficients. Differencing correlates the points through the shared ε(0);
# with that covariance written out and no intercept — Δε(0) = 0 by
# construction — the differenced fit is identical to the five-point fit above,
# which `main` asserts. A free intercept in the differenced model is a redundant
# parameter and widens the slope error by 1.7×.
#
# The abscissa is the absorber setting as recorded, which the original labelled
# "x (μm)". Read as micrometres the slope would be about three times the
# electronic stopping power of Mylar between 4 and 5 MeV, 115–135 keV µm⁻¹
# (ASTAR, Berger et al., NIST Standard Reference Database 124,
# doi:10.18434/T4NC7P), and a 1.85 MeV residual after 8 µm is not compatible
# with the ≈ 29 µm CSDA range of a 4.9 MeV α in Mylar. The setting is more
# likely a foil count; it is labelled neutrally.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, LinearAlgebra, Distributions
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "radiation_matter_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Absorber setting, as recorded."
const X = [0.0, 2.0, 4.0, 6.0, 8.0]
"Residual α energy [keV]."
const ε = [4880.0, 4080.0, 3850.0, 2340.0, 1850.0]
"Uncertainty of each residual energy [keV], as recorded."
const σε = [500.0, 700.0, 970.0, 1100.0, 1740.0]
"""
ASTAR electronic stopping power of α particles in Mylar at 5 and 4 MeV
[keV µm⁻¹]: 813.8 and 948.5 MeV cm² g⁻¹ at ρ = 1.40 g cm⁻³ (doi:10.18434/T4NC7P).
"""
const MYLAR_STOPPING = (114.0, 133.0)
"ASTAR CSDA range of a 5 MeV α in Mylar [µm]: 4.037 × 10⁻³ g cm⁻² at 1.40 g cm⁻³."
const MYLAR_RANGE_5MEV = 28.8

function main()
    n = length(X)
    fit = generalised_least_squares(hcat(ones(n), -X), ε, Diagonal(σε .^ 2))
    ε₀, S = fit.p
    σε₀, σS = sqrt.(diag(fit.cov))

    # the differenced form of the original, with the covariance the shared ε(0)
    # induces and no intercept: the same estimate, the same error, the same χ²
    x, Δε = X[2:end], ε[1] .- ε[2:end]
    C = [i == j ? σε[1]^2 + σε[i + 1]^2 : σε[1]^2 for i in eachindex(x), j in eachindex(x)]
    differenced = generalised_least_squares(reshape(x, :, 1), Δε, C)
    @assert isapprox(differenced.p[1], S; rtol = 1e-9) &&
            isapprox(differenced.cov[1, 1], σS^2; rtol = 1e-9) &&
            isapprox(differenced.χ², fit.χ²; rtol = 1e-7) "differenced fit through the " *
                                                          "origin differs from the five-point fit"
    # ... and with the free intercept the original's model implied
    free = generalised_least_squares(hcat(ones(length(x)), x), Δε, C)
    σS_free = sqrt(free.cov[2, 2])

    @printf("ε(x) = ε₀ − S x, weighted least squares on %d points, %d dof\n", n, fit.dof)
    @printf("  S  = %s keV per setting\n", pm_string(S, σS))
    @printf("  ε₀ = %s keV  (measured at x = 0: %.0f ± %.0f)\n", pm_string(ε₀, σε₀), ε[1],
        σε[1])
    @printf("  χ² = %.2f, p = %.2f\n", fit.χ², ccdf(Chisq(fit.dof), fit.χ²))
    @printf("differenced, free intercept: S = %s, intercept %s keV; slope error ×%.2f\n",
        pm_string(free.p[2], σS_free), pm_string(free.p[1], sqrt(free.cov[1, 1])),
        σS_free / σS)
    @printf("read as µm: %.0f keV/µm against ASTAR %.0f–%.0f keV/µm for Mylar at 4–5 MeV; ",
        S, MYLAR_STOPPING...)
    @printf("CSDA range of a 5 MeV α in Mylar %.1f µm\n", MYLAR_RANGE_5MEV)

    fig = Figure(size = (900, 600))
    ax = Axis(fig[2, 1], xlabel = "Absorber setting",
        ylabel = L"Residual energy $\varepsilon$ [keV]", xticks = 0:2:8,)
    xf = range(-0.3, maximum(X) + 0.3, length = 50)
    l_fit = lines!(ax, xf, ε₀ .- S .* xf, color = PALETTE.blue)
    errorbars_unstroked!(
        ax, X, ε, σε, color = PALETTE.orange, linewidth = GUIDE_WIDTH, whiskerwidth = 10,)
    l_dat = scatter!(ax, X, ε, color = PALETTE.orange)
    text!(ax, 0.97, 0.95;
        text = rich(
            rich(
                it("S"), " = ", pm_string(S, σS), " keV per setting", color = PALETTE.blue,),
            "\n", rich(it("ε"), subscript("0"), " = ", pm_string(ε₀, σε₀), " keV",
                color = PALETTE.blue,),
            "\n", rich(
                it("χ"), superscript("2"), @sprintf("/dof = %.2f/%d", fit.χ², fit.dof),
                color = PALETTE.blue,),),
        space = :relative, align = (:right, :top), justification = :right,
        fontsize = ANNOTATION_SIZE,)
    xlims!(ax, -0.4, 8.4)

    Legend(fig[1, 1], [l_dat, l_fit], ["Measured", "Weighted linear fit"])
    println("wrote ", savefigure(fig, FIGURES, "alpha_attenuation_mylar"))
end

main()
