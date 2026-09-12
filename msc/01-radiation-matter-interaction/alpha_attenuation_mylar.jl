# Energy loss of α particles through stacked Mylar absorber foils. The residual
# energy ε is measured after each absorber thickness and Δε = ε₀ - ε(x) fitted
# linearly, which is the constant-stopping-power (thin-absorber) approximation
# to the electronic stopping power of Mylar.
#
# Ported from Atenuare_alpha.jl. Corrections:
#
#   1. **The fit was unweighted** although a per-point uncertainty was computed
#      and drawn as error bars — so the bars displayed did not enter the χ².
#   2. The uncertainties are correlated: every Δε shares the same reference
#      measurement ε₀, so treating the points as independent understates the
#      slope uncertainty. The covariance is built explicitly below.
#   3. The model has a free intercept, but Δε(0) = 0 by construction. The
#      intercept is now fitted, reported, and tested against zero rather than
#      ignored.
#   4. `stderror` was never called and the fitted coefficients were never
#      printed; the script's only output was a figure.
#
# A caveat on the abscissa, which the original labelled "x (μm)": the fitted
# slope is about 2.5× the tabulated electronic stopping power of Mylar in this
# energy range, and a 1.85 MeV residual after 8 μm is hard to reconcile with the
# ≈30 μm range of a 4.9 MeV α in Mylar. The abscissa is more likely a foil
# count, or an areal thickness, than a length in micrometres. It is left as
# supplied and labelled neutrally.

using Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Absorber setting, as recorded."
const X = [0.0, 2.0, 4.0, 6.0, 8.0]
"Residual α energy in keV."
const ε = [4880.0, 4080.0, 3850.0, 2340.0, 1850.0]
"Peak width / uncertainty in keV."
const σε = [500.0, 700.0, 970.0, 1100.0, 1740.0]

function main()
    x = X[2:end]
    Δε = ε[1] .- ε[2:end]
    # full covariance: every point shares ε₀, so the off-diagonal is σ₀²
    C = [i == j ? σε[1]^2 + σε[i+1]^2 : σε[1]^2 for i in eachindex(x), j in eachindex(x)]

    M = hcat(ones(length(x)), x)
    W = inv(C)
    p = (M' * W * M) \ (M' * W * Δε)
    cov_p = inv(M' * W * M)
    σp = sqrt.(diag(cov_p))
    residuals = Δε .- M * p
    χ² = residuals' * W * residuals

    @printf("intercept = %8.1f ± %.1f keV   (expected 0; that is %.2fσ away)\n",
            p[1], σp[1], abs(p[1]) / σp[1])
    @printf("slope     = %8.1f ± %.1f keV per absorber unit\n", p[2], σp[2])
    @printf("χ² = %.2f on %d dof\n", χ², length(x) - 2)
    @printf("\nif the abscissa were µm this would be %.0f keV/µm; the tabulated\n", p[2])
    @printf("electronic stopping power of Mylar near 4 MeV is ≈150 keV/µm, so the\n")
    @printf("abscissa is probably not a length in micrometres.\n")

    fig = Figure(size = (800, 480))
    ax = Axis(fig[2, 1], xlabel = "Absorber setting", ylabel = L"$\Delta\varepsilon$ [keV]")
    xf = range(0, maximum(x) * 1.08, length = 100)
    l_fit = lines!(ax, xf, p[1] .+ p[2] .* xf, color = PALETTE.blue, linewidth = 1.6)
    errorbars!(ax, x, Δε, sqrt.(diag(C)), color = PALETTE.orange, whiskerwidth = 10)
    l_dat = scatter!(ax, x, Δε, color = PALETTE.orange, markersize = 11)
    scatter!(ax, [0.0], [0.0], color = PALETTE.black, markersize = 11, marker = :diamond)
    text!(ax, 0.03, 0.93;
        text = @sprintf("slope = %.0f ± %.0f keV/unit\nintercept = %.0f ± %.0f keV",
                        p[2], σp[2], p[1], σp[1]),
        space = :relative, align = (:left, :top), fontsize = 15)

    Legend(fig[1, 1], [l_dat, l_fit], ["Measured, correlated errors", "Weighted linear fit"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("wrote ", savefigure(fig, FIGURES, "alpha_attenuation_mylar"))
end

main()
