# Wavelength calibration of a prism spectroscope against the Hg I emission
# spectrum, using the Cauchy dispersion relation
#
#   n(λ) = A + B/λ²
#
# under the working assumption that the drum reading is linear in the refractive
# index, which holds near minimum deviation.
#
# **This is MSc coursework, not BSc.** It was filed under geometrical optics
# until the source file turned up as `CALIBRARE.jl` in the Metode Experimentale
# in Fizica directory of MSc year 1, alongside the submitted reports for the
# interferometry and polarimetry labs of the same course. Its own title line,
# "Etalonarea spectroscopului cu lampa de Hg", is a standard experimental-methods
# lab rather than an optical-design exercise. Moved here.
#
# Ported from Calibrare_Hg.jl. Three corrections.
#
#   1. **A transcription error in the reference data.** The original listed
#      5789.66 Å next to 5790.65 Å at adjacent drum divisions 35 and 36. No Hg I
#      line exists at 5789.66; the second member of the yellow doublet is
#      5769.60 Å, and a 1 Å splitting across one division is impossible at the
#      ~13 Å/division local dispersion. Corrected below.
#   2. The model is **linear in its parameters**, so Levenberg-Marquardt with a
#      finite-difference Jacobian started from p0 = [1, 1] — nine orders of
#      magnitude from B ≈ 2×10⁹ — was gratuitous. Solved here in one least-
#      squares step.
#   3. The original computed A and B and then discarded them: no `display`, the
#      `savefig` commented out, and the fitted values never printed. Running it
#      as a script produced nothing at all. Residuals and standard errors are
#      reported here, which is the actual deliverable of an instrument
#      calibration.
#
# Regressing drum reading on wavelength is the statistically correct direction:
# λ is the reference and the drum reading carries the measurement error.

using Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Hg I air wavelengths in ångström, as read on the spectroscope."
const λ = [7091.99, 7081.88, 6907.16, 6716.17, 5790.65, 5769.60,
           5460.74, 4916.04, 4358.35, 4077.81, 4046.56]
"Drum divisions."
const x = [2.0, 7, 15, 17, 35, 36, 60, 85, 120, 180, 200]

"""
    cauchy_fit(λ, x)

Least-squares fit of x = A + B/λ². Returns the coefficients, their standard
errors, the residuals and the coefficient of determination.
"""
function cauchy_fit(λ, x)
    M = hcat(ones(length(λ)), 1 ./ λ .^ 2)
    p = M \ x
    residuals = x .- M * p
    dof = length(x) - 2
    s² = sum(abs2, residuals) / dof
    covariance = s² * inv(M' * M)
    R² = 1 - sum(abs2, residuals) / sum(abs2, x .- mean(x))
    return p, sqrt.(diag(covariance)), residuals, R²
end

function main()
    p, σ, residuals, R² = cauchy_fit(λ, x)
    @printf("A = %10.4f ± %.4f divisions\n", p[1], σ[1])
    @printf("B = %10.4e ± %.2e division·Å²\n", p[2], σ[2])
    @printf("R² = %.6f,  RMS residual = %.3f divisions\n",
            R², sqrt(mean(abs2, residuals)))
    @printf("largest residual %.3f divisions at λ = %.2f Å\n",
            maximum(abs, residuals), λ[argmax(abs.(residuals))])

    # what the uncorrected 5789.66 value would have done
    λ_original = copy(λ); λ_original[6] = 5789.66
    _, _, res_original, R²_original = cauchy_fit(λ_original, x)
    @printf("with the original value 5789.66 Å: R² = %.6f, RMS residual = %.3f\n",
            R²_original, sqrt(mean(abs2, res_original)))

    fig = Figure(size = (900, 520))

    ax1 = Axis(fig[2, 1], ylabel = "Drum reading [divisions]")
    λf = range(minimum(λ) * 0.97, maximum(λ) * 1.03, length = 400)
    l_fit = lines!(ax1, λf, p[1] .+ p[2] ./ λf .^ 2, color = PALETTE.blue, linewidth = 1.5)
    l_dat = scatter!(ax1, λ, x, color = PALETTE.orange, markersize = 10)
    hidexdecorations!(ax1, grid = false)
    text!(ax1, 0.97, 0.93;
        text = @sprintf("A = %.2f,  B = %.3e,  R² = %.5f", p[1], p[2], R²),
        space = :relative, align = (:right, :top), fontsize = 15)

    ax2 = Axis(fig[3, 1], xlabel = L"Wavelength $\lambda$ [Å]",
        ylabel = "Residual [div]")
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    stem!(ax2, λ, residuals, color = PALETTE.orange,
          stemcolor = PALETTE.orange, markersize = 9)

    linkxaxes!(ax1, ax2)
    Legend(fig[1, 1], [l_dat, l_fit],
        ["Hg I lines", L"Cauchy $x = A + B/\lambda^2$"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24)

    rowsize!(fig.layout, 2, Relative(0.62))
    rowgap!(fig.layout, 8)
    path = savefigure(fig, FIGURES, "hg_spectroscope_calibration")
    println("wrote ", path)
end

main()
