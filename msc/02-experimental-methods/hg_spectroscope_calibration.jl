# Wavelength calibration of a prism spectroscope against the Hg I emission
# spectrum, using the Cauchy dispersion relation
#
#   n(λ) = A + B/λ²
#
# under the working assumption that the drum reading is linear in the refractive
# index, which holds near minimum deviation.
#
# Ported from Calibrare_Hg.jl in Julia-Workflow-FFUB/Single_Files/ on the
# `legacy` branch, which lists 5789.66 Å beside 5790.65 Å at adjacent drum
# divisions 35 and 36. There is no Hg I line at 5789.66 Å: the yellow doublet
# is 5769.60 and 5790.66 Å (NIST Atomic Spectra Database, doi:10.18434/T4W30F),
# and with the drum reading increasing towards shorter wavelengths division 36
# is the 5769.60 Å member. The value is corrected; it moves the fit by less
# than its own scatter. The model is linear in its coefficients, so it is
# solved in one least-squares step rather than with the Levenberg–Marquardt
# iteration from p₀ = [1, 1] of the original, which also never printed A and B
# nor displayed its figure. Regressing the drum reading on the wavelength is
# the right direction: λ is the reference and the reading carries the error.
#
# The fitted dispersion also says how far apart two lines should read. The
# far-red pair 7091.99 and 7081.88 Å sits at divisions 2 and 7, five divisions
# for 10 Å where the fit gives a quarter of one; those readings are kept as
# recorded and the discrepancy is printed.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Hg I air wavelengths [Å] of the lines read on the spectroscope."
const λ = [7091.99, 7081.88, 6907.16, 6716.17, 5790.65, 5769.60,
    5460.74, 4916.04, 4358.35, 4077.81, 4046.56,]
"Drum readings [div], as recorded."
const x = [2.0, 7, 15, 17, 35, 36, 60, 85, 120, 180, 200]
"The wavelength the original listed for division 36 [Å]."
const λ_ORIGINAL_36 = 5789.66

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
    @printf("A = %10.4f ± %.4f div\n", p[1], σ[1])
    @printf("B = %10.4e ± %.2e div Å²\n", p[2], σ[2])
    @printf("R² = %.6f,  RMS residual = %.3f divisions\n",
        R², sqrt(mean(abs2, residuals)))
    @printf("largest residual %.3f divisions at λ = %.2f Å\n",
        maximum(abs, residuals), λ[argmax(abs.(residuals))])

    # what the uncorrected 5789.66 value would have done
    λ_original = copy(λ)
    λ_original[6] = λ_ORIGINAL_36
    _, _, res_original, R²_original = cauchy_fit(λ_original, x)
    @printf("with the original value %.2f Å: R² = %.6f, RMS residual = %.3f\n",
        λ_ORIGINAL_36, R²_original, sqrt(mean(abs2, res_original)))

    # the drum spacing of each adjacent pair against the fitted dispersion
    println("\nadjacent lines: recorded spacing against the fitted dispersion 2B Δλ/λ³")
    for i in 1:(length(λ) - 1)
        predicted = 2 * p[2] * (λ[i] - λ[i + 1]) / ((λ[i] + λ[i + 1]) / 2)^3
        @printf("  %.2f–%.2f Å: %5.1f div recorded, %5.2f div from the fit\n",
            λ[i], λ[i + 1], x[i + 1] - x[i], predicted)
    end

    fig = Figure(size = (900, 760))

    ax1 = Axis(fig[2, 1], ylabel = "Drum reading [div]")
    λf = range(minimum(λ) * 0.97, maximum(λ) * 1.03, length = 400)
    l_fit = lines!(ax1, λf, p[1] .+ p[2] ./ λf .^ 2, color = PALETTE.blue)
    l_dat = scatter!(ax1, λ, x, color = PALETTE.orange)
    hidexdecorations!(ax1, grid = false)
    text!(ax1, 0.97, 0.93;
        text = rich(
            it("A"), replace(@sprintf(" = %.1f ± %.1f div\n", p[1], σ[1]), "-" => "−"),
            it("B"), " = ", rsci(p[2]; digits = 3), " ± ", rsci(σ[2]; digits = 2), " div Å²\n",
            it("R"), superscript("2"), @sprintf(" = %.3f", R²)),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue, justification = :right,)

    ax2 = Axis(fig[3, 1], xlabel = L"Wavelength $\lambda$ [Å]", ylabel = "Residual [div]")
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    stem!(ax2, λ, residuals, color = PALETTE.orange, stemcolor = PALETTE.orange,
        trunkcolor = :transparent, stemwidth = GUIDE_WIDTH,)
    text!(ax2, 0.97, 0.06;
        text = @sprintf("rms %.1f div on a 2–200 div range", sqrt(mean(abs2, residuals))),
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.orange,)

    linkxaxes!(ax1, ax2)
    Legend(fig[1, 1], [l_dat, l_fit], ["Hg I lines", L"Cauchy $x = A + B/\lambda^2$"])
    rowsize!(fig.layout, 2, Relative(0.62))
    path = savefigure(fig, FIGURES, "hg_spectroscope_calibration")
    println("wrote ", path)
end

main()
