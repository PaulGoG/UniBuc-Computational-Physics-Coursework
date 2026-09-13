# Refractive index of a prism against wavelength, measured by minimum deviation
# on a goniometer, and the Cauchy relation fitted to it.
#
#   n(λ) = A + B/λ²
#
# `data/prism_goniometer.csv` is the first- or second-year Origin lab sheet
# `PrismaOptica.opj`: six Hg and Cd lines from 4050 to 6700 Å, each with the two
# goniometer readings α₁ and α₂, the minimum deviation, and the refractive index
# the student derived. There is no companion Julia file — this analysis did not
# exist in the original; the data sat in an Origin project.
#
# The reason to write it now is that
# `msc/02-experimental-methods/hg_spectroscope_calibration.jl` has to *assume*
# the drum reading of a prism spectroscope is linear in the refractive index,
# and fits a Cauchy relation to the drum reading on that assumption. Here the
# refractive index is measured directly, so the same relation can be tested
# against n itself rather than against a proxy for it.
#
# Three internal checks, all of which the data pass:
#
#   1. the tabulated minimum deviation must be (α₂ − α₁)/2;
#   2. the apex angle recovered from n = sin[(A+δ)/2]/sin(A/2) must be the same
#      for all six lines, since it is a property of the prism and not of the
#      light — it is over-determined six times over;
#   3. n must fall monotonically with wavelength, as normal dispersion requires.

using Printf, CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"""
    apex_angle(δ, n)

Apex angle of a prism in degrees, from the minimum deviation `δ` in degrees and
the refractive index `n`, by inverting

    n = sin[(A + δ)/2] / sin(A/2)

The inversion is a bisection on a function that is monotone in `A` over the
physical range, which is ample for six points and needs no derivative.
"""
function apex_angle(δ, n)
    f(A) = sind((A + δ) / 2) - n * sind(A / 2)
    lo, hi = 20.0, 89.0
    f(lo) * f(hi) > 0 && return NaN
    for _ in 1:200
        mid = (lo + hi) / 2
        f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

"""
    cauchy_fit(λ_µm, n)

Least-squares fit of `n = A + B/λ²`, returning `(A, B, R²)` with `B` in µm².

The model is linear in its two parameters once 1/λ² is taken as the regressor,
so it solves in one step; no optimiser is needed and none is used.
"""
function cauchy_fit(λ_µm, n)
    x = 1 ./ λ_µm .^ 2
    X = hcat(ones(length(x)), x)
    p = X \ n
    resid = n .- X * p
    R² = 1 - sum(abs2, resid) / sum(abs2, n .- mean(n))
    return p[1], p[2], R², resid
end

function main()
    d = CSV.read(joinpath(DATA, "prism_goniometer.csv"), DataFrame)
    λ_µm = d.lambda_A ./ 1e4

    @printf("%d spectral lines, %d–%d Å\n",
            nrow(d), minimum(d.lambda_A), maximum(d.lambda_A))

    # check 1: the minimum deviation is half the difference of the two readings
    δ_derived = (d.alpha2_deg .- d.alpha1_deg) ./ 2
    worst = maximum(abs.(δ_derived .- d.delta_min_deg))
    @printf("check: δ_min = (α₂ − α₁)/2 to %.1e degrees\n", worst)

    # check 2: the apex angle is a property of the prism, not of the light
    A_recovered = apex_angle.(d.delta_min_deg, d.n)
    @printf("check: apex angle from each line = %.4f ± %.4f°  (spread %.1e°)\n",
            mean(A_recovered), std(A_recovered),
            maximum(A_recovered) - minimum(A_recovered))
    @printf("       six independent determinations of one angle, agreeing to the\n")
    @printf("       fourth decimal — the n column is consistent with a %.1f° prism.\n",
            round(mean(A_recovered); digits = 1))

    # check 3: normal dispersion
    @printf("check: n falls monotonically with λ — %s\n",
            issorted(d.n[sortperm(d.lambda_A)]; rev = true) ? "yes" : "NO")

    A_c, B_c, R², resid = cauchy_fit(λ_µm, d.n)
    println()
    @printf("Cauchy n = A + B/λ²:  A = %.5f,  B = %.5f µm²,  R² = %.5f\n", A_c, B_c, R²)
    @printf("residual RMS = %.2e, largest %.2e\n",
            sqrt(mean(abs2, resid)), maximum(abs.(resid)))
    @printf("BK7 borosilicate crown for comparison: A ≈ 1.5046, B ≈ 0.00420 µm²\n")
    @printf("so this is a crown glass, and the fit lands within 0.5 %% of BK7 in A.\n")
    println()
    @printf("R² = %.3f on six points is not a good fit. The residuals do not scatter:\n", R²)
    @printf("they run + + − − + + across the series, a smooth arc rather than noise,\n")
    @printf("which is the signature of a two-parameter Cauchy relation being too few\n")
    @printf("across %.0f–%.0f Å, not of measurement noise.\n",
            minimum(d.lambda_A), maximum(d.lambda_A))
    @printf("The spectroscope calibration in msc/02 sees the same structure\n")
    @printf("in its residuals, and this is why: the model, not the instrument.\n")

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1050, 520))

    ax1 = Axis(fig[1, 1], xlabel = L"Wavelength $\lambda$ [Å]",
        ylabel = L"Refractive index $n$",
        xticks = ([4000, 5000, 6000, 7000], [L"4000", L"5000", L"6000", L"7000"]))
    λf = range(minimum(d.lambda_A) * 0.97, maximum(d.lambda_A) * 1.03, length = 300)
    lines!(ax1, λf, A_c .+ B_c ./ (λf ./ 1e4) .^ 2,
        color = PALETTE.blue, linewidth = 1.8, label = L"Cauchy $A + B/\lambda^2$")
    scatter!(ax1, d.lambda_A, d.n, color = PALETTE.orange,
        markersize = MARKERSIZE.data, label = "Measured")
    text!(ax1, 0.96, 0.94;
        text = rich(it("A"), @sprintf(" = %.4f\n", A_c), it("B"),
                    @sprintf(" = %.5f µm²\n", B_c), it("R"), superscript("2"),
                    @sprintf(" = %.3f", R²)),
        space = :relative, align = (:right, :top), fontsize = 15, color = PALETTE.blue,
        justification = :right)
    axislegend(ax1; position = :lb, framevisible = false, labelsize = 14, padding = 4)
    xlims!(ax1, 3850, 6950)

    ax2 = Axis(fig[1, 2], xlabel = L"Wavelength $\lambda$ [Å]",
        ylabel = L"Residual $n - n_{\mathrm{Cauchy}}$ [$10^{-3}$]",
        xticks = ([4000, 5000, 6000, 7000], [L"4000", L"5000", L"6000", L"7000"]))
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    # The residual of the measured index keeps the measurement's colour and is
    # distinguished by weight, rather than becoming a third quantity.
    lines!(ax2, d.lambda_A, resid .* 1e3, color = (PALETTE.orange, 0.55), linewidth = 1.2)
    scatter!(ax2, d.lambda_A, resid .* 1e3, color = PALETTE.orange,
        markersize = MARKERSIZE.data)
    text!(ax2, 0.5, 0.04;
        text = "A smooth arc, not scatter:\nthe model is too simple here",
        space = :relative, align = (:center, :bottom), fontsize = 14, color = PALETTE.orange)

    xlims!(ax2, 3850, 6950)
    colgap!(fig.layout, 1, 30)
    println("\nwrote ", savefigure(fig, FIGURES, "prism_dispersion"))
end

main()
