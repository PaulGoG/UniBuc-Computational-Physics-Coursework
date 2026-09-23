# Refractive index of a prism against wavelength, measured by minimum deviation
# on a goniometer, and the Cauchy relation fitted to it,
#
#   n(λ) = A + B/λ².
#
# `data/prism_goniometer.csv` holds six lines from 4050 to 6700 Å with the two
# goniometer readings α₁ and α₂ at minimum deviation on either side of the
# straight-through direction, δ_min = (α₂ − α₁)/2, and the index
# n = sin[(A_p + δ_min)/2] / sin(A_p/2) computed at the time. The readings are
# from a first-year optics laboratory of mine; the file was assembled from that
# record and no code accompanied it.
#
# Checks made in `main`: δ_min is (α₂ − α₁)/2 to round-off; the n column is
# reproduced from δ_min with an apex angle of 59.9° to 3 × 10⁻⁶, which fixes the
# apex angle the column was computed with and checks its arithmetic — it is not
# an independent measurement of the prism, since n was derived from δ_min and
# A_p in the first place; and n falls with λ.
#
# The two-parameter Cauchy fit is reported with the standard errors of its
# coefficients and its residuals are expressed, through dn/dδ, as the error of
# the minimum-deviation setting they would correspond to. A three-parameter
# fit is made for comparison; it lowers the residuals but returns B < 0, which
# no glass has, so the residual structure is read as measurement error in the
# setting rather than as a failure of the dispersion model. The index is also
# compared with N-BK7 from the Schott Sellmeier coefficients, as a crown glass
# of the same class, not as an identification.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Apex angle [deg] that reproduces the n column of the data file from its δ_min column."
const APEX = 59.9
"Tolerance on the reproduction of the n column, which carries five decimals."
const INDEX_TOLERANCE = 1e-5
"Goniometer reading resolution [deg]."
const READING_RESOLUTION = 0.01
"""
Sellmeier coefficients of Schott N-BK7 for λ in µm,
n² − 1 = Σ Bᵢ λ² / (λ² − Cᵢ) (Schott optical glass data sheet, N-BK7).
"""
const BK7 = (B = (1.03961212, 0.231792344, 1.01046945),
    C = (0.00600069867, 0.0200179144, 103.560653),)

"Refractive index of a prism of apex angle `A` [deg] at minimum deviation `δ` [deg]."
prism_index(A, δ) = sind((A + δ) / 2) / sind(A / 2)

"dn/dδ [deg⁻¹] of `prism_index` at the same arguments."
index_slope(A, δ) = cosd((A + δ) / 2) / (2 * sind(A / 2)) * π / 180

"N-BK7 refractive index at `λ_µm` from the Sellmeier equation."
bk7_index(λ_µm) = sqrt(1 + sum(BK7.B[i] * λ_µm^2 / (λ_µm^2 - BK7.C[i]) for i in 1:3))

"""
    cauchy_fit(λ_µm, n; terms = 2)

Least-squares fit of `n = A + B/λ² (+ C/λ⁴)`, linear in its coefficients.
Returns the coefficients, their standard errors from the residual variance on
`length(n) − terms` degrees of freedom, the residuals and the residual standard
deviation.
"""
function cauchy_fit(λ_µm, n; terms = 2)
    terms in (2, 3) || throw(ArgumentError("terms must be 2 or 3, got $terms"))
    X = hcat(ones(length(n)), 1 ./ λ_µm .^ 2, (1 ./ λ_µm .^ 4)[:, 1:(terms - 2)])
    p = X \ n
    r = n .- X * p
    dof = length(n) - terms
    s² = sum(abs2, r) / dof
    σ = sqrt.(diag(s² * inv(X' * X)))
    return p, σ, r, sqrt(s²)
end

function main()
    d = CSV.read(joinpath(DATA, "prism_goniometer.csv"), DataFrame)
    λ_µm = d.lambda_A ./ 1e4
    @printf("%d lines, %d–%d Å\n", nrow(d), minimum(d.lambda_A), maximum(d.lambda_A))

    δ_derived = (d.alpha2_deg .- d.alpha1_deg) ./ 2
    worst_δ = maximum(abs.(δ_derived .- d.delta_min_deg))
    worst_n = maximum(abs.(prism_index.(APEX, d.delta_min_deg) .- d.n))
    @printf("δ_min = (α₂ − α₁)/2 to %.1e°; n column reproduced with A_p = %.1f° to %.1e\n",
        worst_δ, APEX, worst_n)
    worst_δ < 1e-10 || error("δ_min is not half the difference of the readings")
    worst_n < INDEX_TOLERANCE || error("the n column is not prism_index($APEX°, δ_min)")
    issorted(d.n[sortperm(d.lambda_A)]; rev = true) || error("n does not fall with λ")

    slope = index_slope.(APEX, d.delta_min_deg)
    @printf("dn/dδ = %.4f per degree; a reading resolution of %.2f° in each α is %.1e in n\n",
        mean(slope), READING_RESOLUTION, READING_RESOLUTION / sqrt(2) * mean(slope))

    p2, σ2, r2, s2 = cauchy_fit(λ_µm, d.n)
    p3, σ3, r3, s3 = cauchy_fit(λ_µm, d.n; terms = 3)
    rms2, rms3 = sqrt(mean(abs2, r2)), sqrt(mean(abs2, r3))
    arcmin(residual) = 60 * residual / mean(slope)
    @printf("\nCauchy, two terms:   A = %.4f ± %.4f,  B = %.5f ± %.5f µm² (%.0f %%);  residual rms %.1e = %.1f′ in δ_min, largest %.1f′\n",
        p2[1], σ2[1], p2[2], σ2[2], 100σ2[2] / p2[2], rms2,
        arcmin(rms2), arcmin(maximum(abs.(r2))))
    @printf("Cauchy, three terms: A = %.4f ± %.4f,  B = %.5f ± %.5f µm²,  C = %.6f ± %.6f µm⁴;  residual rms %.1e = %.1f′\n",
        p3[1], σ3[1], p3[2], σ3[2], p3[3], σ3[3], rms3, arcmin(rms3))
    println("residuals of the two-term fit, ×10³: ", join(round.(r2 .* 1e3; digits = 2), ", "))
    Δbk7 = d.n .- bk7_index.(λ_µm)
    @printf("N-BK7 at the same wavelengths: n − n_BK7 from %+.1e to %+.1e; its two-term Cauchy coefficients over this range are A = %.4f, B = %.5f µm²\n",
        minimum(Δbk7), maximum(Δbk7), cauchy_fit(λ_µm, bk7_index.(λ_µm))[1]...)

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1400, 620))
    ticks = ([4000, 5000, 6000, 7000], [L"4000", L"5000", L"6000", L"7000"])
    ax1 = Axis(
        fig[2, 1], xlabel = L"Wavelength $\lambda$ [Å]", ylabel = L"Refractive index $n$",
        xticks = ticks,)
    λf = range(3900, 6900, length = 300)
    l_fit = lines!(ax1, λf, p2[1] .+ p2[2] ./ (λf ./ 1e4) .^ 2, color = PALETTE.blue)
    l_bk7 = lines!(
        ax1, λf, bk7_index.(λf ./ 1e4), color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    l_dat = scatter!(ax1, d.lambda_A, d.n, color = PALETTE.orange)
    text!(ax1, 0.96, 0.94;
        text = rich(it("A"), @sprintf(" = %.4f ± %.4f\n", p2[1], σ2[1]), it("B"),
            @sprintf(" = %.5f ± %.5f µm²", p2[2], σ2[2])),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue, justification = :right,)
    xlims!(ax1, 3850, 6950)

    ax2 = Axis(fig[2, 2], xlabel = L"Wavelength $\lambda$ [Å]",
        ylabel = L"Residual $n - n_{\mathrm{Cauchy}}$ [$10^{-3}$]", xticks = ticks,)
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    band = 1e3 * READING_RESOLUTION / sqrt(2) * mean(slope)
    hspan!(ax2, -band, band, color = (PALETTE.black, 0.08))
    scatterlines!(
        ax2, d.lambda_A, r2 .* 1e3, color = PALETTE.orange, linewidth = GUIDE_WIDTH,)
    text!(ax2, 0.97, 0.04;
        text = rich(
            @sprintf("rms %.1f × 10⁻³, %.0f′ in the deviation setting\n", rms2 * 1e3,
                arcmin(rms2)),
            "Band: reading resolution",),
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,
        justification = :right, color = PALETTE.orange,)
    xlims!(ax2, 3850, 6950)

    Legend(fig[1, 1:2], [l_dat, l_fit, l_bk7],
        ["Measured", L"Cauchy $A + B/\lambda^2$", "N-BK7, Sellmeier"],)
    println("\nwrote ", savefigure(fig, FIGURES, "prism_dispersion"))
end

main()
