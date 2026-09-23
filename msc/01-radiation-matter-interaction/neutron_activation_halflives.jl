# Half-lives of three neutron-activation products from counts accumulated in
# successive equal acquisition intervals,
#
#   ln Nᵢ = a − λ tᵢ,    Var[ln Nᵢ] = 1/Nᵢ  (Poisson),
#
# fitted by weighted least squares on every recorded point, and the NaI(Tl)
# energy calibration that identifies the ¹²⁸I photopeak.
#
#   ²⁷Al(n,γ)²⁸Al,  ²⁶Mg(n,γ)²⁷Mg,  ¹²⁷I(n,γ)¹²⁸I in the NaI(Tl) crystal itself
#
# Ported from ReactiiNeutronice.jl and FitActivareNaITl.jl on the `legacy`
# branch. Both linearised the decay law as ln(N₀/Nᵢ) = λ(tᵢ − t₀) with the first
# count as reference; the first fitted unweighted, the second with weights Nᵢ.
#
# Why ln N and not ln(N₀/N): every ln(N₀/Nᵢ) shares N₀, so the points have the
# covariance Dᵢⱼ + 1/N₀ with D = diag(1/Nᵢ). With a free intercept in the model
# the common term 1/N₀ is absorbed by the intercept, and the slope and its
# variance are exactly those of a weighted fit to N₂ … Nₙ alone — the reference
# count, the most precise of the series, drops out of λ. `main` asserts this
# identity numerically. Fitting ln N keeps all points, with independent errors.
#
# Equal acquisition intervals matter: the factor (1 − e^{−λΔ}) relating counts
# to activity is then common to every point and changes only the intercept.
#
# The calibration is an unweighted straight line through four lines. The
# uncertainty of the peak energy is the full quadratic form [1, c]ᵀ C [1, c];
# intercept and slope are anticorrelated (ρ = −0.96), so dropping the covariance
# term nearly doubles it.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, LinearAlgebra, Distributions
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "radiation_matter_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Energies of the NaI(Tl) calibration lines in keV, as recorded."
const CAL_ENERGY = [511.0, 1173.0, 1274.0, 1332.0]
"Channels of the calibration lines."
const CAL_CHANNEL = [77.0, 172.0, 185.0, 195.0]
"Channel of the photopeak whose decay is followed in the ¹²⁸I series."
const PEAK_CHANNEL = 67.0
"Energy in keV of the 2⁺ → 0⁺ transition of ¹²⁸Xe fed by ¹²⁸I β⁻ decay (doi:10.1016/j.nds.2015.09.002)."
const IODINE_GAMMA = 442.9

"""
The three series: acquisition times in s, counts per interval, and the NUBASE2020
half-life with its uncertainty in min (doi:10.1088/1674-1137/abddae).

The ¹²⁸I reference count is placed one interval before the first time of
`FitActivareNaITl.jl`, which lists 520, 780, 1040 s for the second to fourth
counts.
"""
const SERIES = [
    (name = "28Al", mass = "28", symbol = "Al", t = collect(100.0:100.0:500.0),
        N = [2985.0, 1728.0, 1116.0, 632.0, 402.0], T_ref = 2.245, σT_ref = 0.005,),
    (name = "27Mg", mass = "27", symbol = "Mg", t = collect(100.0:100.0:500.0),
        N = [1459.0, 1018.0, 1175.0, 953.0, 829.0], T_ref = 9.435, σT_ref = 0.027,),
    (name = "128I", mass = "128", symbol = "I", t = collect(260.0:260.0:1040.0),
        N = [3691.0, 3275.0, 2827.0, 2487.0], T_ref = 24.99, σT_ref = 0.02,),
]

"""
    decay_fit(t, N)

Weighted least squares of ln N = a − λt with Poisson weights, on all points.
Returns the [`LinearFit`](@ref) with `p = [a, λ]`.
"""
function decay_fit(t, N)
    all(>(0), N) || throw(ArgumentError("counts must be positive"))
    return generalised_least_squares(hcat(ones(length(t)), -t), log.(N), Diagonal(1 ./ N))
end

"""
    shared_reference_fit(t, N)

The linearisation of the original scripts, ln(N₁/Nᵢ) = b + λtᵢ for i ≥ 2, with
the full covariance 1/Nᵢ δᵢⱼ + 1/N₁ that the shared reference count induces.
Returns the [`LinearFit`](@ref) with `p = [b, λ]`.
"""
function shared_reference_fit(t, N)
    y = log.(N[1] ./ N[2:end])
    C = Diagonal(1 ./ N[2:end]) + fill(1 / N[1], length(y), length(y))
    return generalised_least_squares(hcat(ones(length(y)), t[2:end]), y, C)
end

half_life_min(λ) = log(2) / λ / 60
half_life_sigma_min(λ, σλ) = log(2) / λ^2 / 60 * σλ

function calibrate()
    fit = ordinary_least_squares(hcat(ones(length(CAL_CHANNEL)), CAL_CHANNEL), CAL_ENERGY)
    σa, σb = sqrt.(diag(fit.cov))
    v = [1.0, PEAK_CHANNEL]
    E = dot(v, fit.p)
    σE = prediction_sigma(fit, v)
    σE_diagonal = sqrt(σa^2 + (PEAK_CHANNEL * σb)^2)
    @printf("calibration E = a + b·channel, %d lines, %d dof\n", length(CAL_ENERGY),
        fit.dof)
    @printf("  a = %.1f ± %.1f keV,  b = %.3f ± %.3f keV/channel,  ρ(a,b) = %.3f\n",
        fit.p[1], σa, fit.p[2], σb, fit.cov[1, 2] / (σa * σb))
    @printf("  residual standard deviation = %.1f keV\n", sqrt(fit.χ² / fit.dof))
    @printf("peak at channel %.0f (lowest calibration channel %.0f): E = %.1f ± %.1f keV\n",
        PEAK_CHANNEL, minimum(CAL_CHANNEL), E, σE)
    @printf("  without the covariance term: ± %.1f keV\n", σE_diagonal)
    @printf("  ¹²⁸I γ line %.1f keV: difference %.1f keV = %.2f σ\n\n",
        IODINE_GAMMA, IODINE_GAMMA - E, abs(IODINE_GAMMA - E) / σE)
    return E, σE
end

function main()
    calibrate()

    results = map(SERIES) do s
        fit = decay_fit(s.t, s.N)
        λ, σλ = fit.p[2], sqrt(fit.cov[2, 2])

        # the shared-reference fit with a free intercept is the weighted fit
        # without the first point
        shared = shared_reference_fit(s.t, s.N)
        dropped = decay_fit(s.t[2:end], s.N[2:end])
        @assert isapprox(shared.p[2], dropped.p[2]; rtol = 1e-9) &&
                isapprox(shared.cov[2, 2], dropped.cov[2, 2]; rtol = 1e-9) &&
                isapprox(shared.χ², dropped.χ²; rtol = 1e-7) "shared-reference fit " *
                                                             "differs from the fit without the first point"

        T, σT = half_life_min(λ), half_life_sigma_min(λ, σλ)
        reduced = fit.χ² / fit.dof
        @printf("%-5s λ = %.3e ± %.1e s⁻¹   T½ = %s min   χ²/ν = %.2f/%d = %.2f   p = %.2g\n",
            s.name, λ, σλ, pm_string(T, σT), fit.χ², fit.dof,
            reduced,
            ccdf(Chisq(fit.dof), fit.χ²))
        reduced > 1 &&
            @printf("      σ(T½) scaled by √(χ²/ν): %.2f min\n", σT * sqrt(reduced))
        @printf("      without the first point (shared-reference form): T½ = %s min\n",
            pm_string(half_life_min(dropped.p[2]),
            half_life_sigma_min(dropped.p[2], sqrt(dropped.cov[2, 2])),))
        @printf("      NUBASE2020: %g ± %g min\n", s.T_ref, s.σT_ref)
        (; fit, T, σT, reduced)
    end

    colors = (PALETTE.blue, PALETTE.orange, PALETTE.green)
    markers = (:circle, :utriangle, :diamond)
    nuclide(s) = rich(superscript(s.mass), s.symbol)

    fig = Figure(size = (1300, 660))
    ax1 = Axis(
        fig[2, 1], xlabel = L"Acquisition time $t$ [s]", ylabel = "Counts per interval",
        yscale = log10,
        yticks = (
            [300, 500, 1000, 2000, 5000], [L"300", L"500", L"1000", L"2000", L"5000"],),)
    ax2 = Axis(fig[2, 2], ylabel = L"Half-life $T_{1/2}$ [min]", yscale = log10,
        yticks = ([2, 5, 10, 20, 50], [L"2", L"5", L"10", L"20", L"50"]),
        xticks = (1:length(SERIES), [nuclide(s) for s in SERIES]), xgridvisible = false,)

    for (k, (s, r)) in enumerate(zip(SERIES, results))
        a, λ = r.fit.p
        tf = range(0, maximum(s.t) * 1.06, length = 50)
        lines!(ax1, tf, exp.(a .- λ .* tf), color = colors[k])
        errorbars_unstroked!(
            ax1, s.t, s.N, sqrt.(s.N), color = colors[k], whiskerwidth = 10,)
        scatter!(ax1, s.t, s.N, color = colors[k], marker = markers[k])

        errorbars_unstroked!(ax2, [k], [r.T], [r.σT], color = colors[k], whiskerwidth = 16)
        scatter!(ax2, [k], [r.T], color = colors[k], marker = markers[k],
            markersize = MARKERSIZE.emphasis,)
        scatter!(ax2, [k], [s.T_ref], color = PALETTE.black, marker = :hline,
            markersize = 3 * MARKERSIZE.emphasis,)
        # the ¹²⁸I reference lies above its point, the other two on or below theirs
        below = r.T < s.T_ref
        text!(ax2, k, below ? r.T / 1.22 : r.T * 1.22;
            text = rich(
                pm_string(r.T, r.σT), "\n", it("χ"), superscript("2"), "/", it("ν"),
                @sprintf(" = %.1f", r.reduced)),
            align = (:center, below ? :top : :bottom), color = colors[k],
            fontsize = ANNOTATION_SIZE, justification = :center,)
    end
    xlims!(ax1, 0, 1120)
    ylims!(ax1, 300, 5000)
    xlims!(ax2, 0.4, length(SERIES) + 0.6)
    ylims!(ax2, 1.5, 50)

    product_keys = [[LineElement(color = colors[k]),
                        MarkerElement(color = colors[k], marker = markers[k],
                            markersize = MARKERSIZE.key, strokewidth = 1.5,),]
                    for k in eachindex(SERIES)]
    half_life_keys = [
        MarkerElement(color = :grey40, marker = :circle, markersize = MARKERSIZE.key,
            strokewidth = 1.5,),
        MarkerElement(color = PALETTE.black, marker = :hline, markersize = 2 *
                                                                           MARKERSIZE.key),]
    Legend(fig[1, 1:2], [product_keys, half_life_keys],
        [[nuclide(s) for s in SERIES], ["Fitted", "NUBASE2020"]],
        ["Counts and fit", "Half-life"], titleposition = :left, nbanks = 1,)
    colsize!(fig.layout, 2, Relative(0.32))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_activation_halflives"))
end

main()
