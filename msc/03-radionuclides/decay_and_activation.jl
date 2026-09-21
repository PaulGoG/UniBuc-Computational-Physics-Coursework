# Two time-dependent activity problems.
#
#   (a) Two-member Bateman chain, ²³⁸U → ²³⁴Th:
#         Λ₁(t) = Λ₀e^{-λ₁t}
#         Λ₂(t) = Λ₀ λ₂/(λ₂-λ₁) (e^{-λ₁t} - e^{-λ₂t}),  maximum at
#         t_m = ln(λ₂/λ₁)/(λ₂-λ₁)
#
#   (b) ²⁷Al(n,γ)²⁸Al under a pulsed neutron flux: irradiation and decay
#       alternating with period τ, solved analytically per cycle.
#
# Ported from Radionuclizi_4_1.jl and Radionuclizi_4_2.jl.
#
# Corrections:
#
#   1. **The Bateman figure was degenerate.** The time axis was built from
#      `T_scalare = max(T½₁, T½₂)` — the ²³⁸U half-life, 4.468 × 10⁹ yr — giving
#      a step of 4.5 × 10⁸ yr, while the ²³⁴Th transient peaks at about 2 years.
#      The entire interesting region, the vertical marker and the printed
#      `round(t_m/T_scalare, digits=2) = 0.0` all collapsed onto x = 0. A
#      logarithmic time axis is used here, which shows both scales at once.
#   2. **A wrong half-life.** ²⁸Al was given as `2.3*60` s; the evaluated value
#      is 2.245 min = 134.7 s, a 2.4 % error carried straight into λ.
#   3. `savefig` interpolated a `LaTeXString` into the filename, so the generated
#      name would have contained `$^{238}\textrm{U}$` — broken as soon as it was
#      uncommented.
#   4. The activation code recomputed its cycle sum Cₙ from scratch at every one
#      of 100 time points per cycle, an O(n) inner sum for a quantity constant
#      within the cycle, and located a plotting value by float equality after
#      `round` on a time vector rather than evaluating the closed form it had.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

const YEAR = 3.15576e7                       # s
"²³⁸U and ²³⁴Th half-lives in seconds."
const T_U238 = 4.468e9 * YEAR
const T_TH234 = 24.10 * 86400
"²⁸Al half-life in seconds — the evaluated value, not the original 2.3 min."
const T_AL28 = 2.245 * 60

decay_constant(T½) = log(2) / T½

"Daughter activity of a two-member chain, normalised to the parent's initial activity."
function bateman(t, λ₁, λ₂)
    return λ₂ / (λ₂ - λ₁) * (exp(-λ₁ * t) - exp(-λ₂ * t))
end

"Time of maximum daughter activity."
bateman_peak(λ₁, λ₂) = log(λ₂ / λ₁) / (λ₂ - λ₁)

"""
    pulsed_activity(t, λ, τ, K)

Activity under a flux alternating between on and off every `τ` seconds, with
saturation activity `K`. Solved cycle by cycle in closed form.
"""
function pulsed_activity(t, λ, τ, K)
    cycle = floor(Int, t / τ)
    phase = t - cycle * τ
    A = 0.0
    for c in 0:(cycle - 1)                    # completed cycles
        A = isodd(c) ? A * exp(-λ * τ) : K + (A - K) * exp(-λ * τ)
    end
    return iseven(cycle) ? K + (A - K) * exp(-λ * phase) : A * exp(-λ * phase)
end

function main()
    λ₁, λ₂ = decay_constant(T_U238), decay_constant(T_TH234)
    t_m = bateman_peak(λ₁, λ₂)
    @printf("²³⁸U T½ = %.3e yr, ²³⁴Th T½ = %.2f d\n", T_U238 / YEAR, T_TH234 / 86400)
    @printf("λ₂/λ₁ = %.3e — the daughter equilibrates %.0e times faster than the parent\n",
        λ₂ / λ₁, λ₂ / λ₁)
    @printf("daughter maximum at t = %.3f yr = %.1f d\n", t_m / YEAR, t_m / 86400)
    @printf("  the original time grid had a step of %.2e yr, so this sat at x = 0\n",
        T_U238 / YEAR / 10)
    @printf("secular equilibrium ratio Λ₂/Λ₁ → %.6f\n\n", λ₂ / (λ₂ - λ₁))

    λ_al = decay_constant(T_AL28)
    # the original set τ = 5·T½ and ran three activation/pause cycles
    τ, K = 5 * T_AL28, 2.565e7
    @printf("²⁸Al T½ = %.1f s (the original used %.1f s, %.1f %% high)\n",
        T_AL28, 2.3 * 60, 100 * (2.3 * 60 / T_AL28 - 1))
    @printf("pulsed flux, τ = %.0f s: saturation activity K = %.3e s⁻¹\n", τ, K)
    for n in (1, 3, 5, 10)
        @printf("  end of irradiation %2d: Λ = %.4e s⁻¹ (%.1f %% of saturation)\n",
            n, pulsed_activity((2n - 1) * τ, λ_al, τ, K),
            100 * pulsed_activity((2n - 1) * τ, λ_al, τ, K) / K)
    end

    fig = Figure(size = (1000, 450))

    ax1 = Axis(fig[2, 1], xlabel = L"Time $t$ [d]", ylabel = L"$\Lambda_2/\Lambda_0$",
        xscale = log10, xticks = logticks(-1, 3),)
    td = 10 .^ range(-1, 3, length = 500)
    lines!(ax1, td, bateman.(td .* 86400, λ₁, λ₂), color = PALETTE.blue, linewidth = 1.8)
    v = vlines!(ax1, [t_m / 86400], color = PALETTE.red, linestyle = :dash, linewidth = 1.3)
    # Right-aligned to the left of the line. The maximum falls at the right-hand
    # end of the decade range, so a label placed to its right was sliced in two
    # by the axis frame and read "t_n".
    text!(ax1, t_m / 86400 * 0.85, 0.35;
        text = rich(it("t"), subscript("m"), @sprintf(" = %.0f d", t_m / 86400)),
        color = PALETTE.red, align = (:right, :center), fontsize = 15,)

    # Activity in units of 1e7: every tick otherwise repeated the factor, and
    # Makie wrote the first of them as 1x10^7.
    scale = 1e7
    ax2 = Axis(fig[2, 2], xlabel = L"Time $t$ [s]",
        ylabel = L"Activity $\Lambda$ [$10^{7}$ s$^{-1}$]",)
    ts = range(0, 6τ, length = 3000)
    for c in 0:2
        vspan!(ax2, 2c * τ, (2c + 1) * τ, color = (PALETTE.orange, 0.14))
    end
    lines!(ax2, ts, pulsed_activity.(ts, λ_al, τ, K) ./ scale,
        color = PALETTE.green, linewidth = 1.5,)
    h = hlines!(ax2, [K / scale], color = PALETTE.black, linestyle = :dash, linewidth = 1.1)
    reached = pulsed_activity(τ, λ_al, τ, K) / K
    # headroom above the saturation line, the only band this panel leaves clear
    ylims!(ax2, -0.1, K / scale * 1.28)
    text!(ax2, 0.98, 0.97;
        text = @sprintf("%.1f %% of saturation by the end of the first irradiation",
            100 * reached),
        space = :relative, align = (:right, :top), fontsize = 14,
        color = PALETTE.green,)

    Legend(fig[1, 1:2],
        [v, h, PolyElement(color = (PALETTE.orange, 0.35))],
        [rich(superscript("234"), "Th maximum"), "Saturation activity",
            "Flux on",],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 26,)
    rowsize!(fig.layout, 2, Relative(0.85))
    println("\nwrote ", savefigure(fig, FIGURES, "decay_and_activation"))
end

main()
