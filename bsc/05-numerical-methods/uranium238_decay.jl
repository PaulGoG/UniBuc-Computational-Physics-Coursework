# Radioactive decay of a single nuclide, dN/dt = −N/τ, by forward Euler against
# the analytic N₀e^{−t/τ}, on the grid of Dez_U_28.m on the `legacy` branch:
# Δt = 10⁷ yr, 1000 points.
#
# The original sets
#
#     tau = 4.4e9;   % timpul mediu de viata U238
#
# and uses it as the mean lifetime in both the Euler loop and the analytic curve.
# 4.4 × 10⁹ yr is within 1.5 % of the half-life of ²³⁸U, (4.4683 ± 0.0024) × 10⁹ yr
# (Jaffey et al., Phys. Rev. C 4, 1889 (1971), doi:10.1103/PhysRevC.4.1889), and
# 32 % below the mean lifetime τ = T½/ln 2 = 6.446 × 10⁹ yr, so I take it to be
# the half-life entered where the mean lifetime belongs. The two Euler curves
# below differ only in that number.
#
# The file is named for a decay but models no chain: ²³⁸U reaches ²⁰⁶Pb through
# fourteen intermediate nuclides, none of which appears; one nuclide, one rate.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Half-life of ²³⁸U in years, Jaffey et al. (1971), doi:10.1103/PhysRevC.4.1889."
const T_HALF = 4.4683e9
"Mean lifetime in years, τ = T½/ln 2."
const τ = T_HALF / log(2)
"The value assigned to `tau` in Dez_U_28.m, in years."
const τ_LEGACY = 4.4e9
"Time step of the original, in years."
const Δt_LEGACY = 1e7
"Number of grid points of the original."
const N_POINTS = 1000
"Years per gigayear, for the axis."
const GYR = 1e9

"""
    euler_decay(N₀, τ, Δt, n)

`n` forward-Euler steps of dN/dt = −N/τ from `N₀`. Returns times and populations.
"""
function euler_decay(N₀, τ, Δt, n)
    t = Vector{Float64}(undef, n + 1)
    N = Vector{Float64}(undef, n + 1)
    t[1], N[1] = 0.0, N₀
    for i in 1:n
        N[i + 1] = N[i] - N[i] / τ * Δt
        t[i + 1] = i * Δt
    end
    return t, N
end

"""
    crossing_time(t, N, level)

Time at which the decreasing series `N` crosses `level`, by interpolation of
log N between the two bracketing grid points, which is exact for a geometric
sequence such as the Euler solution.
"""
function crossing_time(t, N, level)
    i = findfirst(<=(level), N)
    i === nothing && throw(ArgumentError("the series never falls to $level"))
    i == 1 && return t[1]
    w = (log(N[i - 1]) - log(level)) / (log(N[i - 1]) - log(N[i]))
    return t[i - 1] + w * (t[i] - t[i - 1])
end

function main()
    N₀ = 1.0
    n = N_POINTS - 1
    t, N_euler = euler_decay(N₀, τ, Δt_LEGACY, n)
    _, N_legacy = euler_decay(N₀, τ_LEGACY, Δt_LEGACY, n)
    analytic = N₀ .* exp.(-t ./ τ)

    x = Δt_LEGACY / τ
    max_relative = maximum(abs.(N_euler ./ analytic .- 1))
    @printf("T½ = %.4e yr, τ = T½/ln 2 = %.4e yr, legacy τ = %.1e yr\n",
        T_HALF, τ, τ_LEGACY)
    @printf("Δt/τ = %.3e; max relative Euler error = %.3e, leading term n(Δt/τ)²/2 = %.3e\n",
        x, max_relative, n * x^2 / 2)
    isapprox(max_relative, n * x^2 / 2; rtol = 0.01) ||
        error("Euler error $(max_relative) is not the first-order n(Δt/τ)²/2")

    # The Euler solution is N₀(1 − Δt/τ)ⁿ, which halves after
    # Δt ln 2 / (−ln(1 − Δt/τ)) = T½ [1 − Δt/(2τ) + …].
    half_euler = crossing_time(t, N_euler, N₀ / 2)
    half_closed = Δt_LEGACY * log(2) / -log1p(-x)
    half_legacy = crossing_time(t, N_legacy, N₀ / 2)
    @printf("half-life of the Euler solution: %.4e yr (closed form %.4e), %.3f %% below T½\n",
        half_euler, half_closed, 100 * (1 - half_euler / T_HALF))
    @printf("half-life of the Euler solution with the legacy τ: %.3e yr\n", half_legacy)
    isapprox(half_euler, half_closed; rtol = 1e-9) ||
        error("interpolated half-life $(half_euler) is off the closed form $(half_closed)")

    at_τ = exp(-1.0)
    legacy_at_τ = exp(-τ / τ_LEGACY)
    @printf("surviving fraction at t = τ: %.4f; with the legacy τ: %.4f (%.0f %% lower)\n",
        at_τ, legacy_at_τ, 100 * (1 - legacy_at_τ / at_τ))
    @printf("decay rate with the legacy τ: %.1f %% too high\n", 100 * (τ / τ_LEGACY - 1))

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (900, 640))
    ax = Axis(
        fig[2, 1], xlabel = L"Time $t$ [Gyr]", ylabel = L"Surviving fraction $N/N_0$",
        xticks = 0:2:10, yticks = 0:0.25:1,)

    hlines!(ax, [0.5], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    vlines!(ax, [T_HALF / GYR], color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    text!(ax, T_HALF / GYR, 1.0;
        text = rich(it("T"), subscript("1/2"), @sprintf(" = %.3f Gyr", T_HALF / GYR)),
        align = (:left, :top), offset = (8, 0), fontsize = ANNOTATION_SIZE,)

    l_analytic = lines!(ax, t ./ GYR, analytic, color = PALETTE.black)
    l_legacy = lines!(ax, t ./ GYR, N_legacy, color = PALETTE.red, linestyle = :dash)
    l_euler = scatter!(ax, t[1:40:end] ./ GYR, N_euler[1:40:end], color = PALETTE.blue)
    xlims!(ax, -0.2, 10.2)
    ylims!(ax, 0, 1.04)

    text!(ax, 0.97, 0.62;
        text = rich(
            rich("Half-life ", @sprintf("%.3f Gyr", half_euler / GYR), color = PALETTE.blue),
            "\n",
            rich("Half-life ", @sprintf("%.2f Gyr", half_legacy / GYR),
                color = PALETTE.red,),),
        space = :relative, align = (:right, :bottom), justification = :right,
        fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1], [l_analytic, l_euler, l_legacy],
        # all three as LaTeX strings, so that τ is one glyph across the legend
        [L"N_0 e^{-t/\tau}",
            latexstring(@sprintf("\\text{Euler, }\\tau = %.3f\\ \\text{Gyr}", τ / GYR)),
            latexstring(@sprintf("\\text{Euler, original }\\tau = %.1f\\ \\text{Gyr}",
                τ_LEGACY / GYR)),],)

    path = savefigure(fig, FIGURES, "uranium238_decay")
    println("wrote ", path)
end

main()
