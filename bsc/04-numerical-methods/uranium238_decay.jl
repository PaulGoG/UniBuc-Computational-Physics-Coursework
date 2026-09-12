# Radioactive decay of a single isotope, dN/dt = -N/τ, integrated with forward
# Euler against the analytic N₀e^{-t/τ}.
#
# Ported from Dez_U_28.m (Octave). Two things to correct.
#
# First, it is **not** a decay chain, despite the name: one nuclide, no
# daughters. ²³⁸U in fact decays through a fourteen-member chain to ²⁰⁶Pb, but
# none of that is modelled and the single-isotope law is what the file solves.
#
# Second, and more consequential, the original wrote
#
#     tau = 4.4e9;   % timpul mediu de viata U238
#
# 4.468 × 10⁹ yr is the **half-life** of ²³⁸U, not its mean lifetime. The mean
# life is T½/ln 2 = 6.446 × 10⁹ yr, so the sample was decayed about 44 per cent
# too fast. Both are shown below.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Half-life of ²³⁸U in years."
const T_HALF = 4.468e9
"Mean lifetime, τ = T½/ln 2."
const τ = T_HALF / log(2)

"""
    euler_decay(N₀, τ, Δt, n)

Forward Euler on dN/dt = -N/τ.
"""
function euler_decay(N₀, τ, Δt, n)
    t = Vector{Float64}(undef, n + 1)
    N = Vector{Float64}(undef, n + 1)
    t[1], N[1] = 0.0, N₀
    for i in 1:n
        N[i+1] = N[i] - N[i] / τ * Δt
        t[i+1] = i * Δt
    end
    return t, N
end

function main()
    N₀ = 1.0
    # the original's grid: dt = 1e7 yr over 1000 points
    Δt = 1e7
    n = 1000

    t, N_correct = euler_decay(N₀, τ, Δt, n)
    _, N_2018 = euler_decay(N₀, T_HALF, Δt, n)   # half-life used as mean life
    analytic = N₀ .* exp.(-t ./ τ)

    @printf("T½ = %.3e yr, mean life τ = T½/ln2 = %.3e yr\n", T_HALF, τ)
    @printf("Euler step Δt/τ = %.2e, max relative error against analytic = %.2e\n",
            Δt / τ, maximum(abs.(N_correct .- analytic) ./ analytic))
    @printf("at t = τ: correct N/N₀ = %.4f, using T½ as τ = %.4f (%.1f%% low)\n",
            N_correct[findfirst(>=(τ), t)], N_2018[findfirst(>=(τ), t)],
            100 * (1 - N_2018[findfirst(>=(τ), t)] / N_correct[findfirst(>=(τ), t)]))
    @printf("half-life recovered from the correct curve: %.3e yr\n",
            t[findfirst(<=(0.5), N_correct)])

    fig = Figure(size = (820, 480))
    ax = Axis(fig[2, 1], xlabel = L"Time $t$ [Gyr]", ylabel = L"$N/N_0$")

    l_a = lines!(ax, t ./ 1e9, analytic, color = PALETTE.black, linewidth = 1.6)
    l_c = scatter!(ax, t[1:40:end] ./ 1e9, N_correct[1:40:end],
        color = PALETTE.blue, markersize = 8)
    l_w = lines!(ax, t ./ 1e9, N_2018, color = PALETTE.red,
        linewidth = 1.5, linestyle = :dash)

    l_h = hlines!(ax, [0.5], color = PALETTE.green, linestyle = :dot, linewidth = 1.2)
    vlines!(ax, [T_HALF / 1e9], color = PALETTE.green, linestyle = :dot, linewidth = 1.2)
    text!(ax, T_HALF / 1e9 + 0.3, 0.55;
        text = L"$T_{1/2} = 4.468$ Gyr", color = PALETTE.green,
        align = (:left, :bottom), fontsize = 16)

    Legend(fig[1, 1], [l_a, l_c, l_w],
        [L"Analytic $N_0 e^{-t/\tau}$", "Forward Euler", L"2018: $T_{1/2}$ used as $\tau$"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24)

    rowsize!(fig.layout, 2, Relative(0.86))
    path = savefigure(fig, FIGURES, "uranium238_decay")
    println("wrote ", path)
end

main()
