# Bayesian updating of the bias of a coin. With a Beta(α, β) prior and a
# Bernoulli likelihood the posterior after h heads in N tosses is conjugate,
#
#   p(θ | data) = Beta(α + h, β + N - h),
#
# so the whole sequence is available in closed form and the posterior
# concentrates as 1/√N about the true bias.
#
# Ported from Statis.jl. The statistics were right; the presentation was the
# problem. The original produced a 1001-frame animation at 1280×900 that came to
# 14.6 MB — larger than every other file in the repository put together — and
# carried a plot title, which this project's figures do not use. The animation
# here is decimated and scaled, and a static figure carries the actual result.
#
# The credible interval and the 1/√N concentration are new: the original showed
# the posterior widening and narrowing but never quantified it.

using Printf, Distributions, StableRNGs, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

const N_TOSSES = 1000
const P_TRUE = 0.5
const SEED = 12

"Posterior after `h` heads in `n` tosses, from a Beta(α, β) prior."
posterior(h, n; α = 1.0, β = 1.0) = Beta(α + h, β + n - h)

function main()
    rng = StableRNG(SEED)
    tosses = rand(rng, Bernoulli(P_TRUE), N_TOSSES)
    heads = cumsum(tosses)

    checkpoints = (1, 5, 20, 100, 1000)
    for n in checkpoints
        d = posterior(heads[n], n)
        lo, hi = quantile(d, 0.025), quantile(d, 0.975)
        @printf("N = %4d  h = %4d  mean %.4f  95%% CI [%.4f, %.4f]  width %.4f\n",
                n, heads[n], mean(d), lo, hi, hi - lo)
    end

    # the width should fall as 1/sqrt(N)
    ns = [10, 30, 100, 300, 1000]
    widths = map(ns) do n
        d = posterior(heads[n], n)
        quantile(d, 0.975) - quantile(d, 0.025)
    end
    slope = (log(widths[end]) - log(widths[1])) / (log(ns[end]) - log(ns[1]))
    @printf("credible-interval width scales as N^%.3f (expected -0.5)\n", slope)

    fig = Figure(size = (980, 440))

    ax1 = Axis(fig[2, 1], xlabel = L"Bias $\theta$", ylabel = "Posterior density")
    θ = range(0, 1, length = 600)
    handles = []
    cols = (PALETTE.sky, PALETTE.green, PALETTE.orange, PALETTE.purple, PALETTE.blue)
    for (k, n) in enumerate(checkpoints)
        d = posterior(heads[n], n)
        push!(handles, lines!(ax1, θ, pdf.(d, θ), color = cols[k], linewidth = 1.6))
    end
    vlines!(ax1, [P_TRUE], color = PALETTE.red, linestyle = :dash, linewidth = 1.2)
    xlims!(ax1, 0.2, 0.8)

    ax2 = Axis(fig[2, 2], xlabel = L"Tosses $N$", ylabel = "95 % credible width",
        xscale = log10, yscale = log10,
        xticks = ([10, 100, 1000], ["10", "100", "1000"]))
    scatterlines!(ax2, ns, widths, color = PALETTE.blue, markersize = 9)
    lines!(ax2, ns, widths[1] .* (ns ./ ns[1]) .^ (-0.5),
        color = PALETTE.black, linestyle = :dash, linewidth = 1.2)
    text!(ax2, 0.96, 0.92; text = L"$\propto N^{-1/2}$",
        space = :relative, align = (:right, :top), fontsize = 16)

    Legend(fig[1, 1:2], handles, [L"N = %$n" for n in checkpoints],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 22)
    rowsize!(fig.layout, 2, Relative(0.85))
    println("wrote ", savefigure(fig, FIGURES, "bayesian_coin_updating"))

    # animation: every fifth toss, 640x420, which keeps it near 1 MB rather
    # than the 14.6 MB of the 2018 version
    mkpath(FIGURES)
    gif_path = joinpath(FIGURES, "bayesian_coin_updating.gif")
    frames = 1:5:N_TOSSES
    fa = Figure(size = (640, 420))
    axa = Axis(fa[1, 1], xlabel = L"Bias $\theta$", ylabel = "Posterior density")
    obs = Observable(pdf.(posterior(heads[1], 1), θ))
    lines!(axa, θ, obs, color = PALETTE.green, linewidth = 2)
    vlines!(axa, [P_TRUE], color = PALETTE.red, linestyle = :dash, linewidth = 1.2)
    label = Observable("N = 1")
    text!(axa, 0.97, 0.93; text = label, space = :relative,
          align = (:right, :top), fontsize = 18)
    xlims!(axa, 0, 1); ylims!(axa, 0, 30)
    record(fa, gif_path, frames; framerate = 25) do n
        obs[] = pdf.(posterior(heads[n], n), θ)
        label[] = "N = $n"
    end
    @printf("wrote %s (%.2f MB)\n", gif_path, filesize(gif_path) / 1024^2)
end

main()
