# Bayesian updating of the bias of a coin. With a Beta(α, β) prior and a
# Bernoulli likelihood the posterior after h heads in N tosses is conjugate,
#
#   p(θ | data) = Beta(α + h, β + N − h),
#
# so the whole sequence is available in closed form. Ported from Statis.jl in
# Julia-Workflow-FFUB/Statistics_symul_L_4/ on the `legacy` branch, whose
# statistics were right and whose animation was 1001 frames at 1280 × 900,
# 14.6 MB.
#
# The width of the central 95 % credible interval is that of the normal
# approximation,
#
#   w_N ≈ 2 z₀.₉₇₅ σ_N,   σ_N² = θ̂(1 − θ̂)/(N + 3),   θ̂ = (h + 1)/(N + 2),
#
# the variance of Beta(a, b) being ab/((a + b)²(a + b + 1)) with a + b = N + 2
# for the flat prior. `main` asserts the approximation at N = 10 and N = 1000.
# The interval therefore narrows as N^{−1/2} only asymptotically: over any
# finite range of N the fitted exponent departs from −1/2 through the +3 of the
# prior and through the wandering of θ̂ along the sequence, and `main`
# separates the two by refitting with θ̂ held at the true bias.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Distributions, StableRNGs, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Number of tosses, true bias and seed of the sequence."
const N_TOSSES = 1000
const P_TRUE = 0.5
const SEED = 12
"Tosses at which the posterior is drawn."
const CHECKPOINTS = (1, 5, 20, 100, 1000)
"Range of N over which the exponent of the credible width is fitted."
const FIT_RANGE = 10:N_TOSSES
"Relative tolerance of the normal approximation to the credible width at N = 10 and at N = 1000."
const WIDTH_TOLERANCE = (0.03, 0.01)

"Posterior after `h` heads in `n` tosses, from a Beta(α, β) prior."
posterior(h, n; α = 1.0, β = 1.0) = Beta(α + h, β + n - h)

"Width of the central 95 % credible interval of `d`."
credible_width(d) = quantile(d, 0.975) - quantile(d, 0.025)

"Normal approximation to the 95 % credible width of `d`, 2 z₀.₉₇₅ σ."
normal_width(d) = 2 * quantile(Normal(), 0.975) * std(d)

"Least-squares slope of `y` against `x`."
function slope(x, y)
    x̄, ȳ = mean(x), mean(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
end

function main()
    rng = StableRNG(SEED)
    tosses = rand(rng, Bernoulli(P_TRUE), N_TOSSES)
    heads = cumsum(tosses)

    for n in CHECKPOINTS
        d = posterior(heads[n], n)
        @printf("N = %4d  h = %4d  mean %.4f  95 %% interval [%.4f, %.4f]  width %.4f  normal approximation %.4f\n",
            n, heads[n], mean(d), quantile(d, 0.025),
            quantile(d, 0.975), credible_width(d),
            normal_width(d))
    end
    for (n, tol) in zip((first(FIT_RANGE), N_TOSSES), WIDTH_TOLERANCE)
        d = posterior(heads[n], n)
        isapprox(credible_width(d), normal_width(d); rtol = tol) ||
            error("normal approximation to the credible width off by more than $tol at N = $n")
    end

    ns = collect(FIT_RANGE)
    posteriors = [posterior(heads[n], n) for n in ns]
    widths = credible_width.(posteriors)
    exponent = slope(log.(ns), log.(widths))
    exponent_normal = slope(log.(ns), log.(normal_width.(posteriors)))
    exponent_fixed = slope(log.(ns),
        log.(2 * quantile(Normal(), 0.975) .* sqrt.(P_TRUE * (1 - P_TRUE) ./ (ns .+ 3))),)
    @printf("credible width over %d ≤ N ≤ %d: fitted exponent %.3f, %.3f from the normal approximation, %.3f with θ̂ held at %.1f, −1/2 asymptotically\n",
        first(ns), last(ns), exponent, exponent_normal, exponent_fixed, P_TRUE)
    abs(exponent - exponent_normal) < 0.01 ||
        error("the exponents of the exact and approximate widths differ by more than 0.01")

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1400, 620))
    ax1 = Axis(fig[2, 1], xlabel = L"Bias $\theta$", ylabel = "Posterior density")
    θ = range(0, 1, length = 600)
    cols = (PALETTE.sky, PALETTE.green, PALETTE.orange, PALETTE.purple, PALETTE.blue)
    handles = [lines!(ax1, θ, pdf.(posterior(heads[n], n), θ), color = cols[k])
               for (k, n) in enumerate(CHECKPOINTS)]
    vlines!(
        ax1, [P_TRUE], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    limits!(ax1, 0, 1, 0, 30)
    text!(ax1, P_TRUE + 0.01, 29.5; text = rich("True bias ", it("θ"), " = 0.5"),
        align = (:left, :top), fontsize = ANNOTATION_SIZE,)

    ax2 = Axis(fig[2, 2], xlabel = L"Tosses $N$", ylabel = "95 % credible width",
        xscale = log10, yscale = log10, xticks = logticks(1, 3),
        yticks = ([0.05, 0.1, 0.2, 0.5], [L"0.05", L"0.1", L"0.2", L"0.5"]),)
    l_w = lines!(ax2, ns, widths, color = PALETTE.red)
    l_g = lines!(ax2, ns, widths[end] .* (ns ./ ns[end]) .^ (-0.5),
        color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    # Lower left, where the curve has not yet come down; kept short enough to
    # end before it does.
    text!(ax2, 0.03, 0.05;
        text = rich(
            rich(
                @sprintf("Fitted exponent −%.3f over %d ≤ ", -exponent, first(ns)), it("N"),
                @sprintf(" ≤ %d", last(ns)), color = PALETTE.red,),
            "\n",
            @sprintf("−%.3f with the mean held at ½; −½ asymptotically", -exponent_fixed)),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1:2], [handles..., l_w, l_g],
        [[L"N = %$n" for n in CHECKPOINTS]..., "Credible width", L"N^{-1/2}"],)
    println("wrote ", savefigure(fig, FIGURES, "bayesian_coin_updating"))

    # --- animation: every tenth toss ------------------------------------------
    mkpath(FIGURES)
    gif_path = joinpath(FIGURES, "bayesian_coin_updating.gif")
    frames = 1:10:N_TOSSES
    fa = Figure(size = (1200, 560))
    axa = Axis(fa[1, 1], xlabel = L"Bias $\theta$", ylabel = "Posterior density")
    density = Observable(pdf.(posterior(heads[1], 1), θ))
    lines!(axa, θ, density, color = PALETTE.blue)
    vlines!(
        axa, [P_TRUE], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    caption = Observable{Any}(rich(it("N"), " = 1, ", it("h"), " = $(heads[1])"))
    text!(axa, 0.97, 0.95; text = caption, space = :relative, align = (:right, :top),
        fontsize = ANNOTATION_SIZE,)
    limits!(axa, 0, 1, 0, 30)
    record(fa, gif_path, frames; framerate = 12, px_per_unit = 1) do n
        density[] = pdf.(posterior(heads[n], n), θ)
        caption[] = rich(it("N"), " = $n, ", it("h"), " = $(heads[n])")
    end
    @printf("wrote %s (%.2f MB)\n", gif_path, filesize(gif_path) / 1024^2)
end

main()
