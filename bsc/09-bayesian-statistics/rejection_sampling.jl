# Von Neumann rejection sampling from an asymmetric bi-Gaussian density
#
#   f(x) = exp(-2x²)   for x ≥ 0        (σ = 1/2)
#          exp(-x²/2)  for x <  0        (σ = 1)
#
# Ported from MonteCarloDistribution.jl in Julia-Workflow-FFUB/Single_Files/ on
# the `legacy` branch. Four corrections.
#
#   1. The original defined `function Distribution(x)` after `using
#      Distributions`, which exports the abstract type `Distribution`. On Julia
#      ≤ 1.11 that was a hard error; on 1.12+ it is a deprecation warning and
#      silently *extends* `Distributions.Distribution` instead of defining a new
#      function in Main. Renamed here.
#   2. `contor` was initialised to 1 rather than 0 and incremented before the
#      acceptance test, so `contor/Np` reported the reciprocal of the efficiency
#      — and as a bare top-level expression it printed nothing when run as a
#      script.
#   3. The proposal support [-3, 3] truncates the tails without renormalising,
#      which discards 0.18 % of the mass, nearly all of it on the σ = 1 side, so
#      the sampled density is the truncated one. Quantified below against the
#      analytic normalisation.
#
# `main` asserts the acceptance rate against Z_trunc/((b − a) M) and the sample
# mean against quadrature, each within four standard errors.
#   4. No seed, so no run was reproducible; and `rand(Uniform(-3,3))` rebuilt the
#      distribution object on each of ~1.8 million iterations.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, StableRNGs, Statistics, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

const A, B = -3.0, 3.0
const SEED = 20210416

"Unnormalised target density. Maximum is f(0) = 1."
target(x) = x >= 0 ? exp(-2x^2) : exp(-x^2 / 2)

"Standard errors allowed between a Monte Carlo estimate and its analytic value."
const SIGMA_MAX = 4.0

"""
    sample(rng, n, a, b)

Draw `n` accepted samples by rejection against a uniform proposal on `[a, b]`
with envelope M = max f = 1. Returns the samples and the number of trials.
"""
function sample(rng, n, a, b)
    out = Vector{Float64}(undef, n)
    accepted = 0
    trials = 0
    while accepted < n
        trials += 1
        x = a + (b - a) * rand(rng)
        if rand(rng) < target(x)
            accepted += 1
            out[accepted] = x
        end
    end
    return out, trials
end

function main()
    rng = StableRNG(SEED)
    n = 1_000_000
    x, trials = sample(rng, n, A, B)
    rate = n / trials

    # analytic normalisation, full line and truncated
    Z_full, _ = quadgk(target, -Inf, Inf)
    Z_trunc, _ = quadgk(target, A, B)
    p_accept = Z_trunc / (B - A)
    σ_rate = sqrt(p_accept * (1 - p_accept) / trials)
    @printf("acceptance rate    = %.4f  (expected %.4f = Z_trunc / (b-a)·M, %.1f standard errors away)\n",
        rate, p_accept, abs(rate - p_accept) / σ_rate)
    abs(rate - p_accept) < SIGMA_MAX * σ_rate ||
        error("acceptance rate $rate is not $p_accept")
    @printf("normalisation      full line %.6f, truncated to [%.0f, %.0f] %.6f\n",
        Z_full, A, B, Z_trunc)
    @printf("mass discarded by truncation = %.4f %%\n", 100 * (1 - Z_trunc / Z_full))
    μ = quadgk(u -> u * target(u), A, B)[1] / Z_trunc
    σ_mean = std(x) / sqrt(n)
    @printf("sample mean %.5f vs analytic %.5f, %.1f standard errors away\n",
        mean(x), μ, abs(mean(x) - μ) / σ_mean)
    abs(mean(x) - μ) < SIGMA_MAX * σ_mean || error("sample mean $(mean(x)) is not $μ")

    fig = Figure(size = (1000, 640))
    ax = Axis(fig[2, 1], xlabel = L"x", ylabel = "Density")
    h_s = hist!(ax, x, bins = range(A, B, length = 120), normalization = :pdf,
        color = (PALETTE.blue, 0.6), strokewidth = 0.5, strokecolor = PALETTE.blue,)
    xf = range(A, B, length = 500)
    l_t = lines!(ax, xf, target.(xf) ./ Z_trunc, color = PALETTE.red)
    vlines!(ax, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    ylims!(ax, -0.015, 0.68)
    text!(ax, 0.02, 0.98;
        text = rich("Acceptance ", @sprintf("%.2f %%", 100rate), " measured\n",
            @sprintf("%.2f %%", 100p_accept), " analytic, ", rsci(float(n); digits = 0), " samples",),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)
    text!(ax, 0.98, 0.98;
        text = rich(
            it("σ"), " = 1 below 0, ", it("σ"), " = ½ above\ntruncated to [−3, 3]: ",
            @sprintf("%.2f %%", 100 * (1 - Z_trunc / Z_full)), " of the mass lost",),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        justification = :right,)

    Legend(fig[1, 1], [h_s, l_t], ["Rejection samples", "Target, truncated and normalised"])

    path = savefigure(fig, FIGURES, "rejection_sampling")
    println("wrote ", path)
end

main()
