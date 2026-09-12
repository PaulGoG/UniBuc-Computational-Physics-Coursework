# Largest Lyapunov exponent of the Lorenz and Rössler attractors, measured from
# the divergence of two neighbouring trajectories, together with the
# distribution of intervals between successive maxima of that separation.
#
# This is the measurement the 2018 project was built around and never obtained.
# `helpers.cpp` read the two trajectory files through an ifstream constructed on
# the same paths as unflushed ofstreams, so every extraction failed and, since
# C++11, set its target to zero: `Distanta.txt` and `TimpiUC.txt` contained only
# zeros, and the Octave scripts that histogrammed them were histogramming
# nothing. It also used absolute Windows paths, compared doubles with `!=` and
# `== 0`, discarded a transient with an off-by-one, and mutated its `int &N`
# argument across a call chain.
#
# Two estimates of the separation are produced here:
#
#   * un-renormalised, as the original attempted, which grows exponentially only
#     until it saturates at the diameter of the attractor;
#   * Benettin renormalisation, in which the perturbation is rescaled back to d₀
#     at fixed intervals and the accumulated logarithms averaged. This is what
#     actually converges to λ₁.

using Printf, StatsBase
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "attractors_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Initial separation of the perturbed trajectory."
const D₀ = 1e-9

"""
    separation_history(f, u₀, p, h, n, d₀; transient)

Integrate a reference orbit and a perturbed copy without renormalising, and
return the time base together with the Euclidean separation.
"""
function separation_history(f, u₀, p, h, n, d₀; transient = 0)
    u = Tuple(float.(u₀))
    for _ in 1:transient
        u = rk4_step(f, u, h, p)
    end
    v = u .+ (d₀ / sqrt(3),) .* (1.0, 1.0, 1.0)
    t = Vector{Float64}(undef, n)
    d = Vector{Float64}(undef, n)
    for i in 1:n
        u = rk4_step(f, u, h, p)
        v = rk4_step(f, v, h, p)
        t[i] = i * h
        d[i] = sqrt(sum(abs2, u .- v))
    end
    return t, d
end

"""
    benettin_lyapunov(f, u₀, p, h, n_renorm, steps_between, d₀; transient)

Largest Lyapunov exponent by the Benettin algorithm: after every
`steps_between` integration steps the perturbation is rescaled back to `d₀` and
the logarithm of its growth accumulated. Returns the running estimate.
"""
function benettin_lyapunov(f, u₀, p, h, n_renorm, steps_between, d₀; transient = 0)
    u = Tuple(float.(u₀))
    for _ in 1:transient
        u = rk4_step(f, u, h, p)
    end
    v = u .+ (d₀ / sqrt(3),) .* (1.0, 1.0, 1.0)
    τ = steps_between * h
    running = Vector{Float64}(undef, n_renorm)
    accumulated = 0.0
    for k in 1:n_renorm
        for _ in 1:steps_between
            u = rk4_step(f, u, h, p)
            v = rk4_step(f, v, h, p)
        end
        Δ = v .- u
        d = sqrt(sum(abs2, Δ))
        accumulated += log(d / d₀)
        running[k] = accumulated / (k * τ)
        v = u .+ Δ .* (d₀ / d)          # rescale back to d₀ along the same direction
    end
    return running
end

"""
    maxima_intervals(t, d)

Intervals between successive interior local maxima of `d`. This is the
"return time" series the 2018 Octave scripts histogrammed.
"""
function maxima_intervals(t, d)
    peaks = Int[]
    for i in 2:(length(d) - 1)
        d[i] > d[i-1] && d[i] > d[i+1] && push!(peaks, i)
    end
    return length(peaks) > 1 ? diff(t[peaks]) : Float64[]
end

const CASES = (
    (lorenz, LORENZ, (1.0, 1.0, 1.0), "Lorenz", PALETTE.blue, 0.9056),
    (rossler, ROSSLER, (1.0, 1.0, 1.0), "Rössler", PALETTE.green, 0.0714),
)

function main()
    h = 0.002
    transient = 20_000

    fig = Figure(size = (1020, 480))
    ax1 = Axis(fig[2, 1], xlabel = L"Time $t$",
        ylabel = L"Separation $|\Delta\mathbf{u}|$", yscale = log10,
        yticks = ([1e-9, 1e-6, 1e-3, 1.0, 1e3],
                  [L"10^{-9}", L"10^{-6}", L"10^{-3}", L"1", L"10^{3}"]))
    ax2 = Axis(fig[2, 2], xlabel = "Interval between separation maxima",
        ylabel = "Count")

    handles = []
    for (f, p, u₀, name, colour, reference) in CASES
        t, d = separation_history(f, u₀, p, h, 40_000, D₀; transient = transient)
        push!(handles, lines!(ax1, t, d, color = colour, linewidth = 1.2))

        running = benettin_lyapunov(f, u₀, p, h, 4000, 50, D₀; transient = transient)
        λ = last(running)
        @printf("%-8s  lambda_1 = %.4f   literature %.4f   ratio %.3f\n",
                name, λ, reference, λ / reference)

        intervals = maxima_intervals(t, d)
        @printf("%-8s  %d separation maxima, mean interval %.4f\n",
                name, length(intervals) + 1, sum(intervals) / max(length(intervals), 1))
        hist!(ax2, intervals, bins = range(0, 4.0, length = 45),
              color = (colour, 0.55), strokewidth = 0.5, strokecolor = colour)

        # growth predicted by the measured exponent, anchored at the start
        lines!(ax1, t, D₀ .* exp.(λ .* t), color = colour,
               linestyle = :dash, linewidth = 1.1)
    end

    ylims!(ax1, 1e-10, 1e3)
    text!(ax1, 0.03, 0.95;
        text = L"Dashed: $d_0 e^{\lambda_1 t}$ from the Benettin estimate",
        space = :relative, align = (:left, :top), fontsize = 15)
    text!(ax1, 0.97, 0.80;
        text = "Saturation at the\nattractor diameter",
        space = :relative, align = (:right, :top), color = PALETTE.blue, fontsize = 15)

    Legend(fig[1, 1:2], handles, [c[4] for c in CASES],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)

    rowsize!(fig.layout, 2, Relative(0.86))
    path = savefigure(fig, FIGURES, "lyapunov_and_return_times")
    println("wrote ", path)
end

main()
