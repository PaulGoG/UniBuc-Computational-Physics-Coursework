# Largest Lyapunov exponent of the Lorenz and Rössler systems by the Benettin
# algorithm; the un-renormalised divergence of two neighbouring trajectories,
# which is what the 2018 program followed; and the intervals between successive
# extrema of that separation as that program defined them.
#
# Benettin: the perturbation is rescaled to d₀ after every interval τ and
# λ₁ = ⟨ln(d/d₀)⟩/τ over the run. Successive increments are correlated, so the
# standard error is the largest blocked estimate of Flyvbjerg and Petersen
# (J. Chem. Phys. 91, 461 (1989), doi:10.1063/1.457480). Reference values from
# Sprott, Chaos and Time-Series Analysis, Oxford University Press (2003):
# 0.9056 for the Lorenz attractor and 0.0714 for the Rössler attractor. The
# measured values are asserted within three standard errors of them.
#
# Ported from Integratori.cpp, Atractori.cpp and helpers.cpp in
# Code_Archive/Old_2018/C_C++/AtractorI/ on the `legacy` branch.
# Integratori.cpp integrates the reference orbit and a copy displaced by 10⁻⁷ in
# each coordinate; helpers.cpp drops 30 time units at each end of the
# separation series, walks it accepting an extremum when it differs from the
# last accepted one by more than ε = 10⁻⁴ in relative terms, and writes the
# min→max intervals in one column and the max→min ones in the other. That
# definition is kept. The program reads its own output back through an
# ifstream opened on the file its ofstream is still writing, which the C++
# standard does not guarantee; compiled with local paths, N = 200 000 steps
# over 200 time units from (1, 1, 1), it does produce 140 001 distances, none
# zero, and 301 interval rows with mean 0.28. Its Rössler branch was
# unreachable: the selector was assigned and never read.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "attractors_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Displacement of the perturbed orbit in each coordinate, as in Integratori.cpp."
const D₀ = 1e-7
"Integration step."
const H = 0.002
"Steps between two Benettin renormalisations; τ = 0.1 time units."
const RENORM_STEPS = 50
"Transient discarded before every measurement, in steps: 40 time units."
const TRANSIENT = 20_000
"Time units dropped at each end of the separation series by helpers.cpp."
const SKIP = 30.0
"Relative change that qualifies an extremum in helpers.cpp."
const ε_EXTREMUM = 1e-4
"Number of standard errors within which λ₁ must agree with the reference."
const AGREEMENT = 3.0

"""
    separation_history(f, u₀, p, h, T, d₀; transient)

Reference orbit and a copy displaced by `d₀` in each coordinate, integrated
without renormalisation for `T` time units after `transient` steps. Returns the
time base and the Euclidean separation.
"""
function separation_history(f, u₀, p, h, T, d₀; transient = 0)
    u = Tuple(float.(u₀))
    for _ in 1:transient
        u = rk4_step(f, u, h, p)
    end
    v = u .+ d₀
    n = round(Int, T / h)
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
    benettin_increments(f, u₀, p, h, T, steps_between, d₀; transient)

The Benettin increments ln(d/d₀)/τ over `T` time units, one per
renormalisation; their mean is λ₁.
"""
function benettin_increments(f, u₀, p, h, T, steps_between, d₀; transient = 0)
    u = Tuple(float.(u₀))
    for _ in 1:transient
        u = rk4_step(f, u, h, p)
    end
    v = u .+ d₀
    τ = steps_between * h
    increments = Vector{Float64}(undef, round(Int, T / τ))
    for k in eachindex(increments)
        for _ in 1:steps_between
            u = rk4_step(f, u, h, p)
            v = rk4_step(f, v, h, p)
        end
        Δ = v .- u
        d = sqrt(sum(abs2, Δ))
        increments[k] = log(d / d₀) / τ
        v = u .+ Δ .* (d₀ / d)          # back to d₀ along the same direction
    end
    return increments
end

"""
    blocking_error(x)

Standard error of the mean of the correlated series `x` by Flyvbjerg–Petersen
blocking: successive pairs are averaged until fewer than eight blocks remain,
and the largest standard error along the way is taken.
"""
function blocking_error(x)
    y = copy(x)
    worst = 0.0
    while length(y) >= 8
        worst = max(worst, std(y) / sqrt(length(y)))
        y = [(y[2i - 1] + y[2i]) / 2 for i in 1:(length(y) ÷ 2)]
    end
    return worst
end

"""
    extrema_intervals(t, d; skip, ε)

The interval series of helpers.cpp: with `skip` time units dropped at each end,
an interior extremum of `d` is accepted when it differs from the last accepted
one by more than `ε` relatively, and the time since the previous accepted
extremum is recorded as a rise (min → max) or a fall (max → min).
"""
function extrema_intervals(t, d; skip, ε)
    keep = findall(τ -> skip <= τ <= last(t) - skip, t)
    tk, dk = t[keep], d[keep]
    rise, fall = Float64[], Float64[]
    last_value, last_time = dk[1], NaN
    for i in 2:(length(dk) - 1)
        is_max = dk[i] > dk[i - 1] && dk[i] > dk[i + 1]
        is_min = dk[i] < dk[i - 1] && dk[i] < dk[i + 1]
        (is_max || is_min) && abs(last_value - dk[i]) / last_value > ε || continue
        isnan(last_time) || push!(is_max ? rise : fall, tk[i] - last_time)
        last_time, last_value = tk[i], dk[i]
    end
    return rise, fall
end

# (field, parameters, name, colour, reference λ₁, horizon of the separation
# history, horizon of the Benettin run) — the Rössler separation saturates only
# after 200 time units, and its exponent needs a run five times longer for a
# standard error of 2 in the third digit
const CASES = (
    (lorenz, LORENZ, "Lorenz", PALETTE.blue, 0.9056, 80.0, 20_000.0),
    (rossler, ROSSLER, "Rössler", PALETTE.green, 0.0714, 300.0, 100_000.0),
)

function main()
    u₀ = (1.0, 1.0, 1.0)
    fig = Figure(size = (1600, 640))
    ax1 = Axis(
        fig[2, 1], xlabel = L"Time $t$", ylabel = L"Separation $|\Delta\mathbf{u}|$",
        yscale = log10, yticks = logticks(-9, 3; step = 3),)
    ax2 = Axis(fig[2, 2], xlabel = L"Run length $T$", ylabel = L"Running $\lambda_1$",
        xscale = log10, xticks = logticks(1, 5),)
    ax3 = Axis(fig[2, 3], xlabel = "Interval between extrema",
        ylabel = "Count",)

    handles = Any[]
    results = NamedTuple[]
    for (f, p, name, colour, reference, T_sep, T_ben) in CASES
        t, d = separation_history(f, u₀, p, H, T_sep, D₀; transient = TRANSIENT)
        push!(handles, lines!(ax1, t, d, color = colour))
        saturation = findfirst(>(1.0), d)
        @printf("%-8s separation from %.0e reaches 1 at t = %.0f, %.1f × %.0f at the end\n",
            name, D₀, saturation === nothing ? NaN : t[saturation], last(d), T_sep)

        increments = benettin_increments(f, u₀, p, H, T_ben, RENORM_STEPS, D₀;
            transient = TRANSIENT,)
        λ, σλ = mean(increments), blocking_error(increments)
        running = cumsum(increments) ./ (1:length(increments))
        τ = RENORM_STEPS * H
        lines!(ax2, τ .* (1:length(increments)), running, color = colour)
        hlines!(
            ax2, [reference], color = colour, linestyle = :dash, linewidth = GUIDE_WIDTH,)
        @printf("%-8s λ₁ = %.4f ± %.4f over T = %.0f (blocked standard error); reference %.4f, %.1f σ away\n",
            name, λ, σλ, T_ben, reference, abs(λ - reference) / σλ)
        abs(λ - reference) <= AGREEMENT * σλ ||
            error("$name: λ₁ = $λ ± $σλ is more than $AGREEMENT σ from $reference")
        push!(results, (; name, colour, λ, σλ, reference))
    end

    # the interval series of helpers.cpp, on its own run: 200 000 steps of
    # 0.001 over 200 time units from (1, 1, 1), no transient beyond its own skip
    t, d = separation_history(lorenz, u₀, LORENZ, 0.001, 200.0, D₀)
    rise, fall = extrema_intervals(t, d; skip = SKIP, ε = ε_EXTREMUM)
    @printf("Lorenz, helpers.cpp definition: %d min→max intervals, mean %.3f ± %.3f; %d max→min, mean %.3f ± %.3f\n",
        length(rise), mean(rise), std(rise) / sqrt(length(rise)),
        length(fall), mean(fall),
        std(fall) / sqrt(length(fall)))
    edges = range(0, 1.2, length = 31)
    hist!(ax3, rise, bins = edges, color = (PALETTE.blue, 0.55), strokewidth = 1.0,
        strokecolor = PALETTE.blue,)
    hist!(ax3, fall, bins = edges, color = (PALETTE.orange, 0.55), strokewidth = 1.0,
        strokecolor = PALETTE.orange,)
    text!(ax3, 0.97, 0.96;
        text = rich(
            rich(@sprintf("min → max: %d, mean %.3f", length(rise), mean(rise)),
                color = PALETTE.blue,),
            "\n",
            rich(@sprintf("max → min: %d, mean %.3f", length(fall), mean(fall)),
                color = PALETTE.orange,), "\nLorenz, 200 time units",),
        space = :relative, align = (:right, :top), justification = :right,
        fontsize = ANNOTATION_SIZE,)

    ylims!(ax1, 1e-9, 1e3)
    text!(ax1, 0.97, 0.04; text = "Saturation at the attractor diameter",
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)
    for (k, r) in enumerate(results)
        text!(ax2, 0.97, 0.96 - 0.09(k - 1);
            text = rich(it("λ"), subscript("1"), @sprintf(" = %.4f ± %.4f", r.λ, r.σλ)),
            space = :relative, align = (:right, :top), color = r.colour,
            fontsize = ANNOTATION_SIZE,)
    end
    text!(ax2, 0.97, 0.74; text = "Dashed: reference values",
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,)
    xlims!(ax2, 8, 1.5e5)
    ylims!(ax2, -0.1, 1.6)

    Legend(fig[1, 1:3], handles, [r.name for r in results])
    println("wrote ", savefigure(fig, FIGURES, "lyapunov_and_return_times"))
end

main()
