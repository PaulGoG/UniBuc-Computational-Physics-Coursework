# Separation of two neighbouring trajectories of the linear oscillator
#
#   ẍ + 2δẋ + ω₀²x = F cos(ωt)
#
# integrated with RK4 from initial conditions a distance d₀ apart in phase space.
# The system is linear, so the difference of two solutions obeys the homogeneous
# equation whatever the driving does. Writing ω₁² = ω₀² - δ², the quadratic form
#
#   Q = (Δv + δΔx)² + ω₁²Δx²
#
# decays as e^{-2δt} exactly, for under-damped, critically damped (ω₁ = 0) and
# undamped motion alike. The Euclidean separation |Δ(x, v)| that the 2018
# programs wrote out oscillates inside that envelope; √Q does not, so the decay
# rate is asserted on √Q and only reported for |Δ|. A straight-line fit to
# ln|Δ| is biased by a few parts in 10⁴ according to where in the oscillation
# the window ends, which its formal standard error does not cover.
#
# Ported from OscilatorArmonique.cpp, OscilatorAmortizat.cpp and
# OscilatorFortat.cpp. Their three parameter sets are integrated as written,
# next to two under-damped ones. The originals read both initial conditions, the
# step count and the horizon from standard input; the values below are mine.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Initial phase-space separation of the two trajectories, in the units of x and v."
const D₀ = 1e-7

"Separations below this multiple of machine epsilon are left out of the rate fits."
const ROUNDOFF_MARGIN = 1e3

"Tolerance on the decay rate of √Q against the input δ, in inverse time units."
const RATE_TOLERANCE = 1e-5

"Tolerance on the relative agreement of the two δ = 0 separations; round-off on a 1e-7 difference."
const COINCIDENCE_TOLERANCE = 1e-6

"""
    rhs(t, x, v, ω₀, δ, F, ω)

Right-hand side of the driven damped oscillator as a first-order system.
"""
rhs(t, x, v, ω₀, δ, F, ω) = (v, F * cos(ω * t) - 2δ * v - ω₀^2 * x)

"""
    rk4_step(t, x, v, Δt, p)

One classical RK4 step. `p = (ω₀, δ, F, ω)`; the driving is evaluated at the
stage times.
"""
function rk4_step(t, x, v, Δt, p)
    k1x, k1v = rhs(t, x, v, p...)
    k2x, k2v = rhs(t + Δt / 2, x + Δt * k1x / 2, v + Δt * k1v / 2, p...)
    k3x, k3v = rhs(t + Δt / 2, x + Δt * k2x / 2, v + Δt * k2v / 2, p...)
    k4x, k4v = rhs(t + Δt, x + Δt * k3x, v + Δt * k3v, p...)
    return (x + Δt / 6 * (k1x + 2k2x + 2k3x + k4x),
        v + Δt / 6 * (k1v + 2k2v + 2k3v + k4v),)
end

"""
    separation(p, x₀, v₀, d₀, Δt, n)

Integrate two trajectories `d₀` apart and return the time base, the Euclidean
separation |Δ(x, v)| and the envelope norm √Q.
"""
function separation(p, x₀, v₀, d₀, Δt, n)
    ω₀, δ = p[1], p[2]
    ω₁² = ω₀^2 - δ^2
    t = collect(range(0, step = Δt, length = n + 1))
    x1, v1 = x₀, v₀
    x2, v2 = x₀ + d₀ / sqrt(2), v₀ + d₀ / sqrt(2)
    d = Vector{Float64}(undef, n + 1)
    q = Vector{Float64}(undef, n + 1)
    for i in 1:(n + 1)
        Δx, Δv = x2 - x1, v2 - v1
        d[i] = hypot(Δx, Δv)
        q[i] = sqrt((Δv + δ * Δx)^2 + ω₁² * Δx^2)
        i > n && break
        x1, v1 = rk4_step(t[i], x1, v1, Δt, p)
        x2, v2 = rk4_step(t[i], x2, v2, Δt, p)
    end
    return t, d, q
end

"""
    decay_rate(t, y)

Least-squares slope of -ln y against t with its standard error, over the range
where `y` exceeds `ROUNDOFF_MARGIN * eps()`.
"""
function decay_rate(t, y)
    keep = y .> ROUNDOFF_MARGIN * eps()
    x, z = t[keep], log.(y[keep])
    n = length(x)
    x̄, z̄ = sum(x) / n, sum(z) / n
    sxx = sum(abs2, x .- x̄)
    slope = sum((x .- x̄) .* (z .- z̄)) / sxx
    residual = z .- z̄ .- slope .* (x .- x̄)
    return -slope, sqrt(sum(abs2, residual) / (n - 2) / sxx)
end

# (legend label, ω₀, δ, F, ω, colour, linestyle, family, name printed)
const CASES = (
    (rich("Undamped"), 1.0, 0.0, 0.0, 0.0, PALETTE.blue, :solid, :legacy, "undamped"),
    (rich("Critical, ", it("δ"), " = ", it("ω"), subscript("0"), " = 1"),
        1.0, 1.0, 0.0, 0.0, PALETTE.red, :solid, :legacy, "critical",),
    (rich("Driven, ", it("δ"), " = 0"), 1.0, 0.0, 10.0,
        3.0, PALETTE.purple, :dash, :legacy, "driven, δ = 0",),
    (rich(it("δ"), " = 0.15"), 1.0, 0.15, 0.0, 0.0,
        PALETTE.orange, :dashdot, :added, "δ = 0.15",),
    (rich("Driven, ", it("δ"), " = 0.3"), 1.0, 0.3, 1.0,
        1.6, PALETTE.green, :dot, :added, "driven, δ = 0.3",),
)

function main()
    Δt, t_end = 0.002, 40.0
    n = round(Int, t_end / Δt)

    fig = Figure(size = (1100, 680))
    ax = Axis(fig[2, 1], xlabel = L"Time $t$",
        ylabel = L"Separation $|\Delta(x, v)|$",
        yscale = log10, yticks = logticks(-25, -5; step = 5),)

    legacy, added = Lines[], Lines[]
    results = Dict{String,Vector{Float64}}()
    println("case                     |Δ| rate (formal fit error)   √Q rate        input δ")
    for (_, ω₀, δ, F, ω, colour, style, family, label) in CASES
        t, d, q = separation((ω₀, δ, F, ω), 1.0, 0.0, D₀, Δt, n)
        results[label] = d
        rate_d, σ_d = decay_rate(t, d)
        rate_q, _ = decay_rate(t, q)
        @printf("%-22s  %8.5f (%.5f)            %11.8f     %.2f\n",
            label, rate_d, σ_d, rate_q, δ)
        @assert abs(rate_q - δ) < RATE_TOLERANCE "$label: √Q decays at $rate_q, not δ = $δ"

        # the undamped curve is drawn broad so that the dashed driven one, which
        # lies on it, leaves it visible through its gaps
        width = label == "undamped" ? 6 : 3
        h = lines!(ax, t, d, color = colour, linestyle = style, linewidth = width)
        push!(family == :legacy ? legacy : added, h)
        δ > 0 && text!(ax, t_end, last(d); text = @sprintf(" %.3f", rate_q),
            align = (:left, :center), color = colour, fontsize = ANNOTATION_SIZE,)
    end

    coincidence = maximum(abs.(results["undamped"] ./ results["driven, δ = 0"] .- 1))
    @printf("largest relative difference between the two δ = 0 separations: %.1e\n",
        coincidence)
    @assert coincidence < COINCIDENCE_TOLERANCE "the driving does not cancel: $coincidence"

    text!(ax, 0.97, 0.95;
        text = L"\text{Envelope decay rate of }\sqrt{Q}\text{ at each curve's end}",
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,)
    text!(ax, 0.03, 0.05;
        text = rich("Driven and undamped agree to ", rsci(coincidence; digits = 1)),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.purple,)
    xlims!(ax, -1, t_end * 1.16)
    ylims!(ax, 2e-26, 5e-5)

    Legend(fig[1, 1], [legacy, added],
        [[c[1] for c in CASES if c[8] == :legacy], [c[1] for c in CASES if c[8] == :added]],
        ["2018 sets", "Under-damped"]; titleposition = :left, titlesize = 22,
        titlegap = 14, groupgap = 36, nbanks = 2,)

    path = savefigure(fig, FIGURES, "oscillator_trajectory_separation")
    println("wrote ", path)
end

main()
