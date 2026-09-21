# Separation of two neighbouring trajectories of a linear oscillator, in its
# undamped, under-damped and driven forms:
#
#   ẍ + 2δẋ + ω₀²x = F cos(ω t)
#
# Two initial conditions a distance d₀ apart in phase space are integrated with
# RK4 and their separation |Δ(x, v)| followed. Because the system is linear, the
# difference of two solutions satisfies the *homogeneous* equation whatever the
# driving does -- the forcing is identical along both trajectories and cancels.
# The separation therefore decays with envelope
#
#   |Δ(t)| ~ e^{-δt}
#
# For δ = 0 and ω₀ = 1 the difference vector rotates rigidly in phase space and
# |Δ| is exactly constant. For δ > 0 the rotation is no longer rigid, so |Δ|
# oscillates within the decaying envelope by a relative amount of order δ/ω₀.
# What never happens is growth, which is why the return-time statistics these
# three programs were built around carry no information for a linear system.
# They only become meaningful in 02-chaos-and-attractors, where the same
# measurement on the Lorenz system grows exponentially and yields a positive
# Lyapunov exponent.
#
# Ported from OscilatorArmonique.cpp, OscilatorAmortizat.cpp and
# OscilatorFortat.cpp (2018). All three opened an ifstream on the same file as
# an unflushed ofstream, so every extraction failed and, since C++11, silently
# set its target to zero: the separation and extremum-interval results those
# programs printed were identically zero and no real number was ever produced.
# The analysis is done here in memory, so the failure mode cannot recur. The
# parameter choices were also degenerate -- the damped case used δ = ω₀ = 1,
# exactly critical damping, which does not oscillate at all, and the driven case
# used δ = 0, so it never reaches a steady state. Both are given physically
# useful values below — but the original parameter sets are integrated here too,
# so the record shows what those three programs were actually configured to do
# and why the analysis built on top of them could not work.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Initial phase-space separation of the two trajectories."
const D₀ = 1e-7

"""
    rhs(t, x, v, ω₀, δ, F, ω)

Right-hand side of the driven damped oscillator written as a first-order system.
"""
rhs(t, x, v, ω₀, δ, F, ω) = (v, F * cos(ω * t) - 2δ * v - ω₀^2 * x)

"""
    rk4_step(t, x, v, Δt, p)

One classical RK4 step of the two-component system. `p` carries `(ω₀, δ, F, ω)`.
The driving term is evaluated at the correct stage times, which is what makes
this valid for the non-autonomous case.
"""
function rk4_step(t, x, v, Δt, p)
    k1x, k1v = rhs(t, x, v, p...)
    k2x, k2v = rhs(t + Δt/2, x + Δt*k1x/2, v + Δt*k1v/2, p...)
    k3x, k3v = rhs(t + Δt/2, x + Δt*k2x/2, v + Δt*k2v/2, p...)
    k4x, k4v = rhs(t + Δt, x + Δt*k3x, v + Δt*k3v, p...)
    return (x + Δt/6 * (k1x + 2k2x + 2k3x + k4x),
        v + Δt/6 * (k1v + 2k2v + 2k3v + k4v),)
end

"""
    separation(p, x₀, v₀, d₀, Δt, n)

Integrate two trajectories differing initially by `d₀` in both components and
return the time base together with their Euclidean phase-space separation.
"""
function separation(p, x₀, v₀, d₀, Δt, n)
    t = range(0, step = Δt, length = n + 1)
    x1, v1 = x₀, v₀
    x2, v2 = x₀ + d₀ / sqrt(2), v₀ + d₀ / sqrt(2)
    d = Vector{Float64}(undef, n + 1)
    d[1] = hypot(x2 - x1, v2 - v1)
    for i in 1:n
        x1, v1 = rk4_step(t[i], x1, v1, Δt, p)
        x2, v2 = rk4_step(t[i], x2, v2, Δt, p)
        d[i + 1] = hypot(x2 - x1, v2 - v1)
    end
    return collect(t), d
end

"""
    fitted_decay_rate(t, d)

Least-squares slope of ln|Δ| against t, over the range where the separation is
still well above round-off. For this system the slope recovers -δ regardless of
the oscillation inside the envelope, which a pointwise comparison against
d₀e^{-δt} would not.
"""
function fitted_decay_rate(t, d)
    keep = d .> 1e3 * eps()
    x, y = t[keep], log.(d[keep])
    n = length(x)
    slope = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(abs2, x) - sum(x)^2)
    return -slope
end

# (label, ω₀, δ, F, ω, colour, linestyle, linewidth)
#
# The first three are the parameter sets of OscilatorArmonique.cpp,
# OscilatorAmortizat.cpp and OscilatorFortat.cpp exactly as written. Two of them
# are degenerate: δ = ω₀ = 1 is critical damping, which does not oscillate, and
# the driven case has δ = 0, so it never reaches a steady state and its damping
# term is dead code. The last two are non-degenerate replacements.
#
# The undamped and driven 2018 cases have *identical* separation to every digit,
# so one curve would otherwise be buried under the other. They are drawn broad
# solid and narrow dashed, in that order, so that both are visible: the blue
# shows through the gaps in the purple. That coincidence is the result, and a
# figure that hid it would be hiding the point.
const CASES = (
    ("Undamped (2018)", 1.0, 0.0, 0.0, 0.0, PALETTE.blue, :solid, 4.0),
    ("Critically damped (2018)", 1.0, 1.0, 0.0, 0.0, PALETTE.red, :solid, 1.6),
    ("Driven, δ = 0 (2018)", 1.0, 0.0, 10.0, 3.0, PALETTE.purple, :dash, 1.8),
    ("Under-damped, δ = 0.15", 1.0, 0.15, 0.0, 0.0, PALETTE.orange, :dash, 1.6),
    ("Driven, δ = 0.3, ω = 1.6", 1.0, 0.3, 1.0, 1.6, PALETTE.green, :dash, 1.6),
)

function main()
    Δt, t_end = 0.002, 60.0
    n = round(Int, t_end / Δt)

    fig = Figure(size = (900, 560))
    ax = Axis(fig[2, 1],
        xlabel = L"Time $t$",
        ylabel = L"Phase-space separation $|\Delta(x, v)|$",
        yscale = log10, yticks = logticks(-30, -5; step = 5),)

    # the reference goes down first and stays thin, so that the two curves
    # sitting exactly on it are not overdrawn by their own guide
    h_ref = hlines!(ax, [D₀], color = (PALETTE.black, 0.45), linewidth = 0.8)

    handles = []
    for (label, ω₀, δ, F, ω, colour, ls, lw) in CASES
        t, d = separation((ω₀, δ, F, ω), 1.0, 0.0, D₀, Δt, n)
        push!(handles, lines!(ax, t, d, color = colour, linewidth = lw, linestyle = ls))

        @printf("%-26s |Δ(%.0f)| = %.3e   fitted decay rate = %.5f   input δ = %.5f\n",
            label, t_end, last(d), fitted_decay_rate(t, d), δ)
    end
    println()
    println("the driving cancels in the difference, so the δ = 0 driven case has")
    println("exactly the same separation as the δ = 0 undamped one — which is why")
    println("OscilatorFortat.cpp could not have learned anything from it even if")
    println("its file handles had worked.")

    text!(ax, 0.02, 0.04;
        text = rich(
            "The difference obeys the homogeneous equation: it decays inside the envelope ",
            it("e"), superscript(rich("−δt", font = :italic)), " and never grows.",),
        space = :relative, align = (:left, :bottom), fontsize = 15,)
    # Headroom above d₀ so that the note sits beside the two curves it is about
    # rather than across the ones it is not.
    ylims!(ax, 1e-33, 1e-4)
    text!(ax, 0.98, 0.99;
        text = "The two 2018 δ = 0 cases coincide to every digit: the driving\nis the same along both trajectories and cancels in their difference",
        space = :relative, align = (:right, :top), fontsize = 15,
        color = PALETTE.purple,)

    Legend(fig[1, 1], [handles; h_ref],
        [[c[1] for c in CASES]; L"Initial separation $d_0 = 10^{-7}$"],
        orientation = :horizontal, framevisible = false, labelsize = 14,
        colgap = 14, nbanks = 2,)

    rowsize!(fig.layout, 2, Relative(0.86))
    path = savefigure(fig, FIGURES, "oscillator_trajectory_separation")
    println("wrote ", path)
end

main()
