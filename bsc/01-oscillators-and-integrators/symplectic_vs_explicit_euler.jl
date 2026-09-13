# Explicit Euler against the two symplectic (Euler-Cromer) variants, on the
# harmonic oscillator ẍ = -Ω²x written as the first-order system
#
#   ẋ = v,   v̇ = -Ω²x
#
# The point is that the three schemes are identical to first order in Δt yet
# behave completely differently over long times: explicit Euler pumps energy in
# without bound, while both symplectic variants conserve a nearby "shadow"
# energy and so keep their orbits closed.
#
# Ported from Spatiul_FazelorEcDiff.cpp and PbEcDiffExamen.cpp (2018). The
# original named its first scheme EulerImplicit, but it advances both components
# from the old values, which is explicit Euler; the misnomer propagated into six
# output filenames, so everything labelled "implicit" in that project was in
# fact about the explicit method. It also worked in single precision on a study
# whose entire subject is accumulated integration error, drew initial conditions
# as integers on a 101x101 lattice with modulo bias, and seeded from the clock.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Angular frequency of the oscillator."
const Ω = 1.0
"Integration horizon used for the phase portraits."
const T_END = 40.0

energy(x, v, Ω) = 0.5 * v^2 + 0.5 * Ω^2 * x^2

"""
    explicit_euler(x, v, Δt, Ω)

Both components advanced from the old state. Amplification factor per step is
√(1 + Ω²Δt²) > 1 for every Δt, so the orbit spirals outward without bound.
"""
explicit_euler(x, v, Δt, Ω) = (x + Δt * v, v - Δt * Ω^2 * x)

"""
    symplectic_position_first(x, v, Δt, Ω)

Position updated first, then velocity using the *new* position. Symplectic.
"""
function symplectic_position_first(x, v, Δt, Ω)
    x_new = x + Δt * v
    return (x_new, v - Δt * Ω^2 * x_new)
end

"""
    symplectic_velocity_first(x, v, Δt, Ω)

Velocity updated first, then position using the *new* velocity. Symplectic.
"""
function symplectic_velocity_first(x, v, Δt, Ω)
    v_new = v - Δt * Ω^2 * x
    return (x + Δt * v_new, v_new)
end

"""
    integrate(step, x₀, v₀, Δt, n, Ω)

Apply `step` `n` times. Returns the time, position, velocity and energy traces.
"""
function integrate(step, x₀, v₀, Δt, n, Ω)
    t = Vector{Float64}(undef, n + 1)
    x = Vector{Float64}(undef, n + 1)
    v = Vector{Float64}(undef, n + 1)
    e = Vector{Float64}(undef, n + 1)
    t[1], x[1], v[1], e[1] = 0.0, x₀, v₀, energy(x₀, v₀, Ω)
    for i in 1:n
        x[i + 1], v[i + 1] = step(x[i], v[i], Δt, Ω)
        t[i + 1] = i * Δt
        e[i + 1] = energy(x[i + 1], v[i + 1], Ω)
    end
    return t, x, v, e
end

const SCHEMES = (
    (explicit_euler, "Explicit Euler", PALETTE.red),
    (symplectic_position_first, "Symplectic, position first", PALETTE.blue),
    (symplectic_velocity_first, "Symplectic, velocity first", PALETTE.green),
)

function main()
    Δt = 0.05
    n = round(Int, T_END / Δt)

    # distinct radii: a ring of equal-energy points would all trace one orbit
    starts = [(r, 0.0) for r in (0.8, 1.6, 2.4, 3.2)]
    push!(starts, (1.0, -1.0))   # the initial condition used by PbEcDiffExamen.cpp

    fig = Figure(size = (1050, 780))

    # Two nested layouts rather than one grid: the phase portraits are square by
    # DataAspect, and if they shared columns with the row below they would be
    # spaced by that row's very different column widths.
    top = GridLayout(fig[2, 1])
    bottom = GridLayout(fig[3, 1])

    for (k, (step, name, colour)) in enumerate(SCHEMES)
        ax = Axis(top[1, k], aspect = DataAspect(),
            xlabel = L"Position $x$", ylabel = k == 1 ? L"Velocity $v$" : "")
        for (x₀, v₀) in starts
            _, x, v, _ = integrate(step, x₀, v₀, Δt, n, Ω)
            lines!(ax, x, v, color = colour, linewidth = 0.9)
        end
        limits!(ax, -6, 6, -6, 6)
        k > 1 && hideydecorations!(ax, grid = false)
    end

    # energy drift along one orbit
    ax_e = Axis(bottom[1, 1],
        xlabel = L"Time $t$", ylabel = L"Energy $E(t)/E(0)$", yscale = log10,
        yticks = ([1, 2, 3, 5, 8], [L"1", L"2", L"3", L"5", L"8"]))
    handles = []
    for (step, name, colour) in SCHEMES
        t, _, _, e = integrate(step, 1.0, -1.0, Δt, n, Ω)
        push!(handles, lines!(ax_e, t, e ./ first(e), color = colour, linewidth = 1.5))
    end
    hlines!(ax_e, [1.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    text!(ax_e, 0.5, 1.075; text = L"Initial energy $E(0)$", space = :data,
        align = (:left, :bottom), fontsize = 15, color = PALETTE.black)

    # growth factor predicted for explicit Euler: (1 + Ω²Δt²)^(N/2) on the amplitude,
    # hence (1 + Ω²Δt²)^N on the energy
    predicted = (1 + Ω^2 * Δt^2)^n
    @printf("explicit Euler energy growth over %.0f time units: predicted %.3f, ", T_END, predicted)
    _, _, _, e_ex = integrate(explicit_euler, 1.0, -1.0, Δt, n, Ω)
    observed = last(e_ex) / first(e_ex)
    @printf("observed %.3f\n", observed)

    text!(ax_e, 0.97, 0.93;
        text = latexstring(@sprintf("(1 + \\Omega^2 \\Delta t^2)^N = %.3f \\text{ predicted,~} \
                                     %.3f \\text{ measured}", predicted, observed)),
        space = :relative, align = (:right, :top),
        fontsize = 15, color = PALETTE.red)

    # energy error against step count, the sweep PbEcDiffExamen.cpp attempted
    ax_n = Axis(bottom[1, 2],
        xlabel = L"Steps $N$ over $t \in [0,\, 10]$",
        ylabel = L"$|E(10)/E(0) - 1|$",
        xscale = log10, yscale = log10,
        xticks = logticks(1, 4; style = :decimal),
        yticks = logticks(-4, 2; step = 2))
    Ns = [2^k for k in 4:16]
    for (step, name, colour) in SCHEMES
        errs = map(Ns) do N
            _, _, _, e = integrate(step, 1.0, -1.0, 10.0 / N, N, Ω)
            abs(last(e) / first(e) - 1)
        end
        scatterlines!(ax_n, Ns, errs, color = colour, linewidth = 1.3,
            markersize = MARKERSIZE.dense)
    end

    Legend(fig[1, 1], handles, [s[2] for s in SCHEMES],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)

    colsize!(bottom, 1, Relative(0.66))
    rowsize!(fig.layout, 2, Relative(0.46))
    rowgap!(fig.layout, 12)
    colgap!(bottom, 30)
    path = savefigure(fig, FIGURES, "symplectic_vs_explicit_euler")
    println("wrote ", path)
end

main()
