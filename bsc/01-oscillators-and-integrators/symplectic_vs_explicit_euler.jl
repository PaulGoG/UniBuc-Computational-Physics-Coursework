# Explicit Euler against the two symplectic (Euler–Cromer) variants on the
# harmonic oscillator ẍ = -Ω²x written as the first-order system
#
#   ẋ = v,   v̇ = -Ω²x
#
# The three schemes agree to first order in Δt and differ completely over long
# times: explicit Euler multiplies the amplitude by √(1 + Ω²Δt²) at every step,
# so its energy grows as (1 + Ω²Δt²)^N without bound, while each symplectic
# variant conserves the modified energy E ∓ ½ΔtΩ²xv exactly and so keeps its
# orbit closed and its energy oscillating about E(0). Both statements are
# asserted below.
#
# Ported from Spatiul_FazelorEcDiff.cpp and PbEcDiffExamen.cpp on the legacy
# branch. The first names its explicit scheme EulerImplicit, works in single
# precision, and draws its initial conditions as integers on a 101 × 101 lattice
# with modulo bias from a clock-seeded rand(); the second sweeps the step count
# for the velocity-first variant from the initial condition (1, -1) used here.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Angular frequency of the oscillator."
const Ω = 1.0
"Integration horizon of the phase portraits and the energy traces."
const T_END = 40.0
"Step of the phase portraits and the energy traces."
const Δt = 0.05
"Relative tolerance on the conserved modified energy of the symplectic schemes."
const SHADOW_TOLERANCE = 1e-12
"Relative tolerance on the energy growth of explicit Euler against (1 + Ω²Δt²)^N."
const GROWTH_TOLERANCE = 1e-12

energy(x, v, Ω) = 0.5 * v^2 + 0.5 * Ω^2 * x^2

"""
    explicit_euler(x, v, Δt, Ω)

Both components advanced from the old state. The amplification factor per step
is √(1 + Ω²Δt²) > 1 for every Δt, so the orbit spirals outward without bound.
"""
explicit_euler(x, v, Δt, Ω) = (x + Δt * v, v - Δt * Ω^2 * x)

"""
    symplectic_position_first(x, v, Δt, Ω)

Position updated first, then velocity from the new position. Conserves
`E + ½ΔtΩ²xv` exactly.
"""
function symplectic_position_first(x, v, Δt, Ω)
    x_new = x + Δt * v
    return (x_new, v - Δt * Ω^2 * x_new)
end

"""
    symplectic_velocity_first(x, v, Δt, Ω)

Velocity updated first, then position from the new velocity. Conserves
`E - ½ΔtΩ²xv` exactly.
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

"""
    shadow_energy_drift(x, v, Δt, Ω, sign)

Largest relative excursion of the modified energy `E + sign·½ΔtΩ²xv` along a
trajectory; zero up to round-off for the symplectic scheme it belongs to.
"""
function shadow_energy_drift(x, v, Δt, Ω, sign)
    h = energy.(x, v, Ω) .+ sign * 0.5 * Δt * Ω^2 .* x .* v
    return maximum(abs.(h ./ h[1] .- 1))
end

# (step, legend label, colour, line style, sign of the conserved xv term)
const SCHEMES = (
    (explicit_euler, "Explicit Euler", PALETTE.blue, :solid, nothing),
    (symplectic_position_first, "Symplectic, position first", PALETTE.orange, :dash, +1),
    (symplectic_velocity_first, "Symplectic, velocity first", PALETTE.green, :dot, -1),
)

function main()
    n = round(Int, T_END / Δt)

    # distinct radii: a ring of equal-energy points would all trace one orbit
    starts = [(r, 0.0) for r in (0.8, 1.6, 2.4, 3.2)]
    push!(starts, (1.0, -1.0))   # the initial condition of PbEcDiffExamen.cpp

    fig = Figure(size = (1200, 1000))

    # Two nested layouts: the phase portraits are square by DataAspect, and if
    # they shared columns with the row below they would be spaced by that row's
    # very different column widths.
    top = GridLayout(fig[2, 1])
    bottom = GridLayout(fig[3, 1])

    for (k, (step, _, colour, style, _)) in enumerate(SCHEMES)
        # ±6 sit on the frame corners and would collide at the panel joints
        ax = Axis(top[1, k], aspect = DataAspect(), xticks = -4:2:4, yticks = -4:2:4,
            xlabel = L"Position $x$", ylabel = k == 1 ? L"Velocity $v$" : "",)
        for (x₀, v₀) in starts
            _, x, v, _ = integrate(step, x₀, v₀, Δt, n, Ω)
            lines!(ax, x, v, color = colour, linestyle = style)
        end
        limits!(ax, -6, 6, -6, 6)
        k > 1 && hideydecorations!(ax, grid = false)
    end

    # energy along the orbit from (1, -1)
    ax_e = Axis(bottom[1, 1],
        xlabel = L"Time $t$", ylabel = L"Energy $E(t)/E(0)$", yscale = log10,
        yticks = ([1, 2, 3, 5, 8], [L"1", L"2", L"3", L"5", L"8"]),)
    handles = Lines[]
    for (step, name, colour, style, sign) in SCHEMES
        t, x, v, e = integrate(step, 1.0, -1.0, Δt, n, Ω)
        push!(handles, lines!(ax_e, t, e ./ first(e), color = colour, linestyle = style))
        sign === nothing && continue
        drift = shadow_energy_drift(x, v, Δt, Ω, sign)
        @printf("%-27s modified energy E %s ½ΔtΩ²xv conserved to %.1e\n",
            name, sign < 0 ? "-" : "+", drift)
        @assert drift < SHADOW_TOLERANCE "$name: modified energy drifts by $drift"
    end
    hlines!(ax_e, [1.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    text!(ax_e, 0.5, 1.09; text = L"Initial energy $E(0)$", space = :data,
        align = (:left, :bottom), fontsize = ANNOTATION_SIZE, color = PALETTE.black,)

    # explicit Euler multiplies the amplitude by √(1 + Ω²Δt²) per step, the energy
    # by (1 + Ω²Δt²) per step
    predicted = (1 + Ω^2 * Δt^2)^n
    _, _, _, e_ex = integrate(explicit_euler, 1.0, -1.0, Δt, n, Ω)
    observed = last(e_ex) / first(e_ex)
    @printf("explicit Euler energy growth over %.0f time units: predicted %.3f, observed %.3f\n",
        T_END, predicted, observed)
    @assert abs(observed / predicted - 1) < GROWTH_TOLERANCE "explicit Euler growth $observed against $predicted"

    text!(ax_e, 0.03, 0.95;
        text = latexstring(@sprintf("(1 + \\Omega^2 \\Delta t^2)^N = %.3f\\text{ predicted, }%.3f\\text{ measured}",
            predicted, observed)),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,)

    # energy error against step count, the sweep PbEcDiffExamen.cpp made
    ax_n = Axis(bottom[1, 2],
        xlabel = L"Steps $N$ over $t \in [0,\, 10]$",
        ylabel = L"$|E(10)/E(0) - 1|$",
        xscale = log10, yscale = log10,
        xticks = logticks(1, 4; style = :decimal),
        yticks = logticks(-4, 2; step = 2),)
    Ns = [2^k for k in 4:16]
    for (step, _, colour, style, _) in SCHEMES
        errs = map(Ns) do N
            _, _, _, e = integrate(step, 1.0, -1.0, 10.0 / N, N, Ω)
            abs(last(e) / first(e) - 1)
        end
        scatterlines!(ax_n, Ns, errs, color = colour, linestyle = style)
    end

    Legend(fig[1, 1], handles, [s[2] for s in SCHEMES])

    colsize!(bottom, 1, Relative(0.64))
    rowsize!(fig.layout, 2, Auto(0.85))
    colgap!(bottom, 24)
    path = savefigure(fig, FIGURES, "symplectic_vs_explicit_euler")
    println("wrote ", path)
end

main()
