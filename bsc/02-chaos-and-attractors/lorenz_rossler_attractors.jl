# The Lorenz and Rössler systems integrated with classical RK4, shown as
# projections of their strange attractors.
#
#   Lorenz    ẋ = σ(y - x),  ẏ = x(ρ - z) - y,  ż = xy - βz
#   Rössler   ẋ = -y - z,    ẏ = x + ay,        ż = b + z(x - c)
#
# Ported from Integratori.cpp and Atractori.cpp (2018). The RK4 stage coupling
# in the original was correct. What was not: β was written as 2.66 rather than
# 8/3, which displaces the fixed points (±√(β(ρ-1)), ±√(β(ρ-1)), ρ-1); the
# Rössler branch was unreachable because the selector variable was assigned and
# never read; both systems wrote to hardcoded Windows paths; and the parameters
# were #define macros, so `a`, `b` and `c` textually replaced any identifier of
# those names anywhere in the translation unit.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "attractors_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"""
    trajectory(f, u₀, h, n, p; transient)

Integrate `f` for `n` steps of size `h`, discarding the first `transient` steps
so that the orbit has settled onto the attractor. Returns three coordinate
vectors.
"""
function trajectory(f, u₀, h, n, p; transient = 0)
    u = Tuple(float.(u₀))
    for _ in 1:transient
        u = rk4_step(f, u, h, p)
    end
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    zs = Vector{Float64}(undef, n)
    for i in 1:n
        u = rk4_step(f, u, h, p)
        xs[i], ys[i], zs[i] = u
    end
    return xs, ys, zs
end


function main()
    h = 0.002
    n = 250_000

    xL, yL, zL = trajectory(lorenz, (1.0, 1.0, 1.0), h, n, LORENZ; transient = 15_000)
    xR, yR, zR = trajectory(rossler, (1.0, 1.0, 1.0), h, n, ROSSLER; transient = 15_000)

    fp = lorenz_fixed_points(LORENZ)
    @printf("Lorenz fixed points at (±%.4f, ±%.4f, %.1f)\n", fp[1][1], fp[1][2], fp[1][3])
    @printf("  with the 2018 value beta = 2.66 they would sit at ±%.4f\n",
            sqrt(2.66 * (LORENZ.ρ - 1)))
    @printf("Lorenz  x range [%.2f, %.2f], z range [%.2f, %.2f]\n",
            minimum(xL), maximum(xL), minimum(zL), maximum(zL))
    @printf("Rossler x range [%.2f, %.2f], z range [%.2f, %.2f]\n",
            minimum(xR), maximum(xR), minimum(zR), maximum(zR))

    fig = Figure(size = (1000, 430))

    ax1 = Axis(fig[1, 1], xlabel = L"x", ylabel = L"z",
        title = L"Lorenz: $\sigma = 10$, $\rho = 28$, $\beta = 8/3$", titlesize = 18)
    lines!(ax1, xL, zL, color = PALETTE.blue, linewidth = 0.25)
    scatter!(ax1, [fp[1][1], fp[2][1]], [fp[1][3], fp[2][3]],
        color = PALETTE.red, markersize = 11, marker = :xcross)
    text!(ax1, 0.5, 0.02; text = L"Fixed points $(\pm\sqrt{\beta(\rho-1)},\ \rho-1)$",
        space = :relative, align = (:center, :bottom), color = PALETTE.red, fontsize = 15)

    ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y",
        title = L"Rössler: $a = b = 0.2$, $c = 5.7$", titlesize = 18)
    lines!(ax2, xR, yR, color = PALETTE.green, linewidth = 0.25)

    path = savefigure(fig, FIGURES, "lorenz_rossler_attractors")
    println("wrote ", path)
end

main()
