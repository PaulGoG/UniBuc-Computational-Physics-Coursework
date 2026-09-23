# The Lorenz and Rössler systems integrated with classical RK4, shown as
# projections of their strange attractors.
#
#   Lorenz    ẋ = σ(y - x),  ẏ = x(ρ - z) - y,  ż = xy - βz
#   Rössler   ẋ = -y - z,    ẏ = x + ay,        ż = b + z(x - c)
#
# Ported from Integratori.cpp and Atractori.cpp in
# Code_Archive/Old_2018/C_C++/AtractorI/ on the `legacy` branch. The RK4 stage
# coupling there is correct. β is written as 2.66 rather than 8/3, which moves
# the fixed points (±√(β(ρ-1)), ±√(β(ρ-1)), ρ-1) from ±8.4853 to ±8.4747; the
# Rössler branch is unreachable because the selector variable is assigned and
# never read; both systems write to hardcoded Windows paths; and the parameters
# are #define macros, so `a`, `b` and `c` textually replace any identifier of
# those names in the translation unit.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "attractors_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"""
Line width of the attractor renderings. A trajectory of 2.5 × 10⁵ points drawn
at the theme's 3 fills each lobe solid; a thin line is the analogue of the
theme's `cloud` marker for a dense series.
"""
const TRACE_WIDTH = 0.8

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

    fig = Figure(size = (1400, 620))

    ax1 = Axis(fig[1, 1], xlabel = L"x", ylabel = L"z")
    lines!(ax1, xL, zL, color = PALETTE.blue, linewidth = TRACE_WIDTH)
    text!(ax1, 0.02, 0.98;
        text = rich("Lorenz: ", it("σ"), " = 10, ", it("ρ"), " = 28, ", it("β"), " = 8/3"),
        space = :relative, align = (:left, :top), color = PALETTE.blue,
        fontsize = ANNOTATION_SIZE,)
    # the fixed points sit inside the lobes; their label goes in the notch
    # between the wings, the one region of this projection no strand enters
    scatter!(ax1, [fp[1][1], fp[2][1]], [fp[1][3], fp[2][3]],
        color = PALETTE.red, markersize = MARKERSIZE.emphasis, marker = :xcross,)
    text!(ax1, 0.5, 0.88;
        text = rich("✕  Fixed points\n(±8.4853, 27)"),
        space = :relative, align = (:center, :top), color = PALETTE.red,
        fontsize = ANNOTATION_SIZE, justification = :center,)
    ylims!(ax1, 0, 53)

    ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y")
    lines!(ax2, xR, yR, color = PALETTE.green, linewidth = TRACE_WIDTH)
    text!(ax2, 0.02, 0.98;
        text = rich("Rössler: ", it("a"), " = ", it("b"), " = 0.2, ", it("c"), " = 5.7"),
        space = :relative, align = (:left, :top), color = PALETTE.green,
        fontsize = ANNOTATION_SIZE,)
    ylims!(ax2, nothing, 12)

    path = savefigure(fig, FIGURES, "lorenz_rossler_attractors")
    println("wrote ", path)
    println("wrote ", animate_attractors(xL, zL, xR, yR, h))
end

"""
    animate_attractors(xL, zL, xR, yR, h)

Trace both attractors out in time, with a head marker on the current state.

The static figure is the invariant set — where the trajectory eventually goes.
The animation is the trajectory itself, and it shows the thing the static plot
cannot: on the Lorenz attractor the state circles one lobe an unpredictable
number of times before crossing to the other, which is the sensitivity that
makes the system chaotic. On the Rössler attractor it spirals outward in a
near-plane and is folded back, which is the simpler mechanism producing the same
kind of set.
"""
function animate_attractors(xL, zL, xR, yR, h)
    n = min(length(xL), length(xR))
    frames = 160
    stride = n ÷ frames
    tail = 25_000            # points kept behind the head, about 50 time units

    headL = Observable(Point2f[])
    headR = Observable(Point2f[])
    dotL = Observable(Point2f[])
    dotR = Observable(Point2f[])
    caption = Observable(rich(it("t"), " = 0.0"))

    fig = Figure(size = (1200, 560))
    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"z")
    lines!(ax1, headL, color = PALETTE.blue, linewidth = TRACE_WIDTH)
    # black for the moving state: red names the fixed points in the static figure
    scatter!(ax1, dotL, color = PALETTE.black, markersize = MARKERSIZE.data)
    xlims!(ax1, minimum(xL) - 2, maximum(xL) + 2)
    ylims!(ax1, minimum(zL) - 2, maximum(zL) + 6)
    text!(ax1, 0.02, 0.98; text = "Lorenz", space = :relative,
        align = (:left, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.blue,)

    ax2 = Axis(fig[2, 2], xlabel = L"x", ylabel = L"y")
    lines!(ax2, headR, color = PALETTE.green, linewidth = TRACE_WIDTH)
    scatter!(ax2, dotR, color = PALETTE.black, markersize = MARKERSIZE.data)
    xlims!(ax2, minimum(xR) - 2, maximum(xR) + 2)
    ylims!(ax2, minimum(yR) - 2, maximum(yR) + 4)
    text!(ax2, 0.02, 0.98; text = "Rössler", space = :relative,
        align = (:left, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.green,)

    Label(fig[1, 1:2], caption, fontsize = 22, tellwidth = false)
    rowgap!(fig.layout, 6)

    path = joinpath(FIGURES, "lorenz_rossler_attractors.gif")
    mkpath(FIGURES)
    record(fig, path, 1:frames; framerate = 14, px_per_unit = 1) do k
        j = k * stride
        i = max(1, j - tail)
        headL[] = Point2f.(view(xL, i:j), view(zL, i:j))
        headR[] = Point2f.(view(xR, i:j), view(yR, i:j))
        dotL[] = [Point2f(xL[j], zL[j])]
        dotR[] = [Point2f(xR[j], yR[j])]
        caption[] = rich(it("t"), @sprintf(" = %.1f", j * h))
    end
    return path
end

main()
