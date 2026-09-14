# The pendulum, linearised and full, and what the original RK4 actually integrated.
#
#   θ̈ = -(g/L) sin θ - q θ̇ + F_D sin(Ω_D t)
#
# Ported from Pendul_simplu.m, Pendul_dampat.m, Pendul_dampat_fortat.m and
# Runge_Kutta_4_pendul_dampat_fortat.m (Octave).
#
# The serious defect is in the RK4 file. For the coupled system θ̇ = ω,
# ω̇ = F(θ, ω, t), the θ-stage increments must use the θ-slopes. The original
# wrote
#
#     k2 = frhs(omega + dt*k1/2, theta + dt*k1/2, t + dt/2)
#     pk1 = omega(step);  pk2 = omega(step) + dt*pk1/2;  ...
#
# so `k1`, which is dω/dt, advanced *both* components, and the θ-stages
# propagated ω as though dω/dt = ω. Expanding the θ update gives
# dt·ω + dt²·ω/2 where RK4 requires dt·ω + dt²·F/2. The scheme is therefore
# globally **first order**, not fourth — it is not the RK4 its filename claims.
# Both versions are integrated below and their observed orders measured.
#
# Two lesser points: Pendul_simplu.m is headed "pendulul matematic" but
# integrates -(g/L)θ, the small-angle approximation, not -(g/L)sin θ; and all
# four files assigned `length = 1` or `length = 9.8`, shadowing the Octave
# builtin.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Giordano's chaotic parameter set for the driven damped pendulum."
const CHAOTIC = (g_over_L = 1.0, q = 0.5, F_D = 1.2, Ω_D = 2/3)

"Angular acceleration of the full nonlinear pendulum."
accel(θ, ω, t, p) = -p.g_over_L * sin(θ) - p.q * ω + p.F_D * sin(p.Ω_D * t)
"Angular acceleration under the small-angle approximation."
accel_linear(θ, ω, t, p) = -p.g_over_L * θ - p.q * ω + p.F_D * sin(p.Ω_D * t)

"""
    rk4_step(θ, ω, t, Δt, p, f)

Correct RK4 for the coupled system: the θ-stages use the ω-values, the ω-stages
use the accelerations.
"""
function rk4_step(θ, ω, t, Δt, p, f)
    k1θ, k1ω = ω, f(θ, ω, t, p)
    k2θ, k2ω = ω + Δt*k1ω/2, f(θ + Δt*k1θ/2, ω + Δt*k1ω/2, t + Δt/2, p)
    k3θ, k3ω = ω + Δt*k2ω/2, f(θ + Δt*k2θ/2, ω + Δt*k2ω/2, t + Δt/2, p)
    k4θ, k4ω = ω + Δt*k3ω, f(θ + Δt*k3θ, ω + Δt*k3ω, t + Δt, p)
    return (θ + Δt/6*(k1θ + 2k2θ + 2k3θ + k4θ),
        ω + Δt/6*(k1ω + 2k2ω + 2k3ω + k4ω),)
end

"""
    rk4_step_original(θ, ω, t, Δt, p, f)

The stage coupling exactly as written in the original, kept so its order can be
measured. `k1` is dω/dt but advances both components, and the θ-stages
propagate ω as though dω/dt = ω.
"""
function rk4_step_original(θ, ω, t, Δt, p, f)
    k1 = f(θ, ω, t, p)
    k2 = f(θ + Δt*k1/2, ω + Δt*k1/2, t + Δt/2, p)
    k3 = f(θ + Δt*k2/2, ω + Δt*k2/2, t + Δt/2, p)
    k4 = f(θ + Δt*k3, ω + Δt*k3, t + Δt, p)
    ω_new = ω + Δt/6*(k1 + 2k2 + 2k3 + k4)

    pk1 = ω
    pk2 = ω + Δt*pk1/2
    pk3 = ω + Δt*pk2/2
    pk4 = ω + Δt*pk3
    θ_new = θ + Δt/6*(pk1 + 2pk2 + 2pk3 + pk4)
    return (θ_new, ω_new)
end

"Wrap an angle to (-π, π]."
wrap(θ) = mod(θ + π, 2π) - π

"""
    wrapped_trace(t, θ)

`θ` wrapped to (-π, π] against `t`, with a `NaN` pair inserted at every wrap, so
that a line plot breaks there instead of drawing a vertical stroke across the
panel. Used for the original stage coupling, which winds rather than oscillates.
"""
function wrapped_trace(t, θ)
    xs, ys = Float64[], Float64[]
    w = wrap.(θ)
    for i in eachindex(w)
        if i > 1 && abs(w[i] - w[i - 1]) > π
            push!(xs, NaN)
            push!(ys, NaN)
        end
        push!(xs, t[i])
        push!(ys, w[i])
    end
    return xs, ys
end

"""
    integrate(step, θ₀, ω₀, Δt, n, p, f)

Apply `step` `n` times, returning time, angle and angular-velocity traces.
"""
function integrate(step, θ₀, ω₀, Δt, n, p, f)
    t = Vector{Float64}(undef, n + 1)
    θ = Vector{Float64}(undef, n + 1)
    ω = Vector{Float64}(undef, n + 1)
    t[1], θ[1], ω[1] = 0.0, θ₀, ω₀
    for i in 1:n
        θ[i + 1], ω[i + 1] = step(θ[i], ω[i], t[i], Δt, p, f)
        t[i + 1] = i * Δt
    end
    return t, θ, ω
end

"""
    observed_order(step, p, f, t_end)

Global error at `t_end` against a finely resolved reference, swept over step
size, reduced to a least-squares slope.
"""
function observed_order(step, p, f, t_end)
    Δt_ref = 1e-5
    _, θr, _ = integrate(rk4_step, 0.2, 0.0, Δt_ref, round(Int, t_end/Δt_ref), p, f)
    reference = last(θr)
    hs = [0.02, 0.01, 0.005, 0.0025, 0.00125]
    errs = map(hs) do h
        _, θ, _ = integrate(step, 0.2, 0.0, h, round(Int, t_end/h), p, f)
        abs(last(θ) - reference)
    end
    lx, ly = log.(hs), log.(errs)
    n = length(lx)
    return (n*sum(lx .* ly) - sum(lx)*sum(ly)) / (n*sum(abs2, lx) - sum(lx)^2), hs, errs
end

function main()
    # a mild, non-chaotic parameter set for the order study
    mild = (g_over_L = 1.0, q = 0.5, F_D = 0.2, Ω_D = 2/3)
    p_correct, hs, e_correct = observed_order(rk4_step, mild, accel, 5.0)
    p_original, _, e_original = observed_order(rk4_step_original, mild, accel, 5.0)
    @printf("observed order, correct RK4 = %.2f\n", p_correct)
    @printf("observed order, original stages = %.2f\n", p_original)

    # small-angle against full nonlinear, undriven and undamped.
    # Pendul_simplu.m used L = 1 m with g = 9.8, so g/L = 9.8 — not the g/L = 1
    # of the driven files, which set L = 9.8 m.
    free = (g_over_L = 9.8, q = 0.0, F_D = 0.0, Ω_D = 0.0)
    Δt = 0.002
    n = round(Int, 12 / Δt)
    t, θ_full, _ = integrate(rk4_step, 2.5, 0.0, Δt, n, free, accel)
    _, θ_lin, _ = integrate(rk4_step, 2.5, 0.0, Δt, n, free, accel_linear)
    @printf("free pendulum from θ₀ = 2.5 rad: full and linearised differ by up to %.3f rad\n",
        maximum(abs.(θ_full .- θ_lin)))

    # chaotic attractor, Poincaré section sampled at the drive period
    Δt_c = 2π / CHAOTIC.Ω_D / 2000
    n_c = 2000 * 3000
    tc, θc, ωc = integrate(rk4_step, 0.2, 0.0, Δt_c, n_c, CHAOTIC, accel)
    stride = 2000
    keep = (300 * stride):stride:n_c          # discard the transient
    @printf("Poincaré section: %d points at the drive period\n", length(keep))

    fig = Figure(size = (1040, 460))

    ax1 = Axis(
        fig[2, 1], xlabel = L"Step size $\Delta t$ [s]", ylabel = "Global error [rad]",
        xscale = log10, yscale = log10,
        xticks = ([0.00125, 0.0025, 0.005, 0.01, 0.02],
            [L"0.00125", L"0.0025", L"0.005", L"0.01", L"0.02"],),
        yticks = logticks(-13, -4; step = 3),)
    ax1.xticklabelrotation = π/4
    # Each panel carries its own legend. One figure-level legend gave blue two
    # meanings -- the correct RK4 of the order study and the full pendulum of
    # the comparison -- and orange two, while leaving the Poincaré section with
    # no entry at all.
    scatterlines!(ax1, hs, e_correct, color = PALETTE.blue,
        markersize = MARKERSIZE.data,
        label = latexstring(@sprintf("\\text{RK4, order } %.2f", p_correct)),)
    scatterlines!(ax1, hs, e_original, color = PALETTE.red,
        markersize = MARKERSIZE.data, marker = :rect,
        label = latexstring(@sprintf("\\text{Original stages, order } %.2f", p_original)),)
    axislegend(ax1, position = :lt, framevisible = false, labelsize = 15, padding = 2)

    ax2 = Axis(fig[2, 2], xlabel = L"Time $t$ [s]", ylabel = L"Angle $\theta$ [rad]")
    lines!(ax2, t, θ_full, color = PALETTE.green, linewidth = 1.5,
        label = "Full pendulum",)
    # the takeaway as a legend entry: the panel is too narrow to carry it as a
    # separate annotation without running into the legend
    lines!(ax2, t, θ_lin, color = PALETTE.purple, linewidth = 1.5, linestyle = :dash,
        label = @sprintf("Small-angle, up to %.1f rad apart",
            maximum(abs.(θ_full .- θ_lin))))
    xlims!(ax2, 0, 12)
    ylims!(ax2, -2.9, 4.1)
    axislegend(ax2, position = :lt, framevisible = false, labelsize = 15, padding = 2)

    ax3 = Axis(fig[2, 3], xlabel = L"\theta \ \mathrm{[rad]}",
        ylabel = L"\omega \ \mathrm{[rad\ s^{-1}]}",
        xticks = ([-π, -π/2, 0, π/2, π],
            [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"],),)
    scatter!(ax3, wrap.(θc[keep]), ωc[keep], color = (PALETTE.orange, 0.6),
        markersize = MARKERSIZE.cloud,)
    # built by hand: a legend entry taken from the plot would inherit the small
    # marker the section needs and be unreadable
    axislegend(ax3,
        [MarkerElement(color = PALETTE.orange, marker = :circle,
            markersize = MARKERSIZE.key,)],
        ["Poincaré section"],
        position = :lt, framevisible = false, labelsize = 15, padding = 2,)

    rowsize!(fig.layout, 2, Relative(0.92))
    path = savefigure(fig, FIGURES, "pendulum_integrators")
    println("wrote ", path)
    println("wrote ", animate_pendulum())
end

"""
    animate_pendulum()

Swing the pendulum under the correct RK4 and under the original stage coupling,
side by side, at a step size where the difference in order is visible.

The order study measures the defect; this shows it. Both schemes are given the
same coarse step, and the one whose θ-stages were fed the ω-slopes drifts in
phase within a few swings — a first-order error accumulating where a
fourth-order one would not.
"""
function animate_pendulum()
    free = (g_over_L = 9.8, q = 0.0, F_D = 0.0, Ω_D = 0.0)
    Δt = 0.05                 # coarse on purpose: at 0.002 the two agree
    n = round(Int, 12 / Δt)
    t, θc, _ = integrate(rk4_step, 2.5, 0.0, Δt, n, free, accel)
    _, θo, _ = integrate(rk4_step_original, 2.5, 0.0, Δt, n, free, accel)

    # Initialised from the first frame: a colour vector cannot be matched
    # against an empty position vector at construction.
    p1 = Point2f(sin(θc[1]), -cos(θc[1]))
    p2 = Point2f(sin(θo[1]), -cos(θo[1]))
    bobs = Observable([p1, p2])
    rods = Observable([Point2f(0, 0), p1, Point2f(0, 0), p2])
    trace_t = Observable([t[1]])
    trace_c = Observable([θc[1]])
    trace_ot = Observable([t[1]])
    trace_o = Observable([wrap(θo[1])])
    caption = Observable{Any}("")   # the frame captions are rich text, not String

    fig = Figure(size = (940, 460))
    axp = Axis(fig[2, 1], aspect = DataAspect())
    hidedecorations!(axp)
    hidespines!(axp)
    # one colour per point: four points make the two rods
    linesegments!(
        axp,
        rods,
        color = [PALETTE.blue, PALETTE.blue, PALETTE.red, PALETTE.red],
        linewidth = 2.5,
    )
    scatter!(axp, bobs, color = [PALETTE.blue, PALETTE.red], markersize = 32)
    scatter!(axp, [Point2f(0, 0)], color = :black, markersize = MARKERSIZE.dense)
    # The bob sits at (sin θ, -cos θ), so it rises above y = 0 for |θ| > π/2 and
    # the swing starts at 2.5 rad: a top limit of 0.25 clipped both bobs over
    # most of every swing, exactly where the two pendulums differ most.
    xlims!(axp, -1.35, 1.35)
    ylims!(axp, -1.35, 1.35)

    axt = Axis(fig[2, 2], xlabel = L"Time $t$ [s]", ylabel = L"Angle $\theta$ [rad]")
    lc = lines!(axt, trace_t, trace_c, color = PALETTE.blue, linewidth = 2)
    lo = lines!(axt, trace_ot, trace_o, color = PALETTE.red, linewidth = 2)
    xlims!(axt, 0, last(t))
    # The original coupling does not oscillate, it winds: θ runs past -70 rad by
    # the end. Wrapping it to (-π, π] keeps it in the panel, so the comparison
    # the animation exists to make is visible for the whole twelve seconds
    # instead of the first 1.7.
    ylims!(axt, -3.4, 3.4)

    Legend(fig[1, 1:2], [lc, lo], ["RK4, correct stages", "Original stage coupling"],
        orientation = :horizontal, framevisible = false, labelsize = 16,)
    Label(fig[3, 1:2], caption, fontsize = 16, tellwidth = false)
    colsize!(fig.layout, 1, Relative(0.38))
    rowgap!(fig.layout, 6)

    path = joinpath(FIGURES, "pendulum_integrators.gif")
    mkpath(FIGURES)
    record(fig, path, 1:2:n; framerate = 20) do k
        pc = Point2f(sin(θc[k]), -cos(θc[k]))
        po = Point2f(sin(θo[k]), -cos(θo[k]))
        bobs[] = [pc, po]
        rods[] = [Point2f(0, 0), pc, Point2f(0, 0), po]
        trace_t[] = t[1:k]
        trace_c[] = θc[1:k]
        ot, ow = wrapped_trace(t[1:k], θo[1:k])
        trace_ot[] = ot
        trace_o[] = ow
        caption[] = rich(it("t"), @sprintf(" = %.2f s,  ", t[k]), "Δ", it("t"),
            @sprintf(" = %.2f s      ", Δt), "phase error ",
            replace(@sprintf("%+.2f", θo[k] - θc[k]), "-" => "−"), " rad",)
    end
    return path
end

main()
