# The pendulum, linearised and full, and what the RK4 of the original files
# integrated:
#
#   θ̈ = −(g/L) sin θ − q θ̇ + F_D sin(Ω_D t)
#
# Ported from Pendul_simplu.m, Pendul_dampat.m, Pendul_dampat_fortat.m and
# Runge_Kutta_4_pendul_dampat_fortat.m on the `legacy` branch.
#
# The defect that matters is in the RK4 file. For the coupled system θ̇ = ω,
# ω̇ = F(θ, ω, t), the θ-stages must be fed the θ-slopes. The original wrote
#
#     k2 = frhs(omega + dt*k1/2, theta + dt*k1/2, t + dt/2)
#     pk1 = omega(step);  pk2 = omega(step) + dt*pk1/2;  ...
#
# so k1, which is dω/dt, advances both components, and the θ-stages propagate
# ω as though dω/dt = ω. The θ update expands to dt·ω + dt²·ω/2 where RK4 needs
# dt·ω + dt²·F/2: the scheme is first order. Both couplings are integrated
# below against a reference in 256-bit arithmetic and their orders measured.
#
# Two lesser points: Pendul_simplu.m is headed "pendulul matematic" but
# integrates −(g/L)θ, the small-angle approximation, not −(g/L) sin θ; and all
# four files assign `length = 1` or `length = 9.8`, shadowing the Octave builtin.
#
# The chaotic parameter set is the one of Giordano and Nakanishi, Computational
# Physics, 2nd ed., Pearson (2006), §3.3.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
using SpecialFunctions: ellipk
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "numerics_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Chaotic parameter set of the driven damped pendulum, g/L in s⁻², q in s⁻¹, Ω_D in rad s⁻¹."
const CHAOTIC = (g_over_L = 1.0, q = 0.5, F_D = 1.2, Ω_D = 2 / 3)
"A weakly driven, non-chaotic set for the order study."
const MILD = (g_over_L = 1.0, q = 0.5, F_D = 0.2, Ω_D = 2 / 3)
"Undamped, undriven pendulum with the L = 1 m, g = 9.8 m s⁻² of Pendul_simplu.m."
const FREE = (g_over_L = 9.8, q = 0.0, F_D = 0.0, Ω_D = 0.0)
"Initial angle of the free-pendulum comparison, in radians."
const θ₀_FREE = 2.5
"Horizon of the order study, in seconds; a power of two times 0.01 so every step divides it."
const T_ORDER = 5.12
"Step counts of the order study, from 0.32 s down to 0.0025 s."
const N_ORDER = [2^k for k in 4:11]
"Tolerance on the observed orders of the two couplings."
const ORDER_TOLERANCE = 0.05
"Relative tolerance on the period of the free pendulum against 4K(k²)/ω₀."
const PERIOD_TOLERANCE = 1e-6

"Angular acceleration of the full pendulum."
accel(θ, ω, t, p) = -p.g_over_L * sin(θ) - p.q * ω + p.F_D * sin(p.Ω_D * t)
"Angular acceleration under the small-angle approximation."
accel_linear(θ, ω, t, p) = -p.g_over_L * θ - p.q * ω + p.F_D * sin(p.Ω_D * t)

"""
    correct_step(t, u, h, p, a)

One classical RK4 step of the pendulum written as the system `(θ, ω)′ =
(ω, a(θ, ω, t, p))`, the θ-stages fed the ω-values.
"""
correct_step(t, u, h, p, a) = rk4_step((s, w) -> (w[2], a(w[1], w[2], s, p)), t, u, h)

"""
    original_step(t, u, h, p, a)

The stage coupling as written in Runge_Kutta_4_pendul_dampat_fortat.m, kept so
that its order can be measured: `k1` is dω/dt but advances both components, and
the θ-stages propagate ω as though dω/dt = ω.
"""
function original_step(t, u, h, p, a)
    θ, ω = u
    k1 = a(θ, ω, t, p)
    k2 = a(θ + h * k1 / 2, ω + h * k1 / 2, t + h / 2, p)
    k3 = a(θ + h * k2 / 2, ω + h * k2 / 2, t + h / 2, p)
    k4 = a(θ + h * k3, ω + h * k3, t + h, p)
    ω_new = ω + h / 6 * (k1 + 2k2 + 2k3 + k4)

    pk1 = ω
    pk2 = ω + h * pk1 / 2
    pk3 = ω + h * pk2 / 2
    pk4 = ω + h * pk3
    return (θ + h / 6 * (pk1 + 2pk2 + 2pk3 + pk4), ω_new)
end

"""
    integrate(step, u₀, h, n, p, a)

Apply `step` `n` times from `u₀ = (θ₀, ω₀)`. Returns the time, angle and
angular-velocity traces.
"""
function integrate(step, u₀, h, n, p, a)
    t = Vector{Float64}(undef, n + 1)
    θ = Vector{Float64}(undef, n + 1)
    ω = Vector{Float64}(undef, n + 1)
    t[1], (θ[1], ω[1]) = 0.0, u₀
    for i in 1:n
        θ[i + 1], ω[i + 1] = step(t[i], (θ[i], ω[i]), h, p, a)
        t[i + 1] = i * h
    end
    return t, θ, ω
end

"""
    endpoint(step, u₀, h, n, p, a)

The angle after `n` steps, without storing the trace; `u₀`, `h` and `p` may be
`BigFloat` for a reference solution.
"""
function endpoint(step, u₀, h, n, p, a)
    u = u₀
    for i in 1:n
        u = step((i - 1) * h, u, h, p, a)
    end
    return u[1]
end

"""
    reference_angle(p, t_end, n)

θ(t_end) by RK4 in 256-bit arithmetic at `n` steps, accepted only if it agrees
with the `n ÷ 2` result to 1e-17, far below double precision.
"""
function reference_angle(p, t_end, n)
    setprecision(BigFloat, 256) do
        pb = map(big, p)
        fine = endpoint(correct_step, (big(0.2), big(0.0)), big(t_end) / n, n, pb, accel)
        coarse = endpoint(
            correct_step, (big(0.2), big(0.0)), 2big(t_end) / n, n ÷ 2, pb, accel,)
        abs(fine - coarse) < big"1e-17" ||
            error("reference not converged: |Δθ| = $(Float64(abs(fine - coarse)))")
        Float64(fine)
    end
end

"""
    period(t, θ, ω)

Period of a free oscillation started at rest: the time of the first return of
ω to zero from above, by linear interpolation between the bracketing samples.
"""
function period(t, θ, ω)
    i = findfirst(k -> ω[k] > 0 && ω[k + 1] <= 0, 2:(length(ω) - 1))
    i === nothing && throw(ArgumentError("no full oscillation in the trace"))
    k = i + 1
    return t[k] + ω[k] / (ω[k] - ω[k + 1]) * (t[k + 1] - t[k])
end

"Wrap an angle to (−π, π]."
wrap(θ) = mod(θ + π, 2π) - π

"""
    wrapped_trace(t, θ)

`θ` wrapped to (−π, π] against `t`, with a `NaN` pair inserted at every wrap so
that a line breaks there instead of striking across the panel. For the original
coupling, which winds rather than oscillates.
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

function main()
    # --- order of the two stage couplings -----------------------------------
    reference = reference_angle(MILD, T_ORDER, 2^16)
    hs = T_ORDER ./ N_ORDER
    e_correct = [abs(endpoint(correct_step, (0.2, 0.0), h, n, MILD, accel) - reference)
                 for (h, n) in zip(hs, N_ORDER)]
    e_original = [abs(endpoint(original_step, (0.2, 0.0), h, n, MILD, accel) - reference)
                  for (h, n) in zip(hs, N_ORDER)]
    o_correct = convergence_order(hs, e_correct)
    o_original = convergence_order(hs, e_original)
    println("order study, θ(t = $T_ORDER) against RK4 in 256-bit arithmetic")
    println("      h      RK4 error   order   original error  order")
    for i in eachindex(hs)
        @printf("%8.5f   %.3e   %5s     %.3e     %5s\n", hs[i], e_correct[i],
            i == 1 ? "" : @sprintf("%.3f", o_correct.pairwise[i - 1]), e_original[i],
            i == 1 ? "" : @sprintf("%.3f", o_original.pairwise[i - 1]))
    end
    p_correct, p_original = last(o_correct.pairwise), last(o_original.pairwise)
    @printf("order between the two finest steps: RK4 %.3f, original coupling %.3f\n",
        p_correct, p_original)
    abs(p_correct - 4) < ORDER_TOLERANCE || error("RK4 is not fourth order: $p_correct")
    abs(p_original - 1) < ORDER_TOLERANCE ||
        error("the original coupling is not first order: $p_original")

    # --- free pendulum, full against small-angle ----------------------------
    Δt = 0.002
    n = round(Int, 12 / Δt)
    t, θ_full, ω_full = integrate(correct_step, (θ₀_FREE, 0.0), Δt, n, FREE, accel)
    _, θ_lin, _ = integrate(correct_step, (θ₀_FREE, 0.0), Δt, n, FREE, accel_linear)
    gap = maximum(abs.(θ_full .- θ_lin))
    ω₀ = sqrt(FREE.g_over_L)
    T_exact = 4 * ellipk(sin(θ₀_FREE / 2)^2) / ω₀
    T_linear = 2π / ω₀
    T_measured = period(t, θ_full, ω_full)
    @printf("\nfree pendulum from θ₀ = %.1f rad: period %.6f s, exact 4K(k²)/ω₀ = %.6f s, small-angle 2π/ω₀ = %.4f s\n",
        θ₀_FREE, T_measured, T_exact, T_linear)
    @printf("full and small-angle solutions differ by up to %.3f rad over 12 s\n", gap)
    isapprox(T_measured, T_exact; rtol = PERIOD_TOLERANCE) ||
        error("measured period $T_measured is off the closed form $T_exact")

    # --- Poincaré section of the chaotic set --------------------------------
    stride = 2000
    periods = 3000
    Δt_c = 2π / CHAOTIC.Ω_D / stride
    _, θc, ωc = integrate(correct_step, (0.2, 0.0), Δt_c, stride * periods, CHAOTIC, accel)
    keep = (300 * stride):stride:(stride * periods)   # discard the transient
    @printf("Poincaré section: %d points, one per drive period after the first 300\n",
        length(keep))

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1600, 680))

    ax1 = Axis(
        fig[2, 1], xlabel = L"Step size $\Delta t$ [s]", ylabel = "Global error [rad]",
        xscale = log10, yscale = log10,
        xticks = ([0.0025, 0.01, 0.04, 0.16], [L"0.0025", L"0.01", L"0.04", L"0.16"]),
        yticks = logticks(-13, -1; step = 3),)
    guide_4 = 4 * last(e_correct) .* (hs ./ last(hs)) .^ 4
    guide_1 = 4 * last(e_original) .* (hs ./ last(hs))
    lines!(
        ax1, hs, guide_4, color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    lines!(
        ax1, hs, guide_1, color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    text!(
        ax1, hs[3], guide_4[3]; text = L"\propto \Delta t^{4}", align = (:right, :bottom),
        offset = (-6, 4), fontsize = ANNOTATION_SIZE,)
    text!(ax1, hs[3], guide_1[3]; text = L"\propto \Delta t", align = (:right, :bottom),
        offset = (-6, 4), fontsize = ANNOTATION_SIZE,)
    l_correct = scatterlines!(ax1, hs, e_correct, color = PALETTE.blue)
    l_original = scatterlines!(ax1, hs, e_original, color = PALETTE.red, marker = :rect)
    text!(ax1, 0.97, 0.04;
        text = rich(
            rich(@sprintf("RK4: order %.2f", p_correct), color = PALETTE.blue), "\n",
            rich(@sprintf("Original stages: order %.2f", p_original), color = PALETTE.red),),
        space = :relative, align = (:right, :bottom), justification = :right,
        fontsize = ANNOTATION_SIZE,)
    xlims!(ax1, 0.0017, 0.45)
    ylims!(ax1, 3e-14, 1)

    ax2 = Axis(fig[2, 2], xlabel = L"Time $t$ [s]", ylabel = L"Angle $\theta$ [rad]",
        xticks = 0:2:12,)
    l_full = lines!(ax2, t, θ_full, color = PALETTE.green)
    l_lin = lines!(ax2, t, θ_lin, color = PALETTE.purple, linestyle = :dash)
    text!(ax2, 0.03, 0.97;
        text = rich(rich(@sprintf("Period %.2f s", T_exact), color = PALETTE.green), "\n",
            rich(@sprintf("Period %.2f s", T_linear), color = PALETTE.purple), "\n",
            @sprintf("Up to %.1f rad apart", gap),),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)
    xlims!(ax2, 0, 12)
    ylims!(ax2, -3.0, 5.2)

    ax3 = Axis(fig[2, 3], xlabel = L"Angle $\theta$ [rad]",
        ylabel = L"Angular velocity $\omega$ [rad s$^{-1}$]",
        xticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]),)
    scatter!(ax3, wrap.(θc[keep]), ωc[keep], color = (PALETTE.orange, 0.6),
        markersize = MARKERSIZE.cloud, strokewidth = 0,)
    text!(ax3, 0.03, 0.04; text = @sprintf("%d drive periods", length(keep)),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.orange,)
    xlims!(ax3, -π, π)

    # the section's own marker is too small to read in a key, so its entry is
    # built by hand
    section_key = MarkerElement(color = PALETTE.orange, marker = :circle,
        markersize = MARKERSIZE.key,)
    Legend(fig[1, 1:3],
        [[l_correct, l_original], [l_full, l_lin], [section_key]],
        [["RK4", "Original stages"], ["Full", "Small-angle"], ["Poincaré section"]],
        ["Order study", "Free pendulum", "Chaotic drive"];
        titleposition = :left, nbanks = 1, groupgap = 36,)

    path = savefigure(fig, FIGURES, "pendulum_integrators")
    println("wrote ", path)
    println("wrote ", animate_pendulum())
end

"""
    animate_pendulum()

Swing the free pendulum under the correct RK4 and under the original stage
coupling side by side at a step coarse enough for the difference in order to
show: the first-order scheme drifts in phase within a few swings.
"""
function animate_pendulum()
    Δt = 0.05                 # coarse on purpose: at 0.002 the two agree
    n = round(Int, 12 / Δt)
    t, θc, _ = integrate(correct_step, (θ₀_FREE, 0.0), Δt, n, FREE, accel)
    _, θo, _ = integrate(original_step, (θ₀_FREE, 0.0), Δt, n, FREE, accel)

    # initialised from the first frame: a colour vector cannot be matched
    # against an empty position vector at construction
    p1 = Point2f(sin(θc[1]), -cos(θc[1]))
    p2 = Point2f(sin(θo[1]), -cos(θo[1]))
    bobs = Observable([p1, p2])
    rods = Observable([Point2f(0, 0), p1, Point2f(0, 0), p2])
    trace_t = Observable([t[1]])
    trace_c = Observable([θc[1]])
    trace_ot = Observable([t[1]])
    trace_o = Observable([wrap(θo[1])])
    frame_caption(k) = rich(it("t"), @sprintf(" = %.2f s,  Δ", t[k]), it("t"),
        @sprintf(" = %.2f s,  phase error ", Δt),
        replace(@sprintf("%+.2f", θo[k] - θc[k]), "-" => "−"), " rad",)
    caption = Observable(frame_caption(1))

    fig = Figure(size = (1200, 560))
    axp = Axis(fig[2, 1], aspect = DataAspect())
    hidedecorations!(axp)
    hidespines!(axp)
    # one colour per point: four points make the two rods
    linesegments!(axp, rods, color = [PALETTE.blue, PALETTE.blue, PALETTE.red, PALETTE.red])
    scatter!(axp, bobs, color = [PALETTE.blue, PALETTE.red], markersize = 32)
    scatter!(axp, [Point2f(0, 0)], color = :black, markersize = MARKERSIZE.dense)
    # the bob sits at (sin θ, −cos θ) and rises above y = 0 for |θ| > π/2
    xlims!(axp, -1.35, 1.35)
    ylims!(axp, -1.35, 1.35)

    axt = Axis(fig[2, 2], xlabel = L"Time $t$ [s]", ylabel = L"Angle $\theta$ [rad]",
        xticks = 0:2:12,)
    lc = lines!(axt, trace_t, trace_c, color = PALETTE.blue)
    lo = lines!(axt, trace_ot, trace_o, color = PALETTE.red, linestyle = :dash)
    xlims!(axt, 0, last(t))
    # the original coupling winds rather than oscillates, past −70 rad by the
    # end; wrapped to (−π, π] it stays in the panel for the whole twelve seconds
    ylims!(axt, -3.4, 3.4)

    Legend(fig[1, 1:2], [lc, lo], ["RK4, correct stages", "Original stage coupling"])
    Label(fig[3, 1:2], caption, fontsize = 22, tellwidth = false)
    colsize!(fig.layout, 1, Relative(0.36))
    rowgap!(fig.layout, 6)

    path = joinpath(FIGURES, "pendulum_integrators.gif")
    mkpath(FIGURES)
    record(fig, path, 1:2:n; framerate = 20, px_per_unit = 1) do k
        pc = Point2f(sin(θc[k]), -cos(θc[k]))
        po = Point2f(sin(θo[k]), -cos(θo[k]))
        bobs[] = [pc, po]
        rods[] = [Point2f(0, 0), pc, Point2f(0, 0), po]
        trace_t[] = t[1:k]
        trace_c[] = θc[1:k]
        ot, ow = wrapped_trace(t[1:k], θo[1:k])
        trace_ot[] = ot
        trace_o[] = ow
        caption[] = frame_caption(k)
    end
    return path
end

main()
