# The focusing nonlinear Schrödinger equation by the method of lines,
#
#   i ∂ₜΨ = −½ ∂ₓₓΨ − |Ψ|²Ψ,     periodic in x,
#
# solved as ∂ₜΨ = i(½ ∂ₓₓΨ + |Ψ|²Ψ) with a second-order centred stencil in space
# and RK4 in time, on the two-soliton initial condition of SolitonicEq_MOL.jl in
# Julia-Workflow-FFUB/Examen_PDF_MN_II_L_4/ on the `legacy` branch.
#
# A soliton η sech(η(x − x₀)) e^{ikx} of this equation travels with velocity
# v = k, so the original's pair — centred at x = ∓5 with k = ∓0.1 — moves apart
# and meets across the periodic boundary. `main` measures the speed of the
# left-hand peak before the encounter and asserts v = −0.1, asserts the norm
# ∫|Ψ|² dx conserved to round-off, and reports the drift of the Hamiltonian
# ∫(½|∂ₓΨ|² − ½|Ψ|⁴) dx, which RK4 does not conserve exactly.
#
# The periodic grid holds one copy of each point: x = −a, …, a − dx. The
# original's time stepping perturbed all three stencil neighbours by the same
# increment, the one belonging to the centre point,
#
#     k2 = Sⁿ(Ψⱼ₋₁ + dt*k1/2, Ψⱼ + dt*k1/2, Ψⱼ₊₁ + dt*k1/2, dx)
#
# where a method-of-lines stage must be formed on the whole solution vector
# before the spatial operator is applied; it was a three-stage scheme named
# RungeKutta3; and its dx = 0.5 put about one grid point across a soliton of
# width 0.5.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"""
    rhs!(dΨ, Ψ, dx)

∂ₜΨ = i(½ ∂ₓₓΨ + |Ψ|²Ψ) with periodic boundaries, applied to the whole vector.
"""
function rhs!(dΨ, Ψ, dx)
    n = length(Ψ)
    @inbounds for j in 1:n
        jm = j == 1 ? n : j - 1
        jp = j == n ? 1 : j + 1
        laplacian = (Ψ[jp] - 2Ψ[j] + Ψ[jm]) / dx^2
        dΨ[j] = im * (0.5 * laplacian + abs2(Ψ[j]) * Ψ[j])
    end
    return dΨ
end

"""
    step_rk4!(Ψ, work, dx, dt)

One RK4 step of the semi-discrete system. Stages are formed on the full vector,
which is the point the 2021 version missed.
"""
function step_rk4!(Ψ, work, dx, dt)
    k1, k2, k3, k4, tmp = work
    rhs!(k1, Ψ, dx)
    @. tmp = Ψ + dt/2 * k1
    rhs!(k2, tmp, dx)
    @. tmp = Ψ + dt/2 * k2
    rhs!(k3, tmp, dx)
    @. tmp = Ψ + dt * k3
    rhs!(k4, tmp, dx)
    @. Ψ += dt/6 * (k1 + 2k2 + 2k3 + k4)
    return Ψ
end

function main()
    # domain and initial condition as in the original; the resolution is not
    a, dx = 10.0, 0.02
    n = round(Int, 2a / dx)
    x = collect(range(-a, a - dx, length = n))   # one copy of each periodic point
    dt = 0.2 * dx^2                      # within the stability limit of the centred stencil
    t_end = 50.0
    n_steps = round(Int, t_end / dt)

    # the original pair: centred at ∓5, phases ∓0.1, so they travel apart
    Ψ = @. 2 * exp(-im * 0.1 * x) / cosh(2 * (x + 5)) +
           2 * exp(+im * 0.1 * x) / cosh(2 * (x - 5))
    work = ntuple(_ -> similar(Ψ), 5)

    hamiltonian(Ψ) = sum(0.5 * abs2(Ψ[j == n ? 1 : j + 1] - Ψ[j]) / dx^2 -
                         0.5 * abs2(Ψ[j])^2
    for j in 1:n) * dx
    norm0, H0 = sum(abs2, Ψ) * dx, hamiltonian(Ψ)
    n_save = 400
    save_every = max(n_steps ÷ n_save, 1)
    frames = Vector{Vector{Float64}}()
    times = Float64[]
    norms = Float64[]

    for s in 1:n_steps
        step_rk4!(Ψ, work, dx, dt)
        if s % save_every == 0
            push!(frames, abs2.(Ψ))
            push!(times, s * dt)
            push!(norms, sum(abs2, Ψ) * dx)
        end
    end

    @printf("grid %d points, dt = %.2e, %d steps to t = %.1f\n", n, dt, n_steps, t_end)
    norm_drift = maximum(abs.(norms .- norm0)) / norm0
    @printf("norm ∫|Ψ|²dx: initial %.6f, largest relative drift %.1e; Hamiltonian drift %.1e\n",
        norm0, norm_drift, abs(hamiltonian(Ψ) - H0) / abs(H0))
    norm_drift < 1e-12 || error("the norm drifts by $norm_drift")

    # speed of the left-hand soliton before the encounter, from the peak of |Ψ|²
    # on x < 0 located by parabolic interpolation, fitted over 2 ≤ t ≤ 20
    function peak_position(ρ_t)
        left = findall(<(0), x)
        k = left[argmax(ρ_t[left])]
        km, kp = mod1(k - 1, n), mod1(k + 1, n)
        denom = ρ_t[km] - 2ρ_t[k] + ρ_t[kp]
        return x[k] - dx / 2 * (ρ_t[kp] - ρ_t[km]) / denom
    end
    window = findall(t -> 2 <= t <= 20, times)
    tw, pw = times[window], peak_position.(frames[window])
    v = sum((tw .- sum(tw) / length(tw)) .* (pw .- sum(pw) / length(pw))) /
        sum(abs2, tw .- sum(tw) / length(tw))
    @printf("left soliton: peak speed %.4f over 2 ≤ t ≤ 20, against v = k = −0.1\n", v)
    isapprox(v, -0.1; rtol = 0.03) || error("the soliton speed $v is not −0.1")
    peak_t = times[argmax(maximum.(frames))]
    @printf("the two peaks coincide at t = %.1f with |Ψ|² = %.2f\n", peak_t,
        maximum(maximum.(frames)))

    ρ = reduce(hcat, frames)

    fig = Figure(size = (1500, 640))
    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"Time $t$")
    # the collision reaches |Ψ|² = 16 for an instant while the solitons sit at
    # 4; the square root keeps both visible, and the colour bar is ticked in
    # the original units
    heatmap!(ax1, x, times, sqrt.(ρ), colormap = :viridis)

    ax2 = Axis(fig[2, 2], xlabel = L"x", ylabel = L"|\Psi|^2")
    shown = (1, length(times) ÷ 2, length(times))
    profiles = [lines!(ax2, x, frames[idx], color = (
                    PALETTE.blue, PALETTE.orange, PALETTE.green,)[k])
                for (k, idx) in enumerate(shown)]
    ylims!(ax2, -0.2, 5.1)
    text!(ax2, 0.03, 0.96;
        text = rich("Left peak speed ", @sprintf("%.3f", v), " against −0.1\n",
            "Norm drift ", rsci(norm_drift; digits = 1),),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)

    ticks = [0, 1, 4, 9, 16]
    Colorbar(fig[2, 3], limits = (sqrt(minimum(ρ)), sqrt(maximum(ρ))),
        colormap = :viridis, label = L"|\Psi|^2",
        ticks = (sqrt.(ticks), [latexstring(string(t)) for t in ticks]),)
    colsize!(fig.layout, 3, Relative(0.03))
    Legend(fig[1, 1:3], profiles,
        [rich(it("t"), @sprintf(" = %.1f", times[idx])) for idx in shown],
        rich("Profile |", it("Ψ"), "|²"); titleposition = :left,)

    path = savefigure(fig, FIGURES, "nonlinear_schrodinger_mol")
    println("wrote ", path)
    println("wrote ", animate_solitons(x, frames, times, norms, norm0))
end

"""
    animate_solitons(x, frames, times, norms, norm0)

Animate the two solitons colliding, with the conserved norm tracked beneath.

The heatmap shows the collision as a static interference pattern; the animation
shows what actually happens — the two pulses pass through one another and come
out unchanged in shape and speed, which is what makes them solitons rather than
merely localised wave packets. The norm panel is the numerical check running
alongside: a method-of-lines scheme that distorted the pulses would not hold
the norm flat.
"""
function animate_solitons(x, frames, times, norms, norm0)
    every = max(length(frames) ÷ 120, 1)
    ks = 1:every:length(frames)

    profile = Observable(frames[1])
    trace_t = Observable([times[1]])
    trace_n = Observable([(norms[1] / norm0 - 1) * 1e15])
    caption = Observable(rich(it("t"), " = 0.0"))

    fig = Figure(size = (1200, 680))
    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"|\Psi|^2")
    lines!(ax1, x, profile, color = PALETTE.blue)
    xlims!(ax1, first(x), last(x))
    ylims!(ax1, 0, maximum(maximum, frames) * 1.08)

    # the drift is of order 1e-14: the trace is the departure from unity in
    # units of 1e-15, where the scheme's behaviour is visible
    drift = (norms ./ norm0 .- 1) .* 1e15
    span = 1.15 * maximum(abs, drift)
    ax2 = Axis(fig[3, 1], xlabel = L"Time $t$", ylabel = L"Norm drift [$10^{-15}$]")
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    lines!(ax2, trace_t, trace_n, color = PALETTE.green)
    text!(ax2, 0.99, 0.95; text = L"$\int|\Psi|^2\mathrm{d}x$ conserved to round-off",
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,)
    xlims!(ax2, 0, last(times))
    ylims!(ax2, -span, span)

    Label(fig[1, 1], caption, fontsize = 22, tellwidth = false)
    rowsize!(fig.layout, 2, Relative(0.62))
    rowgap!(fig.layout, 10)

    path = joinpath(FIGURES, "nonlinear_schrodinger_mol.gif")
    mkpath(FIGURES)
    record(fig, path, ks; framerate = 15, px_per_unit = 1) do k
        profile[] = frames[k]
        trace_t[] = times[1:k]
        trace_n[] = drift[1:k]
        caption[] = rich(it("t"), @sprintf(" = %.1f      ", times[k]), "norm drift ",
            rsci(norms[k] / norm0 - 1; digits = 2),)
    end
    return path
end

main()
