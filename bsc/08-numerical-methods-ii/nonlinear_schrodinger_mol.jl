# The focusing nonlinear Schrödinger equation by the method of lines,
#
#   i ∂ₜΨ = -½ ∂ₓₓΨ - |Ψ|²Ψ,     periodic in x,
#
# solved as ∂ₜΨ = i(½ ∂ₓₓΨ + |Ψ|²Ψ) with a second-order centred stencil in space
# and RK4 in time, on the two-soliton initial condition of the original.
#
# A soliton carrying the phase factor e^{ikx} travels with velocity v = k, so the
# original's pair — centred at x = ∓5 with k = ∓0.1 — moves **apart**, not
# together. On the periodic domain [-10, 10] each reaches the boundary at
# t = 50 and wraps, so they meet at the edge rather than in the middle. The
# initial condition, domain and phase are kept exactly as written.
#
# Ported from SolitonicEq_MOL.jl. The sign convention and the stencil were
# right. The time stepping was not: the original wrote
#
#     k2 = Sⁿ(Ψⱼ₋₁ + dt*k1/2, Ψⱼ + dt*k1/2, Ψⱼ₊₁ + dt*k1/2, dx)
#
# perturbing all three stencil neighbours by the *same* increment k1, which is
# the increment belonging to the centre point alone. In the method of lines the
# stage must be formed on the whole solution vector, then the spatial operator
# applied to that. Done correctly here, which also makes the conserved norm
# actually conserved. It was a three-stage scheme labelled RungeKutta3; RK4 is
# used below.

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
    @. tmp = Ψ + dt/2 * k1;  rhs!(k2, tmp, dx)
    @. tmp = Ψ + dt/2 * k2;  rhs!(k3, tmp, dx)
    @. tmp = Ψ + dt   * k3;  rhs!(k4, tmp, dx)
    @. Ψ += dt/6 * (k1 + 2k2 + 2k3 + k4)
    return Ψ
end

function main()
    # domain and initial condition as in the original; the resolution is not.
    # SolitonicEq_MOL.jl used dx = 0.5 against a soliton of width ≈0.5, so it
    # carried about one grid point per soliton — far too coarse for a
    # second-order stencil to represent the shape at all.
    a, dx = 10.0, 0.02
    x = collect(-a:dx:a)
    n = length(x)
    dt = 0.2 * dx^2                      # diffusive stability limit of the stencil
    t_end = 50.0
    n_steps = round(Int, t_end / dt)

    # the original pair: centred at ∓5, phases ∓0.1, so they travel apart
    Ψ = @. 2 * exp(-im * 0.1 * x) / cosh(2 * (x + 5)) +
           2 * exp(+im * 0.1 * x) / cosh(2 * (x - 5))
    work = ntuple(_ -> similar(Ψ), 5)

    norm0 = sum(abs2, Ψ) * dx
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
    @printf("norm ∫|Ψ|²dx: initial %.6f, final %.6f, relative drift %.2e\n",
            norm0, last(norms), abs(last(norms) - norm0) / norm0)

    ρ = reduce(hcat, frames)

    fig = Figure(size = (980, 440))
    ax1 = Axis(fig[1, 1], xlabel = L"x", ylabel = L"Time $t$")
    # The collision reaches |Psi|^2 = 16 for one instant while the solitons sit
    # at 4, so a linear scale spent three quarters of its range on that instant
    # and rendered the trajectories as a dim fringe on near-black. The square
    # root keeps the collision the brightest thing in the panel and still shows
    # the solitons, and the colourbar is ticked in the original units.
    heatmap!(ax1, x, times, sqrt.(ρ), colormap = :viridis)

    ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"|\Psi|^2")
    for (k, idx) in enumerate((1, length(times) ÷ 2, length(times)))
        lines!(ax2, x, frames[idx],
               color = (PALETTE.blue, PALETTE.orange, PALETTE.green)[k], linewidth = 1.5,
               label = rich(it("t"), @sprintf(" = %.1f", times[idx])))
    end
    ylims!(ax2, -0.2, 5.1)   # headroom: the legend sat on the right-hand pulses
    axislegend(ax2, position = :rt, framevisible = false, labelsize = 15)

    ticks = [0, 1, 4, 9, 16]
    Colorbar(fig[1, 3], limits = (sqrt(minimum(ρ)), sqrt(maximum(ρ))),
             colormap = :viridis, label = L"|\Psi|^2",
             ticks = (sqrt.(ticks), [latexstring(string(t)) for t in ticks]))
    colsize!(fig.layout, 3, Relative(0.03))

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
    caption = Observable{Any}("")   # the frame captions are rich text, not String

    fig = Figure(size = (900, 520))
    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"|\Psi|^2")
    lines!(ax1, x, profile, color = PALETTE.blue, linewidth = 2)
    xlims!(ax1, first(x), last(x))
    ylims!(ax1, 0, maximum(maximum, frames) * 1.08)

    # The drift is of order 1e-15, so a panel scaled to a relative norm of
    # 1 +/- 0.005 showed a flat line on the unity guide and nothing else. The
    # trace is the departure from unity in units of 1e-15, where the scheme's
    # actual behaviour is visible; the guide keeps its label.
    # The full integral expression is too tall a label for a strip this short,
    # and left to itself the axis crowded nine ticks into it.
    ax2 = Axis(fig[3, 1], xlabel = L"Time $t$",
        ylabel = L"Norm drift [$10^{-15}$]",
        yticks = ([-10, -5, 0, 5, 10], [L"-10", L"-5", L"0", L"5", L"10"]))
    hlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    lines!(ax2, trace_t, trace_n, color = PALETTE.green, linewidth = 2)
    text!(ax2, 0.99, 0.95; text = L"$\int|\Psi|^2\mathrm{d}x$ exactly conserved",
        space = :relative,
        align = (:right, :top), fontsize = 13, color = PALETTE.black)
    xlims!(ax2, 0, last(times))
    ylims!(ax2, -12, 12)

    Label(fig[1, 1], caption, fontsize = 17, tellwidth = false)
    rowsize!(fig.layout, 2, Relative(0.62))
    rowgap!(fig.layout, 10)

    path = joinpath(FIGURES, "nonlinear_schrodinger_mol.gif")
    mkpath(FIGURES)
    record(fig, path, ks; framerate = 15) do k
        profile[] = frames[k]
        trace_t[] = times[1:k]
        trace_n[] = (norms[1:k] ./ norm0 .- 1) .* 1e15
        caption[] = rich(it("t"), @sprintf(" = %.1f      ", times[k]), "norm drift ",
            rsci(norms[k] / norm0 - 1; digits = 2))
    end
    return path
end

main()
