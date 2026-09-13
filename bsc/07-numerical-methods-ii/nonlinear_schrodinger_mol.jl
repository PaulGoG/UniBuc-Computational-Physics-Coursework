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

    norm₀ = sum(abs2, Ψ) * dx
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
            norm₀, last(norms), abs(last(norms) - norm₀) / norm₀)

    ρ = reduce(hcat, frames)

    fig = Figure(size = (980, 440))
    ax1 = Axis(fig[1, 1], xlabel = L"x", ylabel = L"Time $t$",
        title = L"|\Psi(x,t)|^2", titlesize = 17)
    heatmap!(ax1, x, times, ρ, colormap = :viridis)

    ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"|\Psi|^2")
    for (k, idx) in enumerate((1, length(times) ÷ 2, length(times)))
        lines!(ax2, x, frames[idx],
               color = (PALETTE.blue, PALETTE.orange, PALETTE.green)[k], linewidth = 1.5,
               label = @sprintf("t = %.1f", times[idx]))
    end
    axislegend(ax2, position = :rt, framevisible = false, labelsize = 15)

    Colorbar(fig[1, 3], limits = (minimum(ρ), maximum(ρ)), colormap = :viridis,
             label = L"|\Psi|^2")
    colsize!(fig.layout, 3, Relative(0.03))

    path = savefigure(fig, FIGURES, "nonlinear_schrodinger_mol")
    println("wrote ", path)
end

main()
