# Two initial-value problems from the fourth-year numerical-methods exam, and a
# well-posedness trap in the first of them.
#
#   (a)  y' = √|sin y|          (b)  y' = z,  z' = -y sin y
#
# Ported from ODE_RK4.jl and ODE_system_RK4.jl.
#
# **(a) is not well posed as its own comment states it.** The file is headed
# `y'(x) = sqrt(abs(sin(y(x)))) & y(0) = 0`, but √|sin y| is not Lipschitz at
# y = 0 — its derivative diverges there — so Picard–Lindelöf does not apply and
# uniqueness fails. y ≡ 0 is a solution, and so is a family of solutions that
# leave the origin after an arbitrary delay. A numerical method cannot choose
# between them: RK4 started exactly at 0 stays at 0 forever, while starting at
# any ε > 0 climbs away. The 2021 code quietly set y[1] = 1, contradicting its
# own comment and stepping around the problem. Both behaviours are shown.
#
# (b) was integrated correctly but called RungeKutta4 twice per step, once for
# each component, doubling the work and discarding half of each result.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Right-hand side of the scalar problem. Not Lipschitz where sin y = 0."
f_scalar(y) = sqrt(abs(sin(y)))

"Right-hand side of the second-order problem, as a first-order system."
f_system(y, z) = (z, -sin(y) * y)

function rk4_scalar(f, y₀, Δx, n)
    y = Vector{Float64}(undef, n + 1)
    y[1] = y₀
    for i in 1:n
        k₁ = f(y[i]);              k₂ = f(y[i] + Δx*k₁/2)
        k₃ = f(y[i] + Δx*k₂/2);    k₄ = f(y[i] + Δx*k₃)
        y[i+1] = y[i] + Δx/6*(k₁ + 2k₂ + 2k₃ + k₄)
    end
    return y
end

"One RK4 step of the coupled system — computed once, not twice as in 2021."
function rk4_system(f, y₀, z₀, Δx, n)
    y = Vector{Float64}(undef, n + 1)
    z = Vector{Float64}(undef, n + 1)
    y[1], z[1] = y₀, z₀
    for i in 1:n
        k1y, k1z = f(y[i], z[i])
        k2y, k2z = f(y[i] + Δx*k1y/2, z[i] + Δx*k1z/2)
        k3y, k3z = f(y[i] + Δx*k2y/2, z[i] + Δx*k2z/2)
        k4y, k4z = f(y[i] + Δx*k3y,   z[i] + Δx*k3z)
        y[i+1] = y[i] + Δx/6*(k1y + 2k2y + 2k3y + k4y)
        z[i+1] = z[i] + Δx/6*(k1z + 2k2z + 2k3z + k4z)
    end
    return y, z
end

function main()
    Δx, n = 0.02, 5000
    x = range(0, step = Δx, length = n + 1)

    starts = (0.0, 1e-12, 1e-6, 1.0)
    solutions = [rk4_scalar(f_scalar, y₀, Δx, n) for y₀ in starts]
    for (y₀, y) in zip(starts, solutions)
        @printf("y(0) = %-8.0e  ->  y(%.0f) = %.6f\n", y₀, last(x), last(y))
    end
    println("all four are legitimate solutions of the same IVP: uniqueness fails at y = 0")

    y, z = rk4_system(f_system, 1.0, 0.0, Δx, n)
    # the system is conservative: E = z²/2 + ∫₀^y u sin u du = z²/2 + sin y - y cos y
    energy(y, z) = z^2/2 + sin(y) - y*cos(y)
    E = energy.(y, z)
    @printf("system invariant z²/2 + sin y - y cos y: drift %.3e over %d steps\n",
            maximum(abs.(E .- E[1])), n)

    fig = Figure(size = (1000, 440))

    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"y(x)")
    handles = []
    # The 10^-12 and 10^-6 solutions agree to four digits at x = 100 and one
    # would simply erase the other at equal width; broad solid under narrow
    # dashed keeps both visible, which is the point -- a perturbation twelve
    # orders down lands in the same place.
    widths = (1.8, 4.5, 1.8, 1.8)
    styles = (:solid, :solid, :dash, :solid)
    for (k, (y₀, ysol)) in enumerate(zip(starts, solutions))
        col = (PALETTE.black, PALETTE.blue, PALETTE.orange, PALETTE.green)[k]
        push!(handles, lines!(ax1, x, ysol, color = col,
              linewidth = widths[k], linestyle = styles[k]))
    end
    xlims!(ax1, 0, 100)
    text!(ax1, 0.03, 0.96;
        text = latexstring(@sprintf("\\text{A perturbation of } 10^{-12} \\text{ moves } \
                                     y(100) \\text{ from } 0 \\text{ to } %.1f", 
                                    last(solutions[2]))),
        space = :relative, align = (:left, :top), fontsize = 15)

    ax2 = Axis(fig[2, 2], xlabel = L"y", ylabel = L"z = y'")
    lines!(ax2, y, z, color = PALETTE.purple, linewidth = 1.2)
    text!(ax2, 0.97, 0.06;
        text = rich("Invariant drift ", rsci(maximum(abs.(E .- E[1])); digits = 1),
                    @sprintf(" over %d RK4 steps", n)),
        space = :relative, align = (:right, :bottom), fontsize = 15)

    Legend(fig[1, 1:2], handles,
        [L"y(0) = 0", L"y(0) = 10^{-12}", L"y(0) = 10^{-6}", L"y(0) = 1"],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 22)

    rowsize!(fig.layout, 2, Relative(0.85))
    path = savefigure(fig, FIGURES, "nonunique_ivp_and_rk4")
    println("wrote ", path)
end

main()
