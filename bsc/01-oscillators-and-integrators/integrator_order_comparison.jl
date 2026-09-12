# Order of convergence of forward Euler against classical Runge-Kutta 4, on a
# linear scalar initial-value problem whose exact solution is available in
# closed form.
#
#   y'(x) = 3e^{-x} - 0.4 y,   y(0) = y₀
#   y(x)  = (y₀ + 5) e^{-0.4x} - 5 e^{-x}
#
# Ported from RKtrial.cpp and DiffEqEuler.cpp (2018). Both tableaux in the
# original were correct, but neither program compared against the exact
# solution or swept the step size, so the O(h) against O(h⁴) claim the files
# were written to demonstrate was never actually shown. RKtrial.cpp printed
# only the endpoint value; DiffEqEuler.cpp used one variable as both the
# initial abscissa and the length of the integration interval, and never
# advanced the abscissa at all.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Right-hand side of the test problem."
rhs(x, y) = 3 * exp(-x) - 0.4 * y

"Exact solution of the test problem for the initial value `y₀` at `x = 0`."
exact(x, y₀) = (y₀ + 5) * exp(-0.4x) - 5 * exp(-x)

"""
    euler(rhs, x₀, y₀, h, n)

`n` steps of forward Euler. Returns the value at `x₀ + n·h`.
"""
function euler(rhs, x₀, y₀, h, n)
    x, y = x₀, y₀
    for _ in 1:n
        y += h * rhs(x, y)
        x += h
    end
    return y
end

"""
    rk4(rhs, x₀, y₀, h, n)

`n` steps of the classical fourth-order Runge-Kutta method.
"""
function rk4(rhs, x₀, y₀, h, n)
    x, y = x₀, y₀
    for _ in 1:n
        k₁ = rhs(x, y)
        k₂ = rhs(x + h / 2, y + h * k₁ / 2)
        k₃ = rhs(x + h / 2, y + h * k₂ / 2)
        k₄ = rhs(x + h, y + h * k₃)
        y += (h / 6) * (k₁ + 2k₂ + 2k₃ + k₄)
        x += h
    end
    return y
end

"""
    observed_order(h, error)

Least-squares slope of log|error| against log h, i.e. the observed order of
convergence. Points at or below round-off are excluded, since they carry no
information about the discretisation error.
"""
function observed_order(h, error)
    keep = error .> 1e-13
    lx, ly = log.(h[keep]), log.(error[keep])
    n = length(lx)
    return (n * sum(lx .* ly) - sum(lx) * sum(ly)) / (n * sum(abs2, lx) - sum(lx)^2)
end

function main()
    x₀, y₀, x_end = 0.0, 1.0, 5.0
    reference = exact(x_end, y₀)

    steps = [2^k for k in 3:18]
    h = (x_end - x₀) ./ steps
    err_euler = [abs(euler(rhs, x₀, y₀, h[i], steps[i]) - reference) for i in eachindex(steps)]
    err_rk4 = [abs(rk4(rhs, x₀, y₀, h[i], steps[i]) - reference) for i in eachindex(steps)]

    p_euler = observed_order(h, err_euler)
    p_rk4 = observed_order(h, err_rk4)
    @printf("exact y(%.1f) = %.12f\n", x_end, reference)
    @printf("observed order, forward Euler = %.3f  (expected 1)\n", p_euler)
    @printf("observed order, RK4           = %.3f  (expected 4)\n", p_rk4)

    fig = Figure(size = (760, 560))
    ax = Axis(fig[2, 1],
        xlabel = L"Step size $h$",
        ylabel = L"Global error $|y_N - y(x_\mathrm{end})|$",
        xscale = log10, yscale = log10)

    l_e = scatterlines!(ax, h, err_euler, color = PALETTE.blue,
        linewidth = 1.4, markersize = 9)
    l_r = scatterlines!(ax, h, err_rk4, color = PALETTE.orange,
        linewidth = 1.4, markersize = 9, marker = :rect)

    # guide lines anchored on the coarsest step
    guide_e = err_euler[1] .* (h ./ h[1]) .^ 1
    guide_r = err_rk4[1] .* (h ./ h[1]) .^ 4
    g_e = lines!(ax, h, guide_e, color = PALETTE.blue, linestyle = :dash, linewidth = 1.0)
    g_r = lines!(ax, h, guide_r, color = PALETTE.orange, linestyle = :dash, linewidth = 1.0)

    Legend(fig[1, 1],
        [l_e, l_r, g_e, g_r],
        [L"Forward Euler, slope $%$(round(p_euler, digits = 2))$",
         L"RK4, slope $%$(round(p_rk4, digits = 2))$",
         L"$\mathcal{O}(h)$ guide", L"$\mathcal{O}(h^4)$ guide"],
        orientation = :horizontal, framevisible = false, nbanks = 2,
        labelsize = 17, colgap = 20)

    ylims!(ax, 1e-17, 5.0)

    text!(ax, 0.97, 0.06;
        text = L"RK4 reaches round-off near $h \approx 10^{-3}$",
        space = :relative, align = (:right, :bottom),
        color = PALETTE.orange, fontsize = 15)

    rowsize!(fig.layout, 2, Relative(0.84))
    path = savefigure(fig, FIGURES, "integrator_order_comparison")
    println("wrote ", path)
end

main()
