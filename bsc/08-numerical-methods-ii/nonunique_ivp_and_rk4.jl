# Two initial-value problems from the fourth-year numerical-methods exam,
#
#   (a)  y' = √|sin y|            (b)  y' = z,  z' = −y sin y,  y(0) = 1, z(0) = 0
#
# ported from ODE_RK4.jl and ODE_system_RK4.jl in
# Julia-Workflow-FFUB/Examen_PDF_MN_II_L_4/ on the `legacy` branch.
#
# (a) is not well posed as an initial-value problem on any of the lines y = kπ:
# √|sin y| vanishes there with infinite slope, so it is not Lipschitz and
# uniqueness fails. From y(0) = 0, which is how the file's comment states the
# problem, y ≡ 0 is a solution and so is every solution that waits at 0 for an
# arbitrary time before leaving; the file itself starts from y(0) = 1, which
# only moves the ambiguity to the first crossing of π. Away from the zeros the
# equation separates, x = ∫ dy/√|sin y|, and the integral converges at every
# zero: the solution that never waits crosses each multiple of π after a
# finite time, T₁ = ∫₀^π dy/√(sin y) = √π Γ(1/4)/Γ(3/4) = 5.2441 per half-turn.
# `main` asserts that closed form by quadrature and uses the inverse of x(y)
# as the reference for the never-waiting solution.
#
# RK4 cannot resolve the non-Lipschitz points. Started exactly at 0 it stays at
# 0; started at any ε > 0, or at 1, it climbs, but at every crossing of kπ the
# scheme loses accuracy, so its value at x = 100 converges only slowly as the
# step is halved. The table printed below shows the digits that survive.
#
# (b) is conservative: E = z²/2 + sin y − y cos y is asserted invariant. The
# original integrated it correctly but called its stepper twice per step, once
# per component, discarding half of each result.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, QuadGK, SpecialFunctions
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Right end of the interval."
const X_END = 100.0
"Step sizes of the convergence table; the first is the original's 100/999 ≈ 0.1, rounded."
const STEPS = (0.1, 0.02, 0.01, 0.005, 0.0025, 0.00125)
"Initial values followed for the scalar problem."
const STARTS = (0.0, 1e-12, 1e-6, 1.0)
"Tolerance on the invariant of the system and on the half-turn time."
const CHECK_TOLERANCE = 1e-6

"Right-hand side of the scalar problem; not Lipschitz where sin y = 0."
f_scalar(y) = sqrt(abs(sin(y)))
"Right-hand side of the second-order problem as a first-order system."
f_system(y, z) = (z, -sin(y) * y)

"RK4 for the scalar autonomous problem; returns the whole trajectory."
function rk4_scalar(f, y₀, h, n)
    y = Vector{Float64}(undef, n + 1)
    y[1] = y₀
    for i in 1:n
        k₁ = f(y[i])
        k₂ = f(y[i] + h * k₁ / 2)
        k₃ = f(y[i] + h * k₂ / 2)
        k₄ = f(y[i] + h * k₃)
        y[i + 1] = y[i] + h / 6 * (k₁ + 2k₂ + 2k₃ + k₄)
    end
    return y
end

"RK4 for the two-component system, each step computed once."
function rk4_system(f, y₀, z₀, h, n)
    y = Vector{Float64}(undef, n + 1)
    z = Vector{Float64}(undef, n + 1)
    y[1], z[1] = y₀, z₀
    for i in 1:n
        k1y, k1z = f(y[i], z[i])
        k2y, k2z = f(y[i] + h * k1y / 2, z[i] + h * k1z / 2)
        k3y, k3z = f(y[i] + h * k2y / 2, z[i] + h * k2z / 2)
        k4y, k4z = f(y[i] + h * k3y, z[i] + h * k3z)
        y[i + 1] = y[i] + h / 6 * (k1y + 2k2y + 2k3y + k4y)
        z[i + 1] = z[i] + h / 6 * (k1z + 2k2z + 2k3z + k4z)
    end
    return y, z
end

"Time for the never-waiting solution of (a) to go from y = 0 to y = π, in closed form."
half_turn_closed() = sqrt(π) * gamma(1 / 4) / gamma(3 / 4)

"""
    partial_half_turn(u)

∫₀ᵘ dv/√(sin v) for 0 ≤ u ≤ π, with the substitution v = w² that removes the
integrable singularity at v = 0.
"""
function partial_half_turn(u)
    u <= 0 && return 0.0
    return quadgk(w -> 2w / sqrt(sin(w^2)), 0, sqrt(u); rtol = 1e-12)[1]
end

"""
    travel_time(y₀, y)

x needed by the never-waiting solution of (a) to go from `y₀` to `y ≥ y₀`:
∫ du/√|sin u| assembled from whole half-turns of length T₁ and the two partial
half-turns at either end, |sin| being π-periodic.
"""
function travel_time(y₀, y)
    T₁ = half_turn_closed()
    k₀, k = floor(Int, y₀ / π), floor(Int, y / π)
    k == k₀ && return partial_half_turn(y - k * π) - partial_half_turn(y₀ - k₀ * π)
    return (T₁ - partial_half_turn(y₀ - k₀ * π)) + (k - k₀ - 1) * T₁ +
           partial_half_turn(y - k * π)
end

"Value at `x` of the never-waiting solution of (a) started at `y₀`, by bisection on `travel_time`."
function reference_solution(y₀, x)
    lo, hi = y₀, y₀ + (x / half_turn_closed() + 2) * π
    for _ in 1:100
        mid = (lo + hi) / 2
        travel_time(y₀, mid) < x ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

function main()
    T₁ = partial_half_turn(π)
    @printf("half-turn time of the never-waiting solution: quadrature %.6f, closed form √π Γ(1/4)/Γ(3/4) = %.6f\n",
        T₁, half_turn_closed())
    isapprox(T₁, half_turn_closed(); rtol = CHECK_TOLERANCE) ||
        error("half-turn time $T₁ differs from the closed form")

    # (a): RK4 from the four starts at the original's step, and the convergence
    # of y(100) with the step for the two starts that leave
    references = Dict(y₀ => reference_solution(y₀, X_END) for y₀ in STARTS[2:end])
    @printf("\ny(%.0f) of the never-waiting solution: %.4f from 10⁻¹², %.4f from 10⁻⁶, %.4f from 1\n",
        X_END, references[1e-12], references[1e-6], references[1.0])
    println("RK4 y(100) against the step:")
    println("       h      from 0    from 10⁻¹²   from 10⁻⁶     from 1")
    table = Dict{Float64,Vector{Float64}}()
    for h in STEPS
        n = round(Int, X_END / h)
        table[h] = [last(rk4_scalar(f_scalar, y₀, h, n)) for y₀ in STARTS]
        @printf("  %7.5f   %8.5f   %9.4f   %9.4f   %9.4f\n", h, table[h]...)
    end
    all(table[h][1] == 0 for h in STEPS) ||
        error("RK4 from y(0) = 0 left the zero solution")
    drift(k) = table[STEPS[end]][k] - table[STEPS[end - 1]][k]
    @printf("change of y(100) between the two finest steps: %.4f from 10⁻¹², %.4f from 1 — the second decimal is not converged\n",
        drift(2), drift(4))

    # (b): the conservative system at the original's resolution
    h, n = 0.02, 5000
    y, z = rk4_system(f_system, 1.0, 0.0, h, n)
    energy(y, z) = z^2 / 2 + sin(y) - y * cos(y)
    E = energy.(y, z)
    E_drift = maximum(abs.(E .- E[1]))
    @printf("\n(b) invariant z²/2 + sin y − y cos y: drift %.2e over %d steps of %.2f\n",
        E_drift, n, h)
    E_drift < CHECK_TOLERANCE || error("the invariant of (b) drifts by $E_drift")

    # --- figure ---------------------------------------------------------------
    h_plot = 0.02
    x = range(0, X_END, length = round(Int, X_END / h_plot) + 1)
    solutions = [rk4_scalar(f_scalar, y₀, h_plot, length(x) - 1) for y₀ in STARTS]
    x_ref = range(0, X_END, length = 400)
    y_ref = [reference_solution(1e-12, xi) for xi in x_ref]

    fig = Figure(size = (1400, 640))
    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"y(x)", xticks = 0:20:100)
    colours = (PALETTE.black, PALETTE.blue, PALETTE.orange, PALETTE.green)
    styles = (:solid, :solid, :dash, :solid)
    handles = [lines!(ax1, x, sol, color = colours[k], linestyle = styles[k])
               for (k, sol) in enumerate(solutions)]
    l_ref = lines!(
        ax1, x_ref, y_ref, color = PALETTE.red, linestyle = :dot, linewidth = GUIDE_WIDTH,)
    xlims!(ax1, 0, X_END)
    ylims!(ax1, -3, 70)
    text!(ax1, 0.03, 0.96;
        text = rich(
            rich(
                @sprintf("Never-waiting solution from 0⁺: y(100) = %.2f",
                    references[1e-12]),
                color = PALETTE.red,),
            "\n",
            rich(
                @sprintf("RK4 from 10⁻¹²: %.2f at h = 0.02, %.2f at h = 0.00125",
                    table[0.02][2], table[0.00125][2]),
                color = PALETTE.blue,), "\n",
            rich(
                @sprintf("RK4 from 1: %.2f at h = 0.02, %.2f at h = 0.00125",
                    table[0.02][4], table[0.00125][4]),
                color = PALETTE.green,),),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)

    ax2 = Axis(fig[2, 2], xlabel = L"y", ylabel = L"z = y'")
    lines!(ax2, y, z, color = PALETTE.purple)
    text!(ax2, 0.97, 0.06;
        text = rich("Invariant drift ", rsci(E_drift; digits = 1), @sprintf(" over %d RK4 steps",
            n)),
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1:2], [handles..., l_ref],
        [L"y(0) = 0", L"y(0) = 10^{-12}", L"y(0) = 10^{-6}",
            L"y(0) = 1", "Quadrature reference",],)
    println("\nwrote ", savefigure(fig, FIGURES, "nonunique_ivp_and_rk4"))
end

main()
