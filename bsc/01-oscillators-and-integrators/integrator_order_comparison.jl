# Observed order of convergence of forward Euler and classical Runge–Kutta 4 on
# two scalar initial-value problems:
#
#   y' = 3e^{-x} - 0.4y,   y(0) = 1,  x ∈ [0, 5]    (RKtrial.cpp)
#   y' = y + sin(0.4y),    y(0) = 1,  x ∈ [0, 1]    (DiffEqEuler.cpp)
#
# The first has the closed-form solution y = (y₀ + 5)e^{-0.4x} - 5e^{-x}. The
# second has none, so its reference is RK4 in 256-bit arithmetic, checked
# against itself at half the step count.
#
# The order is read from pairs of successive step sizes,
#
#   p(h) = ln[e(2h)/e(h)] / ln 2,
#
# which tends to the order of the scheme as h → 0 and is meaningful only while
# e(h) stands clear of accumulated round-off. Each pair carries two error terms:
# the pre-asymptotic drift |p(h) - p(2h)|, which for p = p₀ + c·h equals the
# distance still to go, and the round-off term ε/(e ln 2), with ε the round-off
# floor measured from the sweep itself. A pair is used while the second is below
# the first; the order quoted is that of the finest such pair, with the sum of
# the two terms as its uncertainty. A single straight-line fit through the whole
# sweep mixes the pre-asymptotic points at large h with the floor at small h.
#
# Both schemes are coded as in the 2018 C++ sources. Those programs read N and
# the interval from standard input and printed the endpoint value only; the
# step-size sweep and the comparison with a reference are added here.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"""
A pairwise order further than this from its predecessor marks the end of the
power-law regime: the errors from there on are accumulated round-off.
"""
const ORDER_BREAKDOWN = 1.0

"Smallest uncertainty quoted for an order; the orders are printed to three decimals."
const ORDER_RESOLUTION = 1e-3

"Right-hand side of the linear problem of RKtrial.cpp."
rhs_linear(x, y) = 3 * exp(-x) - 0.4 * y

"Right-hand side of the autonomous nonlinear problem of DiffEqEuler.cpp."
rhs_nonlinear(x, y) = y + sin(0.4 * y)

"Exact solution of the linear problem for the initial value `y₀` at `x = 0`."
exact_linear(x, y₀) = (y₀ + 5) * exp(-0.4x) - 5 * exp(-x)

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

`n` steps of the classical fourth-order Runge–Kutta method.
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
    raw_orders(h, err)

Observed order between every pair of successive step sizes,
`ln(eᵢ/eᵢ₊₁) / ln(hᵢ/hᵢ₊₁)`, with no selection applied.
"""
raw_orders(h, err) = [log(err[i] / err[i + 1]) / log(h[i] / h[i + 1])
                      for i in 1:(length(h) - 1)]

"""
    roundoff_floor(h, err)

Accumulated round-off of the sweep, measured as the largest error after the
power law has broken down (`ORDER_BREAKDOWN`). Round-off grows with the step
count, so this overstates it at the coarser steps the order is read from, which
makes the quoted uncertainty conservative. Zero if the sweep never reaches the
floor.
"""
function roundoff_floor(h, err)
    p = raw_orders(h, err)
    k = findfirst(i -> abs(p[i] - p[i - 1]) > ORDER_BREAKDOWN, 2:length(p))
    return k === nothing ? 0.0 : maximum(err[(k + 2):end])
end

"""
    pairwise_orders(h, err)

Usable pairwise orders: those whose round-off term `ε / (e ln 2)` is below their
drift `|p(h) - p(2h)|`. Returns the geometric-mean step, the order, the drift
and the round-off term of each.
"""
function pairwise_orders(h, err)
    p = raw_orders(h, err)
    ε = roundoff_floor(h, err)
    h_pair, order, drift, roundoff = Float64[], Float64[], Float64[], Float64[]
    for i in 2:length(p)
        δ = abs(p[i] - p[i - 1])
        r = ε / (err[i + 1] * log(2))
        (δ < ORDER_BREAKDOWN && r < δ) || (ε > 0 && break)
        r < δ || continue
        push!(h_pair, sqrt(h[i] * h[i + 1]))
        push!(order, p[i])
        push!(drift, δ)
        push!(roundoff, r)
    end
    return (; h = h_pair, p = order, drift, roundoff)
end

"""
    asymptotic_order(orders)

Order at the finest usable pair, with drift plus round-off as its uncertainty.
Both terms are bounds on a bias, not variances, so they add linearly. Never
quoted below `ORDER_RESOLUTION`.
"""
function asymptotic_order(orders)
    isempty(orders.p) &&
        throw(ArgumentError("no pairwise order clears the round-off floor"))
    return orders.p[end], max(orders.drift[end] + orders.roundoff[end], ORDER_RESOLUTION)
end

"""
    high_precision_reference(rhs, x₀, y₀, x_end, n)

RK4 in 256-bit arithmetic at `n` steps, accepted only if it agrees with the
`n ÷ 2` result to better than 1e-17, which puts its own error (a further factor
16 smaller) far below double precision.
"""
function high_precision_reference(rhs, x₀, y₀, x_end, n)
    setprecision(BigFloat, 256) do
        L = big(x_end) - big(x₀)
        fine = rk4(rhs, big(x₀), big(y₀), L / n, n)
        coarse = rk4(rhs, big(x₀), big(y₀), 2L / n, n ÷ 2)
        abs(fine - coarse) < big"1e-17" ||
            error("high-precision reference not converged: |Δ| = $(Float64(abs(fine - coarse)))")
        Float64(fine)
    end
end

"""
    sweep(rhs, x₀, y₀, x_end, reference, exponents)

Global error of both schemes at `x_end` for `2^k` steps, `k` in `exponents`.
"""
function sweep(rhs, x₀, y₀, x_end, reference, exponents)
    steps = [2^k for k in exponents]
    h = (x_end - x₀) ./ steps
    err_euler = [abs(euler(rhs, x₀, y₀, hᵢ, n) - reference) for (hᵢ, n) in zip(h, steps)]
    err_rk4 = [abs(rk4(rhs, x₀, y₀, hᵢ, n) - reference) for (hᵢ, n) in zip(h, steps)]
    return (; h, err_euler, err_rk4)
end

function report(name, s)
    println(name)
    println("       h        Euler error   p      RK4 error    p")
    p_euler, p_rk4 = raw_orders(s.h, s.err_euler), raw_orders(s.h, s.err_rk4)
    for i in eachindex(s.h)
        order(p) = i > 1 ? @sprintf("%6.3f", p[i - 1]) : "     –"
        @printf("  %.3e    %.3e  %s    %.3e  %s\n",
            s.h[i], s.err_euler[i], order(p_euler), s.err_rk4[i], order(p_rk4))
    end
    for (scheme, err) in (("forward Euler", s.err_euler), ("RK4", s.err_rk4))
        o = pairwise_orders(s.h, err)
        @printf("  %-13s round-off floor %.1e; finest usable pair h = %.2e: drift %.4f, round-off %.4f\n",
            scheme, roundoff_floor(s.h, err), o.h[end], o.drift[end], o.roundoff[end])
    end
    p_e, σ_e = asymptotic_order(pairwise_orders(s.h, s.err_euler))
    p_r, σ_r = asymptotic_order(pairwise_orders(s.h, s.err_rk4))
    @printf("  order: forward Euler %.3f ± %.3f, RK4 %.3f ± %.3f\n\n", p_e, σ_e, p_r, σ_r)
    return (p_e, σ_e, p_r, σ_r)
end

"Draw one problem: the error panel with its guides, and the order strip beneath."
function draw_problem!(ax, strip, s, orders, equation)
    p_e, σ_e, p_r, σ_r = orders
    floor_rk4 = roundoff_floor(s.h, s.err_rk4)
    band = hspan!(ax, 1e-17, floor_rk4, color = (:grey, 0.15))
    text!(ax, 0.8, floor_rk4 / 3; text = "Round-off floor",
        space = :data, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        color = :grey35,)

    g_e = lines!(ax, s.h, s.err_euler[1] .* (s.h ./ s.h[1]) .^ 1,
        color = PALETTE.blue, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    g_r = lines!(ax, s.h, s.err_rk4[1] .* (s.h ./ s.h[1]) .^ 4,
        color = PALETTE.orange, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    l_e = scatterlines!(ax, s.h, s.err_euler, color = PALETTE.blue, marker = :circle)
    l_r = scatterlines!(ax, s.h, s.err_rk4, color = PALETTE.orange, marker = :rect)

    text!(ax, 0.03, 0.97; text = equation, space = :relative, align = (:left, :top),
        fontsize = ANNOTATION_SIZE,)
    text!(ax, 0.03, 0.56; text = @sprintf("Euler order %.3f ± %.3f", p_e, σ_e),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.blue,)
    text!(ax, 0.03, 0.48; text = @sprintf("RK4 order %.3f ± %.3f", p_r, σ_r),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.orange,)

    hlines!(strip, [1, 4], color = [PALETTE.blue, PALETTE.orange], linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    o_e, o_r = pairwise_orders(s.h, s.err_euler), pairwise_orders(s.h, s.err_rk4)
    scatterlines!(strip, o_e.h, o_e.p, color = PALETTE.blue, marker = :circle)
    scatterlines!(strip, o_r.h, o_r.p, color = PALETTE.orange, marker = :rect)
    errorbars!(strip, o_r.h, o_r.p, o_r.drift .+ o_r.roundoff, color = PALETTE.orange,
        whiskerwidth = 8,)
    return (l_e, l_r, g_e, g_r, band)
end

function main()
    linear = sweep(rhs_linear, 0.0, 1.0, 5.0, exact_linear(5.0, 1.0), 3:18)
    reference = high_precision_reference(rhs_nonlinear, 0.0, 1.0, 1.0, 2^16)
    nonlinear = sweep(rhs_nonlinear, 0.0, 1.0, 1.0, reference, 3:16)

    @printf("linear problem:    exact y(5) = %.15f\n", exact_linear(5.0, 1.0))
    @printf("nonlinear problem: reference y(1) = %.15f\n\n", reference)
    orders_linear = report("y' = 3e^{-x} - 0.4y on [0, 5]", linear)
    orders_nonlinear = report("y' = y + sin(0.4y) on [0, 1]", nonlinear)

    # each scheme must reach its nominal order within the quoted uncertainty
    for (p_e, σ_e, p_r, σ_r) in (orders_linear, orders_nonlinear)
        @assert abs(p_e - 1) <= 2σ_e "forward Euler order $p_e is not 1 within $(2σ_e)"
        @assert abs(p_r - 4) <= 2σ_r "RK4 order $p_r is not 4 within $(2σ_r)"
    end

    fig = Figure(size = (1200, 880))
    error_label = L"Global error $|y_N - y(x_\mathrm{end})|$"
    ax_l = Axis(fig[2, 1], ylabel = error_label, xscale = log10, yscale = log10,
        xticks = logticks(-4, 0), yticks = logticks(-15, 0; step = 5),)
    ax_n = Axis(fig[2, 2], xscale = log10, yscale = log10,
        xticks = logticks(-4, 0), yticks = logticks(-15, 0; step = 5),)
    strip_l = Axis(fig[3, 1], xlabel = L"Step size $h$", ylabel = L"Order $p(h)$",
        xscale = log10, xticks = logticks(-4, 0), yticks = 1:4,)
    strip_n = Axis(fig[3, 2], xlabel = L"Step size $h$",
        xscale = log10, xticks = logticks(-4, 0), yticks = 1:4,)

    handles = draw_problem!(ax_l, strip_l, linear, orders_linear,
        L"y' = 3e^{-x} - 0.4y",)
    draw_problem!(ax_n, strip_n, nonlinear, orders_nonlinear, L"y' = y + \sin(0.4y)")

    linkxaxes!(ax_l, strip_l, ax_n, strip_n)
    linkyaxes!(ax_l, ax_n)
    linkyaxes!(strip_l, strip_n)
    xlims!(ax_l, 1.2e-5, 1.3)
    ylims!(ax_l, 1e-17, 5.0)
    ylims!(strip_l, 0.4, 4.8)
    hidexdecorations!(ax_l, grid = false, ticks = false)
    hidexdecorations!(ax_n, grid = false, ticks = false)
    hideydecorations!(ax_n, grid = false, ticks = false)
    hideydecorations!(strip_n, grid = false, ticks = false)

    l_e, l_r, g_e, g_r, _ = handles
    Legend(fig[1, 1:2], [[l_e, l_r], [g_e, g_r]],
        [["Forward Euler", "RK4"], [L"\mathcal{O}(h)", L"\mathcal{O}(h^4)"]],
        ["Global error", "Guides"]; titleposition = :left, nbanks = 1,
        titlesize = 22, titlegap = 14, groupgap = 40,)

    rowsize!(fig.layout, 3, Auto(0.32))
    path = savefigure(fig, FIGURES, "integrator_order_comparison")
    println("wrote ", path)
end

main()
