# Four-step Adams-Bashforth predictor with three-step Adams-Moulton corrector,
# in PECE form, on
#
#   y'(x) = -y + 2 sin x,   y(0) = 0,   exact  y = sin x - cos x + e^{-x}
#
# Starting values come from classical RK4.
#
# Ported from Adams_Predictor_Corector.m (Octave, not MATLAB). The Adams
# coefficients in the original were correct — h/24·(55, -59, 37, -9) for the
# predictor and h/24·(9, 19, -5, 1) for the corrector. What it discarded was the
# free local-error estimate: the difference between predictor and corrector is
# proportional to the local truncation error, and is the whole reason to run a
# predictor-corrector pair rather than a corrector alone. It is restored here and
# compared against the true local error. The original also grew `y` element by
# element inside the loop, never checked that the grid held at least four
# starting points, and hardcoded the problem at the call site.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

rhs(x, y) = -y + 2 * sin(x)
exact(x) = sin(x) - cos(x) + exp(-x)

"Predictor coefficients, four-step Adams-Bashforth."
const AB4 = (55, -59, 37, -9) ./ 24
"Corrector coefficients, three-step Adams-Moulton."
const AM3 = (9, 19, -5, 1) ./ 24

"""
    rk4_start(rhs, x, y₀, h, n)

`n` classical RK4 steps, used to generate the starting values the multistep
method needs.
"""
function rk4_start(rhs, x, y₀, h, n)
    ys = [y₀]
    y = y₀
    for i in 1:n
        xi = x + (i - 1) * h
        k₁ = rhs(xi, y)
        k₂ = rhs(xi + h/2, y + h*k₁/2)
        k₃ = rhs(xi + h/2, y + h*k₂/2)
        k₄ = rhs(xi + h, y + h*k₃)
        y += h/6 * (k₁ + 2k₂ + 2k₃ + k₄)
        push!(ys, y)
    end
    return ys
end

"""
    adams_pece(rhs, xspan, y₀, h)

Adams-Bashforth-Moulton in PECE form. Returns the grid, the solution, and the
predictor-corrector difference at each step, which estimates the local
truncation error.
"""
function adams_pece(rhs, xspan, y₀, h)
    x = collect(xspan[1]:h:xspan[2])
    length(x) >= 4 ||
        throw(ArgumentError("need at least four grid points, got $(length(x))"))

    y = Vector{Float64}(undef, length(x))
    y[1:4] .= rk4_start(rhs, x[1], y₀, h, 3)
    pc_difference = zeros(length(x))

    f = [rhs(x[i], y[i]) for i in 1:4]
    for n in 4:(length(x) - 1)
        # predict
        p = y[n] + h * (AB4[1]*f[4] + AB4[2]*f[3] + AB4[3]*f[2] + AB4[4]*f[1])
        # evaluate, correct, evaluate
        fp = rhs(x[n + 1], p)
        y[n + 1] = y[n] + h * (AM3[1]*fp + AM3[2]*f[4] + AM3[3]*f[3] + AM3[4]*f[2])
        pc_difference[n + 1] = abs(y[n + 1] - p)
        f = (f[2], f[3], f[4], rhs(x[n + 1], y[n + 1]))
    end
    return x, y, pc_difference
end

function main()
    h = 0.1
    x, y, pc = adams_pece(rhs, (0.0, 6.0), 0.0, h)
    err = abs.(y .- exact.(x))

    @printf("h = %.2f, %d points\n", h, length(x))
    @printf("max global error          = %.3e\n", maximum(err))
    @printf("max predictor-corrector Δ = %.3e\n", maximum(pc))

    # order check
    orders = Float64[]
    for hh in (0.2, 0.1, 0.05, 0.025, 0.0125)
        xx, yy, _ = adams_pece(rhs, (0.0, 6.0), 0.0, hh)
        push!(orders, abs(last(yy) - exact(last(xx))))
    end
    hs = [0.2, 0.1, 0.05, 0.025, 0.0125]
    slope = (log(orders[end]) - log(orders[1])) / (log(hs[end]) - log(hs[1]))
    @printf("observed order = %.2f (Adams-Bashforth-Moulton PECE is 4)\n", slope)

    fig = Figure(size = (960, 440))

    # Both panels are the same abscissa, so they carry the same limits and the
    # same ticks; the error panel had been left to pick its own.
    xticks = ([0, 1, 2, 3, 4, 5, 6], [L"0", L"1", L"2", L"3", L"4", L"5", L"6"])

    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"y(x)", xticks = xticks)
    l_num = scatter!(ax1, x, y, color = PALETTE.blue, markersize = MARKERSIZE.dense)
    xf = range(0, 6, length = 400)
    l_ex = lines!(ax1, xf, exact.(xf), color = PALETTE.black, linewidth = 1.2)
    xlims!(ax1, -0.25, 6.25)

    ax2 = Axis(fig[2, 2], xlabel = L"x", ylabel = "Error", yscale = log10,
        xticks = xticks, yticks = logticks(-8, -6),)
    l_g = lines!(ax2, x[5:end], err[5:end], color = PALETTE.orange, linewidth = 1.5)
    l_p = lines!(ax2, x[5:end], pc[5:end], color = PALETTE.green, linewidth = 1.5)
    xlims!(ax2, -0.25, 6.25)

    Legend(fig[1, 1:2], [l_num, l_ex, l_g, l_p],
        ["Adams PECE", "Exact",
            # the takeaway as a legend entry: the error panel has no band clear of
            # both curves wide enough to hold it
            latexstring(@sprintf("\\text{Global error, observed order } %.2f", slope)),
            "Predictor–corrector Δ",],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24,)

    rowsize!(fig.layout, 2, Relative(0.85))
    path = savefigure(fig, FIGURES, "adams_predictor_corrector")
    println("wrote ", path)
end

main()
