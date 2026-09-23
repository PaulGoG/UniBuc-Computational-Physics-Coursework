# Four-step Adams–Bashforth predictor with three-step Adams–Moulton corrector,
# in PECE form, on
#
#   y′(x) = −y + 2 sin x,   y(0) = 0,   exact  y = sin x − cos x + e⁻ˣ
#
# with starting values from classical RK4, as in Adams_Predictor_Corector.m on
# the `legacy` branch. The coefficients there are correct: h/24·(55, −59, 37, −9)
# for the predictor and h/24·(9, 19, −5, 1) for the corrector.
#
# What the original did not use is the error estimate a predictor–corrector pair
# gives for free. With y the exact solution through the back values, the two
# formulae have local errors
#
#   y(xₙ₊₁) − p = +251/720 h⁵ y⁽⁵⁾,     y(xₙ₊₁) − c = −19/720 h⁵ y⁽⁵⁾,
#
# so c − y(xₙ₊₁) = 19/270 (c − p) to leading order: Milne's device. That is an
# estimate of the one-step error, O(h⁵), and is compared here with the one-step
# error itself, measured by taking a step from exact back values. The global
# error is a different quantity, O(h⁴), and is measured separately.
#
# Reference for the error constants and the device: Hairer, Nørsett and Wanner,
# Solving Ordinary Differential Equations I, 2nd ed., Springer 1993, §III.1–III.2,
# doi:10.1007/978-3-540-78862-1.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "numerics_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

rhs(x, y) = -y + 2 * sin(x)
exact(x) = sin(x) - cos(x) + exp(-x)

"Predictor coefficients, four-step Adams–Bashforth, newest value first."
const AB4 = (55, -59, 37, -9) ./ 24
"Corrector coefficients, three-step Adams–Moulton, newest value first."
const AM3 = (9, 19, -5, 1) ./ 24
"Milne's factor: corrector error constant 19/720 over the sum 251/720 + 19/720."
const MILNE = 19 / 270
"Integration interval of the original."
const XSPAN = (0.0, 6.0)
"Step size of the original."
const H_LEGACY = 0.1
"Step sizes of the convergence sweep, each half the one before."
const H_SWEEP = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625]

"""
    pece_step(rhs, x_next, y_n, h, f)

One predict–evaluate–correct step from `y_n` to `x_next`, with `f` the four most
recent derivative values, oldest first. Returns the predicted and the corrected
value.
"""
function pece_step(rhs, x_next, y_n, h, f)
    p = y_n + h * (AB4[1] * f[4] + AB4[2] * f[3] + AB4[3] * f[2] + AB4[4] * f[1])
    c = y_n + h * (AM3[1] * rhs(x_next, p) + AM3[2] * f[4] + AM3[3] * f[3] + AM3[4] * f[2])
    return p, c
end

"""
    adams_pece(rhs, xspan, y₀, h)

Adams–Bashforth–Moulton in PECE form over `xspan` with step `h`, started by
three RK4 steps. Returns the grid, the solution, and the signed
corrector-minus-predictor difference at each multistep point (zero at the four
starting points).
"""
function adams_pece(rhs, xspan, y₀, h)
    x = collect(range(xspan[1], xspan[2]; step = h))
    length(x) >= 5 ||
        throw(ArgumentError("need at least five grid points, got $(length(x))"))

    y = Vector{Float64}(undef, length(x))
    y[1] = y₀
    for i in 1:3
        y[i + 1] = rk4_step(rhs, x[i], y[i], h)
    end
    difference = zeros(length(x))

    f = ntuple(i -> rhs(x[i], y[i]), 4)
    for n in 4:(length(x) - 1)
        p, y[n + 1] = pece_step(rhs, x[n + 1], y[n], h, f)
        difference[n + 1] = y[n + 1] - p
        f = (f[2], f[3], f[4], rhs(x[n + 1], y[n + 1]))
    end
    return x, y, difference
end

"""
    one_step_errors(rhs, exact, xspan, h)

The one-step error of the PECE pair and Milne's estimate of it, at every
multistep point of the grid. Each step starts from exact back values, so what is
measured is the local error alone, with nothing accumulated. Returns the grid
points, `c − y(x)` and `MILNE · (c − p)`.
"""
function one_step_errors(rhs, exact, xspan, h)
    x = collect(range(xspan[1], xspan[2]; step = h))
    points = 5:length(x)
    local_error = Vector{Float64}(undef, length(points))
    estimate = Vector{Float64}(undef, length(points))
    for (k, m) in enumerate(points)
        f = ntuple(i -> rhs(x[m - 5 + i], exact(x[m - 5 + i])), 4)
        p, c = pece_step(rhs, x[m], exact(x[m - 1]), h, f)
        local_error[k] = c - exact(x[m])
        estimate[k] = MILNE * (c - p)
    end
    return x[points], local_error, estimate
end

"Largest absolute entry and the index where it occurs."
peak(v) = findmax(abs, v)

function main()
    h = H_LEGACY
    x, y, _ = adams_pece(rhs, XSPAN, 0.0, h)
    global_error = abs.(y .- exact.(x))
    x_local, local_error, estimate = one_step_errors(rhs, exact, XSPAN, h)

    @printf("h = %.2f, %d points\n", h, length(x))
    @printf("max global error          = %.3e\n", maximum(global_error))
    @printf("max one-step error        = %.3e\n", peak(local_error)[1])
    @printf("max Milne estimate        = %.3e\n", peak(estimate)[1])

    # convergence sweep: the three quantities in the maximum norm
    sweep = map(H_SWEEP) do hh
        xx, yy, _ = adams_pece(rhs, XSPAN, 0.0, hh)
        _, le, es = one_step_errors(rhs, exact, XSPAN, hh)
        k = peak(le)[2]
        (global_error = maximum(abs.(yy .- exact.(xx))), local_error = abs(le[k]),
            estimate = peak(es)[1], ratio = es[k] / le[k],)
    end
    e_global = [s.global_error for s in sweep]
    e_local = [s.local_error for s in sweep]
    e_milne = [s.estimate for s in sweep]
    ratio = [s.ratio for s in sweep]
    o_global = convergence_order(H_SWEEP, e_global)
    o_local = convergence_order(H_SWEEP, e_local)

    println("\n      h    global error  order   one-step error  order   Milne/true")
    for i in eachindex(H_SWEEP)
        @printf("%7.5f   %.3e   %5s    %.3e    %5s    %.4f\n", H_SWEEP[i], e_global[i],
            i == 1 ? "" : @sprintf("%.2f", o_global.pairwise[i - 1]), e_local[i],
            i == 1 ? "" : @sprintf("%.2f", o_local.pairwise[i - 1]), ratio[i])
    end
    @printf("order between the two finest steps: global %.2f, one-step %.2f\n",
        last(o_global.pairwise), last(o_local.pairwise))

    # The corrector is evaluated at the predicted value, which adds
    # h·(9/24)·∂f/∂y·(p − c) to its error: a relative correction of
    # (9/24)(270/19) h |∂f/∂y| = 5.33 h for this equation, where ∂f/∂y = −1.
    pece_correction = AM3[1] / MILNE
    predicted_ratio = 1 / (1 + pece_correction * last(H_SWEEP))
    @printf("Milne/true at h = %.5f: %.4f, predicted 1/(1 + %.2f h) = %.4f\n",
        last(H_SWEEP), last(ratio), pece_correction, predicted_ratio)

    abs(last(o_global.pairwise) - 4) < 0.1 ||
        error("global error is not fourth order: $(last(o_global.pairwise))")
    abs(last(o_local.pairwise) - 5) < 0.1 ||
        error("one-step error is not fifth order: $(last(o_local.pairwise))")
    isapprox(last(ratio), predicted_ratio; rtol = 0.01) ||
        error("Milne estimate is off the one-step error: ratio $(last(ratio))")

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1500, 660))
    xticks = 0:1:6

    ax1 = Axis(fig[2, 1], xlabel = L"x", ylabel = L"y(x)", xticks = xticks)
    xf = range(XSPAN..., length = 400)
    l_exact = lines!(ax1, xf, exact.(xf), color = PALETTE.black)
    l_pece = scatter!(ax1, x, y, color = PALETTE.blue, markersize = MARKERSIZE.dense)
    xlims!(ax1, -0.2, 6.2)

    decade_floor(v) = floor(Int, log10(minimum(v)))
    decade_ceil(v) = ceil(Int, log10(maximum(v)))
    shown = vcat(global_error[5:end], abs.(local_error), abs.(estimate))
    ax2 = Axis(fig[2, 2], xlabel = L"x", ylabel = "Absolute error", yscale = log10,
        xticks = xticks,
        yticks = logticks(decade_floor(shown), decade_ceil(shown); style = :power),)
    l_global = lines!(ax2, x[5:end], global_error[5:end], color = PALETTE.orange)
    l_local = lines!(ax2, x_local, abs.(local_error), color = PALETTE.green)
    l_milne = lines!(ax2, x_local, abs.(estimate), color = PALETTE.green,
        linestyle = :dash,)
    xlims!(ax2, -0.2, 6.2)
    ylims!(ax2, 10.0^decade_floor(shown), 10.0^decade_ceil(shown))
    text!(ax2, 0.04, 0.04;
        text = rich("Milne / true = ", @sprintf("%.2f", ratio[findfirst(==(h), H_SWEEP)]),
            " at ", it("h"), @sprintf(" = %.1f", h)),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,
        color = PALETTE.green,)

    # Guides of slope 4 and 5 through the finest point of each series, drawn a
    # factor of six clear of the data so that the two can be told apart.
    guide_gap = 6
    guide_4 = guide_gap * last(e_global) .* (H_SWEEP ./ last(H_SWEEP)) .^ 4
    guide_5 = last(e_local) / guide_gap .* (H_SWEEP ./ last(H_SWEEP)) .^ 5
    swept = vcat(e_global, e_local, e_milne, guide_4, guide_5)
    ax3 = Axis(fig[2, 3], xlabel = L"Step size $h$", ylabel = "Maximum error",
        xscale = log10, yscale = log10,
        xticks = (
            [0.01, 0.02, 0.05, 0.1, 0.2], [L"0.01", L"0.02", L"0.05", L"0.1", L"0.2"],),
        yticks = logticks(decade_floor(swept), decade_ceil(swept); step = 2),)
    lines!(ax3, H_SWEEP, guide_4, color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    text!(ax3, H_SWEEP[end - 1], guide_4[end - 1]; text = L"\propto h^{4}",
        align = (:right, :bottom), offset = (-6, 4), fontsize = ANNOTATION_SIZE,)
    lines!(ax3, H_SWEEP, guide_5, color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    text!(ax3, H_SWEEP[2], guide_5[2]; text = L"\propto h^{5}",
        align = (:left, :top), offset = (6, -4), fontsize = ANNOTATION_SIZE,)
    scatterlines!(ax3, H_SWEEP, e_global, color = PALETTE.orange)
    scatterlines!(ax3, H_SWEEP, e_local, color = PALETTE.green, marker = :rect)
    scatterlines!(ax3, H_SWEEP, e_milne, color = PALETTE.green, marker = :diamond,
        linestyle = :dash,)
    xlims!(ax3, 0.0045, 0.28)
    ylims!(ax3, 10.0^decade_floor(swept), 10.0^decade_ceil(swept))
    text!(ax3, 0.97, 0.04;
        text = rich(
            rich(@sprintf("Global: order %.2f", last(o_global.pairwise)),
                color = PALETTE.orange,),
            "\n",
            rich(@sprintf("One-step: order %.2f", last(o_local.pairwise)),
                color = PALETTE.green,),),
        space = :relative, align = (:right, :bottom), justification = :right,
        fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1:3],
        [[l_pece, l_exact], [l_global, l_local, l_milne]],
        [["Adams PECE", "Exact"], ["Global", "One-step", "Milne estimate"]],
        ["Solution", "Error"]; titleposition = :left, nbanks = 1, groupgap = 40,)

    path = savefigure(fig, FIGURES, "adams_predictor_corrector")
    println("\nwrote ", path)
end

main()
