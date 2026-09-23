# Classical RK4 on a three-component linear system with a closed-form solution:
#
#   y₁′ = y₂
#   y₂′ = −y₁ − 2eᵗ + 1
#   y₃′ = −y₁ − eᵗ + 1          on t ∈ [1, 8],  y(1) = (1, 2, 3)
#
# Components 1 and 2 form a driven linear oscillator, y₁″ + y₁ = 1 − 2eᵗ, whose
# solution is y₁ = C₁cos t + C₂sin t + 1 − eᵗ. Component 3 is a quadrature of
# it, y₃′ = −C₁cos t − C₂sin t.
#
# The system, the interval, the initial values and the 100 steps are those of
# Sys_ODE_RK4.m on the `legacy` branch, whose tableau is correct. Its orientation
# check reads
#
#     m = size(alpha);  if m==1  alpha = alpha';  endif
#
# where `size` returns [1 3], so `m==1` is [1 0] and the `if`, which requires
# every element non-zero, never fires; the routine works because assigning a row
# into `w(:,1)` reorients it. The routine also calls its right-hand side by name
# instead of taking it as an argument, prints the 101×4 result twice through
# missing semicolons, unpacks the state by column-major linear indexing, and is
# not compared with the closed form above.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "numerics_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Right-hand side of the system."
f(t, y) = (y[2], -y[1] - 2exp(t) + 1, -y[1] - exp(t) + 1)

"Integration interval and initial values of the original."
const TSPAN = (1.0, 8.0)
const Y₀ = (1.0, 2.0, 3.0)
"Number of steps of the original."
const N_LEGACY = 100
"Step counts of the convergence sweep, each double the one before."
const N_SWEEP = [25, 50, 100, 200, 400, 800]

# integration constants fixed by y(1) = (1, 2, 3)
const C₁ = exp(1) * cos(1) - (2 + exp(1)) * sin(1)
const C₂ = exp(1) * sin(1) + (2 + exp(1)) * cos(1)
const C₃ = 3 + C₁ * sin(1) - C₂ * cos(1)

"Exact solution, all three components."
exact(t) = (C₁ * cos(t) + C₂ * sin(t) + 1 - exp(t),
    -C₁ * sin(t) + C₂ * cos(t) - exp(t),
    -C₁ * sin(t) + C₂ * cos(t) + C₃,)

"""
    rk4_system(f, tspan, y₀, n)

`n` RK4 steps over `tspan` on a system of any dimension, with the right-hand
side `f(t, y)` an argument. Returns the grid and a `length(y₀) × (n + 1)` matrix
of states.
"""
function rk4_system(f, tspan, y₀::NTuple{N,Float64}, n::Integer) where {N}
    n > 0 || throw(ArgumentError("need at least one step, got $n"))
    h = (tspan[2] - tspan[1]) / n
    t = collect(range(tspan[1], tspan[2], length = n + 1))
    y = Matrix{Float64}(undef, N, n + 1)
    u = y₀
    y[:, 1] .= u
    for i in 1:n
        u = rk4_step(f, t[i], u, h)
        y[:, i + 1] .= u
    end
    return t, y
end

function main()
    t, y = rk4_system(f, TSPAN, Y₀, N_LEGACY)
    ref = reduce(hcat, [collect(exact(ti)) for ti in t])

    maximum(abs.(collect(exact(TSPAN[1])) .- collect(Y₀))) < 1e-14 ||
        error("the closed form does not satisfy the initial values")

    err = abs.(y .- ref)
    @printf("%d steps, h = %.2f\n", N_LEGACY, (TSPAN[2] - TSPAN[1]) / N_LEGACY)
    for k in 1:3
        @printf("y%d: max |error| = %.2e, relative to max |y%d| = %.1e;  y%d(8) = %.6f, exact %.6f\n",
            k, maximum(err[k, :]), k, maximum(err[k, :]) / maximum(abs, ref[k, :]),
            k, y[k, end], ref[k, end])
    end

    hs = (TSPAN[2] - TSPAN[1]) ./ N_SWEEP
    e_sweep = map(N_SWEEP) do n
        tt, yy = rk4_system(f, TSPAN, Y₀, n)
        maximum(abs.(yy .- reduce(hcat, [collect(exact(ti)) for ti in tt])))
    end
    order = convergence_order(hs, e_sweep)
    println("\n      h     max error   order")
    for i in eachindex(hs)
        @printf("%7.5f   %.3e   %s\n", hs[i], e_sweep[i],
            i == 1 ? "" : @sprintf("%.3f", order.pairwise[i - 1]))
    end
    @printf("order between the two finest steps: %.3f\n", last(order.pairwise))
    abs(last(order.pairwise) - 4) < 0.05 ||
        error("RK4 is not converging at fourth order: $(last(order.pairwise))")

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1300, 840))
    colours = (PALETTE.blue, PALETTE.orange, PALETTE.green)
    markers = (:circle, :rect, :diamond)
    # y₁ and y₂ coincide at the scale of their panel, so their RK4 points are
    # drawn on alternate grid nodes rather than on top of each other
    shown = (1:10:(N_LEGACY + 1), 6:10:(N_LEGACY + 1), 1:5:(N_LEGACY + 1))

    # y₁ and y₂ both go as −eᵗ and reach −3000, y₃ is a bounded oscillation, so
    # the two scales get a panel each over a shared abscissa.
    left = GridLayout(fig[2, 1])
    ax1 = Axis(left[1, 1], ylabel = L"y_1,\; y_2", xticks = 1:8)
    ax3 = Axis(left[2, 1], xlabel = L"t", ylabel = L"y_3", xticks = 1:8)
    for (k, ax) in zip(1:3, (ax1, ax1, ax3))
        lines!(ax, t, ref[k, :], color = colours[k], linestyle = k == 2 ? :dash : :solid)
        scatter!(ax, t[shown[k]], y[k, shown[k]], color = colours[k],
            marker = markers[k], markersize = MARKERSIZE.dense,)
    end
    hidexdecorations!(ax1, grid = false, ticks = false)
    linkxaxes!(ax1, ax3)
    xlims!(ax3, 0.8, 8.2)
    rowsize!(left, 2, Relative(0.4))

    right = GridLayout(fig[2, 2])
    ax_e = Axis(right[1, 1], xlabel = L"t", ylabel = "Absolute error",
        yscale = log10, xticks = 1:8, yticks = logticks(-9, -3; step = 2),)
    # the errors of y₂ and y₃ coincide, so the dashed one goes on top
    for k in (1, 3, 2)
        lines!(ax_e, t[2:end], err[k, 2:end], color = colours[k],
            linestyle = k == 2 ? :dash : :solid,)
    end
    xlims!(ax_e, 0.8, 8.2)
    ylims!(ax_e, 1e-9, 2e-3)

    guide = 4 * last(e_sweep) .* (hs ./ last(hs)) .^ 4
    ax_c = Axis(right[2, 1], xlabel = L"Step size $h$", ylabel = "Maximum error",
        xscale = log10, yscale = log10,
        xticks = (
            [0.01, 0.02, 0.05, 0.1, 0.2], [L"0.01", L"0.02", L"0.05", L"0.1", L"0.2"],),
        yticks = logticks(-7, -1; step = 2),)
    lines!(ax_c, hs, guide, color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    text!(ax_c, hs[end - 1], guide[end - 1]; text = L"\propto h^{4}",
        align = (:right, :bottom), offset = (-6, 4), fontsize = ANNOTATION_SIZE,)
    scatterlines!(ax_c, hs, e_sweep, color = PALETTE.black)
    text!(ax_c, 0.97, 0.06; text = @sprintf("Order %.3f", last(order.pairwise)),
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)
    xlims!(ax_c, 0.007, 0.36)

    component(k) = [LineElement(color = colours[k], linestyle = k == 2 ? :dash : :solid),
        MarkerElement(color = colours[k], marker = markers[k],
            markersize = MARKERSIZE.dense, strokewidth = 1.5, strokecolor = :black,),]
    Legend(fig[1, 1:2],
        [[component(k) for k in 1:3],
            [LineElement(color = PALETTE.black),
                MarkerElement(color = PALETTE.black, marker = :circle,
                    markersize = MARKERSIZE.dense,),],],
        [[L"y_1", L"y_2", L"y_3"], ["Exact", "RK4"]],
        ["Component", "Solution"]; titleposition = :left, nbanks = 1, groupgap = 40,)

    path = savefigure(fig, FIGURES, "ode_system_rk4")
    println("\nwrote ", path)
end

main()
