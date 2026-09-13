# Classical RK4 on a three-component linear system with a closed-form solution:
#
#   y₁' = y₂
#   y₂' = -y₁ - 2eᵗ + 1
#   y₃' = -y₁ - eᵗ + 1          on t ∈ [1, 8],  y(1) = (1, 2, 3)
#
# Components 1 and 2 form a driven linear oscillator, y₁'' + y₁ = 1 - 2eᵗ, whose
# solution is y₁ = C₁cos t + C₂sin t + 1 - eᵗ. Component 3 is a quadrature of
# it, y₃' = -C₁cos t - C₂sin t.
#
# Ported from Sys_ODE_RK4.m (Octave). The Butcher tableau was correct. The
# problems were structural: an orientation check written as
#
#     m = size(alpha);  if m==1  alpha = alpha';  endif
#
# where `size` returns [1 3], so `m==1` is [1 0] and the `if` — which requires
# every element non-zero — never fired; the code worked only because assigning
# into a column silently reoriented. The routine also called its right-hand side
# by name rather than taking it as an argument, so the "generic system solver"
# was welded to one problem, dumped a 101×4 matrix to the console twice through
# missing semicolons, demultiplexed its state by column-major linear indexing,
# and never compared against the exact solution that exists for this system.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Right-hand side of the system."
f(t, y) = (y[2], -y[1] - 2exp(t) + 1, -y[1] - exp(t) + 1)

# integration constants fixed by y(1) = (1, 2, 3)
const C₁ = exp(1) * cos(1) - (2 + exp(1)) * sin(1)
const C₂ = exp(1) * sin(1) + (2 + exp(1)) * cos(1)
const C₃ = 3 + C₁ * sin(1) - C₂ * cos(1)

"Exact solution, all three components."
exact(t) = (C₁*cos(t) + C₂*sin(t) + 1 - exp(t),
            -C₁*sin(t) + C₂*cos(t) - exp(t),
            -C₁*sin(t) + C₂*cos(t) + C₃)

"""
    rk4_system(f, tspan, y₀, n)

`n` RK4 steps on a system of any dimension. The right-hand side is an argument,
not a hardcoded name.
"""
function rk4_system(f, tspan, y₀, n)
    h = (tspan[2] - tspan[1]) / n
    t = range(tspan[1], tspan[2], length = n + 1)
    y = Matrix{Float64}(undef, length(y₀), n + 1)
    y[:, 1] .= y₀
    for i in 1:n
        u = Tuple(@view y[:, i])
        k₁ = f(t[i], u)
        k₂ = f(t[i] + h/2, u .+ h .* k₁ ./ 2)
        k₃ = f(t[i] + h/2, u .+ h .* k₂ ./ 2)
        k₄ = f(t[i] + h,   u .+ h .* k₃)
        y[:, i+1] .= u .+ (h/6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
    end
    return collect(t), y
end

function main()
    tspan, y₀ = (1.0, 8.0), (1.0, 2.0, 3.0)
    t, y = rk4_system(f, tspan, y₀, 100)

    ref = reduce(hcat, [collect(exact(ti)) for ti in t])
    err = abs.(y .- ref)
    for k in 1:3
        @printf("component y%d: max |error| = %.3e,  final %.6f vs exact %.6f\n",
                k, maximum(err[k, :]), y[k, end], ref[k, end])
    end

    orders = map((25, 50, 100, 200, 400)) do n
        _, yy = rk4_system(f, tspan, y₀, n)
        maximum(abs.(yy[:, end] .- collect(exact(tspan[2]))))
    end
    ns = [25, 50, 100, 200, 400]
    slope = (log(orders[1]) - log(orders[end])) / (log(ns[end]) - log(ns[1]))
    @printf("observed order = %.2f\n", slope)

    fig = Figure(size = (940, 440))
    colours = (PALETTE.blue, PALETTE.orange, PALETTE.green)

    # y₁ and y₂ are both dominated by -e^t and reach -3000; y₃ is a bounded
    # oscillation. On one axis the first two lie on top of each other and the
    # third is a flat line on zero, so the panel showed one curve where it
    # claimed three. Each gets the scale it needs, sharing the abscissa.
    top = GridLayout(fig[2, 1])
    ax1 = Axis(top[1, 1], ylabel = L"y_{1,2}(t)")
    handles = []
    # broad solid under narrow dashed: y₁ and y₂ agree to 0.15 % of their range
    # here, so drawn at equal width the second simply erases the first
    push!(handles, lines!(ax1, t, ref[1, :], color = colours[1], linewidth = 4.0))
    push!(handles, lines!(ax1, t, ref[2, :], color = colours[2], linewidth = 1.6,
        linestyle = :dash))
    hidexdecorations!(ax1, grid = false)
    text!(ax1, 0.04, 0.06;
        text = rich(it("y"), subscript("1"), " and ", it("y"), subscript("2"),
                    " are indistinguishable here:\nboth go as −", it("e"),
                    superscript(rich("t", font = :italic)),
                    ", and differ by the oscillatory part alone"),
        space = :relative, align = (:left, :bottom), fontsize = 14)

    ax3 = Axis(top[2, 1], xlabel = L"t", ylabel = L"y_3(t)")
    lines!(ax3, t, ref[3, :], color = colours[3], linewidth = 1.4)
    push!(handles, scatter!(ax3, t[1:5:end], y[3, 1:5:end],
          color = colours[3], markersize = MARKERSIZE.dense))
    linkxaxes!(ax1, ax3)
    rowsize!(top, 1, Relative(0.62))
    rowgap!(top, 6)

    ax2 = Axis(fig[2, 2], xlabel = L"t", ylabel = "Absolute error",
        yscale = log10, yticks = logticks(-8, -3; step = 2))
    for k in 1:3
        lines!(ax2, t[2:end], err[k, 2:end], color = colours[k], linewidth = 1.4)
    end
    text!(ax2, 0.04, 0.96;
        text = latexstring(@sprintf("\\text{Observed order } %.2f", slope)),
        space = :relative, align = (:left, :top), fontsize = 15)

    Legend(fig[1, 1:2], handles,
        [L"y_1", L"y_2", L"$y_3$, with the RK4 points"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)

    rowsize!(fig.layout, 2, Relative(0.85))
    path = savefigure(fig, FIGURES, "ode_system_rk4")
    println("wrote ", path)
end

main()
