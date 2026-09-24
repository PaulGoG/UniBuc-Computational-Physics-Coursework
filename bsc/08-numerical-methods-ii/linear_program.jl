# A two-variable linear program, solved by vertex enumeration.
#
#   maximise  x₁ + x₂
#   subject to  x₁ + x₂ ≤ 3,  -x₁ + 3x₂ ≤ 1,  x₂ ≤ 3,  x₁, x₂ ≥ 0
#
# Ported from MaximizeLinearSystemEq.jl in Julia-Workflow-FFUB/Single_Files/ on the
# `legacy` branch, which called JuMP with GLPK through
# `Model(with_optimizer(GLPK.Optimizer))`. `with_optimizer` was deprecated in
# JuMP 0.21 and removed in 0.22 (2021), so the file raises UndefVarError on any
# current JuMP.
#
# It is replaced rather than repaired. For two variables and three constraints
# the feasible region has a handful of vertices and enumerating them is exact,
# dependency-free, and shows the structure the solver hides — which matters
# here, because:
#
#   * **the problem is degenerate.** The objective x₁ + x₂ is exactly parallel
#     to the first constraint x₁ + x₂ ≤ 3, so the optimal face is the whole
#     segment from (3, 0) to (2, 1), not a point. The original printed
#     `value(x1)` and `value(x2)` as though the answer were unique; which vertex
#     a simplex implementation returns is arbitrary.
#   * **the third constraint is redundant.** x₂ ≤ 3 is implied by x₁ + x₂ ≤ 3
#     together with x₁ ≥ 0, and never binds.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Constraints as (a₁, a₂, b) meaning a₁x₁ + a₂x₂ ≤ b, including the sign bounds."
const CONSTRAINTS = [(1.0, 1.0, 3.0), (-1.0, 3.0, 1.0), (0.0, 1.0, 3.0),
    (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0),]
objective(x₁, x₂) = x₁ + x₂

"Is the point inside every constraint, to within `tol`?"
feasible(p; tol = 1e-9) = all(c -> c[1]*p[1] + c[2]*p[2] <= c[3] + tol, CONSTRAINTS)

"""
    vertices()

All feasible intersections of constraint pairs — the vertices of the polytope.
"""
function vertices()
    pts = Tuple{Float64,Float64}[]
    for i in 1:length(CONSTRAINTS), j in (i + 1):length(CONSTRAINTS)

        (a₁, a₂, b₁) = CONSTRAINTS[i]
        (c₁, c₂, b₂) = CONSTRAINTS[j]
        det = a₁*c₂ - a₂*c₁
        abs(det) < 1e-12 && continue
        p = ((b₁*c₂ - a₂*b₂) / det, (a₁*b₂ - b₁*c₁) / det)
        feasible(p) &&
            !any(
                q -> isapprox(q[1], p[1]; atol = 1e-9) &&
                     isapprox(q[2], p[2]; atol = 1e-9), pts,) && push!(pts, p)
    end
    return pts
end

function main()
    V = vertices()
    values = [objective(p...) for p in V]
    best = maximum(values)
    optimal = V[findall(v -> isapprox(v, best; atol = 1e-9), values)]

    println("feasible vertices and objective:")
    for (p, v) in zip(V, values)
        @printf("  (%.4f, %.4f)  ->  %.4f%s\n", p[1], p[2], v,
            isapprox(v, best; atol = 1e-9) ? "   <- optimal" : "")
    end
    @printf("optimum = %.4f, attained at %d vertices\n", best, length(optimal))
    # the optimum is the bound of the first constraint, attained on its whole face
    length(optimal) == 2 && isapprox(best, CONSTRAINTS[1][3]; atol = 1e-12) ||
        error("expected the optimum $(CONSTRAINTS[1][3]) on a two-vertex face, got $best at $(length(optimal)) vertices")
    length(optimal) > 1 &&
        println("the optimal face is a segment: the solution is NOT unique")

    # is constraint 3 ever active at a feasible vertex?
    active3 = count(p -> isapprox(p[2], 3.0; atol = 1e-9), V)
    @printf("constraint x₂ ≤ 3 active at %d vertices -> redundant\n", active3)

    fig = Figure(size = (820, 820))
    ax = Axis(fig[2, 1], xlabel = L"x_1", ylabel = L"x_2", aspect = DataAspect())

    order = sortperm([atan(p[2] - 1, p[1] - 1) for p in V])
    l_feas = poly!(ax, Point2f[V[order]...], color = (PALETTE.sky, 0.35),
        strokewidth = 1.5, strokecolor = PALETTE.blue,)
    l_vert = scatter!(ax, first.(V), last.(V), color = PALETTE.blue)

    seg = sort(optimal, by = first)
    l_opt = lines!(ax, [seg[1][1], seg[end][1]], [seg[1][2], seg[end][2]],
        color = PALETTE.red, linewidth = 6,)

    # the objective contours are parallel to the binding constraint, which is
    # why the optimum is a face and not a vertex; drawn up to the optimum
    local l_obj
    for c in 0.5:0.5:3
        l_obj = lines!(ax, [0, 3.6], [c - 0, c - 3.6], color = (PALETTE.black, 0.35),
            linestyle = :dot, linewidth = GUIDE_WIDTH,)
    end

    # x₂ ≤ 3 is redundant; the axis reaches it so that can be seen
    l_red = hlines!(ax, [3.0], color = PALETTE.purple, linestyle = :dashdot,
        linewidth = GUIDE_WIDTH,)

    limits!(ax, -0.2, 3.6, -0.2, 3.4)
    text!(ax, 3.55, 2.94; text = rich(it("x"), subscript("2"), " ≤ 3, redundant"),
        align = (:right, :top), color = PALETTE.purple, fontsize = ANNOTATION_SIZE,)
    text!(ax, 0.97, 0.62;
        text = @sprintf("Optimum %.0f on the whole segment\n(3, 0) – (2, 1)", best),
        space = :relative, align = (:right, :top), color = PALETTE.red,
        fontsize = ANNOTATION_SIZE, justification = :right,)

    Legend(fig[1, 1], [l_feas, l_vert, l_opt, l_obj],
        ["Feasible region", "Vertices", "Optimal face",
            rich("Objective contours ", it("x"), subscript("1"),
                " + ", it("x"), subscript("2"),),];
        nbanks = 2,)
    path = savefigure(fig, FIGURES, "linear_program")
    println("wrote ", path)
end

main()
