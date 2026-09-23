# Paraxial (first-order, Gaussian) ray-transfer analysis of two classical
# optical designs: the Erfle wide-field eyepiece and the Rudolph Tessar
# photographic objective.
#
# The ray state is the column vector [n·u ; y] — reduced angle and height — and
# the system matrix is accumulated from the last surface backwards, so that
# S = M_N T_{N-1} M_{N-1} … T_1 M_1 acts as v_out = S·v_in.
#
# Ported from Erlfe.m and Tessar.m in Code_Archive/Optics/ on the `legacy`
# branch. The matrix convention, the accumulation order and the index
# bookkeeping there are correct. The principal-plane separation is written
# `sum(d) - abs(zH1) - abs(zH2)`, right only when both principal planes fall on
# one side; the files' own diagrams place H1 at -zH1 from the first vertex and
# H2 at sum(d)+zH2 from the last, so the separation is sum(d) + zH1 + zH2, which
# agrees with the original formula for these two prescriptions and would not
# for a design with a principal plane on the other side. `f = -1/S(1,2)` is
# labelled "Convergenta sistemului" there; that is the focal length, and the
# convergence is P = 1/f = -S(1,2). "Erlfe" is Erfle, the 1921 wide-field
# eyepiece.
#
# `main` asserts the engine on a thin lens against the lensmaker formula, the
# unit determinant of each system matrix (air on both sides), and Newton's
# conjugate relation (z₁ − z_F1)(z₂ − z_F2) = f², object distances counted
# to the left of V₁ and image distances to the right of V_N, for one object
# position of
# each design.
#
# Both designs are strictly paraxial here, as they were originally: no real ray
# trace, no aberrations, and a single refractive index per glass. The last point
# matters, because the cemented doublet in the Tessar and the cemented doublet
# and triplet in the Erfle exist precisely to achromatise, and without dispersion
# data that cannot be evaluated at all.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Relative tolerance of the closed-form checks."
const CHECK_TOLERANCE = 1e-9

"Refraction at a single surface of radius `R` between media `n` and `n′`."
refraction(n, n′, R) = [1.0 (n - n′)/R; 0.0 1.0]

"Transfer over reduced distance `d/n`."
transfer(d, n) = [1.0 0.0; d/n 1.0]

"""
    system_matrix(R, n, d)

Ray-transfer matrix of a centred system of `length(R)` surfaces, `length(R)+1`
media and `length(R)-1` axial separations, accumulated from the last surface
backwards.
"""
function system_matrix(R, n, d)
    length(n) == length(R) + 1 ||
        throw(DimensionMismatch("$(length(R)) surfaces need $(length(R)+1) media, got $(length(n))"))
    length(d) == length(R) - 1 ||
        throw(DimensionMismatch("$(length(R)) surfaces need $(length(R)-1) gaps, got $(length(d))"))

    S = Matrix{Float64}(I, 2, 2)
    for i in eachindex(R)
        k = length(R) - i + 1
        S = S * refraction(n[k], n[k + 1], R[k])
        k > 1 && (S = S * transfer(d[k - 1], n[k]))
    end
    return S
end

"""
    cardinal(S, d)

Cardinal elements of a system with matrix `S` and axial thicknesses `d`.

Returns focal length, front and back focal distances measured from the first and
last vertex, the two principal-plane positions, and their separation. The
separation is `sum(d) + zH1 + zH2`, consistent with placing H₁ at `-zH1` from V₁
and H₂ at `sum(d)+zH2` from V₂.
"""
function cardinal(S, d)
    A, B, C, D = S[1, 1], S[1, 2], S[2, 1], S[2, 2]
    f = -1 / B
    return (f = f,
        power = -B,
        zf1 = -A / B,
        zf2 = -D / B,
        zH1 = (1 - A) / B,
        zH2 = (1 - D) / B,
        interstice = sum(d) + (1 - A)/B + (1 - D)/B,)
end

"""
    conjugate(S, z₁)

Image position and transverse magnification for an object at `z₁` before the
first vertex.
"""
function conjugate(S, z₁)
    A, B, C, D = S[1, 1], S[1, 2], S[2, 1], S[2, 2]
    denom = A + z₁ * B
    abs(denom) < 1e-12 && return (z₂ = Inf, magnification = Inf)
    z₂ = -(z₁ * D + C) / denom
    return (z₂ = z₂, magnification = z₂ * B + D)
end

# ---- the two prescriptions, distances in mm --------------------------------

"Erfle Type II wide-field eyepiece: six elements in three groups, nine surfaces."
const ERFLE = (
    name = "Erfle eyepiece",
    R = [-62.89, 42.74, -42.74, 76.92, -81.30, 36.36, -36.36, 32.26, Inf],
    n = [1.0, 1.720, 1.638, 1.0, 1.638, 1.0, 1.638, 1.649, 1.638, 1.0],
    d = [5.4, 12.0, 0.5, 8.0, 0.5, 11.5, 1.8, 5.5],
)

"Rudolph Tessar photographic objective, 1902: four elements in three groups."
const TESSAR = (
    name = "Tessar objective",
    R = [16.28, -275.70, -34.57, 15.82, Inf, 19.20, -24.00],
    n = [1.0, 1.6116, 1.0, 1.6053, 1.0, 1.5123, 1.6116, 1.0],
    d = [3.57, 1.89, 0.81, 3.25, 2.17, 3.96],
)

"""
    check_thin_lens()

Sanity check of the engine: a single thin lens of power P must give
A = D = 1, B = -P, C = 0, hence zf1 = zf2 = f and both principal planes at the
vertex.
"""
function check_thin_lens()
    n_glass, R1, R2 = 1.5, 100.0, -100.0
    P = (n_glass - 1) * (1/R1 - 1/R2)          # lensmaker, zero thickness
    S = system_matrix([R1, R2], [1.0, n_glass, 1.0], [0.0])
    c = cardinal(S, [0.0])
    @printf("thin-lens check: f = %.6f mm, lensmaker 1/P = %.6f mm, zH1 = %.2e, zH2 = %.2e\n",
        c.f, 1/P, c.zH1, c.zH2)
    return isapprox(c.f, 1/P; rtol = 1e-10) && abs(c.zH1) < 1e-9 && abs(c.zH2) < 1e-9
end

function main()
    check_thin_lens() || error("paraxial engine failed its thin-lens check")

    results = map((ERFLE, TESSAR)) do sys
        S = system_matrix(sys.R, sys.n, sys.d)
        c = cardinal(S, sys.d)
        abs(det(S) - 1) < CHECK_TOLERANCE ||
            error("$(sys.name): det S = $(det(S)), not 1 with air on both sides")
        z₁ = 2 * c.f                       # an object at twice the focal length
        z₂ = conjugate(S, z₁).z₂
        isapprox((z₁ - c.zf1) * (z₂ - c.zf2), c.f^2; rtol = CHECK_TOLERANCE) ||
            error("$(sys.name): Newton's relation fails at z₁ = $z₁")
        @printf("\n%s  (%d surfaces, Σd = %.2f mm)\n", sys.name, length(sys.R), sum(sys.d))
        @printf("  focal length f       = %+9.4f mm\n", c.f)
        @printf("  power      P = 1/f   = %+9.6f mm⁻¹\n", c.power)
        @printf("  front focal zf1      = %+9.4f mm from V1\n", c.zf1)
        @printf("  back  focal zf2      = %+9.4f mm from V%d\n", c.zf2, length(sys.R))
        @printf("  principal   zH1, zH2 = %+9.4f, %+9.4f mm\n", c.zH1, c.zH2)
        @printf("  interstice H1H2      = %+9.4f mm   (original formula gave %+9.4f)\n",
            c.interstice, sum(sys.d) - abs(c.zH1) - abs(c.zH2))
        @printf("  det S = %.12f;  Newton (z₁ − z_F1)(z₂ − z_F2) = f² holds at z₁ = 2f\n",
            det(S))
        (sys = sys, S = S, c = c)
    end

    fig = Figure(size = (1200, 720))
    for (row, r) in enumerate(results)
        Σd = sum(r.sys.d)
        span = maximum(abs, [r.c.zf1, r.c.zf2, r.c.zH1, r.c.zH2, Σd]) * 1.25

        # no ordinate: only the axial coordinate means anything
        ax = Axis(fig[row, 1], xlabel = row == 2 ? "Axial position [mm]" : "",
            ylabel = r.sys.name, yticksvisible = false, yticklabelsvisible = false,
            ygridvisible = false,)
        hlines!(ax, [0.0], color = PALETTE.black, linewidth = GUIDE_WIDTH)

        # the glass block, from the first vertex to the last
        poly!(ax, Point2f[(0, -0.5), (Σd, -0.5), (Σd, 0.5), (0, 0.5)],
            color = (PALETTE.sky, 0.25),)
        text!(ax, Σd / 2, -0.62;
            text = rich(
                "Glass, ", it("V"), subscript("1"), " to ", it("V"), subscript("2"),),
            color = PALETTE.sky, align = (:center, :top), fontsize = ANNOTATION_SIZE,)

        # H₁ and H₂ can sit within a millimetre of each other: their labels are
        # staggered vertically
        cardinal(sym, k) = rich(it(sym), subscript(k))
        for (x, sym, k, col, dy) in ((0.0, "V", "1", PALETTE.black, 0.62),
            (Σd, "V", "2", PALETTE.black, 0.62),
            (-r.c.zf1, "F", "1", PALETTE.blue, 0.62),
            (Σd + r.c.zf2, "F", "2", PALETTE.blue, 0.62),
            (-r.c.zH1, "H", "1", PALETTE.green, 1.1),
            (Σd + r.c.zH2, "H", "2", PALETTE.green, 0.62))
            scatter!(ax, [x], [0.0], color = col, markersize = MARKERSIZE.emphasis)
            text!(ax, x, dy; text = cardinal(sym, k), color = col,
                align = (:center, :bottom), fontsize = ANNOTATION_SIZE,)
        end
        xlims!(ax, -span, Σd + span)
        ylims!(ax, -1.45, 1.75)
        text!(ax, 0.99, 0.05;
            text = rich(it("f"), @sprintf(" = %.2f mm,   ", r.c.f), it("H"),
                subscript("1"), it("H"), subscript("2"),
                @sprintf(" = %.2f mm", r.c.interstice)),
            space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)
    end

    path = savefigure(fig, FIGURES, "paraxial_systems")
    println("\nwrote ", path)
end

main()
