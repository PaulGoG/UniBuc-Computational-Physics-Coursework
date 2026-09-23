# Coherent interference of two Breit–Wigner resonances sharing the same quantum
# numbers,
#
#   B_k(E) = (Γ_k/2) / (E_k − E − iΓ_k/2),
#   dσ/dε ∝ |C₁B₁ + C₂B₂ e^{iφ}|²,
#
# for relative phases φ = 30° and 45°, each amplitude normalised to unit
# integral over the window so that the two enter with equal weight. The
# amplitudes add before squaring, so the observed line is neither resonance nor
# their sum, and its apparent position and width depend on φ.
#
# Ported from Breit_Wigner.jl in Julia-Workflow-FFUB/FPECA_M_2/ on the `legacy`
# branch, whose resonance parameters, 2300/150 and 2340/320 MeV, are kept. The
# physics there was right. The width was measured by scanning for samples with
#
#     minimum(y)*0.005 >= abs(maximum(y)/2 - i)
#
# an absolute tolerance keyed to the minimum of the array, which finds nothing
# if no sample lands within it and otherwise returns the first and last
# near-misses; linear interpolation between the bracketing samples replaces it.
# `main` asserts each normalisation against the closed form of ∫|B_k|² over the
# window, (Γ/2)[arctan(2(b − E_k)/Γ) − arctan(2(a − E_k)/Γ)], and the measured
# width of each isolated resonance against its Γ.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Resonance position and width, MeV."
const RES1 = (name = "f₂(2300)", E = 2300.0, Γ = 150.0)
const RES2 = (name = "f₂(2340)", E = 2340.0, Γ = 320.0)
"Lower and upper end of the energy window, MeV."
const E_MIN = 1600.0
const E_MAX = 3000.0
"Grid step, MeV."
const ΔE = 0.05
"Relative phases compared."
const PHASES = (π / 6, π / 4)
"Relative tolerance of the closed-form normalisation check."
const CHECK_TOLERANCE = 1e-6

"Breit–Wigner amplitude."
amplitude(E, r) = (0.5 * r.Γ) / (r.E - E - 0.5im * r.Γ)

"∫ₐᵇ |B|² dE in closed form."
lorentzian_integral(r, a, b) = 0.5 * r.Γ *
                               (atan(2 * (b - r.E) / r.Γ) - atan(2 * (a - r.E) / r.Γ))

"Trapezoidal integral of `y` on the grid `x`."
trapezoid(x, y) = sum((y[1:(end - 1)] .+ y[2:end]) ./ 2 .* diff(x))

"""
    fwhm(x, y)

Full width at half maximum by linear interpolation between the samples that
bracket each half-maximum crossing. Returns the width and the peak position.
"""
function fwhm(x, y)
    peak = argmax(y)
    half = y[peak] / 2
    left = findlast(i -> y[i] < half, 1:peak)
    right = findfirst(i -> y[i] < half, peak:length(y))
    (left === nothing || right === nothing) && return (NaN, x[peak])
    right += peak - 1
    xl = x[left] + (half - y[left]) * (x[left + 1] - x[left]) / (y[left + 1] - y[left])
    xr = x[right - 1] +
         (half - y[right - 1]) * (x[right] - x[right - 1]) / (y[right] - y[right - 1])
    return (xr - xl, x[peak])
end

"""
    interference_shape(E, r1, r2, φ)

Normalised |C₁B₁ + C₂B₂e^{iφ}|² on the grid `E`, each amplitude normalised to
unit integral over the grid.
"""
function interference_shape(E, r1, r2, φ)
    a1, a2 = amplitude.(E, Ref(r1)), amplitude.(E, Ref(r2))
    c1, c2 = 1 / sqrt(trapezoid(E, abs2.(a1))), 1 / sqrt(trapezoid(E, abs2.(a2)))
    y = abs2.(c1 .* a1 .+ c2 .* a2 .* exp(im * φ))
    return y ./ trapezoid(E, y)
end

function main()
    E = collect(E_MIN:ΔE:E_MAX)

    isolated = map((RES1, RES2)) do r
        y = abs2.(amplitude.(E, Ref(r)))
        Z = trapezoid(E, y)
        isapprox(Z, lorentzian_integral(r, E_MIN, E_MAX); rtol = CHECK_TOLERANCE) ||
            error("$(r.name): ∫|B|² = $Z against the closed form $(lorentzian_integral(r, E_MIN, E_MAX))")
        w, pos = fwhm(E, y ./ Z)
        (abs(pos - r.E) <= ΔE && abs(w - r.Γ) <= 2ΔE) ||
            error("$(r.name): measured peak $pos and FWHM $w against E = $(r.E), Γ = $(r.Γ)")
        @printf("%s: E = %.0f, Γ = %.0f MeV; ∫|B|² over the window %.6f, closed form %.6f; measured peak %.2f, FWHM %.2f MeV\n",
            r.name, r.E, r.Γ, Z, lorentzian_integral(r, E_MIN, E_MAX), pos, w)
        y ./ Z
    end

    println()
    interference = map(PHASES) do φ
        y = interference_shape(E, RES1, RES2, φ)
        w, pos = fwhm(E, y)
        @printf("φ = %2.0f°: interference peak %.1f MeV, FWHM %.1f MeV\n", rad2deg(φ), pos,
            w)
        (; φ, y, w, pos)
    end
    @printf("\nbetween the two phases the peak moves by %.1f MeV and the width by %.1f MeV; the line sits at neither %.0f nor %.0f MeV\n",
        abs(interference[2].pos - interference[1].pos), abs(interference[2].w -
                                                            interference[1].w),
        RES1.E, RES2.E)

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1200, 660))
    ax = Axis(fig[2, 1], xlabel = L"$E$ [MeV]", ylabel = L"$P(E)$ [MeV$^{-1}$]")
    # isolated resonances dashed, their interference solid
    l1 = lines!(ax, E, isolated[1], color = PALETTE.sky, linestyle = :dash)
    l2 = lines!(ax, E, isolated[2], color = PALETTE.orange, linestyle = :dash)
    colours = (PALETTE.blue, PALETTE.red)
    li = [lines!(ax, E, r.y, color = colours[k]) for (k, r) in enumerate(interference)]
    top = maximum(maximum, (isolated..., (r.y for r in interference)...))
    limits!(ax, 1900, 2800, 0, 1.3 * top)
    for (k, r) in enumerate(interference)
        text!(ax, 0.97, 0.96 - 0.08 * (k - 1);
            text = rich(it("φ"), @sprintf(" = %.0f°: peak %.0f MeV, FWHM %.0f MeV",
                rad2deg(r.φ), r.pos, r.w)),
            space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
            color = colours[k],)
    end

    label(k, r) = rich(
        "|", it("C"), subscript(k), it("B"), subscript(k), "|", superscript("2"),
        ": ", it("E"), subscript(k), @sprintf(" = %.0f MeV, Γ", r.E), subscript(k),
        @sprintf(" = %.0f MeV", r.Γ))
    Legend(fig[1, 1], [l1, l2, li...],
        [label("1", RES1), label("2", RES2), rich("Interference, ", it("φ"), " = 30°"),
            rich("Interference, ", it("φ"), " = 45°"),]; nbanks = 2,)
    println("wrote ", savefigure(fig, FIGURES, "breit_wigner_interference"))
end

main()
