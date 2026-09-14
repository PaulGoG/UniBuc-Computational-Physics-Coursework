# Coherent interference of two Breit-Wigner resonances sharing the same quantum
# numbers,
#
#   B_k(E) = (Γ_k/2) / (E_k - E - iΓ_k/2),
#   dσ/dε ∝ |C₁B₁ + C₂B₂ e^{iφ}|²,
#
# for relative phases φ = 30° and 45°. Because the quantum numbers coincide the
# amplitudes add before squaring, so the observed line shape is neither of the
# two resonances nor their sum, and its apparent position and width depend on φ.
#
# The two resonances, the two phases and the equal weights all come from the
# assignment (Fenomenologia particulelor elementare, Tema 1), which sets them
# per student. The row for this one gives f₂(2300) at 2297 ± 60 MeV with
# Γ = 149 ± 40 MeV and f₂(2340) at 2339 ± 60 MeV with Γ = 319 ± 70 MeV.
#
# `Breit_Wigner.jl` used 2300/150 and 2340/320 instead — the resonances' names
# and round numbers rather than the measured values it was given. The difference
# is small against widths of 149 and 319 MeV, and it moves the interference peak
# by a couple of MeV, but the uncertainties are the point: at ±60 MeV on each
# position and ±40 to ±70 MeV on each width, the apparent peak position this
# calculation produces is not determined to the 0.01 MeV its arithmetic
# suggests. The assigned values are used here and the uncertainties are carried
# so that they can be said out loud.
#
# Ported from Breit_Wigner.jl. The physics was right. The measurement of the
# resulting peak was not: the full width at half maximum was located by scanning
# for points satisfying
#
#     minimum(y)*0.005 >= abs(maximum(y)/2 - i)
#
# an absolute tolerance keyed to the minimum of the array, which finds nothing
# at all if the sampling happens to straddle the half-maximum and silently
# returns whatever the first and last near-misses were. It is replaced here by
# linear interpolation between the bracketing samples, which is exact to the
# grid spacing and cannot fail to find a crossing that exists.
#
# The struct holding the distributions also had four untyped fields, making
# every access `Any`, and the annotation coordinates were hardcoded in data
# units.

using Printf, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"""
Resonance parameters in MeV: position, width, and the quoted uncertainty on
each.

These are the values the assignment sets for this student, not the rounded
figures implied by the resonances' names. f₂(2300) and f₂(2340) are both 2⁺, so
the amplitudes add coherently, and the assignment specifies equal generation
weights — which is where C₁ = C₂ comes from.
"""
const RES1 = (E = 2297.0, σE = 60.0, Γ = 149.0, σΓ = 40.0)
const RES2 = (E = 2339.0, σE = 60.0, Γ = 319.0, σΓ = 70.0)

"Breit-Wigner amplitude."
amplitude(E, r) = (0.5 * r.Γ) / (r.E - E - 0.5im * r.Γ)

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
    left === nothing && return (NaN, x[peak])
    right === nothing && return (NaN, x[peak])
    right += peak - 1
    xl = x[left] + (half - y[left]) * (x[left+1] - x[left]) / (y[left+1] - y[left])
    xr = x[right-1] + (half - y[right-1]) * (x[right] - x[right-1]) / (y[right] - y[right-1])
    return (xr - xl, x[peak])
end

"Normalise a sampled distribution to unit integral."
normalise(x, y) = y ./ (sum((y[1:end-1] .+ y[2:end]) ./ 2 .* diff(x)))

"""
    interference_shape(E, r1, r2, φ)

Normalised `|C₁B₁ + C₂B₂e^{iφ}|²` on the grid `E`, with each amplitude
separately normalised so that the two enter with equal weight, as the
assignment specifies.
"""
function interference_shape(E, r1, r2, φ)
    a1 = amplitude.(E, Ref(r1)); a2 = amplitude.(E, Ref(r2))
    c1 = 1 / sqrt(sum((abs2.(a1)[1:end-1] .+ abs2.(a1)[2:end]) ./ 2 .* diff(E)))
    c2 = 1 / sqrt(sum((abs2.(a2)[1:end-1] .+ abs2.(a2)[2:end]) ./ 2 .* diff(E)))
    return normalise(E, abs2.(c1 .* a1 .+ c2 .* a2 .* exp(im * φ)))
end

"""
    peak_uncertainty(E, φ)

Spread of the interference peak position and width when each of the four
resonance parameters is moved by its quoted uncertainty in turn.

The arithmetic reports a peak to 0.01 MeV from inputs known to ±60 MeV, so the
sensitivity is worth measuring rather than leaving implied. Each parameter is
moved by ±σ on its own and the **largest** single excursion is returned, which
is an envelope rather than a propagated error: the assignment quotes no
covariances, and with only four parameters the envelope is the honest summary.
"""
function peak_uncertainty(E, φ)
    base_w, base_pos = fwhm(E, interference_shape(E, RES1, RES2, φ))
    δpos = 0.0; δw = 0.0
    for (field, σfield) in ((:E, :σE), (:Γ, :σΓ)), which in 1:2
        for sign in (+1, -1)
            r1, r2 = RES1, RES2
            r = which == 1 ? r1 : r2
            shifted = merge(r, NamedTuple{(field,)}((getfield(r, field) + sign * getfield(r, σfield),)))
            which == 1 ? (r1 = shifted) : (r2 = shifted)
            w, pos = fwhm(E, interference_shape(E, r1, r2, φ))
            δpos = max(δpos, abs(pos - base_pos))
            δw = max(δw, abs(w - base_w))
        end
    end
    return δpos, δw
end

function main()
    E = collect(1600.0:0.05:3000.0)

    B1 = normalise(E, abs2.(amplitude.(E, Ref(RES1))))
    B2 = normalise(E, abs2.(amplitude.(E, Ref(RES2))))
    for (name, r, y) in (("resonance 1", RES1, B1), ("resonance 2", RES2, B2))
        w, pos = fwhm(E, y)
        @printf("%-12s  input E = %6.1f Γ = %5.1f   measured peak %7.2f  FWHM %6.2f\n",
                name, r.E, r.Γ, pos, w)
    end

    println()
    interference = map((π/6, π/4)) do φ
        y = interference_shape(E, RES1, RES2, φ)
        w, pos = fwhm(E, y)
        δpos, δw = peak_uncertainty(E, φ)
        @printf("φ = %3.0f°   interference peak %7.2f ± %.0f MeV   FWHM %6.2f ± %.0f MeV\n",
                rad2deg(φ), pos, δpos, w, δw)
        (φ = φ, y = y, w = w, pos = pos, δpos = δpos, δw = δw)
    end
    @printf("\nthe two phases move the peak by %.1f MeV, well inside the %.0f MeV that\n",
            abs(interference[2].pos - interference[1].pos), interference[1].δpos)
    @printf("the quoted resonance parameters allow: the phase dependence is real but\n")
    @printf("this data cannot resolve it.\n")
    @printf("\nneither resonance sits at %.0f or %.0f: the apparent position and width\n",
            RES1.E, RES2.E)
    @printf("of the observed line depend on the relative phase.\n")

    fig = Figure(size = (900, 480))
    ax = Axis(fig[2, 1], xlabel = L"$E$ [MeV]", ylabel = L"$P(E)$ [MeV$^{-1}$]")
    # Isolated resonances dashed, their interference solid: the two are not
    # variants of one another, and drawn in one hue family each pair read as if
    # they were.
    l1 = lines!(ax, E, B1, color = PALETTE.sky, linewidth = 1.6, linestyle = :dash)
    l2 = lines!(ax, E, B2, color = PALETTE.orange, linewidth = 1.6, linestyle = :dash)
    handles = [l1, l2]
    for (k, r) in enumerate(interference)
        push!(handles, lines!(ax, E, r.y,
              color = (PALETTE.blue, PALETTE.red)[k], linewidth = 2.0))
    end
    xlims!(ax, 1900, 2800)
    # Each uncertainty stays with the quantity it belongs to: a single "± 61 MeV"
    # on its own line read as the uncertainty of the FWHM above it.
    for (k, r) in enumerate(interference)
        text!(ax, 0.97, 0.95 - 0.075 * (k - 1);
            text = rich(it("φ"), @sprintf(" = %.0f°: peak %.0f ± %.0f MeV,  FWHM %.0f ± %.0f MeV",
                                          rad2deg(r.φ), r.pos, r.δpos, r.w, r.δw)),
            space = :relative, align = (:right, :top), fontsize = 15,
            color = (PALETTE.blue, PALETTE.red)[k])
    end
    text!(ax, 0.97, 0.80;
        text = "the quoted resonance parameters allow that much movement,\nso the phase dependence is real but not resolvable here",
        space = :relative, align = (:right, :top), fontsize = 14,
        justification = :right)

    Legend(fig[1, 1], handles,
        [rich("|", it("C"), subscript("1"), it("B"), subscript("1"), "|",
              superscript("2"), ": ", it("f"), subscript("2"), "(2300), ",
              it("E"), subscript("1"), " = 2297 MeV, Γ", subscript("1"), " = 149 MeV"),
         rich("|", it("C"), subscript("2"), it("B"), subscript("2"), "|",
              superscript("2"), ": ", it("f"), subscript("2"), "(2340), ",
              it("E"), subscript("2"), " = 2339 MeV, Γ", subscript("2"), " = 319 MeV"),
         rich("Interference, ", it("φ"), " = 30°"),
         rich("Interference, ", it("φ"), " = 45°")],
        orientation = :horizontal, framevisible = false, labelsize = 15,
        nbanks = 2, colgap = 18)
    rowsize!(fig.layout, 2, Relative(0.82))
    println("wrote ", savefigure(fig, FIGURES, "breit_wigner_interference"))
end

main()
