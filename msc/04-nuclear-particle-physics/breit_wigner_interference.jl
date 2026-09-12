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

"Resonance parameters in MeV: (position, width)."
const RES1 = (E = 2300.0, Γ = 150.0)
const RES2 = (E = 2340.0, Γ = 320.0)

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
        a1 = amplitude.(E, Ref(RES1)); a2 = amplitude.(E, Ref(RES2))
        c1 = 1 / sqrt(sum((abs2.(a1)[1:end-1] .+ abs2.(a1)[2:end]) ./ 2 .* diff(E)))
        c2 = 1 / sqrt(sum((abs2.(a2)[1:end-1] .+ abs2.(a2)[2:end]) ./ 2 .* diff(E)))
        y = normalise(E, abs2.(c1 .* a1 .+ c2 .* a2 .* exp(im * φ)))
        w, pos = fwhm(E, y)
        @printf("φ = %3.0f°   interference peak %7.2f MeV   FWHM %6.2f MeV\n",
                rad2deg(φ), pos, w)
        (φ = φ, y = y, w = w, pos = pos)
    end
    @printf("\nneither resonance sits at %.0f or %.0f: the apparent position and width\n",
            RES1.E, RES2.E)
    @printf("of the observed line depend on the relative phase.\n")

    fig = Figure(size = (900, 480))
    ax = Axis(fig[2, 1], xlabel = L"$E$ [MeV]", ylabel = L"$P(E)$ [MeV$^{-1}$]")
    l1 = lines!(ax, E, B1, color = PALETTE.sky, linewidth = 1.4)
    l2 = lines!(ax, E, B2, color = PALETTE.orange, linewidth = 1.4)
    handles = [l1, l2]
    for (k, r) in enumerate(interference)
        push!(handles, lines!(ax, E, r.y,
              color = (PALETTE.blue, PALETTE.red)[k], linewidth = 1.8))
    end
    xlims!(ax, 1900, 2800)
    text!(ax, 0.97, 0.93;
        text = @sprintf("φ = 30°: peak %.0f, FWHM %.0f\nφ = 45°: peak %.0f, FWHM %.0f",
                        interference[1].pos, interference[1].w,
                        interference[2].pos, interference[2].w),
        space = :relative, align = (:right, :top), fontsize = 15)

    Legend(fig[1, 1], handles,
        [L"$|C_1B_1|^2$, $E_1 = 2300$, $\Gamma_1 = 150$",
         L"$|C_2B_2|^2$, $E_2 = 2340$, $\Gamma_2 = 320$",
         L"interference, $\varphi = 30°$", L"interference, $\varphi = 45°$"],
        orientation = :horizontal, framevisible = false, labelsize = 15,
        nbanks = 2, colgap = 18)
    rowsize!(fig.layout, 2, Relative(0.82))
    println("wrote ", savefigure(fig, FIGURES, "breit_wigner_interference"))
end

main()
