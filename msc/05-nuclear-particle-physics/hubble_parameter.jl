# Hubble parameter from the redshifts of three galaxies, measured on a
# laboratory spectrograph against five reference lines.
#
#   z = (λ - λ₀)/λ₀,   v = cz,   v = H₀ d
#
# Distances are obtained from the apparent size s of each galaxy relative to a
# standard of known distance: d = d_std · s_std / s.
#
# Ported from Hubble.jl. The correction that matters is a one-character typo:
#
#     Suma_σ² =+ σ_z[i]^2
#
# which Julia parses as `Suma_σ² = +σ_z[i]^2`, an assignment, not `+=`. The
# running sum of variances was therefore overwritten on every iteration and kept
# only the last term, so every averaged redshift uncertainty was wrong — too
# small by roughly √n. The weights in the H₀ fit derive from those, so the fit
# was affected too.
#
# Also corrected: the mean was divided by the count of valid entries while the
# variance sum was divided by that count *outside* the square root, and no
# uncertainty was reported on the individual line redshifts.

using Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

const C_KM_S = 2.99792458e5

"Rest wavelengths in Å: Ca-K, Ca-H, H-α, H-β, H-δ."
const λ₀ = [3933.7, 3968.5, 6562.8, 4861.3, 4101.7]
"Wavelength reading uncertainty in Å."
const σλ = 3.3 / 2

"(name, observed wavelengths, apparent size in cm, size uncertainty)"
const GALAXIES = [
    ("NGC 3627", [3941.67, 3980.0, 6578.34, 4873.33, 4113.3], 5.0, 0.75),
    ("NGC 3368", [3944.995, 3980.0, NaN, NaN, 4116.33], 3.5, 0.80),
    ("NGC 3147", [3970.0, 4006.66, 6620.0, 4916.66, 4136.66], 1.8, 0.50),
]
"Reference distance for NGC 3627, in Mpc."
const D_REF = 10.73

"Redshift and its uncertainty for one line."
redshift(λ, λ₀, σλ) = ((λ - λ₀) / λ₀, sqrt(σλ^2 * (λ^2 + λ₀^2) / λ₀^4))

"""
    mean_redshift(λs, λ₀s, σλ)

Inverse-variance mean over the usable lines, with the variance summed correctly
— the defect the original's `=+` introduced.
"""
function mean_redshift(λs, λ₀s, σλ)
    zs = Float64[]
    σs = Float64[]
    for (λ, λr) in zip(λs, λ₀s)
        isnan(λ) && continue
        z, σz = redshift(λ, λr, σλ)
        push!(zs, z)
        push!(σs, σz)
    end
    w = 1 ./ σs .^ 2
    return sum(w .* zs) / sum(w), sqrt(1 / sum(w))
end

function main()
    d = Float64[]
    σd = Float64[]
    v = Float64[]
    σv = Float64[]
    names = String[]
    s_ref, σs_ref = GALAXIES[1][3], GALAXIES[1][4]

    for (name, λs, s, σs) in GALAXIES
        z, σz = mean_redshift(λs, λ₀, σλ)
        dist = D_REF * s_ref / s
        σdist = dist * sqrt((σs_ref / s_ref)^2 + (σs / s)^2)
        push!(names, name)
        push!(d, dist)
        push!(σd, σdist)
        push!(v, C_KM_S * z)
        push!(σv, C_KM_S * σz)
        @printf("%-9s  z = %.5f ± %.5f   v = %7.1f ± %5.1f km/s   d = %5.2f ± %4.2f Mpc\n",
            name, z, σz, C_KM_S * z, C_KM_S * σz, dist, σdist)
    end

    # Weighted fit through the origin, v = H₀d. Both coordinates carry error,
    # so the effective variance of each residual is σ_v² + H₀²σ_d² — which
    # depends on H₀ and is therefore iterated. Here the distance term dominates
    # and omitting it would understate σ_H₀ several-fold.
    H₀ = sum(v .* d) / sum(abs2, d)
    local σH₀ = 0.0
    for _ in 1:50
        w = 1 ./ (σv .^ 2 .+ (H₀ .* σd) .^ 2)
        H₀ = sum(w .* d .* v) / sum(w .* d .^ 2)
        σH₀ = sqrt(1 / sum(w .* d .^ 2))
    end
    @printf("\nH₀ = %.1f ± %.1f km s⁻¹ Mpc⁻¹\n", H₀, σH₀)
    @printf("Planck 2018: 67.4 ± 0.5   SH0ES: 73.0 ± 1.0\n")
    @printf("Hubble time 1/H₀ = %.2f Gyr\n", 977.8 / H₀)

    fig = Figure(size = (780, 500))
    ax = Axis(fig[2, 1], xlabel = "Distance [Mpc]", ylabel = L"Recession velocity [km s$^{-1}$]")
    xf = range(0, maximum(d) * 1.15, length = 50)
    l_fit = lines!(ax, xf, H₀ .* xf, color = PALETTE.blue, linewidth = 1.6)
    band!(ax, xf, (H₀ - σH₀) .* xf, (H₀ + σH₀) .* xf, color = (PALETTE.blue, 0.18))
    errorbars!(ax, d, v, σv, color = PALETTE.orange, whiskerwidth = 10)
    errorbars!(ax, d, v, σd, direction = :x, color = PALETTE.orange, whiskerwidth = 10)
    l_dat = scatter!(ax, d, v, color = PALETTE.orange, markersize = MARKERSIZE.data)
    for (n, x, y) in zip(names, d, v)
        text!(ax, x, y - 3 * maximum(σv) / 4; text = n,
            align = (:center, :top), fontsize = 14,)
    end
    text!(ax, 0.04, 0.93;
        # rich rather than the Unicode subscript: MathTeXEngine draws U+2080 as a
        # baseline letter o, so this read "Ho = 73 ± 11"
        text = rich(it("H"), subscript("0"),
            @sprintf(" = %.1f ± %.1f km s", H₀, σH₀), superscript("−1"),
            " Mpc", superscript("−1"),),
        space = :relative, align = (:left, :top), color = PALETTE.blue, fontsize = 16,)

    Legend(fig[1, 1], [l_dat, l_fit], ["Galaxies", L"$v = H_0 d$"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24,)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("wrote ", savefigure(fig, FIGURES, "hubble_parameter"))
end

main()
