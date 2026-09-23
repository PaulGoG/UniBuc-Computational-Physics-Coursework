# Hubble parameter from the redshifts of three galaxies, measured on a
# laboratory spectrograph against five reference lines.
#
#   z = (λ - λ₀)/λ₀,   v = cz,   v = H₀ d
#
# Distances are obtained from the apparent size s of each galaxy relative to a
# standard of known distance: d = d_std · s_std / s.
#
# Ported from Hubble.jl in Julia-Workflow-FFUB/FPECA_M_2/ on the `legacy`
# branch. The correction that matters is a one-character typo:
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
#
# The fit weights each galaxy by the effective variance σ_v² + H₀²σ_d² of its
# residual, iterated in H₀. The distance of the reference galaxy is the input
# D_REF and carries no error; its apparent size s_ref scales the other two
# distances by the same factor, so the error of s_ref is not independent
# between them and enters as a common-mode term, obtained by refitting with
# s_ref moved by ±σ, rather than being folded into each σ_d.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

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
"Published values compared with, km s⁻¹ Mpc⁻¹: Planck 2018 (doi 10.1051/0004-6361/201833910) and SH0ES 2022 (doi 10.3847/2041-8213/ac5c5b)."
const H0_PLANCK = (67.4, 0.5)
const H0_SHOES = (73.04, 1.04)
"Iterations of the effective-variance weights."
const FIT_ITERATIONS = 50

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

"""
    fit_hubble(d, v, σd, σv)

Weighted fit of v = H₀d through the origin with the effective variance
σ_v² + H₀²σ_d² of each residual, iterated. Returns H₀ and its formal error.
"""
function fit_hubble(d, v, σd, σv)
    H₀ = sum(v .* d) / sum(abs2, d)
    σH₀ = 0.0
    for _ in 1:FIT_ITERATIONS
        w = 1 ./ (σv .^ 2 .+ (H₀ .* σd) .^ 2)
        H₀ = sum(w .* d .* v) / sum(w .* d .^ 2)
        σH₀ = sqrt(1 / sum(w .* d .^ 2))
    end
    return H₀, σH₀
end

"""
    distances(s_ref)

Distances of the three galaxies for a reference apparent size `s_ref`, with the
error from each galaxy's own size only; the reference galaxy sits at D_REF by
definition and its own size cancels out of its distance.
"""
function distances(s_ref)
    d = [D_REF * s_ref / g[3] for g in GALAXIES]
    σd = [k == 1 ? 0.0 : d[k] * g[4] / g[3] for (k, g) in enumerate(GALAXIES)]
    return d, σd
end

function main()
    names = [g[1] for g in GALAXIES]
    v = Float64[]
    σv = Float64[]
    for (name, λs, _, _) in GALAXIES
        z, σz = mean_redshift(λs, λ₀, σλ)
        push!(v, C_KM_S * z)
        push!(σv, C_KM_S * σz)
    end
    s_ref, σs_ref = GALAXIES[1][3], GALAXIES[1][4]
    d, σd = distances(s_ref)
    for k in eachindex(GALAXIES)
        @printf("%-9s  v = %7.1f ± %5.1f km/s   d = %5.2f ± %4.2f Mpc\n",
            names[k], v[k], σv[k], d[k], σd[k])
    end

    H₀, σ_fit = fit_hubble(d, v, σd, σv)
    d_up, σd_up = distances(s_ref + σs_ref)
    d_down, σd_down = distances(s_ref - σs_ref)
    H_up, _ = fit_hubble(d_up, v, σd_up, σv)
    H_down, _ = fit_hubble(d_down, v, σd_down, σv)
    σ_common = abs(H_up - H_down) / 2
    σH₀ = hypot(σ_fit, σ_common)
    H_velocity, σ_velocity = fit_hubble(d, v, zero(σd), σv)
    @printf("\nH₀ = %.1f ± %.1f km s⁻¹ Mpc⁻¹: fit %.1f, reference size ±%.0f %% common to two galaxies %.1f\n",
        H₀, σH₀, σ_fit, 100σs_ref / s_ref, σ_common)
    @printf("velocity errors alone would give %.1f ± %.1f\n", H_velocity, σ_velocity)
    @printf("Planck 2018: %.1f ± %.1f   SH0ES 2022: %.2f ± %.2f\n", H0_PLANCK...,
        H0_SHOES...)
    @printf("Hubble time 1/H₀ = %.1f Gyr\n", 977.8 / H₀)
    # the two distant galaxies carry the fit, so the common-mode term is the
    # scale error H₀ σ_s,ref/s_ref itself, to within the asymmetry of 1/s; the
    # reference galaxy's velocity error is too large to anchor it
    isapprox(σ_common, H₀ * σs_ref / s_ref; rtol = 0.05) ||
        error("common-mode term $σ_common is not the scale error $(H₀ * σs_ref / s_ref)")

    fig = Figure(size = (1000, 660))
    ax = Axis(fig[2, 1], xlabel = "Distance [Mpc]", ylabel = L"Recession velocity [km s$^{-1}$]")
    xf = range(0, maximum(d .+ σd) * 1.04, length = 50)
    band!(ax, xf, (H₀ - σH₀) .* xf, (H₀ + σH₀) .* xf, color = (PALETTE.blue, 0.2))
    l_fit = lines!(ax, xf, H₀ .* xf, color = PALETTE.blue, linestyle = :dash)
    errorbars_unstroked!(ax, d, v, σv, color = PALETTE.orange, whiskerwidth = 10,
        linewidth = GUIDE_WIDTH,)
    errorbars_unstroked!(
        ax, d, v, σd, direction = :x, color = PALETTE.orange, whiskerwidth = 10,
        linewidth = GUIDE_WIDTH,)
    l_dat = scatter!(ax, d, v, color = PALETTE.orange)
    for (n, x, y) in zip(names, d, v)
        text!(ax, x + 0.8, y - 1.2 * maximum(σv); text = n, align = (:left, :top),
            fontsize = ANNOTATION_SIZE,)
    end
    text!(ax, 0.04, 0.95;
        text = rich(
            rich(it("H"), subscript("0"), @sprintf(" = %.1f ± %.1f km s", H₀, σH₀),
                superscript("−1"), " Mpc", superscript("−1"), color = PALETTE.blue,),
            @sprintf("\nfit %.1f, common reference size %.1f", σ_fit, σ_common),
            "\nhorizontal bars: each galaxy's own size error",),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1], [l_dat, l_fit], ["Galaxies", L"$v = H_0 d$, 1σ band"])
    println("wrote ", savefigure(fig, FIGURES, "hubble_parameter"))
end

main()
