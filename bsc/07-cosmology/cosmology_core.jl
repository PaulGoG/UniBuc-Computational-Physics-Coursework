# Flat Friedmann–Lemaître cosmology shared by the two scripts in this directory:
#
#   H²(z) = H₀² [Ω_m (1+z)³ + Ω_Λ],      Ω_Λ = 1 − Ω_m,
#
# its luminosity distance, the magnitude–redshift fit built on it, and the
# closed-form expansion history the numerical one is checked against.

using LsqFit, QuadGK, Statistics

"Speed of light [km s⁻¹], exact by the SI definition of the metre."
const C_KM_S = 299_792.458

"Kilometres in a megaparsec: 1 pc = 648000/π au (IAU 2015 B2), 1 au = 149 597 870.7 km (IAU 2012 B2)."
const KM_PER_MPC = 648_000 / π * 149_597_870.7 * 1e6

"Seconds in 10⁹ Julian years of 365.25 days."
const S_PER_GYR = 365.25 * 86_400 * 1e9

"""
Planck 2018 base-ΛCDM parameters (TT,TE,EE+lowE+lensing), H₀ in km s⁻¹ Mpc⁻¹:
Planck Collaboration VI, A&A 641, A6 (2020), doi:10.1051/0004-6361/201833910.
"""
const PLANCK2018 = (H₀ = 67.4, σH₀ = 0.5, Ωm = 0.315, σΩm = 0.007)

"""
Physical range of Ω_m in a flat universe of matter and a non-negative Λ. A fit
that ends on either edge has not determined Ω_m.
"""
const Ωm_BOUNDS = (0.0, 1.0)

"`colour` mixed 35 % towards black: the stroke of a marker filled with `colour`."
function darker(colour)
    c = CairoMakie.Makie.to_color(colour)
    return RGBAf(0.65 * c.r, 0.65 * c.g, 0.65 * c.b, c.alpha)
end

"""
    unstroked_errorbars!(ax, args...; kwargs...)

`errorbars!` with unstroked whisker caps. Makie draws the caps as a scatter plot,
which otherwise takes the marker stroke of the shared theme and turns black.
"""
function unstroked_errorbars!(ax, args...; kwargs...)
    bars = errorbars!(ax, args...; linewidth = GUIDE_WIDTH, whiskerwidth = 7, kwargs...)
    for child in bars.plots
        child isa Scatter && (child.strokewidth = 0)
    end
    return bars
end

"Hubble time 1/H₀ in Gyr, for H₀ in km s⁻¹ Mpc⁻¹; 13.97 Gyr at H₀ = 70."
hubble_time(H₀) = KM_PER_MPC / H₀ / S_PER_GYR

"Dimensionless expansion rate E(z) = H(z)/H₀ of a flat matter + Λ universe."
expansion_rate(z, Ωm) = sqrt(Ωm * (1 + z)^3 + (1 - Ωm))

"""
    luminosity_distance(z, Ωm)

Luminosity distance of a flat universe in units of the Hubble distance c/H₀,
`(1+z) ∫₀^z dz′/E(z′)`, by adaptive Gauss–Kronrod quadrature.
"""
function luminosity_distance(z, Ωm)
    z > 0 || throw(ArgumentError("redshift must be positive, got $z"))
    integral, _ = quadgk(x -> 1 / expansion_rate(x, Ωm), 0, z; rtol = 1e-10)
    return (1 + z) * integral
end

"""
    magnitude(z, Ωm, offset)

`5 log₁₀` of the dimensionless luminosity distance, plus `offset`. For apparent
magnitudes of a standard candle the offset is M + 25 + 5 log₁₀(c/H₀ / Mpc), in
which M and H₀ cannot be separated; for distance moduli M = 0 and the offset
fixes H₀ (see [`hubble_constant`](@ref)).
"""
magnitude(z, Ωm, offset) = 5 * log10(luminosity_distance(z, Ωm)) + offset

"H₀ in km s⁻¹ Mpc⁻¹ from the offset of a distance-modulus fit, μ = 5 log₁₀(d_L/Mpc) + 25."
hubble_constant(offset) = C_KM_S / 10^((offset - 25) / 5)

"""
    fit_flat_lcdm(z, m[, σ]; Ωm = nothing)

Least-squares fit of `magnitude(z, Ωm, offset)` to the magnitudes or distance
moduli `m`. With `σ` the fit is weighted by 1/σ² and the standard errors are
those of the weighted problem, `√diag (JᵀWJ)⁻¹`; without it the errors are
scaled by the residual variance. Passing `Ωm` holds it fixed, which leaves the
offset as a weighted mean.

Returns `(; Ωm, offset, σΩm, σoffset, χ², dof, rms, at_bound)`. `χ²` is the
weighted sum of squared residuals, or the plain sum without `σ`.

Ω_m is confined to [`Ωm_BOUNDS`](@ref). If it ends on an edge, `at_bound` is
set, the offset is recomputed with Ω_m held there — the bounded
Levenberg–Marquardt iteration stops short of that conditional optimum — and
`σΩm` is `NaN`, a standard error having no meaning on a boundary.
"""
function fit_flat_lcdm(z, m, σ = nothing; Ωm = nothing)
    length(z) == length(m) || throw(DimensionMismatch("z and m differ in length"))
    w = σ === nothing ? ones(length(z)) : 1 ./ σ .^ 2
    Ωm === nothing || return _fit_offset(z, m, w, Ωm, σ === nothing)

    model(zz, p) = [magnitude(zi, p[2], p[1]) for zi in zz]
    p₀ = [mean(m .- magnitude.(z, 0.3, 0.0)), 0.3]
    bounds = (lower = [-Inf, Ωm_BOUNDS[1]], upper = [Inf, Ωm_BOUNDS[2]])
    fit = σ === nothing ? curve_fit(model, z, m, p₀; bounds...) :
          curve_fit(model, z, m, w, p₀; bounds...)
    Ωm_fit = fit.param[2]
    edge = findfirst(b -> abs(Ωm_fit - b) < 1e-6, Ωm_BOUNDS)
    if edge !== nothing
        held = _fit_offset(z, m, w, Ωm_BOUNDS[edge], σ === nothing)
        return (; held..., dof = length(z) - 2, at_bound = true)
    end
    σ_offset, σ_Ωm = stderror(fit)
    r = m .- model(z, fit.param)
    return (; Ωm = Ωm_fit, offset = fit.param[1], σΩm = σ_Ωm, σoffset = σ_offset,
        χ² = sum(w .* r .^ 2), dof = length(z) - 2, rms = sqrt(mean(abs2, r)),
        at_bound = false,)
end

function _fit_offset(z, m, w, Ωm, scale_by_residuals)
    shape = magnitude.(z, Ωm, 0.0)
    offset = sum(w .* (m .- shape)) / sum(w)
    r = m .- shape .- offset
    χ² = sum(w .* r .^ 2)
    dof = length(z) - 1
    σ_offset = sqrt((scale_by_residuals ? χ² / dof : 1.0) / sum(w))
    return (; Ωm = float(Ωm), offset, σΩm = NaN, σoffset = σ_offset, χ², dof,
        rms = sqrt(mean(abs2, r)), at_bound = false,)
end

"""
    age(Ωm)

Age of a flat matter + Λ universe in units of 1/H₀:
`(2 / 3√Ω_Λ) asinh √(Ω_Λ/Ω_m)`, which tends to the Einstein–de Sitter 2/3 as
Ω_Λ → 0.
"""
function age(Ωm)
    0 < Ωm <= 1 || throw(ArgumentError("Ωm must lie in (0, 1], got $Ωm"))
    ΩΛ = 1 - Ωm
    ΩΛ == 0 && return 2 / 3
    return 2 / (3 * sqrt(ΩΛ)) * asinh(sqrt(ΩΛ / Ωm))
end

"""
    cosmic_time(a, Ωm)

Time since the Big Bang at scale factor `a` in units of 1/H₀, the inverse of
[`scale_factor`](@ref): `(2 / 3√Ω_Λ) asinh(√(Ω_Λ/Ω_m) a^{3/2})`, and
`(2/3) a^{3/2}` for Ω_Λ = 0. `cosmic_time(1, Ωm)` is [`age`](@ref).
"""
function cosmic_time(a, Ωm)
    a >= 0 || throw(ArgumentError("scale factor must be non-negative, got $a"))
    ΩΛ = 1 - Ωm
    ΩΛ == 0 && return 2 / 3 * a^1.5
    return 2 / (3 * sqrt(ΩΛ)) * asinh(sqrt(ΩΛ / Ωm) * a^1.5)
end

"""
    scale_factor(τ, Ωm)

Closed-form scale factor of a flat matter + Λ universe a time `τ ≥ 0` after the
Big Bang, in units of 1/H₀, normalised to a = 1 at `τ = age(Ωm)`:

    a(τ) = (Ω_m/Ω_Λ)^{1/3} sinh^{2/3}(3√Ω_Λ τ / 2),     a(τ) = (3τ/2)^{2/3} for Ω_Λ = 0.
"""
function scale_factor(τ, Ωm)
    τ >= 0 || throw(ArgumentError("time since the Big Bang must be non-negative, got $τ"))
    ΩΛ = 1 - Ωm
    ΩΛ == 0 && return (3τ / 2)^(2 / 3)
    return cbrt(Ωm / ΩΛ) * sinh(3 * sqrt(ΩΛ) * τ / 2)^(2 / 3)
end
