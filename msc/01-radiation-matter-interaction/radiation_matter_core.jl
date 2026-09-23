# Linear least squares with a general data covariance, shared by
# `neutron_activation_halflives.jl` and `alpha_attenuation_mylar.jl`.

using LinearAlgebra, Printf

"""
    LinearFit

Result of [`generalised_least_squares`](@ref): the coefficients `p`, their
covariance matrix `cov`, the minimum `χ²` and the degrees of freedom `dof`.
"""
struct LinearFit
    p::Vector{Float64}
    cov::Matrix{Float64}
    χ²::Float64
    dof::Int
end

"""
    generalised_least_squares(M, y, C)

Minimise (y − Mp)ᵀ C⁻¹ (y − Mp) for the design matrix `M`, the data `y` and the
data covariance `C` (any `AbstractMatrix`; pass a `Diagonal` of variances for
independent errors). The coefficient covariance is (Mᵀ C⁻¹ M)⁻¹.

# Example

```julia
M = hcat(ones(3), [1.0, 2.0, 3.0])
fit = generalised_least_squares(M, [1.1, 1.9, 3.2], Diagonal(fill(0.01, 3)))
fit.p[2], sqrt(fit.cov[2, 2])      # slope and its standard error
```
"""
function generalised_least_squares(M::AbstractMatrix, y::AbstractVector, C::AbstractMatrix)
    size(M, 1) == length(y) ||
        throw(DimensionMismatch("M has $(size(M, 1)) rows, y has $(length(y)) entries"))
    size(C) == (length(y), length(y)) ||
        throw(DimensionMismatch("C must be $(length(y))×$(length(y)), got $(size(C))"))
    size(M, 1) > size(M, 2) ||
        throw(ArgumentError("need more data points than parameters"))
    W = inv(Matrix(C))
    cov = inv(Symmetric(M' * W * M))
    p = cov * (M' * W * y)
    r = y - M * p
    return LinearFit(p, Matrix(cov), dot(r, W * r), size(M, 1) - size(M, 2))
end

"""
    ordinary_least_squares(M, y)

Unweighted fit for data without individual uncertainties. The coefficient
covariance is scaled by the residual variance s² = Σr²/dof, which is returned as
`χ²/dof` of the unit-variance fit, so `fit.χ² / fit.dof` is s².
"""
function ordinary_least_squares(M::AbstractMatrix, y::AbstractVector)
    fit = generalised_least_squares(M, y, Diagonal(ones(length(y))))
    s² = fit.χ² / fit.dof
    return LinearFit(fit.p, s² .* fit.cov, fit.χ², fit.dof)
end

"""
    prediction_sigma(fit, v)

Standard error of the linear combination vᵀp, from the full quadratic form
vᵀ cov v — variances and covariances of the coefficients together.
"""
prediction_sigma(fit::LinearFit, v::AbstractVector) = sqrt(dot(v, fit.cov * v))

"""
    pm_string(x, σ; sig = 2)

`x ± σ` with σ rounded to `sig` significant digits and `x` to the same decimal
place, so a result never carries more digits than its uncertainty supports:
`pm_string(9.3017, 0.7207)` gives `"9.30 ± 0.72"`, `pm_string(375.19, 148.5)`
gives `"380 ± 150"`. The minus sign is the typographic `−`.
"""
function pm_string(x::Real, σ::Real; sig::Integer = 2)
    (isfinite(x) && isfinite(σ) && σ > 0) ||
        throw(ArgumentError("expected finite x and positive finite σ, got $x ± $σ"))
    decimals = sig - 1 - floor(Int, log10(σ))
    xr, σr = round(x; digits = decimals), round(σ; digits = decimals)
    fmt = Printf.Format("%." * string(max(decimals, 0)) * "f")
    return replace(Printf.format(fmt, xr) * " ± " * Printf.format(fmt, σr), "-" => "−")
end
