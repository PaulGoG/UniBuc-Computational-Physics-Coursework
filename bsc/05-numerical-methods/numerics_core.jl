# Helpers shared by the scripts of this module: one classical Runge–Kutta step
# and the reduction of a step-size sweep to an observed order of convergence.

"""
    rk4_step(f, t, u, h)

One classical fourth-order Runge–Kutta step of size `h` for `u′ = f(t, u)`,
returning `u(t + h)`. The state `u` is a number or a tuple of numbers, and `f`
returns the same shape.

# Example

```julia
rk4_step((t, y) -> -y, 0.0, 1.0, 0.1)            # ≈ exp(-0.1)
rk4_step((t, u) -> (u[2], -u[1]), 0.0, (1.0, 0.0), 0.1)
```
"""
function rk4_step(f, t, u, h)
    k₁ = f(t, u)
    k₂ = f(t + h / 2, u .+ (h / 2) .* k₁)
    k₃ = f(t + h / 2, u .+ (h / 2) .* k₂)
    k₄ = f(t + h, u .+ h .* k₃)
    return u .+ (h / 6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

"""
    convergence_order(h, err) -> (; order, stderr, pairwise)

Observed order of convergence from errors `err` measured at step sizes `h`.

`order` is the least-squares slope of `log err` against `log h` and `stderr` its
standard error from the residual variance on `length(h) - 2` degrees of freedom.
The residuals of a convergence sweep are not noise — they are the higher-order
terms of the error expansion — so `stderr` measures how far the sweep is from a
pure power law rather than a statistical uncertainty. `pairwise[i]` is the order
between `h[i]` and `h[i + 1]`, which is what shows whether the asymptotic regime
has been reached.

Throws an `ArgumentError` for fewer than three points, mismatched lengths, or a
non-positive entry, none of which has a logarithm worth fitting.
"""
function convergence_order(h::AbstractVector{<:Real}, err::AbstractVector{<:Real})
    length(h) == length(err) ||
        throw(ArgumentError("$(length(h)) step sizes against $(length(err)) errors"))
    length(h) >= 3 ||
        throw(ArgumentError("need at least three step sizes, got $(length(h))"))
    all(>(0), h) && all(>(0), err) ||
        throw(ArgumentError("step sizes and errors must be positive"))

    x, y = log.(h), log.(err)
    n = length(x)
    x̄, ȳ = sum(x) / n, sum(y) / n
    Sxx = sum(abs2, x .- x̄)
    order = sum((x .- x̄) .* (y .- ȳ)) / Sxx
    residual = y .- ȳ .- order .* (x .- x̄)
    stderr = sqrt(sum(abs2, residual) / (n - 2) / Sxx)
    pairwise = [(y[i] - y[i + 1]) / (x[i] - x[i + 1]) for i in 1:(n - 1)]
    return (; order, stderr, pairwise)
end
