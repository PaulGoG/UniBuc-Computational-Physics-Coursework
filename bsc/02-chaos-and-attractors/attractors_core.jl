# Vector fields, parameters and the integrator shared by the two scripts in this
# directory.
#
#   Lorenz    ẋ = σ(y - x),  ẏ = x(ρ - z) - y,  ż = xy - βz
#   Rössler   ẋ = -y - z,    ẏ = x + ay,        ż = b + z(x - c)

"Lorenz parameters, the canonical chaotic set. The 2018 code wrote β = 2.66."
const LORENZ = (σ = 10.0, ρ = 28.0, β = 8 / 3)
"Rössler parameters."
const ROSSLER = (a = 0.2, b = 0.2, c = 5.7)

lorenz(u, p) = (p.σ * (u[2] - u[1]),
    u[1] * (p.ρ - u[3]) - u[2],
    u[1] * u[2] - p.β * u[3],)

rossler(u, p) = (-u[2] - u[3],
    u[1] + p.a * u[2],
    p.b + u[3] * (u[1] - p.c),)

"""
    rk4_step(f, u, h, p)

One classical RK4 step of the autonomous three-component system `f`.
"""
function rk4_step(f, u, h, p)
    k₁ = f(u, p)
    k₂ = f(u .+ h .* k₁ ./ 2, p)
    k₃ = f(u .+ h .* k₂ ./ 2, p)
    k₄ = f(u .+ h .* k₃, p)
    return u .+ (h / 6) .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

"Non-trivial fixed points of the Lorenz system."
function lorenz_fixed_points(p)
    s = sqrt(p.β * (p.ρ - 1))
    return ((s, s, p.ρ - 1), (-s, -s, p.ρ - 1))
end
