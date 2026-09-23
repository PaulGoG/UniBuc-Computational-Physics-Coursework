# Annual frequency of an aircraft crash onto an installation of radius R sited a
# perpendicular distance y₀ from a straight air route, compared with the 10⁻⁵
# and 10⁻⁷ yr⁻¹ screening levels.
#
#   F = N · P · A · ∫ f(r(x), φ(x)) dx                           [yr⁻¹]
#
#   N   flights per year on the route                            [yr⁻¹]
#   P   probability of loss of control per km flown              [km⁻¹]
#   A   target area, πR²                                         [km²]
#   f   probability per unit area that an aircraft which lost
#       control at the point x of the route comes down at the
#       site, a distance r = √(x² + y₀²) away and at an angle φ
#       from its heading                                         [km⁻²]
#
# The site is small against every length in f, so f is taken constant over it.
#
# Ported from Frecventa_accident_aviatic.jl, which sets P = 10⁻⁹ km⁻¹,
# N = 7·10⁴ yr⁻¹, g = 0.23 km⁻¹, R = 50 m, y₀ from 5 to 50 km and a route length
# x₀ = 200 km, and integrates two kernels along the route from 0 to x₀:
#
#     integrandRadial          = (x/s) e^{-gs}            s = √(x² + y₀²)
#     integrandRadialUnghiular = (x y₀/s²) e^{-gs}        times g/2
#
# With the aircraft at x > 0 flying towards the foot of the perpendicular,
# x/s = cos φ and y₀/s = sin φ. Both kernels carry the forward factor cos φ and
# both are set to zero for x ≤ 0 by an explicit guard: an aircraft comes down
# ahead of the point where control was lost, never behind it. The integration
# over one side of the route is therefore part of the model, and it is kept.
#
# What the original got wrong is dimensional. The first integral is a length
# and the second a pure number, so after the common prefactor P·N·πR² neither
# was a frequency and the two differed by a length, yet both were drawn on one
# axis against thresholds in yr⁻¹. Each kernel is made a probability density
# per unit area here, changing as little as that requires:
#
#   * radial–angular kernel. The g and the 1/2 of the original are the
#     normalisations of an exponential density g e^{-gr} in r and of cos φ / 2
#     over the forward half-plane. The density per unit area with those two
#     marginals is f = g e^{-gr} cos φ / (2r); the original has sin φ = y₀/r
#     where 1/r belongs, so its value is y₀ times a frequency. This is the
#     primary result.
#   * radial kernel. cos φ e^{-gr}, normalised over the forward half-plane as a
#     density per unit area, is f = (g²/2) cos φ e^{-gr}. For a long route it
#     integrates to (g/2) e^{-g y₀}, the double-exponential airway formula.
#
# The isotropic kernel f = g e^{-gr} / (2πr), integrated over both sides of the
# route, is kept as an alternative: it is a different model, not a correction
# of either original.
#
# Every route integral has a closed form and is asserted against it.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, QuadGK, SpecialFunctions
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Probability of loss of control per km of flight [km⁻¹]."
const P_LOSS = 1e-9
"Flights per year on the route [yr⁻¹]."
const N_FLIGHTS = 7e4
"Decay constant g of the crash-distance density [km⁻¹]."
const G = 0.23
"Length x₀ of the approach leg integrated over [km]."
const ROUTE_LENGTH = 200.0
"Radius of the installation in the original [m]."
const SITE_RADIUS = 50.0
"Range of perpendicular distances y₀ between route and site [km]."
const Y0_RANGE = (5.0, 50.0)
"Screening levels [yr⁻¹]."
const THRESHOLDS = (1e-5, 1e-7)
"Relative tolerance of the comparisons between quadrature and closed form."
const RTOL = 1e-8

"""
    forward_cosine_density(r, cosφ)

Crash-location probability per unit area [km⁻²]: exponential density `g e^{-gr}`
in the distance `r` [km], angular density `cos φ / 2` about the heading over the
forward half-plane, zero behind the aircraft.
"""
forward_cosine_density(r, cosφ) = cosφ <= 0 ? 0.0 : G * exp(-G * r) * cosφ / (2r)

"""
    forward_exponential_density(r, cosφ)

Crash-location probability per unit area [km⁻²] of shape `cos φ e^{-gr}`,
normalised over the forward half-plane: `(g²/2) cos φ e^{-gr}`.
"""
forward_exponential_density(r, cosφ) = cosφ <= 0 ? 0.0 : G^2 / 2 * cosφ * exp(-G * r)

"Isotropic crash-location probability per unit area [km⁻²], `g e^{-gr} / (2πr)`."
isotropic_density(r) = G * exp(-G * r) / (2π * r)

"""
    route_integral(density, y₀; both_sides = false)

`∫ f dx` along the route [km⁻¹] for a site at perpendicular distance `y₀` [km]:
over the approach leg `0 < x < x₀` for a forward kernel, over `−x₀ < x < x₀`
when `both_sides` is set. `density` takes `(r, cos φ)`.
"""
function route_integral(density, y₀; both_sides = false)
    integrand(x) = (r = hypot(x, y₀); density(r, x / r))
    a = both_sides ? -ROUTE_LENGTH : 0.0
    I, _ = quadgk(integrand, a, 0.0, ROUTE_LENGTH; rtol = 1e-10)
    return I
end

"Largest distance between the approach leg and the site [km]."
farthest(y₀) = hypot(ROUTE_LENGTH, y₀)

"Closed form of the forward-cosine route integral, `(g/2)[E₁(g y₀) − E₁(g s_max)]`."
forward_cosine_closed(y₀) = G / 2 * (expint(G * y₀) - expint(G * farthest(y₀)))

"Closed form of the forward-exponential route integral, `(g/2)[e^{-g y₀} − e^{-g s_max}]`."
forward_exponential_closed(y₀) = G / 2 * (exp(-G * y₀) - exp(-G * farthest(y₀)))

"Closed form of the isotropic integral over an unbounded route, `(g/π) K₀(g y₀)`."
isotropic_closed(y₀) = G / π * besselk(0, G * y₀)

"Target area [km²] of a disc of radius `R_m` metres; no skid or shadow term."
target_area(R_m) = π * (R_m * 1e-3)^2

"Annual impact frequency [yr⁻¹] from a route integral [km⁻¹]."
frequency(route, R_m) = N_FLIGHTS * P_LOSS * target_area(R_m) * route

"The two estimators exactly as `Frecventa_accident_aviatic.jl` computes them."
function original_estimators(y₀, R_m)
    radial(x) = x > 0 ? x / hypot(x, y₀) * exp(-G * hypot(x, y₀)) : 0.0
    angular(x) = x > 0 ? x * y₀ / (x^2 + y₀^2) * exp(-G * hypot(x, y₀)) : 0.0
    prefactor = P_LOSS * N_FLIGHTS * π * (R_m * 1e-3)^2
    I_r = first(quadgk(radial, 0.0, ROUTE_LENGTH; rtol = 1e-10))
    I_a = G / 2 * first(quadgk(angular, 0.0, ROUTE_LENGTH; rtol = 1e-10))
    return prefactor * I_r, prefactor * I_a
end

function check(name, value, reference)
    isapprox(value, reference; rtol = RTOL) ||
        throw(ErrorException("$name: $value differs from $reference by more than $RTOL"))
    return nothing
end

"Assert normalisation, closed forms and the relation to the original estimators."
function validate()
    # each kernel integrates to one over the plane; the radial and angular parts separate
    radial_cos = first(quadgk(r -> G * exp(-G * r), 0, Inf))
    radial_exp = first(quadgk(r -> G^2 * r * exp(-G * r), 0, Inf))
    angular = first(quadgk(φ -> cos(φ) / 2, -π / 2, π / 2))
    check("normalisation, forward cosine", radial_cos * angular, 1.0)
    check("normalisation, forward exponential", radial_exp * angular, 1.0)
    check("normalisation, isotropic", radial_cos, 1.0)

    for y₀ in range(Y0_RANGE..., length = 10)
        I_cos = route_integral(forward_cosine_density, y₀)
        I_exp = route_integral(forward_exponential_density, y₀)
        I_iso = route_integral((r, _) -> isotropic_density(r), y₀; both_sides = true)
        check("forward cosine at y₀ = $y₀", I_cos, forward_cosine_closed(y₀))
        check("forward exponential at y₀ = $y₀", I_exp, forward_exponential_closed(y₀))
        check("isotropic at y₀ = $y₀", I_iso, isotropic_closed(y₀))

        f_r, f_a = original_estimators(y₀, SITE_RADIUS)
        check("original radial–angular estimator / y₀", f_a / y₀,
            frequency(I_cos, SITE_RADIUS),)
        check("original radial estimator × g²/2", f_r * G^2 / 2,
            frequency(I_exp, SITE_RADIUS),)
    end
    return nothing
end

function main()
    validate()
    println("normalisation, closed forms and the relation to the original: asserted to ",
        RTOL,)
    @printf("N·P = %.2e km⁻¹ yr⁻¹,  g = %.2f km⁻¹,  x₀ = %.0f km\n\n", N_FLIGHTS * P_LOSS,
        G, ROUTE_LENGTH)

    println("F [yr⁻¹] at R = $(Int(SITE_RADIUS)) m")
    println("  y₀ [km]   forward cosine   forward exponential   isotropic     original f_R(r,θ)   original f_R(r)")
    for y₀ in (5.0, 10.0, 25.0, 50.0)
        f_r, f_a = original_estimators(y₀, SITE_RADIUS)
        @printf("  %5.1f     %.3e        %.3e             %.3e     %.3e           %.3e\n",
            y₀,
            frequency(forward_cosine_closed(y₀), SITE_RADIUS),
            frequency(forward_exponential_closed(y₀), SITE_RADIUS),
            frequency(isotropic_closed(y₀), SITE_RADIUS), f_a, f_r)
    end
    println("  (the two original columns are not frequencies: [km yr⁻¹] and [km² yr⁻¹])")

    y_min = Y0_RANGE[1]
    println("\nforward cosine kernel at y₀ = $(Int(y_min)) km:")
    for R in (50.0, 200.0)
        @printf("  R = %3.0f m   F = %.3e yr⁻¹\n", R,
            frequency(forward_cosine_closed(y_min), R))
    end
    for τ in THRESHOLDS
        R_crit = 1e3 * sqrt(τ / (π * N_FLIGHTS * P_LOSS * forward_cosine_closed(y_min)))
        @printf("  F = %.0e yr⁻¹ is reached at R = %.0f m\n", τ, R_crit)
    end

    y₀s = range(Y0_RANGE..., length = 181)
    radii = [25.0, 50.0, 100.0, 200.0, 400.0]
    colours = (PALETTE.sky, PALETTE.blue, PALETTE.green, PALETTE.orange, PALETTE.purple)
    styles = (:dot, :solid, :dash, :dashdot, :dashdotdot)

    # right padding: the last x tick label is centred on the frame's edge
    fig = Figure(size = (900, 800), figure_padding = (10, 28, 10, 10))
    ax = Axis(fig[2, 1],
        ylabel = rich("Impact frequency ", it("F"), " [yr", superscript("−1"), "]"),
        yscale = log10, yticks = logticks(-13, -5; step = 2), xticklabelsvisible = false,)
    curves = [lines!(ax, y₀s, frequency.(forward_cosine_closed.(y₀s), R);
                  color = colours[k], linestyle = styles[k],)
              for (k, R) in enumerate(radii)]
    hlines!(ax, collect(THRESHOLDS); color = PALETTE.red, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    for (τ, exponent) in zip(THRESHOLDS, ("−5", "−7"))
        text!(ax, Y0_RANGE[2] - 0.5, τ * 1.35;
            text = rich("Screening level 10", superscript(exponent), " yr", superscript("−1")),
            align = (:right, :bottom), color = PALETTE.red, fontsize = ANNOTATION_SIZE,)
    end
    ylims!(ax, 8e-15, 1.5e-4)

    ax_ratio = Axis(fig[3, 1],
        xlabel = rich("Distance from the route ", it("y"), subscript("0"), " [km]"),
        ylabel = "Ratio", yscale = log10, yticks = ([1, 2, 5, 10], ["1", "2", "5", "10"]),)
    primary = forward_cosine_closed.(y₀s)
    l_exp = lines!(ax_ratio, y₀s, forward_exponential_closed.(y₀s) ./ primary;
        color = PALETTE.black, linestyle = :dash,)
    l_iso = lines!(ax_ratio, y₀s, isotropic_closed.(y₀s) ./ primary;
        color = PALETTE.black, linestyle = :dot,)
    hlines!(
        ax_ratio, [1.0]; color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH,)
    text!(ax_ratio, Y0_RANGE[2] - 0.5, 1.05;
        text = "Forward cosine", align = (:right, :bottom),
        fontsize = ANNOTATION_SIZE,)
    ylims!(ax_ratio, 0.8, 18)
    linkxaxes!(ax, ax_ratio)
    xlims!(ax_ratio, Y0_RANGE...)
    rowsize!(fig.layout, 3, Relative(0.26))

    Legend(fig[1, 1], [curves, [l_exp, l_iso]],
        [[rich("$(Int(R)) m") for R in radii], ["Forward exponential", "Isotropic"]],
        [rich("Site radius ", it("R")), "Kernel, relative to forward cosine"];
        nbanks = 2, titleposition = :top,)
    println("\nwrote ", savefigure(fig, FIGURES, "aircraft_crash_frequency"))
end

main()
