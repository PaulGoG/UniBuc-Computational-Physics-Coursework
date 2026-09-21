# Annual frequency of an aircraft crash onto a nuclear installation sited at
# perpendicular distance y₀ from an air corridor — an external-hazard screening
# calculation, to be compared against the 10⁻⁵ and 10⁻⁷ yr⁻¹ design-basis
# thresholds.
#
# **Provenance.** This was not a course assignment. It is "modelarea matematică
# (cu scop pedagogic) a unei analize de risc prin metoda statistică, în speță
# studiul statistic al riscului prăbușirii unei aeronave de pasageri asupra unui
# amplasament nuclear folosind distribuții Poisson" — a deliberately pedagogical
# Poisson risk model, written in late 2023 as professional rather than academic
# work.
#
# That matters for how the result should be read. The quantity computed is a
# Poisson rate, so the probability of at least one impact in t years is
# 1 − e^{−Ft}; for the frequencies here, F t ≪ 1 and the two coincide. The
# values of N, P, g and y₀ are not recorded in anything that survives, so they
# remain as the original file set them and the absolute frequency should be read as
# an exercise, not as a screening result for any real site.
#
# The formulation below follows the standard structure of DOE-STD-3014 and
# NUREG-0800 §3.5.1.6, which is how a calculation of this shape is normally
# presented; the original did not cite either:
#
#   F = N · P · ∫ f(s(x)) dx · A_eff        [yr⁻¹]
#
#   N      flights per year on the route                      [yr⁻¹]
#   P      probability of loss of control per km of flight    [km⁻¹]
#   f(s)   crash-location probability density per unit AREA   [km⁻²]
#   A_eff  effective target area of the site                  [km²]
#
# so that yr⁻¹ · km⁻¹ · km⁻² · km · km² = yr⁻¹. Every factor carries its units
# and they cancel to a frequency by construction.
#
# The crash-location density is the key object. If the aircraft, having lost
# control at a point on the route, comes down at radial distance r with the
# exponential density p(r) = g e^{-gr} per unit r, and the direction is
# isotropic, then the density per unit *area* at that distance is
#
#   f(r) = p(r) / (2πr) = g e^{-gr} / (2πr)   [km⁻²]
#
# It is the 1/(2πr) that converts a radial density into an areal one, and its
# absence is what made the original version dimensionally unsound.
#
# Ported from Frecventa_accident_aviatic.jl, which computed two estimators:
#
#     integrandRadial          = (x/s) e^{-gs}
#     integrandRadialUnghiular = (x y₀/s²) e^{-gs},  scaled by g/2
#
# Neither is an areal density. The first integrates to a length, the second to a
# dimensionless number, so after multiplying both by the same prefactor
# P·N·π(R·10⁻³)² the two results differed by one power of length and *neither*
# was a frequency — yet both were plotted on the same axis and compared against
# the same yr⁻¹ thresholds. The first also dropped the normalisation g of the
# exponential density that the second included, a further factor 1/g = 4.35
# between them.
#
# Two further corrections. The route was integrated from 0 to x₀ only, while the
# comment called x₀ the route *length*; a site beside the middle of a corridor
# receives crashes from both directions, so the integral now runs over the whole
# route and the one-sided result is low by a factor approaching two. And the
# effective target area was the bare geometric disc πR²; DOE-STD-3014 adds skid
# and shadow contributions, so the geometric value is a lower bound and is
# labelled as such.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Probability of loss of control per km of flight."
const P_LOSS = 1e-9
"Flights per year on the route."
const N_FLIGHTS = 7e4
"Decay constant of the radial crash-distance density, km⁻¹."
const G = 0.23
"Half-length of the route segment integrated over, km."
const ROUTE_HALF_LENGTH = 100.0
"Design-basis screening thresholds, yr⁻¹."
const THRESHOLDS = (1e-5, 1e-7)

"Crash-location probability density per unit area at radial distance r [km⁻²]."
areal_density(r) = r <= 0 ? 0.0 : G * exp(-G * r) / (2π * r)

"""
    areal_rate(y₀)

Crash rate per unit area at a site a perpendicular distance `y₀` from the route,
in yr⁻¹ km⁻², integrating along the whole route in both directions.
"""
function areal_rate(y₀)
    integrand(x) = areal_density(hypot(x, y₀))
    I, _ = quadgk(integrand, -ROUTE_HALF_LENGTH, ROUTE_HALF_LENGTH; rtol = 1e-8)
    return N_FLIGHTS * P_LOSS * I
end

"""
    effective_area(R_m)

Effective target area in km² for a circular site of radius `R_m` metres. This is
the geometric footprint only — a full assessment per DOE-STD-3014 adds skid and
shadow terms, so this is a lower bound on the target presented to the aircraft.
"""
effective_area(R_m) = π * (R_m * 1e-3)^2

"Annual impact frequency [yr⁻¹]."
impact_frequency(y₀, R_m) = areal_rate(y₀) * effective_area(R_m)

function main()
    # dimensional audit, printed so the reader can check the cancellation
    println("units:  N [yr⁻¹] · P [km⁻¹] · ∫f dx [km⁻² · km] · A [km²]  =  yr⁻¹")
    @printf("N·P = %.2e km⁻¹ yr⁻¹\n\n", N_FLIGHTS * P_LOSS)

    for y₀ in (5.0, 10.0, 25.0, 50.0)
        Φ = areal_rate(y₀)
        for R in (50.0, 200.0)
            F = impact_frequency(y₀, R)
            verdict = F > THRESHOLDS[1] ? "EXCEEDS 1e-5" :
                      F > THRESHOLDS[2] ? "between 1e-7 and 1e-5" : "below 1e-7"
            @printf("y₀ = %5.1f km, R = %5.0f m:  Φ = %.3e yr⁻¹km⁻²,  F = %.3e yr⁻¹   (%s)\n",
                y₀, R, Φ, F, verdict)
        end
    end

    # the site radius at which each threshold is reached, nearest approach
    y_min = 5.0
    Φ_min = areal_rate(y_min)
    for τ in THRESHOLDS
        R_crit = sqrt(τ / (Φ_min * π)) * 1e3
        @printf("\nat y₀ = %.0f km, the %.0e yr⁻¹ threshold needs a site radius of %.0f m\n",
            y_min, τ, R_crit)
    end

    y₀s = range(5, 50, length = 120)
    Rs = [25.0, 50.0, 100.0, 200.0, 400.0]

    fig = Figure(size = (900, 500))
    # Ticks every two decades from the lowest curve up through both thresholds,
    # so that each threshold falls on a labelled decade and can be read off.
    ax = Axis(fig[2, 1], xlabel = L"Distance from the route $y_0$ [km]",
        ylabel = L"Annual impact frequency $F$ [yr$^{-1}$]", yscale = log10,
        yticks = logticks(-13, -5; step = 2),)

    handles = []
    for (k, R) in enumerate(Rs)
        F = impact_frequency.(y₀s, R)
        push!(handles,
            lines!(ax, y₀s, F,
                color = (PALETTE.blue, PALETTE.sky, PALETTE.green,
                    PALETTE.orange, PALETTE.purple,)[k], linewidth = 1.6,),)
    end
    t1 = hlines!(
        ax, [THRESHOLDS[1]], color = PALETTE.red, linestyle = :dash, linewidth = 1.4,)
    t2 = hlines!(
        ax, [THRESHOLDS[2]], color = PALETTE.red, linestyle = :dot, linewidth = 1.4,)
    # Each threshold labelled on its own line rather than by one note four to
    # seven decades below both of them, which the curves ran through.
    for (thr, lab) in ((THRESHOLDS[1], "10⁻⁵"), (THRESHOLDS[2], "10⁻⁷"))
        text!(ax, maximum(y₀s), thr * 1.5;
            text = rich("Screening threshold ", lab, " yr", superscript("−1")),
            align = (:right, :bottom), color = PALETTE.red, fontsize = 15,)
    end
    ylims!(ax, nothing, THRESHOLDS[1] * 60)

    # the unit is upright: set in maths italic, "m" reads as a variable
    Legend(fig[1, 1], handles,
        [rich(it("R"), " = $(Int(R)) m") for R in Rs],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 20,)
    rowsize!(fig.layout, 2, Relative(0.85))
    println("\nwrote ", savefigure(fig, FIGURES, "aircraft_crash_frequency"))
end

main()
