# Monte Carlo chain of successive Compton scatterings of a single photon, from
# an initial energy down to an absorption threshold.
#
#   E' = E / (1 + (E/mₑc²)(1 - cos θ))            Compton relation
#   tan φ = E' sin θ / (E - E' cos θ)              recoil-electron angle
#
# Ported from EfectulCompton.cpp (2018). The three physics formulae in the
# original were individually correct. The sampling was not: the scattering angle
# was drawn uniformly on [0°, 180°], which is neither the Klein-Nishina
# distribution nor even isotropic scattering — isotropy is uniform in cos θ, not
# in θ. Every energy distribution the original produced was therefore
# unphysical. The angle is sampled here from the Klein-Nishina differential
# cross-section by rejection.
#
# Two further defects, both removed: the original computed the scattered energy
# a second time from the recoil angle by the sine rule and averaged the two,
# but the second value is an algebraic identity of the first and carries no
# independent information; and its "standard deviation" was a relative deviation
# of one sample from a running mean, cast to Int and taken mod 100, with no sum
# of squares anywhere.

using Printf, StableRNGs, StatsBase
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Electron rest energy."
const Mₑc² = 510.999
"Initial photon energy."
const E₀ = 1000.0
"Photons are considered absorbed below this energy."
const E_THRESHOLD = 20.0
const SEED = 20180507

"""
    compton_energy(E, cosθ)

Scattered photon energy for incident energy `E` and scattering angle `θ`.
"""
compton_energy(E, cosθ) = E / (1 + (E / Mₑc²) * (1 - cosθ))

"""
    klein_nishina(E, cosθ)

Differential cross-section dσ/d(cos θ) in units of πrₑ², for incident energy
`E`. The Thomson limit (1 + cos²θ) is recovered as E/mₑc² → 0.
"""
function klein_nishina(E, cosθ)
    P = 1 / (1 + (E / Mₑc²) * (1 - cosθ))
    return P^2 * (P + 1 / P - (1 - cosθ^2))
end

"""
    sample_cosθ(rng, E)

Draw cos θ from the Klein-Nishina distribution by rejection against its maximum
over the interval, which is located by a coarse scan.
"""
function sample_cosθ(rng, E)
    grid = range(-1, 1, length = 512)
    fmax = maximum(klein_nishina(E, c) for c in grid) * 1.05
    while true
        c = 2 * rand(rng) - 1
        rand(rng) * fmax <= klein_nishina(E, c) && return c
    end
end

"""
    scatter_chain(rng, E₀, threshold)

Follow one photon until its energy falls below `threshold`. Returns the energy
after each scattering, starting from `E₀`.
"""
function scatter_chain(rng, E₀, threshold)
    energies = [E₀]
    E = E₀
    while E > threshold
        E = compton_energy(E, sample_cosθ(rng, E))
        push!(energies, E)
        length(energies) > 10_000 && break
    end
    return energies
end

function main()
    rng = StableRNG(SEED)
    n_photons = 20_000

    chains = [scatter_chain(rng, E₀, E_THRESHOLD) for _ in 1:n_photons]
    n_scatters = length.(chains) .- 1

    @printf("initial energy %.0f keV, threshold %.0f keV\n", E₀, E_THRESHOLD)
    @printf("scatterings to threshold: mean %.2f, median %d, range %d-%d\n",
            mean(n_scatters), median(n_scatters), minimum(n_scatters), maximum(n_scatters))

    # first-scattering angular distribution, Klein-Nishina against the uniform-in-θ
    # sampling the 2018 code used
    cosθ_kn = [sample_cosθ(rng, E₀) for _ in 1:200_000]
    θ_kn = acosd.(cosθ_kn)
    θ_uniform = 180 .* rand(rng, 200_000)
    @printf("mean scattering angle: Klein-Nishina %.1f deg, uniform-in-theta %.1f deg\n",
            mean(θ_kn), mean(θ_uniform))

    fig = Figure(size = (1040, 460))

    ax1 = Axis(fig[2, 1], xlabel = "Scattering number",
        ylabel = L"Photon energy $E$ [keV]", yscale = log10,
        yticks = ([20, 50, 100, 200, 500, 1000],
                  ["20", "50", "100", "200", "500", "1000"]))
    for c in chains[1:40]
        lines!(ax1, 0:(length(c) - 1), c, color = (PALETTE.blue, 0.35), linewidth = 1.0)
    end
    h_thr = hlines!(ax1, [E_THRESHOLD], color = PALETTE.black,
        linestyle = :dash, linewidth = 1.0)
    text!(ax1, 0.97, 0.06; text = "40 of $(n_photons) histories",
        space = :relative, align = (:right, :bottom), fontsize = 15, color = PALETTE.blue)

    ax2 = Axis(fig[2, 2], xlabel = "Scatterings to threshold", ylabel = "Photons")
    hist!(ax2, n_scatters, bins = range(-0.5, maximum(n_scatters) + 0.5,
          length = maximum(n_scatters) + 2), color = (PALETTE.blue, 0.7),
          strokewidth = 0.5, strokecolor = PALETTE.blue)
    vlines!(ax2, [mean(n_scatters)], color = PALETTE.red, linestyle = :dash, linewidth = 1.2)
    text!(ax2, 0.96, 0.9; text = @sprintf("Mean %.2f", mean(n_scatters)),
        space = :relative, align = (:right, :top), color = PALETTE.red, fontsize = 15)

    ax3 = Axis(fig[2, 3], xlabel = L"Scattering angle $\theta$ [deg]",
        ylabel = "Normalised density",
        xticks = ([0, 45, 90, 135, 180], ["0", "45", "90", "135", "180"]))
    d_kn = StatsBase.normalize(fit(Histogram, θ_kn, range(0, 180, length = 60)), mode = :pdf)
    d_un = StatsBase.normalize(fit(Histogram, θ_uniform, range(0, 180, length = 60)), mode = :pdf)
    l_kn = stairs!(ax3, midpoints(d_kn.edges[1]), d_kn.weights,
        color = PALETTE.green, linewidth = 1.8, step = :center)
    l_un = stairs!(ax3, midpoints(d_un.edges[1]), d_un.weights,
        color = PALETTE.orange, linewidth = 1.8, step = :center)

    Legend(fig[1, 1:3], [h_thr, l_kn, l_un],
        ["Absorption threshold", L"Klein-Nishina at $E_0$", "Uniform in θ (2018)"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)

    rowsize!(fig.layout, 2, Relative(0.86))
    path = savefigure(fig, FIGURES, "compton_scattering_chain")
    println("wrote ", path)
end

main()
