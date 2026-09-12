# Two-dimensional Ising model on a square lattice, driven to its ground state by
# Metropolis single-spin-flip dynamics under a geometrically cooled temperature.
#
# Ported from IsingFinal.cpp (2018). The original allocated n columns per row but
# indexed n+2 of them, bordered its ghost cells once before the Monte Carlo loop
# and never refreshed them, double-counted every bond in the Hamiltonian so that
# the acceptance test ran at an effective temperature of T/2, cooled once per
# attempted flip rather than per sweep so the temperature underflowed to zero
# within a few thousand attempts, and seeded from the clock so no run could be
# reproduced. Periodic boundaries are handled here by modular indexing, which
# removes the ghost cells entirely.

using Random, StableRNGs, Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Lattice side length."
const L = 64
"Exchange coupling. Positive is ferromagnetic."
const J = 1.0
"Initial temperature, in units of J/k_B. Boltzmann's constant is taken as 1."
const T_INITIAL = 5.0
"Final temperature of the annealing schedule, same units."
const T_FINAL = 0.05
"Monte Carlo sweeps; one sweep is L² attempted flips."
const SWEEPS = 3000
const SEED = 20180101

@enum Boundary periodic open

"""
    neighbour_sum(s, i, j, boundary)

Sum of the four nearest-neighbour spins of site `(i, j)`. Under `periodic` the
indices wrap; under `open` sites beyond the edge contribute nothing.
"""
function neighbour_sum(s::Matrix{Int8}, i::Int, j::Int, boundary::Boundary)
    n = size(s, 1)
    if boundary == periodic
        up    = s[mod1(i - 1, n), j]
        down  = s[mod1(i + 1, n), j]
        left  = s[i, mod1(j - 1, n)]
        right = s[i, mod1(j + 1, n)]
    else
        up    = i > 1 ? s[i - 1, j] : zero(Int8)
        down  = i < n ? s[i + 1, j] : zero(Int8)
        left  = j > 1 ? s[i, j - 1] : zero(Int8)
        right = j < n ? s[i, j + 1] : zero(Int8)
    end
    return Int(up) + Int(down) + Int(left) + Int(right)
end

"""
    energy(s, J, boundary)

Total energy of configuration `s` for H = -J Σ⟨ij⟩ sᵢsⱼ, each bond counted once.
"""
function energy(s::Matrix{Int8}, J::Real, boundary::Boundary)
    total = 0
    for j in axes(s, 2), i in axes(s, 1)
        total += Int(s[i, j]) * neighbour_sum(s, i, j, boundary)
    end
    # every bond was visited from both of its ends
    return -J * total / 2
end

"""
    flip_energy(s, i, j, J, boundary)

Energy change on flipping the spin at `(i, j)`: ΔE = 2 J sᵢⱼ Σ_nb. The original
used a factor of four here, which halves the effective temperature.
"""
flip_energy(s::Matrix{Int8}, i::Int, j::Int, J::Real, boundary::Boundary) =
    2 * J * Int(s[i, j]) * neighbour_sum(s, i, j, boundary)

"""
    anneal(rng, L, J, T_initial, T_final, sweeps, boundary; snapshot_at)

Metropolis sweeps with the temperature cooled geometrically from `T_initial` to
`T_final` once per sweep. Returns the temperature, energy-per-spin and
magnetisation traces, plus a copy of the configuration taken as the schedule
first drops below each temperature in `snapshot_at`.
"""
function anneal(rng, L::Int, J::Real, T_initial::Real, T_final::Real,
                sweeps::Int, boundary::Boundary;
                snapshot_at::Vector{Float64} = Float64[])
    s = rand(rng, Int8[-1, 1], L, L)
    cooling = (T_final / T_initial)^(1 / (sweeps - 1))

    temperatures = Vector{Float64}(undef, sweeps)
    energies = Vector{Float64}(undef, sweeps)
    magnetisations = Vector{Float64}(undef, sweeps)
    snapshots = Pair{Float64,Matrix{Int8}}[]
    pending = sort(snapshot_at; rev = true)
    T = float(T_initial)

    for sweep in 1:sweeps
        for _ in 1:(L * L)
            i, j = rand(rng, 1:L), rand(rng, 1:L)
            ΔE = flip_energy(s, i, j, J, boundary)
            if ΔE <= 0 || rand(rng) < exp(-ΔE / T)
                s[i, j] = -s[i, j]
            end
        end
        temperatures[sweep] = T
        energies[sweep] = energy(s, J, boundary) / (L * L)
        magnetisations[sweep] = sum(Int, s) / (L * L)
        if !isempty(pending) && T <= first(pending)
            push!(snapshots, popfirst!(pending) => copy(s))
        end
        T *= cooling
    end
    return s, temperatures, energies, magnetisations, snapshots
end

"Exact Onsager critical temperature of the 2-D square-lattice Ising model."
const T_CRITICAL = 2 / log1p(sqrt(2))

function main()
    rng = StableRNG(SEED)
    shots = [4.0, T_CRITICAL, 0.3]
    s, temperatures, energies, magnetisations, snapshots =
        anneal(rng, L, J, T_INITIAL, T_FINAL, SWEEPS, periodic; snapshot_at = shots)

    @printf("final energy per spin  E/N = %.4f J  (ground state %.4f J)\n",
            last(energies), -2.0 * J)
    @printf("final magnetisation    m   = %+.4f\n", last(magnetisations))
    @printf("Onsager T_c = %.4f J/k_B\n", T_CRITICAL)

    fig = Figure(size = (1000, 740))

    ax = Axis(fig[2, 1:3],
        xlabel = L"Temperature $T$ [$J/k_\mathrm{B}$]",
        ylabel = L"Energy per spin $E/N$ [$J$]",
        xreversed = true, xscale = log10,
        xticks = ([0.05, 0.1, 0.5, 1, 5], ["0.05", "0.1", "0.5", "1", "5"]),
        yticks = -2:0.5:0)
    l_energy = lines!(ax, temperatures, energies, color = PALETTE.blue, linewidth = 1.6)
    l_ground = hlines!(ax, [-2.0 * J], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    l_tc = vlines!(ax, [T_CRITICAL], color = PALETTE.red, linestyle = :dot, linewidth = 1.2)
    ylims!(ax, -2.25, -0.15)

    axm = Axis(fig[2, 1:3],
        ylabel = L"Magnetisation $|m|$",
        yaxisposition = :right, xreversed = true, xscale = log10,
        ygridvisible = false, xgridvisible = false,
        yticklabelcolor = PALETTE.orange, ylabelcolor = PALETTE.orange,
        ytickcolor = PALETTE.orange)
    hidexdecorations!(axm)
    l_mag = lines!(axm, temperatures, abs.(magnetisations),
                   color = PALETTE.orange, linewidth = 1.6)
    linkxaxes!(ax, axm)
    ylims!(axm, -0.03, 1.08)

    Legend(fig[1, 1:3],
        [l_energy, l_mag, l_tc, l_ground],
        [L"Energy per spin $E/N$", L"Magnetisation $|m|$",
         L"Onsager $T_\mathrm{c} = 2.269\ J/k_\mathrm{B}$", L"Ground state $E/N = -2J$"],
        orientation = :horizontal, framevisible = false,
        padding = (0, 0, 0, 0), labelsize = 17, colgap = 22)

    for (k, (T, config)) in enumerate(snapshots)
        axk = Axis(fig[3, k], aspect = DataAspect(),
            title = L"T = %$(round(T, digits = 2))\ J/k_\mathrm{B}",
            titlesize = 18)
        heatmap!(axk, config', colormap = [PALETTE.orange, PALETTE.blue],
                 colorrange = (-1, 1))
        hidedecorations!(axk)
        hidespines!(axk)
    end

    rowsize!(fig.layout, 2, Relative(0.58))
    rowgap!(fig.layout, 14)
    path = savefigure(fig, FIGURES, "ising_annealing")
    println("wrote ", path)
end

main()
