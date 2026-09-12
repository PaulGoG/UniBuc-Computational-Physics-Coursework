# Half-lives of three neutron-activation products, from the linearised decay law
#
#   ln(N₀/Nᵢ) = λ (tᵢ - t₀),
#
# fitted to counts accumulated in successive equal acquisition intervals.
#
#   ²⁷Al(n,γ)²⁸Al   T½ = 2.245 min
#   ²⁶Mg(n,γ)²⁷Mg   T½ = 9.458 min
#   ¹²⁷I(n,γ)¹²⁸I   T½ = 24.99 min, in the NaI(Tl) crystal itself
#
# Ported from ReactiiNeutronice.jl and FitActivareNaITl.jl, which are the same
# measurement done twice with different conventions — one weighted, one not.
#
# Corrections:
#
#   1. **`ReactiiNeutronice.jl` fitted with no weights at all** and never called
#      `stderror`, so it reported half-lives with no uncertainty, while its
#      sibling weighted the identical kind of data. Both are weighted here.
#   2. The weights that were used, `wt = N`, imply Var[ln(N₀/N)] = 1/N and drop
#      the 1/N₀ term. Poisson propagation gives 1/N₀ + 1/N, and since every point
#      shares N₀ the residuals are correlated — the full covariance is used here.
#   3. `N_Mg = log.(N_Mg[1]./N_Mg)` overwrote the raw counts with their
#      logarithms, destroying the Poisson information and forcing the first point
#      to exactly zero with no assigned error.
#   4. The NaI calibration extrapolated below its own range — the lowest
#      calibration point is channel 77 and the studied peak sits at channel 67 —
#      and reported the resulting energy with no uncertainty, although it is the
#      identification of the nuclide.
#
# One claim from the earlier code review does **not** hold and is not repeated:
# using interval start times rather than midpoints was said to bias λ. For
# equal-length acquisition intervals the (1 - e^{-λΔ}) factor is common to every
# point and cancels in the ratio, so the slope recovers λ exactly. Only the
# intercept shifts.

using Printf, Statistics, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"NaI(Tl) energy calibration: γ lines of known energy against channel."
const CAL_ENERGY = [511.0, 1173.0, 1274.0, 1332.0]     # keV
const CAL_CHANNEL = [77.0, 172.0, 185.0, 195.0]

"(label, acquisition start times in s, counts, literature T½ in min)"
const SERIES = [
    ("²⁸Al", collect(100.0:100.0:500.0), [2985.0, 1728, 1116, 632, 402], 2.245),
    ("²⁷Mg", collect(100.0:100.0:500.0), [1459.0, 1018, 1175, 953, 829], 9.458),
    ("¹²⁸I", collect(520.0:260.0:1300.0), [3691.0, 3275, 2827, 2487], 24.99),
]

"""
    decay_fit(t, N)

Weighted least squares of ln(N₀/Nᵢ) against tᵢ with the full Poisson covariance,
including the correlation induced by the shared reference count N₀. Returns λ,
its standard error, and the fitted values.
"""
function decay_fit(t, N)
    y = log.(N[1] ./ N[2:end])
    x = t[2:end]
    C = [i == j ? 1/N[1] + 1/N[i+1] : 1/N[1] for i in eachindex(x), j in eachindex(x)]
    M = hcat(ones(length(x)), x)
    W = inv(C)
    p = (M' * W * M) \ (M' * W * y)
    σ = sqrt.(diag(inv(M' * W * M)))
    return p[2], σ[2], x, y, p
end

function main()
    # energy calibration
    Mc = hcat(ones(length(CAL_CHANNEL)), CAL_CHANNEL)
    pc = Mc \ CAL_ENERGY
    resid = CAL_ENERGY .- Mc * pc
    s² = sum(abs2, resid) / (length(CAL_ENERGY) - 2)
    σc = sqrt.(diag(s² * inv(Mc' * Mc)))
    peak_channel = 67.0
    E_peak = pc[1] + pc[2] * peak_channel
    σ_peak = sqrt(σc[1]^2 + (peak_channel * σc[2])^2)
    @printf("NaI(Tl) calibration: E = %.2f + %.4f × channel\n", pc[1], pc[2])
    @printf("studied peak at channel %.0f -> %.1f ± %.1f keV\n", peak_channel, E_peak, σ_peak)
    @printf("  (an extrapolation: the lowest calibration point is channel %.0f)\n",
            minimum(CAL_CHANNEL))
    @printf("  ¹²⁸I emits a γ at 442.9 keV\n\n")

    fig = Figure(size = (1000, 450))
    ax1 = Axis(fig[2, 1], xlabel = L"Acquisition start $t$ [s]", ylabel = L"\ln(N_0/N)")
    ax2 = Axis(fig[2, 2], xlabel = "", ylabel = L"$T_{1/2}$ [min]",
        xticks = (1:length(SERIES), [s[1] for s in SERIES]))

    handles = []
    measured = Float64[]; errors = Float64[]; reference = Float64[]
    for (k, (label, t, N, T_lit)) in enumerate(SERIES)
        λ, σλ, x, y, p = decay_fit(t, N)
        T½ = log(2) / λ / 60
        σT = log(2) / λ^2 / 60 * σλ
        push!(measured, T½); push!(errors, σT); push!(reference, T_lit)
        @printf("%-5s  λ = %.3e ± %.1e s⁻¹   T½ = %6.3f ± %.3f min   literature %6.3f   (%.1fσ)\n",
                label, λ, σλ, T½, σT, T_lit, abs(T½ - T_lit) / σT)

        col = (PALETTE.blue, PALETTE.orange, PALETTE.green)[k]
        push!(handles, scatter!(ax1, x, y, color = col, markersize = 10))
        xf = range(0, maximum(x) * 1.05, length = 50)
        lines!(ax1, xf, p[1] .+ λ .* xf, color = col, linewidth = 1.3)
    end

    errorbars!(ax2, 1:length(SERIES), measured, errors,
        color = PALETTE.blue, whiskerwidth = 12)
    scatter!(ax2, 1:length(SERIES), measured, color = PALETTE.blue, markersize = 12)
    scatter!(ax2, 1:length(SERIES), reference, color = PALETTE.red,
        markersize = 14, marker = :hline)
    text!(ax2, 0.5, 0.95; text = "red: literature", space = :relative,
        align = (:center, :top), color = PALETTE.red, fontsize = 15)

    Legend(fig[1, 1:2], handles, [s[1] for s in SERIES],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)
    rowsize!(fig.layout, 2, Relative(0.85))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_activation_halflives"))
end

main()
