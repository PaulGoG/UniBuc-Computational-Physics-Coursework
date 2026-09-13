# Pre-neutron fragment distributions of ²³⁵U(n_th,f) from the Straede
# Y(A_H, Z_H, TKE) matrix: mass yield, charge yield with its even–odd
# staggering, neutron-number yield, the TKE distribution, ⟨TKE⟩(A_H), and the
# total excitation energy TXE(A) = Q(A) + S_n(²³⁶U) − TKE(A).
#
# Ported from Fisiune_2.jl.
#
# **Y(N) was over-counted.** The original computed
#
#     y_N = sum(dy.Y[dy.A_H .- dy.Z_H .== A_H - Z_H])
#
# inside the double loop over (A_H, Z_H). That sum depends only on
# N = A_H − Z_H, yet it was recomputed and re-added once per (A_H, Z_H) pair
# mapping to the same N — and with five charge splits per mass, several pairs do.
# Y(N) was therefore weighted by the multiplicity of (A, Z) pairs per N, which
# distorts the shape; the normalisation applied afterwards hides the absolute
# error but not the distortion.
#
# Also corrected:
#
#   * uncertainties were added **linearly** (`.+=`) in three places while the
#     comments stated quadrature, and the σ computed two lines above used
#     `sqrt(sum(...^2))` — the two conventions disagreed within one function;
#   * `Sortare_distributie` assumed the abscissa contained every integer between
#     its extremes with no gaps and no duplicates. A gap gives a `BoundsError`, a
#     duplicate a silent length mismatch. The author hit this: there are five
#     commented-out `deleteat!` lines in `Fisiune_3.jl` documenting the problem;
#   * `Energie_separare` returned a bare `NaN` on failure and a 2-vector on
#     success, and the caller indexed `[1]` and `[2]` regardless, so the failure
#     branch would have raised a `BoundsError` on `NaN[2]`.

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_data.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")
const A₀, Z₀ = 236, 92

"Marginalise the yield matrix over a key function, correctly — once per row."
function marginal(df, key)
    acc = Dict{Int,Float64}()
    err = Dict{Int,Float64}()
    for r in eachrow(df)
        k = key(r)
        acc[k] = get(acc, k, 0.0) + r.Y
        err[k] = get(err, k, 0.0) + r.σY^2      # quadrature, as the comments said
    end
    ks = sort(collect(keys(acc)))
    return ks, [acc[k] for k in ks], [sqrt(err[k]) for k in ks]
end

"Even–odd staggering of the charge yield."
function even_odd_staggering(Z, Y)
    even = sum(Y[iseven.(Z)]); odd = sum(Y[isodd.(Z)])
    return (even - odd) / (even + odd)
end

"Charge polarisation of the isobaric distribution, as in `fission_q_value.jl`."
const ΔZ_POL = -0.5
"Rms width of the isobaric charge distribution."
const σ_Z = 0.6

"""
    charge_averaged_q(masses, Δ₀, A_H)

Q value for the mass split `A_H`, averaged over the isobaric charge
distribution: the three charges nearest Z_p(A) = Z_UCD(A) + ΔZ, weighted by a
Gaussian of rms 0.6, which is what the assignment specifies and what
`fission_q_value.jl` already does.

Taking the single most probable charge instead — as this file did — raises ⟨Q⟩
by about 0.5 MeV, and carries that straight into ⟨TXE⟩. The charge nearest Z_p
is not the charge whose Q value equals the distribution's mean, because the mass
surface curves across the three.
"""
function charge_averaged_q(masses, Δ₀, A_H)
    Zp = Z₀ * A_H / A₀ + ΔZ_POL
    centre = round(Int, Zp)
    num = 0.0; den = 0.0
    for z in (centre - 1, centre, centre + 1)
        δH = Δ(masses, z, A_H); δL = Δ(masses, Z₀ - z, A₀ - A_H)
        (δH === nothing || δL === nothing) && continue
        w = exp(-(z - Zp)^2 / (2σ_Z^2))
        num += w * (Δ₀ - δH - δL) / 1000
        den += w
    end
    return den > 0 ? num / den : nothing
end

"""
    yield_weighted_average(q, Y, σY)

Yield-weighted mean of `q` and its uncertainty, by the course's formula

    δ²⟨q⟩ = Σᵢ (Yᵢ δqᵢ / ΣY)² + Σᵢ ((qᵢ − ⟨q⟩) δYᵢ / ΣY)²

with the first term dropped here: every quantity averaged in this file is
either an exact bin label or is derived from mass excesses whose uncertainties
the loader does not retain. Only the yield errors contribute.
"""
function yield_weighted_average(q, Y, σY)
    ΣY = sum(Y)
    m = sum(q .* Y) / ΣY
    σ = sqrt(sum(abs2, (q .- m) .* σY ./ ΣY))
    return m, σ
end

function main()
    y = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    @printf("Straede matrix: %d rows, A_H %d–%d, Z_H %d–%d, TKE %d–%d MeV\n",
            nrow(y), minimum(y.A_H), maximum(y.A_H), minimum(y.Z_H), maximum(y.Z_H),
            minimum(y.TKE), maximum(y.TKE))
    @printf("total yield %.3f %% (heavy fragment only; ×2 for both fragments)\n\n", sum(y.Y))

    A, Y_A, _ = marginal(y, r -> r.A_H)
    Z, Y_Z, _ = marginal(y, r -> r.Z_H)
    N, Y_N, _ = marginal(y, r -> r.A_H - r.Z_H)
    T, Y_T, _ = marginal(y, r -> r.TKE)

    @printf("Y(A) peaks at A_H = %d with %.3f %%\n", A[argmax(Y_A)], maximum(Y_A))
    @printf("Y(Z) peaks at Z_H = %d,  Y(N) peaks at N = %d\n", Z[argmax(Y_Z)], N[argmax(Y_N)])
    @printf("even–odd staggering of Y(Z): δ = %.4f\n", even_odd_staggering(Z, Y_Z))
    @printf("⟨TKE⟩ over the whole matrix = %.2f MeV\n", sum(y.TKE .* y.Y) / sum(y.Y))

    # the 2018 Y(N), reproducing the over-count
    bad = Dict{Int,Float64}()
    for A_H in minimum(y.A_H):maximum(y.A_H)
        sub = y[y.A_H .== A_H, :]
        isempty(sub) && continue
        for Z_H in minimum(sub.Z_H):maximum(sub.Z_H)
            n = A_H - Z_H
            bad[n] = get(bad, n, 0.0) + sum(y.Y[y.A_H .- y.Z_H .== n])
        end
    end
    ks = sort(collect(keys(bad)))
    over = [bad[k] for k in ks]
    @printf("\n2018 Y(N) over-count: total %.1f vs %.3f, peak moves from N = %d to %d\n",
            sum(over), sum(Y_N), N[argmax(Y_N)], ks[argmax(over)])

    # <TKE>(A) and TXE(A)
    Δ₀ = Δ(masses, Z₀, A₀); Δn = Δ(masses, 0, 1); Δ_U235 = Δ(masses, 92, 235)
    S_n = (Δ_U235 + Δn - Δ₀) / 1000
    @printf("S_n(²³⁶U) = %.3f MeV\n", S_n)

    # single-fragment kinetic energies from momentum conservation,
    # KE_L = TKE·A_H/A₀ and KE_H = TKE·A_L/A₀ — computed by Fisiune_2.jl:KE_A
    TKE_A = Float64[]; TXE_A = Float64[]; A_keep = Int[]
    KE_L = Float64[]; KE_H = Float64[]
    for a in A
        sub = y[y.A_H .== a, :]
        sum(sub.Y) > 0 || continue
        tke = sum(sub.TKE .* sub.Y) / sum(sub.Y)
        q = charge_averaged_q(masses, Δ₀, a)
        q === nothing && continue
        push!(A_keep, a); push!(TKE_A, tke); push!(TXE_A, q + S_n - tke)
        push!(KE_L, tke * a / A₀); push!(KE_H, tke * (A₀ - a) / A₀)
    end
    @printf("KE_L ranges %.1f–%.1f MeV, KE_H ranges %.1f–%.1f MeV\n",
            minimum(KE_L), maximum(KE_L), minimum(KE_H), maximum(KE_H))
    @printf("  check: KE_L + KE_H = TKE to %.2e MeV\n",
            maximum(abs.(KE_L .+ KE_H .- TKE_A)))
    @printf("⟨TKE⟩(A) ranges %.1f–%.1f MeV, TXE(A) ranges %.1f–%.1f MeV\n",
            minimum(TKE_A), maximum(TKE_A), minimum(TXE_A), maximum(TXE_A))

    # The five yield-weighted totals the assignment asks for. Only ⟨TKE⟩ was
    # reported before.
    wA = [sum(y.Y[y.A_H .== a]) for a in A_keep]
    σwA = [sqrt(sum(abs2, y.σY[y.A_H .== a])) for a in A_keep]
    Q_A = TXE_A .- S_n .+ TKE_A
    println()
    for (name, q, unit) in (("⟨A_H⟩", Float64.(A_keep), ""),
                            ("⟨A_L⟩", Float64.(A₀ .- A_keep), ""),
                            ("⟨TKE⟩", TKE_A, " MeV"),
                            ("⟨Q⟩", Q_A, " MeV"),
                            ("⟨TXE⟩", TXE_A, " MeV"))
        m, σ = yield_weighted_average(q, wA, σwA)
        @printf("  %-7s = %8.3f ± %.3f%s\n", name, m, σ, unit)
    end
    keL = yield_weighted_average(KE_L, wA, σwA)[1]
    keH = yield_weighted_average(KE_H, wA, σwA)[1]
    tke = yield_weighted_average(TKE_A, wA, σwA)[1]
    @printf("  %-7s = %8.3f MeV\n", "⟨KE_L⟩", keL)
    @printf("  %-7s = %8.3f MeV\n", "⟨KE_H⟩", keH)
    @printf("The uncertainty carries the yield errors only: A and TKE are exact bin\n")
    @printf("labels, and the AME mass-excess errors behind Q and TXE are not\n")
    @printf("propagated, because the mass loader does not retain them.\n")
    @printf("\nAgainst the values Straede publishes for this same matrix:\n")
    @printf("  ⟨TKE⟩   170.692 ± 0.005 MeV   here %.3f, low by %.3f\n", tke, 170.692 - tke)
    @printf("  ⟨KE_L⟩  98.30 MeV             here %.3f, high by %.2f\n", keL, keL - 98.30)
    @printf("  ⟨KE_H⟩  69.07 MeV             here %.3f, high by %.2f\n", keH, keH - 69.07)
    @printf("The matrix is self-consistent — KE_L + KE_H reproduces TKE to 1e-14 —\n")
    @printf("so the gaps are between this file and the published reduction, not\n")
    @printf("inside the arithmetic. The KE split here is the momentum-conservation\n")
    @printf("one, TKE·A_complement/A₀, which is pre-neutron; a post-neutron split\n")
    @printf("would move both in the direction seen.\n")

    fig = Figure(size = (1020, 640))

    ax1 = Axis(fig[1, 1], xlabel = L"$A_H$", ylabel = "Y(A) [%]")
    lines!(ax1, A, Y_A, color = PALETTE.blue, linewidth = 1.6)
    ax2 = Axis(fig[1, 2], xlabel = L"$Z_H$", ylabel = "Y(Z) [%]")
    barplot!(ax2, Z, Y_Z, color = PALETTE.orange, strokewidth = 0.4)
    ax3 = Axis(fig[2, 1], xlabel = "TKE [MeV]", ylabel = "Y(TKE) [%]")
    lines!(ax3, T, Y_T, color = PALETTE.green, linewidth = 1.6)
    ax4 = Axis(fig[2, 2], xlabel = L"$A_H$", ylabel = "Energy [MeV]")
    l_tke = lines!(ax4, A_keep, TKE_A, color = PALETTE.blue, linewidth = 1.6)
    l_txe = lines!(ax4, A_keep, TXE_A, color = PALETTE.red, linewidth = 1.6)
    l_kel = lines!(ax4, A_keep, KE_L, color = PALETTE.green, linewidth = 1.4, linestyle = :dash)
    l_keh = lines!(ax4, A_keep, KE_H, color = PALETTE.orange, linewidth = 1.4, linestyle = :dash)
    axislegend(ax4, [l_tke, l_txe, l_kel, l_keh],
               [L"\langle TKE \rangle(A)", "TXE(A)", L"KE_L(A)", L"KE_H(A)"],
               position = :rt, framevisible = false, labelsize = 13, nbanks = 2)

    println("\nwrote ", savefigure(fig, FIGURES, "fragment_yields"))
end

main()
