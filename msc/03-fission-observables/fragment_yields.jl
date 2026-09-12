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
        Zp = round(Int, Z₀ * a / A₀ - 0.5)
        δH = Δ(masses, Zp, a); δL = Δ(masses, Z₀ - Zp, A₀ - a)
        (δH === nothing || δL === nothing) && continue
        q = (Δ₀ - δH - δL) / 1000
        push!(A_keep, a); push!(TKE_A, tke); push!(TXE_A, q + S_n - tke)
        push!(KE_L, tke * a / A₀); push!(KE_H, tke * (A₀ - a) / A₀)
    end
    @printf("KE_L ranges %.1f–%.1f MeV, KE_H ranges %.1f–%.1f MeV\n",
            minimum(KE_L), maximum(KE_L), minimum(KE_H), maximum(KE_H))
    @printf("  check: KE_L + KE_H = TKE to %.2e MeV\n",
            maximum(abs.(KE_L .+ KE_H .- TKE_A)))
    @printf("⟨TKE⟩(A) ranges %.1f–%.1f MeV, TXE(A) ranges %.1f–%.1f MeV\n",
            minimum(TKE_A), maximum(TKE_A), minimum(TXE_A), maximum(TXE_A))

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
