# Pre-neutron fragment distributions of ²³⁵U(n_th,f) from the Straede
# Y(A_H, Z_H, TKE) matrix: mass, charge, neutron-number and TKE yields, ⟨TKE⟩(A),
# the total excitation energy TXE(A) = Q(A) + Sₙ(²³⁶U) − TKE(A), the
# single-fragment kinetic energies, and the yield-weighted totals with their
# uncertainties.
#
# What the matrix measures is Y(A_H, TKE): the five charge rows of each
# (A_H, TKE) cell split one measured yield by fractions that are identical
# across TKE for a given mass, to 3 × 10⁻⁶ (asserted below), so Y(Z), Y(N) and
# the even–odd staggering of Y(Z) describe the charge model imposed on the
# matrix, not a measurement. The uncertainties are propagated at the level of
# the measured cells: the rows of one cell share one relative error and add
# linearly, distinct cells add in quadrature, and every total is a weighted
# mean over cells of a quantity fixed per cell, so its variance is
# Σ_c [(q_c − ⟨q⟩) σY_c / ΣY]². Averaging first over TKE and propagating at
# the mass level, as an earlier version of this script did, drops the
# within-mass spread of TKE and halves σ⟨TKE⟩.
#
# Ported from Fisiune_2.jl on the `legacy` branch (Julia-Workflow-FFUB/
# Fisiune_M_2/). Its Y(N) summed over all rows with a given N inside the
# double loop over (A_H, Z_H), once per pair mapping to that N — the total
# came to 805 % instead of 100 %, reproduced below; three of its uncertainty
# accumulations added linearly where the comments said quadrature; its
# `Sortare_distributie` assumed an abscissa without gaps or duplicates; and its
# `Energie_separare` returned a bare `NaN` on failure that the caller indexed.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_core.jl"))

"Tolerance on the charge fractions of a mass being the same in every TKE cell."
const CHARGE_MODEL_TOLERANCE = 1e-5
"Yield below which a charge or TKE bin is not drawn [%]."
const DRAW_FLOOR = 0.05

"""
    charge_fraction_spread(rows)

Largest difference, over all masses, between the charge fractions Y(Z | A, TKE)
of two TKE cells of the same mass. Zero if the charge split is a fixed model.
"""
function charge_fraction_spread(rows)
    by_cell = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
    for r in rows
        d = get!(by_cell, (r.A_H, r.TKE), Dict{Int,Float64}())
        d[r.Z_H] = get(d, r.Z_H, 0.0) + r.Y
    end
    reference = Dict{Int,Dict{Int,Float64}}()
    spread = 0.0
    for ((A, _), d) in by_cell
        total = sum(values(d))
        total > 0 || continue
        f = Dict(Z => y / total for (Z, y) in d)
        if haskey(reference, A)
            for (Z, v) in f
                spread = max(spread, abs(v - get(reference[A], Z, NaN)))
            end
        else
            reference[A] = f
        end
    end
    return spread
end

"Even–odd staggering δ = (ΣY_even − ΣY_odd)/(ΣY_even + ΣY_odd) of a charge yield."
function even_odd_staggering(Z, Y)
    even, odd = sum(Y[iseven.(Z)]), sum(Y[isodd.(Z)])
    return (even - odd) / (even + odd)
end

"The Y(N) of Fisiune_2.jl, summed once per (A_H, Z_H) pair mapping to each N."
function overcounted_neutron_yield(rows)
    by_N = Dict{Int,Float64}()
    for r in rows
        by_N[r.A_H - r.Z_H] = get(by_N, r.A_H - r.Z_H, 0.0) + r.Y
    end
    pairs = unique((r.A_H, r.Z_H) for r in rows)
    over = Dict{Int,Float64}()
    for (A, Z) in pairs
        over[A - Z] = get(over, A - Z, 0.0) + by_N[A - Z]
    end
    ks = sort!(collect(keys(over)))
    return ks, [over[k] for k in ks]
end

function main()
    rows = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    cells = yield_cells(rows)
    @printf("Straede matrix: %d rows in %d (A_H, TKE) cells; A_H %d–%d, Z_H %d–%d, TKE %d–%d MeV\n",
        length(rows), length(cells), extrema(r.A_H for r in rows)...,
        extrema(r.Z_H for r in rows)...,
        extrema(r.TKE for r in rows)...)
    @printf("total yield %.3f %% (heavy fragment only)\n", sum(c.Y for c in cells))

    spread = charge_fraction_spread(rows)
    @printf("charge fractions Y(Z | A, TKE) differ across TKE by at most %.1e: an imposed charge split\n",
        spread)
    spread < CHARGE_MODEL_TOLERANCE ||
        error("the charge fractions vary across TKE by $spread; they are not a fixed model")

    A, Y_A, σY_A = marginal(rows, r -> r.A_H)
    Z, Y_Z, _ = marginal(rows, r -> r.Z_H)
    N, Y_N, _ = marginal(rows, r -> r.A_H - r.Z_H)
    T, Y_T, _ = marginal(rows, r -> r.TKE)
    δ = even_odd_staggering(Z, Y_Z)
    @printf("\nY(A) peaks at A_H = %d with %.3f ± %.3f %%\n", A[argmax(Y_A)], maximum(Y_A),
        σY_A[argmax(Y_A)])
    @printf("Y(Z) peaks at Z_H = %d, Y(N) at N = %d; even–odd staggering of the charge model δ = %.4f\n",
        Z[argmax(Y_Z)], N[argmax(Y_N)], δ)
    N_over, Y_over = overcounted_neutron_yield(rows)
    @printf("Y(N) of Fisiune_2.jl: total %.1f %% instead of %.1f %%, peak at N = %d instead of %d\n",
        sum(Y_over), sum(Y_N), N_over[argmax(Y_over)], N[argmax(Y_N)])

    # per-mass quantities: ⟨TKE⟩(A) from the cells, Q(A) charge-averaged as in
    # fission_q_value.jl, TXE and the momentum-conservation split of TKE
    bins = tke_by_mass(cells)
    S_n = separation_energy(masses, Z₀, A₀)[1]
    @printf("\nSₙ(²³⁶U) = %.3f MeV\n", S_n)
    A_keep = Int[]
    TKE_A, TXE_A, Q_A = Float64[], Float64[], Float64[]
    for a in sort!(collect(keys(bins)))
        q = charge_average(Zh -> q_value(masses, Zh, a), a, nearest_three_charges(a))
        q === nothing && continue
        push!(A_keep, a)
        push!(TKE_A, bins[a].TKE)
        push!(Q_A, q[1])
        push!(TXE_A, q[1] + S_n - bins[a].TKE)
    end
    KE_L = TKE_A .* A_keep ./ A₀
    KE_H = TKE_A .* (A₀ .- A_keep) ./ A₀
    @printf("⟨TKE⟩(A) %.1f–%.1f MeV, TXE(A) %.1f–%.1f MeV, KE_L %.1f–%.1f MeV, KE_H %.1f–%.1f MeV\n",
        extrema(TKE_A)..., extrema(TXE_A)..., extrema(KE_L)..., extrema(KE_H)...)

    # the five totals, each a weighted mean over cells of a quantity fixed per
    # cell: A_H, A_L and TKE are the cell's labels, Q and TXE take the mass's Q(A)
    Q_of = Dict(zip(A_keep, Q_A))
    used = [c for c in cells if haskey(Q_of, c.A_H)]
    Y_c = [c.Y for c in used]
    σ_c = [c.σY for c in used]
    totals = (
        ("⟨A_H⟩", [float(c.A_H) for c in used], ""),
        ("⟨A_L⟩", [float(A₀ - c.A_H) for c in used], ""),
        ("⟨TKE⟩", [float(c.TKE) for c in used], " MeV"),
        ("⟨Q⟩", [Q_of[c.A_H] for c in used], " MeV"),
        ("⟨TXE⟩", [Q_of[c.A_H] + S_n - c.TKE for c in used], " MeV"),
    )
    println("\nyield-weighted totals, uncertainty from the cell yields:")
    means = Dict{String,Float64}()
    for (name, q, unit) in totals
        m, σ = yield_weighted_mean(q, Y_c, σ_c)
        means[name] = m
        @printf("  %-6s = %9.4f ± %.4f%s\n", name, m, σ, unit)
    end
    assert_close("⟨A_H⟩ + ⟨A_L⟩", means["⟨A_H⟩"] + means["⟨A_L⟩"], A₀; atol = 1e-9)
    assert_close("⟨TXE⟩ − ⟨Q⟩ − Sₙ + ⟨TKE⟩",
        means["⟨TXE⟩"] - means["⟨Q⟩"] - S_n + means["⟨TKE⟩"], 0.0;
        atol = 1e-9,)
    ke_L = yield_weighted_mean(KE_L, [bins[a].Y for a in A_keep], zeros(length(A_keep)))[1]
    ke_H = yield_weighted_mean(KE_H, [bins[a].Y for a in A_keep], zeros(length(A_keep)))[1]
    @printf("  ⟨KE_L⟩ = %.3f MeV, ⟨KE_H⟩ = %.3f MeV, sum %.3f\n", ke_L, ke_H, ke_L + ke_H)
    # at the mass level the within-mass spread of TKE is lost
    tke_mass = yield_weighted_mean(TKE_A, [bins[a].Y for a in A_keep],
        [sqrt(sum(abs2, c.σY for c in used if c.A_H == a)) for a in A_keep],)
    @printf("  σ⟨TKE⟩ propagated at the mass level instead: %.4f MeV\n", tke_mass[2])

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1400, 1080))
    yield_label(v) = rich(it("Y"), "(", v, ") [%]")

    ax1 = Axis(fig[2, 1], xlabel = rich("Heavy-fragment mass ", it("A"), subscript("H")),
        ylabel = yield_label(rich(it("A"), subscript("H"))),)
    errorbars_unstroked!(ax1, A, Y_A, σY_A, color = PALETTE.blue, linewidth = GUIDE_WIDTH,
        whiskerwidth = 6,)
    scatterlines!(ax1, A, Y_A, color = PALETTE.blue, markersize = MARKERSIZE.dense)
    text!(ax1, 0.97, 0.95;
        text = rich("Peak ", it("A"), subscript("H"), @sprintf(" = %d, %.2f %%",
            A[argmax(Y_A)], maximum(Y_A))),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.blue,)

    ax2 = Axis(fig[2, 2], xlabel = rich("Heavy-fragment charge ", it("Z"), subscript("H")),
        ylabel = yield_label(rich(it("Z"), subscript("H"))),)
    barplot!(
        ax2, Z, Y_Z, color = PALETTE.blue, strokewidth = 1.5, strokecolor = PALETTE.black,)
    drawn = Y_Z .> DRAW_FLOOR
    xlims!(ax2, minimum(Z[drawn]) - 1, maximum(Z[drawn]) + 1)
    text!(ax2, 0.03, 0.95;
        text = rich("Imposed charge split\n", it("δ"), @sprintf(" = %.3f", δ)),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.blue,)

    ax3 = Axis(fig[3, 1], xlabel = "TKE [MeV]", ylabel = yield_label("TKE"))
    scatterlines!(ax3, T, Y_T, color = PALETTE.blue, markersize = MARKERSIZE.dense)
    drawn_T = Y_T .> 1e-3
    xlims!(ax3, minimum(T[drawn_T]) - 3, maximum(T[drawn_T]) + 3)
    vlines!(ax3, [means["⟨TKE⟩"]], color = PALETTE.black,
        linestyle = :dash, linewidth = GUIDE_WIDTH,)
    text!(ax3, 0.97, 0.5; text = rich("⟨TKE⟩ = ", @sprintf("%.2f MeV", means["⟨TKE⟩"])),
        space = :relative, align = (:right, :center), fontsize = ANNOTATION_SIZE,)

    ax4 = Axis(fig[3, 2], xlabel = rich("Heavy-fragment mass ", it("A"), subscript("H")),
        ylabel = "Energy [MeV]",)
    l_tke = lines!(ax4, A_keep, TKE_A, color = PALETTE.blue)
    l_txe = lines!(ax4, A_keep, TXE_A, color = PALETTE.red)
    l_kel = lines!(ax4, A_keep, KE_L, color = PALETTE.green, linestyle = :dash)
    l_keh = lines!(ax4, A_keep, KE_H, color = PALETTE.purple, linestyle = :dashdot)
    ylims!(ax4, 0, maximum(TKE_A) * 1.12)

    ke = rich("KE")
    Legend(fig[1, 1:2], [l_tke, l_txe, l_kel, l_keh],
        [rich("⟨TKE⟩(", it("A"), subscript("H"), ")"),
            rich("TXE(", it("A"), subscript("H"), ")"),
            rich(ke, subscript("L"), "(", it("A"), subscript("H"), ")"),
            rich(ke, subscript("H"), "(", it("A"), subscript("H"), ")"),],
        "Energies"; titleposition = :left,)
    println("\nwrote ", savefigure(fig, FIGURES, "fragment_yields"))
end

main()
