# Q value of ²³⁵U(n_th,f) from the AME2020 mass excesses,
#
#   Q(A_H, Z_H) = Δ(²³⁶U) − Δ(A_H, Z_H) − Δ(A₀ − A_H, Z₀ − Z_H),
#
# for heavy-fragment masses 118 to 160, averaged at each mass over the three
# charges nearest the most probable one, Z_p(A_H) = Z_UCD(A_H) − 0.5, with the
# weights of a Gaussian isobaric charge distribution of rms width 0.6:
#
#   Q(A_H) = Σ_Z p(Z) Q(A_H, Z) / Σ_Z p(Z),   σ_Q(A_H) = √Σ_Z [p(Z) σ_Q(A_H, Z)]² / Σ_Z p(Z).
#
# The polarisation of 0.5 and the width of 0.6 are rounded averages over the
# actinides of Wahl's Z_p systematics (At. Data Nucl. Data Tables 39, 1 (1988),
# doi:10.1016/0092-640X(88)90016-2).
#
# Ported from Fisiune_1.jl, whose physics is unchanged, the mass-excess
# uncertainties included. The original parsed the mass table inside the
# calculation, filtered the whole table four times per charge, and looped over
# a range of `Float64` charges.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_core.jl"))

function main()
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    bins = tke_by_mass(yield_cells(load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))))
    @printf("Δ(²³⁶U) = %.1f ± %.1f keV, %d nuclides tabulated\n", masses[(Z₀, A₀)]...,
        length(masses))

    A = Int[]
    Q = Float64[]
    σQ = Float64[]
    A_split = Int[]                       # the individual charge splits
    Q_split = Float64[]
    for A_H in A_H_RANGE
        charges = nearest_three_charges(A_H)
        mean_q = charge_average(Z -> q_value(masses, Z, A_H), A_H, charges)
        mean_q === nothing && continue
        push!(A, A_H)
        push!(Q, mean_q[1])
        push!(σQ, mean_q[2])
        for Z in charges
            q = q_value(masses, Z, A_H)
            q === nothing && continue
            push!(A_split, A_H)
            push!(Q_split, q[1])
        end
    end

    # The mean that means something physically weights each mass by its yield. The
    # matrix holds no yield above A_H = 158, so 159 and 160 carry zero weight.
    Y = [haskey(bins, a) ? bins[a].Y : 0.0 for a in A]
    Q_mean = sum(Q .* Y) / sum(Y)
    σ_mass = sqrt(sum(abs2, σQ .* Y)) / sum(Y)

    @printf("Q(A_H) for %d masses, A_H = %d–%d: %.2f to %.2f MeV, σ_Q from %.3f to %.3f MeV\n",
        length(A), first(A), last(A), minimum(Q), maximum(Q), minimum(σQ), maximum(σQ))
    @printf("maximum at A_H = %d: Q = %.2f ± %.2f MeV\n", A[argmax(Q)], maximum(Q),
        σQ[argmax(Q)])
    @printf("symmetric split A_H = %d: Q = %.2f ± %.2f MeV\n", A[1], Q[1], σQ[1])
    @printf("yield-weighted mean ⟨Q⟩ = %.3f ± %.3f MeV (mass-excess errors), %d masses with a yield\n",
        Q_mean, σ_mass, count(>(0), Y))
    @printf("arithmetic mean over the %d mass points = %.2f MeV\n", length(A), mean(Q))

    fig = Figure(size = (900, 600))
    ax = Axis(fig[2, 1], xlabel = rich("Heavy-fragment mass ", it("A"), subscript("H")),
        ylabel = rich(it("Q"), " [MeV]"),)
    s = scatter!(ax, A_split, Q_split; color = (PALETTE.sky, 0.8), marker = :diamond,
        markersize = MARKERSIZE.dense, strokecolor = PALETTE.blue,)
    # σ_Q stays below 0.09 MeV, a third of the line width on this scale: no band is drawn
    l = lines!(ax, A, Q; color = PALETTE.blue)
    hlines!(ax, [Q_mean]; color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    text!(ax, last(A), Q_mean;
        text = rich("Yield-weighted mean ", @sprintf("%.1f MeV", Q_mean)),
        align = (:right, :bottom), offset = (-6, 4), fontsize = ANNOTATION_SIZE,)
    xlims!(ax, first(A) - 1, last(A) + 1)

    Legend(fig[1, 1], [s, l],
        [rich(it("Q"), "(", it("A"), subscript("H"), ", ", it("Z"), subscript("H"), ")"),
            rich(it("Q"), "(", it("A"), subscript("H"), "), charge-averaged"),],)
    println("wrote ", savefigure(fig, FIGURES, "fission_q_value"))
end

main()
