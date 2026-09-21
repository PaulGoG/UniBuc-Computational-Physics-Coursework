# Mean binding energy per nucleon across the chart of nuclides, and the neutron,
# proton, deuteron and α separation energies, from the AME mass evaluations.
#
# Ported from Radionuclizi_1.jl and Radionuclizi_2_1.jl.
#
# Corrections:
#
#   1. **A broken existence check.** `Radionuclizi_1.jl` tested
#      `isassigned(df2.Z[df2.A .== i], j - Z_min + 1)` where `Z_min` came from
#      the *other* library — so it verified only that library 2 had at least that
#      many isobars at mass i, never that the nuclide (i, j) was present, and
#      then indexed it unguarded. It survived only because AME95 is a strict
#      subset of AME2021 for the call order used; swapping the arguments crashes
#      it. Verified here: 0 of 2931 AME95 nuclides are absent from AME2021.
#   2. **Dimensionally inconsistent error propagation.** For ε = 100|D₁−D₂|/D₁
#      the partials are 100·D₂/D₁² and −100/D₁, but the file wrote
#      `sqrt((1 + D2/D1^2)^2 σ₁² + (1 - 1/D1)^2 σ₂²)`, adding a dimensionless 1
#      to quantities carrying keV⁻¹. The result was never plotted — the `yerr`
#      line is commented out — so it was dead weight as well as wrong.
#   3. `Grafic_fitare_simplu_neutron(librarie, scalare)` ignored its own argument
#      and read the global `audi95`, while labelling the plot with `librarie`;
#      passing AME2021 would have produced an AME95 fit labelled AME2021.
#   4. The CSV was re-parsed on every call — six full parses of the same file.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "mass_tables.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

function main()
    ame95 = load_masses(joinpath(DATA, "AUDI95.csv"))
    ame21 = load_masses(joinpath(DATA, "AUDI2021.csv"))
    @printf("AME1995 %d nuclides, AME2021 %d nuclides\n", length(ame95), length(ame21))

    missing_from_21 = count(k -> !haskey(ame21, k), keys(ame95))
    @printf("AME95 nuclides absent from AME2021: %d — the subset relation the original\n",
        missing_from_21)
    @printf("existence check silently depended on\n\n")

    # binding energy per nucleon
    B = [(n.A, binding_energy(ame21, n.Z, n.A) / 1000)
         for n in values(ame21)
         if binding_energy(ame21, n.Z, n.A) !== nothing && n.A >= 2]
    As = first.(B)
    Bs = last.(B)
    peak = argmax(Bs)
    @printf("maximum B/A = %.4f MeV at A = %d\n", Bs[peak], As[peak])
    fe56 = binding_energy(ame21, 26, 56) / 1000
    ni62 = binding_energy(ame21, 28, 62) / 1000
    @printf("  ⁵⁶Fe %.4f MeV,  ⁶²Ni %.4f MeV  (⁶²Ni is the true maximum)\n\n", fe56, ni62)

    # separation energies
    channels = [("neutron", 0, 1), ("proton", 1, 1), ("deuteron", 1, 2), ("α", 2, 4)]
    series = map(channels) do (name, Zx, Ax)
        pts = [(n.A - n.Z, separation_energy(ame21, n.Z, n.A, Zx, Ax) / 1000)
               for n in values(ame21)
               if separation_energy(ame21, n.Z, n.A, Zx, Ax) !== nothing && n.A >= 4]
        @printf("S_%-9s %5d nuclides, median %6.3f MeV, range %7.3f to %6.3f\n",
            name, length(pts), median(last.(pts)),
            minimum(last.(pts)), maximum(last.(pts)))
        (name = name, N = first.(pts), S = last.(pts))
    end

    fig = Figure(size = (1020, 450))
    ax1 = Axis(fig[2, 1], xlabel = L"Mass number $A$",
        ylabel = L"$B/A$ [MeV]",)
    scatter!(ax1, As, Bs, color = (PALETTE.blue, 0.35), markersize = MARKERSIZE.cloud)
    scatter!(ax1, [56, 62], [fe56, ni62], color = PALETTE.red,
        markersize = MARKERSIZE.emphasis,)
    # Both markers are labelled and both labels are above the curve: the single
    # label sat half an MeV below its marker, inside the data cloud, and did not
    # say which of the two nuclides it named.
    text!(ax1, 82, ni62 + 0.72;
        text = rich(superscript("56"), "Fe, ",
            @sprintf("%.4f MeV", fe56)),
        color = PALETTE.red, align = (:left, :bottom), fontsize = 15,)
    text!(ax1, 82, ni62 + 0.26;
        text = rich(superscript("62"), "Ni, ",
            @sprintf("%.4f MeV", ni62)),
        color = PALETTE.red, align = (:left, :bottom), fontsize = 15,)
    lines!(ax1, [80, 56], [ni62 + 0.80, fe56 + 0.10], color = PALETTE.red, linewidth = 1.0)
    lines!(ax1, [80, 62], [ni62 + 0.34, ni62 + 0.10], color = PALETTE.red, linewidth = 1.0)
    # The whole of this curve lies between 7 and 8.8 MeV, and from zero its
    # structure -- the maximum included -- occupied the top tenth of the panel.
    ylims!(ax1, 6.9, 9.9)

    ax2 = Axis(fig[2, 2], xlabel = L"Neutron number $N$",
        ylabel = L"Separation energy $S$ [MeV]",)
    cols = (PALETTE.blue, PALETTE.orange, PALETTE.green, PALETTE.purple)
    for (k, s) in enumerate(series)
        scatter!(ax2, s.N, s.S, markersize = MARKERSIZE.cloud - 2,
            color = (cols[k], 0.35),)
        # The running median: four overlapping clouds of three thousand points
        # hide the shell steps this panel exists to show.
        Ns = sort(unique(s.N))
        med = [median(s.S[s.N .== n]) for n in Ns]
        lines!(ax2, Ns, med, color = cols[k], linewidth = 2.0)
    end
    ylims!(ax2, -20, 42)
    for magic in (28, 50, 82, 126)
        vlines!(ax2, [magic], color = PALETTE.black, linestyle = :dot, linewidth = 1.1)
        text!(ax2, magic + 2, 41; text = string(magic),
            align = (:left, :top), fontsize = 14, color = PALETTE.black,)
    end
    text!(ax2, 0.98, 0.03;
        text = rich("Dotted: neutron shell closures\n",
            "Lines are the running median at each ", it("N"),),
        space = :relative, align = (:right, :bottom), fontsize = 14,
        justification = :right,)

    # Built by hand: a legend entry taken from the scatter would inherit the
    # three-pixel, 35 %-opacity marker the clouds need and be unreadable.
    Legend(fig[1, 1:2],
        [[MarkerElement(color = c, marker = :circle, markersize = MARKERSIZE.key),
             LineElement(color = c, linewidth = 2.0),] for c in cols],
        [s.name for s in series],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 24,)
    rowsize!(fig.layout, 2, Relative(0.85))
    println("\nwrote ", savefigure(fig, FIGURES, "binding_and_separation"))
end

main()
