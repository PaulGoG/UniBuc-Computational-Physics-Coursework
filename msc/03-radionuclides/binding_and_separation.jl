# Mean binding energy per nucleon across the chart of nuclides, and the neutron,
# proton, deuteron and α separation energies, from the AME mass evaluations.
#
# Ported from Radionuclizi_1.jl and Radionuclizi_2_1.jl on the `legacy` branch
# (Julia-Workflow-FFUB/Radionuclizi_M_1/). Radionuclizi_1.jl tested the presence
# of a nuclide in the second table with
# `isassigned(df2.Z[df2.A .== i], j - Z_min + 1)`, Z_min taken from the first
# table, which checks only that the second table has that many isobars at mass
# i; it worked because every AME1995 nuclide is in AME2020, which `main`
# verifies. Both files re-parsed the CSV on every call.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "mass_tables.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

function main()
    ame95 = load_masses(joinpath(DATA, "AUDI95.csv"))
    ame21 = load_masses(joinpath(DATA, "AUDI2021.csv"))
    @printf("AME1995 %d nuclides, AME2020 %d nuclides\n", length(ame95), length(ame21))

    missing_from_21 = count(k -> !haskey(ame21, k), keys(ame95))
    @printf("AME1995 nuclides absent from AME2020: %d\n\n", missing_from_21)
    missing_from_21 == 0 ||
        error("$missing_from_21 AME1995 nuclides are missing from AME2020")

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

    fig = Figure(size = (1400, 640))
    ax1 = Axis(fig[2, 1], xlabel = L"Mass number $A$", ylabel = L"$B/A$ [MeV]")
    scatter!(ax1, As, Bs, color = (PALETTE.blue, 0.35), markersize = MARKERSIZE.cloud,
        strokewidth = 0,)
    scatter!(ax1, [56, 62], [fe56, ni62], color = PALETTE.red,
        markersize = MARKERSIZE.emphasis,)
    text!(ax1, 82, ni62 + 0.62;
        text = rich(superscript("56"), "Fe  ", @sprintf("%.4f MeV", fe56)),
        color = PALETTE.red, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)
    text!(ax1, 82, ni62 + 0.20;
        text = rich(superscript("62"), "Ni  ", @sprintf("%.4f MeV", ni62)),
        color = PALETTE.red, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)
    lines!(ax1, [80, 56], [ni62 + 0.70, fe56 + 0.10], color = PALETTE.red,
        linewidth = GUIDE_WIDTH,)
    lines!(ax1, [80, 62], [ni62 + 0.28, ni62 + 0.10], color = PALETTE.red,
        linewidth = GUIDE_WIDTH,)
    # the curve lies between 7 and 8.8 MeV; from zero its structure would fill
    # only the top tenth of the panel
    ylims!(ax1, 6.9, 9.9)

    ax2 = Axis(fig[2, 2], xlabel = L"Neutron number $N$",
        ylabel = L"Separation energy $S$ [MeV]",)
    cols = (PALETTE.blue, PALETTE.orange, PALETTE.green, PALETTE.purple)
    for (k, s) in enumerate(series)
        scatter!(ax2, s.N, s.S, markersize = MARKERSIZE.cloud, strokewidth = 0,
            color = (cols[k], 0.35),)
        # the running median: four overlapping clouds hide the shell steps
        Ns = sort(unique(s.N))
        med = [median(s.S[s.N .== n]) for n in Ns]
        lines!(ax2, Ns, med, color = cols[k])
    end
    ylims!(ax2, -20, 42)
    for magic in (28, 50, 82, 126)
        vlines!(
            ax2, [magic], color = PALETTE.black, linestyle = :dot, linewidth = GUIDE_WIDTH,)
        # the first label goes left of its line so that it clears the N = 50 one
        left = magic == 28
        text!(ax2, magic + (left ? -2 : 2), 41; text = latexstring("N = $magic"),
            align = (left ? :right : :left, :top), fontsize = ANNOTATION_SIZE,
            color = PALETTE.black,)
    end
    text!(ax2, 0.98, 0.03; text = rich("Lines: running median at each ", it("N")),
        space = :relative, align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)

    # built by hand: an entry taken from the scatter would inherit the cloud marker
    Legend(fig[1, 1:2],
        [
            [[MarkerElement(color = c, marker = :circle, markersize = MARKERSIZE.key),
                 LineElement(color = c),] for c in cols],
            [LineElement(color = PALETTE.black, linestyle = :dot, linewidth = GUIDE_WIDTH)],],
        [["Neutron", "Proton", "Deuteron", "α"], ["Neutron shell closure"]],
        ["Separation energy", "Guide"]; titleposition = :left, nbanks = 1, groupgap = 36,)
    println("\nwrote ", savefigure(fig, FIGURES, "binding_and_separation"))
end

main()
