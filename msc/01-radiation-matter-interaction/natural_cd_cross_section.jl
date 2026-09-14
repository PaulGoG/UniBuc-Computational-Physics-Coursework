# Neutron capture cross-section of natural cadmium as the abundance-weighted sum
# over its eight stable isotopes,
#
#   σ_nat = Σ aᵢ σᵢ,
#
# dominated by the ¹¹³Cd resonance at 0.178 eV.
#
# Ported from SigmaCdNatural.jl. The method was right. Three things were not.
#
#   1. **The result was never printed.** `σCd = sum(...)` as a bare top-level
#      expression means `julia SigmaCdNatural.jl` produces no output whatsoever;
#      the value existed only if the file was pasted into a REPL.
#   2. **The energy was not carried into the result.** The cross-sections are an
#      IAEA set evaluated at 0.25 eV, and the answer is ≈2.20×10³ b. Any reader
#      will assume the conventional thermal value, which is 2520 b at the
#      2200 m s⁻¹ point (0.0253 eV, where σ(¹¹³Cd) = 20 600 b). The energy is now
#      part of the variable name and of the printed result.
#   3. The arrays were positional only — nothing tied σ = 18 000 b to ¹¹³Cd — and
#      `abundente_izotopice./100 # Transformam in procente` has its comment
#      backwards: dividing by 100 converts *from* percent. It was also a
#      destructive rebinding, so re-evaluating that one line in a REPL divided
#      by 100 again.

using Printf
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Energy at which the IAEA cross-section set is evaluated."
const E_EVAL_EV = 0.25

"Mass number, natural abundance in per cent, capture cross-section in barn."
const CADMIUM = [
    (A = 106, abundance = 1.25, σ = 0.30),
    (A = 108, abundance = 0.89, σ = 0.34),
    (A = 110, abundance = 12.47, σ = 3.50),
    (A = 111, abundance = 12.80, σ = 6.80),
    (A = 112, abundance = 24.11, σ = 0.60),
    (A = 113, abundance = 12.23, σ = 18000.0),
    (A = 114, abundance = 28.75, σ = 0.05),
    (A = 116, abundance = 7.51, σ = 0.015),
]

function main()
    total_abundance = sum(i.abundance for i in CADMIUM)
    @printf("abundances sum to %.2f %% (tabulated values, not renormalised)\n",
        total_abundance)

    σ_nat = sum(i.abundance / 100 * i.σ for i in CADMIUM)
    @printf("\nσ_capture(natural Cd) at %.2f eV = %.1f b\n", E_EVAL_EV, σ_nat)
    @printf("for contrast, the 2200 m/s thermal value (0.0253 eV) is 2520 b\n\n")

    println("contribution by isotope:")
    for i in CADMIUM
        c = i.abundance / 100 * i.σ
        @printf("  Cd-%3d  a = %5.2f %%  σ = %9.3f b  ->  %8.3f b  (%5.2f %% of total)\n",
            i.A, i.abundance, i.σ, c, 100c / σ_nat)
    end
    dominant = argmax([i.abundance / 100 * i.σ for i in CADMIUM])
    @printf("\nCd-%d alone supplies %.2f %% of the natural cross-section\n",
        CADMIUM[dominant].A,
        100 * CADMIUM[dominant].abundance / 100 * CADMIUM[dominant].σ / σ_nat)

    fig = Figure(size = (820, 460))
    ax = Axis(fig[2, 1],
        xlabel = "Mass number A", ylabel = "Contribution to σ [b]",
        yscale = log10,
        xticks = ([i.A for i in CADMIUM], [latexstring(string(i.A)) for i in CADMIUM]),
        yticks = logticks(-3, 4; step = 2),)
    contributions = [i.abundance / 100 * i.σ for i in CADMIUM]
    barplot!(ax, [i.A for i in CADMIUM], contributions, color = PALETTE.blue,
        strokewidth = 0.5, strokecolor = PALETTE.black,)
    hlines!(ax, [σ_nat], color = PALETTE.red, linestyle = :dash, linewidth = 1.3)
    # Headroom above the total, with the note in it: centred at 93 % of the
    # panel its last word was drawn over the 113-Cd bar, orange on saturated
    # blue.
    ylims!(ax, nothing, σ_nat * 30)
    text!(ax, 0.98, 0.97;
        text = rich(
            it("σ"), @sprintf("(natural Cd) = %.0f b at %.2f eV\n", σ_nat, E_EVAL_EV),
            @sprintf("%.2f %% of it from ",
                100 * CADMIUM[dominant].abundance / 100 *
                CADMIUM[dominant].σ / σ_nat),
            superscript("113"), "Cd alone",),
        space = :relative, align = (:right, :top), color = PALETTE.red, fontsize = 16,
        justification = :right,)

    Legend(fig[1, 1],
        [PolyElement(color = PALETTE.blue),
            LineElement(color = PALETTE.red, linestyle = :dash, linewidth = 1.3),],
        ["Isotopic contribution, abundance × σ", "Total, natural Cd"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 24,)
    rowsize!(fig.layout, 2, Relative(0.86))
    println("wrote ", savefigure(fig, FIGURES, "natural_cd_cross_section"))
end

main()
