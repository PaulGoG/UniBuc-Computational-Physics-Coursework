# Two microscopic corrections to the liquid-drop picture, from the AME1995
# masses: the nucleonic pairing gap, and the shell correction against the
# Pearson liquid-drop formula benchmarked to Möller–Nix FRDM.
#
# Pairing gap, three-point (Guttormsen) indicator:
#     Δ = ¼|S(A+1) − 2S(A) + S(A−1)|
#
# Shell correction:
#     δW = W_LDM − W_exp,
#     W_LDM = a_v A − a_s A^{2/3} − a_c Z²A^{−1/3} − A(a_sym − a_ss A^{−1/3}) I²
#
# Ported from Radionuclizi_3.jl and Radionuclizi_5.jl.
#
# Corrections:
#
#   1. **The hottest loop in the archive.** `Radionuclizi_3.jl` read a
#      non-`const` global `df` inside a function containing six full boolean
#      scans of a 2931-row table, called six times per (A, Z) over ~10⁵ pairs —
#      of order 2×10⁹ type-unstable comparisons. Replaced by a dictionary lookup.
#   2. **A type-unstable return.** `Energie_separare` returned `[true, S]`
#      (Vector{Float64}) or `[false, 0]` (Vector{Int64}), then indexed `[1]` as a
#      boolean and `[2]` as an energy. Now `nothing` or a `Float64`.
#   3. **A zig-zag reference curve.** `pairing.A` is built in an order that is
#      not monotonic, and the 12/√A guide was plotted against it directly, so
#      `plot!` joined consecutive points into a scribble. Plotted against a
#      sorted range here.
#   4. `Radionuclizi_5.jl` computed `Calcul_W_exp` three times, `Calcul_W_LDM`
#      three times and `Calcul_P_a` twice per nuclide through a chain of
#      functions each re-deriving what its caller already had.
#   5. `σᴾ` was pushed as a hardcoded 1 MeV placeholder because the real
#      propagation was commented out, and was never used in any plot.

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "mass_tables.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"""
Pearson liquid-drop coefficients in MeV.

The course portfolio names the parameterisation — "corecțiile de pături bazate
pe parametrizarea Pearson a modelului picătură de lichid" — which is where this
choice comes from; the coefficients themselves are not tabulated anywhere that
survives, so they stand as the original file set them.
"""
const LDM = (a_v = 15.65, a_s = 17.63, a_sym = 27.72, a_ss = 25.60,
             a_c = 0.864 / 1.233)      # (3/5)e²/r₀ with e² = 1.44 MeV·fm, r₀ = 1.233 fm

"Liquid-drop binding energy without pairing, in MeV."
function W_ldm(Z, A)
    I = (A - 2Z) / A
    return LDM.a_v * A - LDM.a_s * A^(2/3) - LDM.a_c * Z^2 * A^(-1/3) -
           A * (LDM.a_sym - LDM.a_ss * A^(-1/3)) * I^2
end

"Experimental binding energy in MeV from the mass excesses."
function W_exp(t, Z, A)
    b = binding_energy(t, Z, A)
    return b === nothing ? nothing : b * A / 1000
end

"""
    pairing_gap(t, Z, A, kind, formula)

Pairing indicator from the separation energies of neighbouring nuclides, in MeV.
`kind` is `:neutron` or `:proton`. Both indicators the original offered are
available: `:guttormsen`, the three-point form |¼[S(A+1) − 2S(A) + S(A−1)]|, and
`:vladuca`, the two-point |S(A) − S(A−1)|.
"""
function pairing_gap(t, Z, A, kind, formula)
    Zx, Ax = kind === :neutron ? (0, 1) : (1, 1)
    dZ, dA = kind === :neutron ? (0, 1) : (1, 1)
    s = map((-1, 0, 1)) do k
        separation_energy(t, Z + k*dZ, A + k*dA, Zx, Ax)
    end
    any(isnothing, s) && return nothing
    if formula === :guttormsen
        return abs(s[3] - 2s[2] + s[1]) / 4 / 1000
    elseif formula === :vladuca
        return abs(s[2] - s[1]) / 1000
    else
        throw(ArgumentError("formula must be :guttormsen or :vladuca, got $formula"))
    end
end

function main()
    t = load_masses(joinpath(DATA, "AUDI95.csv"))
    moller = load_moller(joinpath(DATA, "MOLLER.csv"))
    @printf("AME1995 %d nuclides, Möller–Nix %d entries\n", length(t), length(moller))

    gap_sets = Dict{Symbol,Vector{Tuple{Int,Float64}}}()
    for formula in (:guttormsen, :vladuca)
        g_list = Tuple{Int,Float64}[]
        for n in values(t)
            n.A < 4 && continue
            g = pairing_gap(t, n.Z, n.A, :neutron, formula)
            g !== nothing && 0 < g < 5 && push!(g_list, (n.A, g))
        end
        gap_sets[formula] = g_list
        resid = [g - 12/sqrt(A) for (A, g) in g_list]
        @printf("neutron pairing gap, %-11s %d nuclides, median %.3f MeV, ",
                string(formula) * ":", length(g_list), median(last.(g_list)))
        @printf("RMS residual against 12/√A %.3f MeV\n", sqrt(mean(abs2, resid)))
    end
    gaps_n = gap_sets[:guttormsen]

    δW = Tuple{Int,Float64,Float64}[]   # (A, ours, Möller–Nix)
    for n in values(t)
        n.A < 16 && continue
        w = W_exp(t, n.Z, n.A)
        w === nothing && continue
        haskey(moller, (n.Z, n.A)) || continue
        push!(δW, (n.A, W_ldm(n.Z, n.A) - w, moller[(n.Z, n.A)]))
    end
    ours = [d[2] for d in δW]; mn = [d[3] for d in δW]
    @printf("\nshell correction on %d nuclides in common\n", length(δW))
    @printf("  this LDM:     mean %+.3f MeV, range %+.2f to %+.2f\n",
            mean(ours), minimum(ours), maximum(ours))
    @printf("  Möller–Nix:   mean %+.3f MeV, range %+.2f to %+.2f\n",
            mean(mn), minimum(mn), maximum(mn))
    @printf("  correlation coefficient %.3f\n", cor(ours, mn))
    for (Z, A, name) in ((82, 208, "²⁰⁸Pb"), (50, 132, "¹³²Sn"), (28, 78, "⁷⁸Ni"))
        haskey(moller, (Z, A)) &&
            @printf("  %s: Möller–Nix δW = %+.2f MeV\n", name, moller[(Z, A)])
    end

    fig = Figure(size = (1000, 450))

    ax1 = Axis(fig[2, 1], xlabel = L"Mass number $A$",
        ylabel = L"Neutron pairing gap $\Delta_n$ [MeV]")
    scatter!(ax1, first.(gap_sets[:vladuca]), last.(gap_sets[:vladuca]),
        color = (PALETTE.orange, 0.18), markersize = 3)
    scatter!(ax1, first.(gaps_n), last.(gaps_n),
        color = (PALETTE.blue, 0.3), markersize = 3)
    text!(ax1, 0.03, 0.93; text = "blue: Guttormsen 3-point\norange: Vlăduca 2-point",
        space = :relative, align = (:left, :top), fontsize = 14)
    Af = sort(unique(first.(gaps_n)))
    l_emp = lines!(ax1, Af, 12 ./ sqrt.(Af), color = PALETTE.red, linewidth = 2)
    ylims!(ax1, 0, 3.5)
    text!(ax1, 0.97, 0.93; text = L"$12/\sqrt{A}$", space = :relative,
        align = (:right, :top), color = PALETTE.red, fontsize = 16)

    ax2 = Axis(fig[2, 2], xlabel = L"Möller–Nix $\delta W$ [MeV]",
        ylabel = L"This LDM $\delta W$ [MeV]")
    scatter!(ax2, mn, ours, color = (PALETTE.green, 0.3), markersize = 3)
    lo, hi = extrema(vcat(mn, ours))
    lines!(ax2, [lo, hi], [lo, hi], color = PALETTE.black,
        linestyle = :dash, linewidth = 1.2)
    text!(ax2, 0.04, 0.93; text = @sprintf("r = %.3f", cor(ours, mn)),
        space = :relative, align = (:left, :top), fontsize = 16)

    rowsize!(fig.layout, 2, Relative(0.88))
    println("\nwrote ", savefigure(fig, FIGURES, "pairing_and_shell_corrections"))
end

main()
