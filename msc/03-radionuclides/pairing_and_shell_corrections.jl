# Two microscopic corrections to the liquid-drop picture from the AME1995
# masses: the nucleonic pairing gap from differences of separation energies,
# and the shell correction against a liquid-drop binding energy, compared with
# the microscopic correction of the finite-range droplet model.
#
# Pairing indicators, on S = Sₙ or Sₚ of neighbouring nuclides:
#
#   three-point in S (Guttormsen):  Δ  = ¼ |S(A+1) − 2S(A) + S(A−1)|
#   two-point in S (Vlăduca):       Δ′ = |S(A) − S(A−1)|
#
# With S(A) = B(A) − B(A−1) the second is the second difference of the binding
# energy, |B(A) − 2B(A−1) + B(A−2)| = 2Δ⁽³⁾, and the first is the third
# difference over four, Δ⁽⁴⁾. For a pure staggering B = B̄ ± Δ/2 the first
# returns Δ and the second 2Δ, so the two-point indicator is read against twice
# the guide of the three-point one. `main` asserts Δ′ = |B(A) − 2B(A−1) + B(A−2)|
# numerically. Nuclides are classed as in the 2021 file: odd Z with even N → Δₙ;
# even Z with odd N → Δₚ; even–even → Δₙ + Δₚ; against the guides 12/√A and
# 24/√A (three-point) and 24/√A and 48/√A (two-point).
#
# Shell correction, with Δ(Z, A) the mass excess:
#
#   δW₀ = W_LDM − W_exp
#   W_LDM = a_v A − a_s A^{2/3} − a_c Z² A^{−1/3} − A (a_sym − a_ss A^{−1/3}) I²
#   P_a = ½ [Δ(Z+1, A+2) − 2Δ(Z, A) + Δ(Z−1, A−2)],       δW = δW₀ + ½ P_a
#
# as in Radionuclizi_5.jl on the `legacy` branch (Julia-Workflow-FFUB/
# Radionuclizi_M_1/), which names the coefficients Pearson's parametrisation
# and takes every nuclide whose two P_a neighbours exist. δW is compared with
# the microscopic correction of the FRDM (Möller, Nix, Myers and Swiatecki,
# At. Data Nucl. Data Tables 59, 185 (1995), doi:10.1006/adnd.1995.1002) on the
# nuclides the two have in common.
#
# Radionuclizi_3.jl read a non-`const` global table inside a function holding
# six full boolean scans, called six times per (A, Z); its `Energie_separare`
# returned `[true, S]` or `[false, 0]`, two element types from one function;
# and its guides were drawn against an unsorted A, as a scribble.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "mass_tables.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"""
Liquid-drop coefficients [MeV]: the set of von Egidy and Bucurescu, Phys. Rev.
C 72, 044311 (2005), Eq. (9), taken there from Pearson, Hyperfine Interact.
132, 59 (2001); a_c = (3/5) e²/r₀ with e² = 1.44 MeV fm and r₀ = 1.233 fm.
Radionuclizi_5.jl set the same values, naming Pearson without the reference.
"""
const LDM = (a_v = 15.65, a_s = 17.63, a_sym = 27.72, a_ss = 25.60, a_c = 0.864 / 1.233)

"Pairing gap of the semi-empirical mass formula [MeV], 12/√A."
const PAIRING_CONSTANT = 12.0

"Upper limit of the pairing panels [MeV]; points above it are counted, not drawn."
const PAIRING_FRAME = (guttormsen = 5.0, vladuca = 10.0)

"Relative tolerance of the two-point identity."
const IDENTITY_TOLERANCE = 1e-9

"Liquid-drop binding energy without pairing [MeV]."
function W_ldm(Z, A)
    I = (A - 2Z) / A
    return LDM.a_v * A - LDM.a_s * A^(2 / 3) - LDM.a_c * Z^2 * A^(-1 / 3) -
           A * (LDM.a_sym - LDM.a_ss * A^(-1 / 3)) * I^2
end

"Experimental binding energy [MeV] from the mass excesses, or `nothing`."
function W_exp(t, Z, A)
    b = binding_energy(t, Z, A)
    return b === nothing ? nothing : b * A / 1000
end

"""
    pairing_term(t, Z, A)

½[Δ(Z+1, A+2) − 2Δ(Z, A) + Δ(Z−1, A−2)] in MeV, the second difference of the
mass excess along the (Z+1, N+1) direction that Radionuclizi_5.jl adds, halved,
to the shell correction; `nothing` if a neighbour is absent.
"""
function pairing_term(t, Z, A)
    (haskey(t, (Z + 1, A + 2)) && haskey(t, (Z, A)) && haskey(t, (Z - 1, A - 2))) ||
        return nothing
    return 0.5 * (t[(Z + 1, A + 2)].Δ - 2t[(Z, A)].Δ + t[(Z - 1, A - 2)].Δ) / 1000
end

"""
    pairing_gap(t, Z, A, kind, formula)

Pairing indicator in MeV from the separation energies of three neighbouring
nuclides along the neutron (`kind = :neutron`) or proton (`:proton`) line, or
`nothing` if any of the three is absent. `formula` is `:guttormsen`,
¼|S(A+1) − 2S(A) + S(A−1)|, or `:vladuca`, |S(A) − S(A−1)|; both require the
same three nuclides so that they are evaluated on one set.
"""
function pairing_gap(t, Z, A, kind, formula)
    formula in (:guttormsen, :vladuca) ||
        throw(ArgumentError("formula must be :guttormsen or :vladuca, got $formula"))
    Zx, Ax = kind === :neutron ? (0, 1) : (1, 1)
    s = map((-1, 0, 1)) do k
        separation_energy(t, Z + k * Zx, A + k * Ax, Zx, Ax)
    end
    any(isnothing, s) && return nothing
    formula === :guttormsen && return abs(s[3] - 2s[2] + s[1]) / 4 / 1000
    return abs(s[2] - s[1]) / 1000
end

"""
    pairing_by_class(t, formula)

The three parity classes of the 2021 file: `:odd_Z` (odd Z, even N) carries Δₙ,
`:odd_N` (even Z, odd N) carries Δₚ, `:even_even` carries Δₙ + Δₚ. Returns a
`Dict` of `(A, value)` vectors.
"""
function pairing_by_class(t, formula)
    classes = Dict(:even_even => Tuple{Int,Float64}[], :odd_Z => Tuple{Int,Float64}[],
        :odd_N => Tuple{Int,Float64}[],)
    for n in values(t)
        N = n.A - n.Z
        Δn = pairing_gap(t, n.Z, n.A, :neutron, formula)
        Δp = pairing_gap(t, n.Z, n.A, :proton, formula)
        if isodd(n.Z) && iseven(N)
            Δn === nothing || push!(classes[:odd_Z], (n.A, Δn))
        elseif iseven(n.Z) && isodd(N)
            Δp === nothing || push!(classes[:odd_N], (n.A, Δp))
        elseif iseven(n.Z) && iseven(N)
            (Δn === nothing || Δp === nothing) || push!(classes[:even_even], (n.A, Δn + Δp))
        end
    end
    return classes
end

"""
    two_point_identity(t)

Largest relative departure of |Sₙ(A) − Sₙ(A−1)| from |W(A) − 2W(A−1) + W(A−2)|
over the table, W the binding energy: zero up to round-off.
"""
function two_point_identity(t)
    worst = 0.0
    for n in values(t)
        two_point = pairing_gap(t, n.Z, n.A, :neutron, :vladuca)
        two_point === nothing && continue
        w = (W_exp(t, n.Z, n.A), W_exp(t, n.Z, n.A - 1), W_exp(t, n.Z, n.A - 2))
        any(isnothing, w) && continue
        second_difference = abs(w[1] - 2w[2] + w[3])
        worst = max(worst, abs(two_point - second_difference) /
                           max(second_difference, 1e-6))
    end
    return worst
end

"Class label, colour and marker of each parity class."
const CLASSES = (
    (key = :even_even, label = "Even–even", colour = PALETTE.blue, marker = :circle),
    (key = :odd_Z, label = rich("Odd ", it("Z"), ", even ", it("N")),
        colour = PALETTE.orange,
        marker = :utriangle,),
    (key = :odd_N, label = rich("Even ", it("Z"), ", odd ", it("N")),
        colour = PALETTE.green,
        marker = :diamond,),
    (key = :odd_odd, label = "Odd–odd", colour = PALETTE.purple, marker = :rect),
)

function main()
    t = load_masses(joinpath(DATA, "AUDI95.csv"))
    moller = load_moller(joinpath(DATA, "MOLLER.csv"))
    @printf("AME1995 %d nuclides, FRDM microscopic corrections %d entries\n\n",
        length(t), length(moller))

    # --- pairing ---------------------------------------------------------------
    identity = two_point_identity(t)
    @printf("two-point indicator against the second difference of W: largest relative departure %.1e\n",
        identity)
    identity < IDENTITY_TOLERANCE || error("the two-point identity fails: $identity")

    gaps = Dict(f => pairing_by_class(t, f) for f in (:guttormsen, :vladuca))
    println("\npairing gap: median of Δ√A per class [MeV], guide 12 (odd A) and 24 (even–even)")
    for f in (:guttormsen, :vladuca), c in CLASSES[1:3]

        pts = gaps[f][c.key]
        A, Δ = first.(pts), last.(pts)
        @printf("  %-11s %-22s n = %4d   median Δ√A = %5.2f   median Δ = %.3f\n",
            string(f), c.key, length(pts), median(Δ .* sqrt.(A)), median(Δ))
    end

    # --- shell corrections -------------------------------------------------------
    rows = NamedTuple{(:Z, :A, :δW₀, :P_a, :δW),NTuple{5,Float64}}[]
    for n in values(t)
        n.A < 3 && continue
        p = pairing_term(t, n.Z, n.A)
        w = W_exp(t, n.Z, n.A)
        (p === nothing || w === nothing) && continue
        δW₀ = W_ldm(n.Z, n.A) - w
        push!(rows, (Z = n.Z, A = n.A, δW₀ = δW₀, P_a = p, δW = δW₀ + 0.5p))
    end
    common = [r for r in rows if haskey(moller, (r.Z, r.A))]
    mn = [moller[(r.Z, r.A)] for r in common]
    δW = [r.δW for r in common]
    δW₀ = [r.δW₀ for r in common]
    @printf("\nshell correction on %d nuclides with both P_a neighbours, %d in common with the FRDM\n",
        length(rows), length(common))
    @printf("  δW₀ = W_LDM − W_exp:   mean %+.3f MeV, r against FRDM %.3f, rms difference %.2f MeV\n",
        mean(δW₀), cor(δW₀, mn), sqrt(mean(abs2, δW₀ .- mn)))
    @printf("  δW = δW₀ + ½P_a:       mean %+.3f MeV, r against FRDM %.3f, rms difference %.2f MeV\n",
        mean(δW), cor(δW, mn), sqrt(mean(abs2, δW .- mn)))
    @printf("  FRDM E_mic:            mean %+.3f MeV, range %+.2f to %+.2f\n",
        mean(mn), minimum(mn), maximum(mn))
    println("  P_a by parity class, median [MeV]:")
    parity(r) = iseven(r.Z) ? (iseven(r.A - r.Z) ? :even_even : :odd_N) :
                (iseven(r.A - r.Z) ? :odd_Z : :odd_odd)
    for c in CLASSES
        sel = [r.P_a for r in rows if parity(r) == c.key]
        @printf("    %-10s n = %4d   median %+.3f\n", c.key, length(sel), median(sel))
    end
    magic = ((82, 208, "²⁰⁸Pb"), (50, 132, "¹³²Sn"), (28, 78, "⁷⁸Ni"))
    for (Z, A, name) in magic
        k = findfirst(r -> r.Z == Z && r.A == A, common)
        k === nothing && continue
        @printf("  %s: δW = %+.2f MeV, FRDM %+.2f MeV\n", name, common[k].δW, mn[k])
    end

    # --- figure -------------------------------------------------------------------
    fig = Figure(size = (1500, 1180))
    pairing_axes = (
        Axis(fig[2, 1], ylabel = L"Pairing indicator $\Delta$ [MeV]",
            xticks = 0:50:250,),
        Axis(fig[2, 2], xticks = 0:50:250),)
    for (ax, f, factor) in zip(pairing_axes, (:guttormsen, :vladuca), (1, 2))
        frame = PAIRING_FRAME[f]
        above = 0
        for c in CLASSES[1:3]
            pts = gaps[f][c.key]
            A, Δ = first.(pts), last.(pts)
            above += count(>(frame), Δ)
            scatter!(ax, A, Δ, color = (c.colour, 0.4), marker = c.marker,
                markersize = MARKERSIZE.cloud, strokewidth = 0,)
        end
        Af = 2:270
        for (m, tag) in ((1, "odd \$A\$"), (2, "even--even"))
            g = factor * m * PAIRING_CONSTANT ./ sqrt.(Af)
            lines!(ax, Af, g, color = PALETTE.black, linestyle = m == 1 ? :dash : :dashdot,
                linewidth = GUIDE_WIDTH,)
            text!(ax, 268, last(g);
                text = latexstring("\$$(factor * m * 12)/\\sqrt{A}\$, ", tag),
                align = (:right, :bottom), offset = (0, 4), fontsize = ANNOTATION_SIZE,)
        end
        text!(ax, 0.97, 0.96;
            text = rich(f === :guttormsen ? "Three-point in " : "Two-point in ", it("S"),
                @sprintf("\n%d points above the frame", above)),
            space = :relative, align = (:right, :top), justification = :right,
            fontsize = ANNOTATION_SIZE,)
        xlims!(ax, 0, 275)
        ylims!(ax, 0, frame)
    end
    pairing_axes[1].xlabel = L"Mass number $A$"
    pairing_axes[2].xlabel = L"Mass number $A$"

    ax3 = Axis(fig[3, 1], xlabel = L"Mass number $A$", ylabel = L"$P_a$ [MeV]",
        xticks = 0:50:250,)
    for c in CLASSES
        sel = [r for r in rows if parity(r) == c.key]
        scatter!(ax3, [r.A for r in sel], [r.P_a for r in sel], color = (c.colour, 0.4),
            marker = c.marker, markersize = MARKERSIZE.cloud, strokewidth = 0,)
    end
    hlines!(ax3, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    xlims!(ax3, 0, 275)
    ylims!(ax3, -8.5, 8.5)

    ax4 = Axis(fig[3, 2], xlabel = L"FRDM microscopic correction $E_\mathrm{mic}$ [MeV]",
        ylabel = L"$\delta W = \delta W_0 + P_a/2$ [MeV]",)
    scatter!(ax4, mn, δW, color = (PALETTE.sky, 0.4), markersize = MARKERSIZE.cloud,
        strokewidth = 0,)
    lo, hi = floor(min(minimum(mn), minimum(δW))), ceil(max(maximum(mn), maximum(δW)))
    lines!(ax4, [lo, hi], [lo, hi], color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    text!(ax4, hi - 0.5, hi - 0.5; text = "Equality", align = (:right, :top),
        fontsize = ANNOTATION_SIZE,)
    # ²⁰⁸Pb and ¹³²Sn lie within 1.5 MeV of each other: one label each side
    placements = (
        ("²⁰⁸Pb", (:center, :bottom), (0, 12)), ("¹³²Sn", (:left, :top), (10, -8)),
        ("⁷⁸Ni", (:left, :bottom), (10, 8)),)
    for ((Z, A, name), (_, align, offset)) in zip(magic, placements)
        k = findfirst(r -> r.Z == Z && r.A == A, common)
        k === nothing && continue
        scatter!(ax4, [mn[k]], [δW[k]], color = PALETTE.red, marker = :xcross,
            markersize = MARKERSIZE.emphasis,)
        text!(ax4, mn[k], δW[k]; text = name, align = align, offset = offset,
            fontsize = ANNOTATION_SIZE, color = PALETTE.red,)
    end
    text!(ax4, 0.04, 0.96;
        text = rich(it("r"), @sprintf(" = %.3f on %d nuclides\n", cor(δW, mn), length(mn)),
            @sprintf("rms difference %.2f MeV\n", sqrt(mean(abs2, δW .- mn))),
            rich("✕ doubly magic", color = PALETTE.red),),
        space = :relative, align = (:left, :top), fontsize = ANNOTATION_SIZE,)
    limits!(ax4, lo, hi, lo, hi)

    Legend(fig[1, 1:2],
        [MarkerElement(color = c.colour, marker = c.marker, markersize = MARKERSIZE.key)
         for c in CLASSES],
        [c.label for c in CLASSES], "Parity class"; titleposition = :left,)
    println("\nwrote ", savefigure(fig, FIGURES, "pairing_and_shell_corrections"))
end

main()
