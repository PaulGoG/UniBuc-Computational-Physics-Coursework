# Expansion history of a flat homogeneous and isotropic universe from the
# Friedmann equation, and flat ΛCDM fitted to three datasets.
#
#   H²(z) = H₀² [Ω_m(1+z)³ + 1 − Ω_m]
#
# Joint work with Alexandru Crăciun at the Laboratory of Information
# Technologies, JINR Dubna, July 2019. Ported from `Friedmann_Eq.jl` and
# `Friedmann_Magnitude.jl` (`Julia-Workflow-FFUB/Dubna_2019_Cosmology/` on the
# `legacy` branch).
#
#   1. H(z): the twelve rows of `data/hubble_parameter.csv`, weighted by their
#      uncertainties. The z = 0 row is a local value of H₀, not a measurement of
#      H(z), so the fit is also made without it.
#   2. a(t): the Friedmann equation integrated with the scale factor as the
#      independent variable, dt/da = 1/(a H), which is regular at a = 0 and so
#      reaches the Big Bang instead of stepping through it. Checked against the
#      closed form a ∝ sinh^{2/3}(3√Ω_Λ H₀t/2) and the ages it implies.
#   3. μ(z): 2376 Type Ia distance moduli reduced from NED-D, unweighted because
#      the reduction kept no uncertainties, refitted above a rising redshift
#      floor to show how far Ω_m depends on the sample.
#   4. The host-galaxy magnitudes that `Friedmann_Magnitude.jl` fitted, which
#      fix an offset and nothing else.
#
# `Friedmann_Eq.jl` fitted a straight line to H(z) and, run top to bottom,
# stops at its first plot, which uses `fit.param` before `fit` exists; its
# Friedmann right-hand side mixes G in SI units with H in km s⁻¹ Mpc⁻¹.
# `Friedmann_Magnitude.jl` assigns `u[3] = …` inside the right-hand side of its
# ODE, writing a magnitude into the state vector instead of a derivative, and
# reads its data from an absolute Windows path.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames, LsqFit, OrdinaryDiffEq, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "cosmology_core.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Lower redshift cuts applied to the distance-modulus sample."
const REDSHIFT_FLOORS = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12]

"Matter densities of the three expansion histories drawn: Planck 2018, Einstein–de Sitter, Λ-dominated."
const HISTORIES = (PLANCK2018.Ωm, 1.0, 0.05)

"Largest scale factor the expansion histories are integrated to."
const A_MAX = 4.0
"Span of future time drawn [Gyr]."
const T_FUTURE = 15.0

"Tolerance on the integrated τ(a) against its closed form, in units of 1/H₀."
const TIME_TOLERANCE = 1e-8

"Planck 2018 age of the universe [Gyr], same reference as `PLANCK2018`."
const PLANCK2018_AGE = (t₀ = 13.797, σ = 0.023)

hubble_rate(z, p) = p[1] .* expansion_rate.(z, p[2])

"""
    fit_hubble_rate(z, H, σH)

Weighted fit of H₀ and Ω_m to measurements of H(z). Returns the parameters,
their standard errors, χ² and the degrees of freedom.
"""
function fit_hubble_rate(z, H, σH)
    fit = curve_fit(hubble_rate, z, H, 1 ./ σH .^ 2, [70.0, 0.3])
    σ = stderror(fit)
    return (; H₀ = fit.param[1], Ωm = fit.param[2], σH₀ = σ[1], σΩm = σ[2],
        χ² = sum(abs2, fit.resid), dof = length(z) - 2,)
end

"""
    expansion_history(Ωm; a_max = A_MAX, n = 400)

Integrate dτ/da = 1/(a E) = √a / √(Ω_m + Ω_Λ a³), τ = H₀t, from today (a = 1,
τ = 0) down to a = 0 and up to `a_max`. Returns `(τ, a)` in increasing order;
`-τ[1]` is the age in units of 1/H₀.
"""
function expansion_history(Ωm; a_max = A_MAX, n = 400)
    rhs(τ, p, a) = sqrt(a) / sqrt(Ωm + (1 - Ωm) * a^3)
    grid_past = range(1.0, 0.0, length = n)
    grid_future = range(1.0, a_max, length = n)
    past = solve(ODEProblem(rhs, 0.0, (1.0, 0.0)), Tsit5();
        saveat = grid_past, reltol = 1e-10, abstol = 1e-12,)
    future = solve(ODEProblem(rhs, 0.0, (1.0, a_max)), Tsit5();
        saveat = grid_future, reltol = 1e-10, abstol = 1e-12,)
    τ = vcat(reverse(past.u), future.u[2:end])
    a = vcat(reverse(past.t), future.t[2:end])
    return τ, a
end

"Median μ in logarithmic redshift bins holding at least `minimum_count` points."
function binned_medians(z, μ; nbins = 20, minimum_count = 15)
    edges = 10 .^ range(log10(minimum(z)), log10(maximum(z)), length = nbins + 1)
    bz, bμ = Float64[], Float64[]
    for i in 1:nbins
        sel = (z .>= edges[i]) .& (z .< edges[i + 1])
        count(sel) < minimum_count && continue
        push!(bz, median(z[sel]))
        push!(bμ, median(μ[sel]))
    end
    return bz, bμ
end

"Rows of the SAI catalogue extract with a numeric host magnitude and 0.001 < z < 2."
function load_host_galaxy_magnitudes(path)
    raw = CSV.read(path, DataFrame; header = ["name", "mag", "z"],
        skipto = 3, silencewarnings = true,)
    zs, ms = Float64[], Float64[]
    for r in eachrow(raw)
        (ismissing(r.mag) || ismissing(r.z)) && continue
        zv = tryparse(Float64, strip(string(r.z)))
        mv = tryparse(Float64, strip(string(r.mag)))
        (zv === nothing || mv === nothing) && continue
        (0.001 < zv < 2.0 && mv > 0) || continue
        push!(zs, zv)
        push!(ms, mv)
    end
    return zs, ms
end

function main()
    # --- H(z) -----------------------------------------------------------------
    hz = CSV.read(joinpath(DATA, "hubble_parameter.csv"), DataFrame;
        header = ["z", "H", "σH"],)
    measured = hz.z .> 0
    all_rows = fit_hubble_rate(hz.z, hz.H, hz.σH)
    no_prior = fit_hubble_rate(hz.z[measured], hz.H[measured], hz.σH[measured])
    # rows z = 0.48 and 0.88 as printed in Table 2 of Stern et al. (2010),
    # JCAP 02 (2010) 008, doi:10.1088/1475-7516/2010/02/008: the file has their
    # H values interchanged; see the README
    stern = copy(hz)
    stern.H[stern.z .== 0.48] .= 97.0
    stern.σH[stern.z .== 0.48] .= 60.0
    stern.H[stern.z .== 0.88] .= 90.0
    stern.σH[stern.z .== 0.88] .= 40.0
    restored = fit_hubble_rate(stern.z, stern.H, stern.σH)
    println("flat ΛCDM to H(z), H₀ in km/s/Mpc")
    for (name, f) in (("12 rows", all_rows), ("11 rows, z > 0", no_prior),
        ("12 rows, 0.48 and 0.88 as printed", restored))
        @printf("  %-32s H₀ = %.1f ± %.1f   Ω_m = %.3f ± %.3f   χ²/dof = %.2f/%d\n",
            name, f.H₀, f.σH₀, f.Ωm, f.σΩm, f.χ², f.dof)
    end

    # --- a(t) -----------------------------------------------------------------
    t_H = hubble_time(PLANCK2018.H₀)
    @printf("\n1/H₀ = %.3f Gyr at H₀ = %.1f km/s/Mpc\n", t_H, PLANCK2018.H₀)
    # The integrated quantity is τ(a), so the closed form is compared in τ: near
    # a = 0 the map τ ↦ a has infinite slope and would magnify any error in τ.
    histories = map(HISTORIES) do Ωm
        τ, a = expansion_history(Ωm)
        age_numeric, age_exact = -τ[1], age(Ωm)
        deviation = maximum(abs.(τ .+ age_exact .- cosmic_time.(a, Ωm)))
        @printf("  Ω_m = %.3f: age %.3f Gyr, numerical − closed form %.1e /H₀, max |Δτ| = %.1e /H₀\n",
            Ωm, age_numeric * t_H, age_numeric - age_exact, deviation)
        deviation < TIME_TOLERANCE ||
            error("τ(a) for Ω_m = $Ωm departs from the closed form by $deviation")
        (; Ωm, t = τ .* t_H, a, age = age_numeric * t_H)
    end
    @printf("  Planck 2018 age: %.3f ± %.3f Gyr\n", PLANCK2018_AGE.t₀, PLANCK2018_AGE.σ)
    abs(histories[1].age - PLANCK2018_AGE.t₀) < 3 * PLANCK2018_AGE.σ ||
        error("age at the Planck parameters is $(histories[1].age) Gyr")

    # --- SN Ia distance moduli ------------------------------------------------
    sn = CSV.read(joinpath(DATA, "sn_ia_distance_moduli.csv"), DataFrame;
        header = ["mu", "z"], skipto = 2,)
    z, μ = sn.z, sn.mu
    sweep = [fit_flat_lcdm(z[z .> c], μ[z .> c]) for c in REDSHIFT_FLOORS]
    H₀(f) = hubble_constant(f.offset)
    σH₀(f) = H₀(f) * log(10) / 5 * f.σoffset
    @printf("\n%d SN Ia distance moduli, z = %.4f–%.3f; flat ΛCDM above a redshift floor\n",
        length(z), minimum(z), maximum(z))
    for (c, f) in zip(REDSHIFT_FLOORS, sweep)
        @printf("  z > %-5.3f  N = %4d   H₀ = %5.2f ± %.2f   Ω_m = %s   RMS = %.3f mag\n",
            c, f.dof + 2, H₀(f), σH₀(f),
            f.at_bound ? "1 (upper bound)" : @sprintf("%.3f ± %.3f", f.Ωm, f.σΩm), f.rms)
    end
    full = first(sweep)
    n_bound = count(f -> f.at_bound, sweep)
    Ωm_range = maximum(f.Ωm for f in sweep) - minimum(f.Ωm for f in sweep)
    @printf("  Ω_m range %.2f, H₀ range %.1f; at the bound in %d of %d\n", Ωm_range,
        maximum(H₀.(sweep)) - minimum(H₀.(sweep)), n_bound, length(sweep))
    @printf("  Planck 2018: H₀ = %.1f ± %.1f, Ω_m = %.3f ± %.3f\n",
        PLANCK2018.H₀, PLANCK2018.σH₀, PLANCK2018.Ωm, PLANCK2018.σΩm)

    # --- host-galaxy magnitudes -------------------------------------------------
    zg, mg = load_host_galaxy_magnitudes(joinpath(DATA, "supernova_magnitudes.csv"))
    host = fit_flat_lcdm(zg, mg)
    M_eff = host.offset - 25 - 5 * log10(C_KM_S / all_rows.H₀)
    @printf("\n%d host-galaxy magnitudes: Ω_m = %s, M = %.2f mag at the H(z) value of H₀, RMS = %.2f mag\n",
        length(zg), host.at_bound ? @sprintf("%g (bound)", host.Ωm) :
                    @sprintf("%.3f ± %.3f", host.Ωm, host.σΩm),
        M_eff, host.rms)

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1300, 1040))
    Ho(f) = rich(it("H"), subscript("0"), @sprintf(" = %.1f ± %.1f", f[1], f[2]))
    Om(v, e) = rich("Ω", subscript("m"), @sprintf(" = %.2f ± %.2f", v, e))
    note!(ax, x, y, text; colour = PALETTE.blue, align = (:left, :top)) = text!(ax, x, y;
        text, space = :relative, align, color = colour, fontsize = ANNOTATION_SIZE,)

    ax1 = Axis(fig[2, 1], xlabel = L"Redshift $z$",
        ylabel = L"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]",)
    zf = range(0, maximum(hz.z) * 1.04, length = 200)
    l_lcdm = lines!(ax1, zf, hubble_rate(zf, [all_rows.H₀, all_rows.Ωm]),
        color = PALETTE.blue,)
    unstroked_errorbars!(ax1, hz.z, hz.H, hz.σH, color = PALETTE.orange)
    l_hz = scatter!(ax1, hz.z, hz.H, color = PALETTE.orange,
        strokecolor = darker(PALETTE.orange),)
    note!(ax1, 0.04, 0.96, Ho((all_rows.H₀, all_rows.σH₀)))
    note!(ax1, 0.04, 0.87, Om(all_rows.Ωm, all_rows.σΩm))
    note!(ax1,
        0.04,
        0.78,
        rich(it("χ"), superscript("2"),
            @sprintf("/dof = %.1f/%d", all_rows.χ², all_rows.dof)),)

    ax2 = Axis(fig[2, 2], xlabel = L"Time from today $t$ [Gyr]",
        ylabel = L"Scale factor $a$",)
    styles = ((PALETTE.blue, :solid), (PALETTE.red, :dash), (PALETTE.green, :dashdot))
    l_hist = map(histories, styles) do h, (colour, style)
        shown = h.t .<= T_FUTURE
        lines!(ax2, h.t[shown], h.a[shown], color = colour, linestyle = style)
    end
    vlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    note!(ax2, 0.04, 0.96, "Age [Gyr]"; colour = PALETTE.black)
    for (k, (h, (colour, _))) in enumerate(zip(histories, styles))
        note!(ax2, 0.04, 0.96 - 0.09k,
            rich("Ω", subscript("m"), @sprintf(" = %g:  %.1f", h.Ωm, h.age)); colour,)
    end
    text!(ax2, 0.4, 0.05; text = "Today", align = (:left, :bottom),
        fontsize = ANNOTATION_SIZE,)
    xlims!(ax2, -1.04 * maximum(h.age for h in histories), T_FUTURE + 0.8)
    ylims!(ax2, 0, nothing)

    ax3 = Axis(fig[3, 1], xlabel = L"Redshift $z$",
        ylabel = L"Distance modulus $\mu$ [mag]", xscale = log10,
        xticks = logticks(-3, 0),)
    scatter!(ax3, z, μ, color = (PALETTE.orange, 0.2), markersize = MARKERSIZE.cloud,
        strokewidth = 0,)
    zfit = 10 .^ range(log10(minimum(z)), log10(maximum(z)), length = 250)
    lines!(ax3, zfit, magnitude.(zfit, full.Ωm, full.offset), color = PALETTE.blue)
    l_med = scatter!(ax3, binned_medians(z, μ)..., color = PALETTE.black,
        markersize = MARKERSIZE.dense,)
    note!(ax3, 0.04, 0.96, Ho((H₀(full), σH₀(full))))
    note!(ax3, 0.04, 0.87, Om(full.Ωm, full.σΩm))
    note!(ax3, 0.04, 0.78, @sprintf("RMS = %.2f mag", full.rms))

    ax4 = Axis(fig[3, 2], xlabel = L"Redshift floor $z_{\mathrm{min}}$",
        ylabel = rich("Fitted Ω", subscript("m")),)
    free = [!f.at_bound for f in sweep]
    Ωs = [f.Ωm for f in sweep]
    hlines!(ax4, [Ωm_BOUNDS[2]], color = PALETTE.black, linestyle = :dot,
        linewidth = GUIDE_WIDTH,)
    hlines!(ax4, [PLANCK2018.Ωm], color = PALETTE.black, linestyle = :dash,
        linewidth = GUIDE_WIDTH,)
    lines!(ax4, REDSHIFT_FLOORS, Ωs, color = PALETTE.blue)
    unstroked_errorbars!(ax4, REDSHIFT_FLOORS[free], Ωs[free],
        [f.σΩm for f in sweep[free]], color = PALETTE.blue,)
    scatter!(ax4, REDSHIFT_FLOORS[free], Ωs[free], color = PALETTE.blue,
        strokecolor = darker(PALETTE.blue),)
    scatter!(ax4, REDSHIFT_FLOORS[.!free], Ωs[.!free], color = PALETTE.red,
        marker = :xcross, markersize = MARKERSIZE.emphasis, strokecolor = darker(PALETTE.red),)
    text!(ax4, 0.125, PLANCK2018.Ωm + 0.015; text = "Planck 2018",
        align = (:right, :bottom), fontsize = ANNOTATION_SIZE,)
    text!(ax4, 0.0, Ωm_BOUNDS[2] + 0.015; text = rich("Bound Ω", subscript("m"), " = 1"),
        align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)
    text!(ax4, 0.125, Ωm_BOUNDS[2] - 0.04;
        text = @sprintf("× At the bound: %d of %d", n_bound, length(sweep)),
        align = (:right, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.red,)
    text!(ax4, 0.125, 0.62; text = @sprintf("Range %.2f", Ωm_range),
        align = (:right, :top), fontsize = ANNOTATION_SIZE, color = PALETTE.blue,)
    ylims!(ax4, 0, 1.14)

    cloud = MarkerElement(color = (PALETTE.orange, 0.6), marker = :circle,
        markersize = MARKERSIZE.dense, strokewidth = 0,)
    Legend(fig[1, 1:2], [[l_hz, cloud, l_med], [l_lcdm, l_hist[2], l_hist[3]]],
        [[L"H(z)", "SN Ia moduli", "Binned medians"],
            ["ΛCDM", "Einstein–de Sitter", rich("Ω", subscript("m"), " = 0.05")],],
        ["Data", "Flat models"]; titleposition = :left, titlesize = 22,)

    colgap!(fig.layout, 24)
    rowgap!(fig.layout, 2, 18)
    println("\nwrote ", savefigure(fig, FIGURES, "friedmann_expansion"))
end

main()
