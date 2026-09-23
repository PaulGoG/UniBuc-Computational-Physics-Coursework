# Prompt-fission neutron spectrum of ²³⁵U(n_th,f): the Los Alamos (Madland–Nix)
# model averaged over the mass yield, and Maxwellian fits to four measured
# spectra.
#
# Madland–Nix with a constant compound-nucleus cross-section and a triangular
# residual-temperature distribution (Madland and Nix, Nucl. Sci. Eng. 81, 213
# (1982), doi:10.13182/NSE82-5):
#
#   N(E) = 1/(3√(E_f T_m)) [u₂^{3/2}E₁(u₂) − u₁^{3/2}E₁(u₁) + γ(3/2,u₂) − γ(3/2,u₁)]
#   u₁ = (√E − √E_f)²/T_m,   u₂ = (√E + √E_f)²/T_m
#
# which integrates to one over E ∈ [0, ∞) for every E_f and T_m; `main` asserts
# this for every fragment. Fisiune_4.jl on the `legacy` branch
# (Julia-Workflow-FFUB/Fisiune_M_2/) wrote the prefactor as
# `(1/3*sqrt(E_F*T_MAX))`, which Julia parses as (1/3)√(E_f T_m): it multiplied
# where the formula divides, and its spectrum integrates to E_f T_m instead of
# one. E_f and T_m both vary with fragment mass, so the error reweights the
# mass average and does not cancel under the renormalisation applied
# afterwards. Both forms are evaluated below.
#
# Per fragment mass, as in Fisiune_4.jl: ⟨TKE⟩(A) from the yield cells, TXE =
# Q(A, Z_p) + Sₙ(²³⁶U) − TKE(A) at the single most probable charge
# Z_p = round(Z_UCD − 0.5), T_m = √(C·TXE/A₀) with C = 10 MeV, E_f per nucleon
# from momentum conservation, and the light- and heavy-fragment spectra
# averaged and weighted by Y(A).
#
# The Maxwellian fits compare data and model both normalised to unit area over
# the measured window, the data by the trapezoid rule; χ² is reduced by
# ν = N − 1 for the one fitted parameter, where Fisiune_5.jl divided by N. That
# file also bounded its optimiser at T_M = 0, where T_M^{-3/2} is infinite, and
# both files re-integrated E₁, γ(3/2, x) and the window integral by quadrature
# at every evaluation, where SpecialFunctions has them in closed form.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, SpecialFunctions, Optim, Statistics, QuadGK
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_core.jl"))

"Level-density constant of T_m = √(C·TXE/A₀) [MeV], as set in Fisiune_4.jl."
const C_LEVEL_DENSITY = 10.0
"Neutron-energy grid of the mass-averaged spectrum [MeV]."
const E_GRID = 10 .^ range(-1.3, log10(20), length = 260)
"Bounds of the Maxwellian temperature search [MeV]."
const T_BOUNDS = (0.3, 3.0)
"Relative tolerance on the unit normalisation of the Madland–Nix spectrum."
const NORMALISATION_TOLERANCE = 1e-6

"Exponential integral E₁."
E₁(z) = expint(z)
"Lower incomplete gamma γ(a, x)."
γ_lower(a, x) = gamma(a) * gamma_inc(a, x)[1]

"""
    madland_nix(E, E_f, T_m; correct = true)

Madland–Nix spectrum at neutron energy `E` for the fragment kinetic energy per
nucleon `E_f` and maximum residual temperature `T_m`, unit-normalised over
[0, ∞). `correct = false` uses the prefactor of Fisiune_4.jl, (1/3)√(E_f T_m).
"""
function madland_nix(E, E_f, T_m; correct = true)
    u₁ = (sqrt(E) - sqrt(E_f))^2 / T_m
    u₂ = (sqrt(E) + sqrt(E_f))^2 / T_m
    bracket = u₂^1.5 * E₁(u₂) - u₁^1.5 * E₁(u₁) + γ_lower(1.5, u₂) - γ_lower(1.5, u₁)
    prefactor = correct ? 1 / (3 * sqrt(E_f * T_m)) : sqrt(E_f * T_m) / 3
    return prefactor * bracket
end

"Maxwellian spectrum of temperature `T`, unit-normalised over [0, ∞)."
maxwellian(E, T) = 2 / sqrt(π) * T^(-1.5) * sqrt(E) * exp(-E / T)

"Fraction of the Maxwellian of temperature `T` below `E`: erf(√(E/T)) − (2/√π)√(E/T) e^{−E/T}."
function maxwellian_cdf(E, T)
    x = sqrt(E / T)
    return erf(x) - 2 / sqrt(π) * x * exp(-x^2)
end

"The Maxwellian normalised to unit integral over the window `[E_min, E_max]`."
windowed_maxwellian(E, T, E_min, E_max) = maxwellian(E, T) / (maxwellian_cdf(E_max, T) -
                                           maxwellian_cdf(E_min, T))

trapezoid_area(E, N) = sum((N[1:(end - 1)] .+ N[2:end]) ./ 2 .* diff(E))

"""
    fit_maxwellian(E, N, σN) -> (T, χ²/ν)

Least-squares Maxwellian temperature, data and model both unit-normalised over
the measured window, so that shape rather than scale is fitted.
"""
function fit_maxwellian(E, N, σN)
    area = trapezoid_area(E, N)
    n, s = N ./ area, σN ./ area
    E_min, E_max = extrema(E)
    χ²(T) = sum(((n .- windowed_maxwellian.(E, T, E_min, E_max)) ./ s) .^ 2)
    T = Optim.minimizer(optimize(χ², T_BOUNDS..., Brent()))
    return T, χ²(T) / (length(E) - 1)
end

"""
    mass_averaged_spectrum(E_grid, fragments; correct = true)

N(E) averaged over the mass yield: for each fragment mass the light- and
heavy-fragment spectra are averaged and weighted by Y(A).
"""
function mass_averaged_spectrum(E_grid, fragments; correct = true)
    return map(E_grid) do E
        num = sum(f.Y * (madland_nix(E, f.E_f_L, f.T_m; correct) +
                   madland_nix(E, f.E_f_H, f.T_m; correct)) / 2 for f in fragments)
        num / sum(f.Y for f in fragments)
    end
end

"Trapezoid weights of a non-uniform grid."
trapezoid_weights(x) = [i == 1 ? (x[2] - x[1]) / 2 :
                        i == length(x) ? (x[end] - x[end - 1]) / 2 :
                        (x[i + 1] - x[i - 1]) / 2 for i in eachindex(x)]

function main()
    cells = yield_cells(load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR")))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    bins = tke_by_mass(cells)
    S_n = separation_energy(masses, Z₀, A₀)[1]

    fragments = NamedTuple[]
    for a in sort!(collect(keys(bins)))
        Z_p = round(Int, most_probable_charge(a))
        q = q_value(masses, Z_p, a)
        q === nothing && continue
        TXE = q[1] + S_n - bins[a].TKE
        TXE > 0 || continue
        T_m = sqrt(C_LEVEL_DENSITY * TXE / A₀)
        A_L = A₀ - a
        push!(fragments,
            (A_H = a, Y = bins[a].Y, TXE = TXE, T_m = T_m,
                E_f_L = a / A_L * bins[a].TKE / A₀, E_f_H = A_L / a * bins[a].TKE / A₀,),)
    end
    @printf("mass average over %d fragment masses, A_H %d–%d; T_m %.3f–%.3f MeV, E_f %.3f–%.3f MeV\n",
        length(fragments), fragments[1].A_H, fragments[end].A_H,
        extrema(f.T_m for f in fragments)...,
        minimum(min(f.E_f_L, f.E_f_H) for f in fragments),
        maximum(max(f.E_f_L, f.E_f_H) for f in fragments))

    # the closed form integrates to one for every fragment; the original
    # prefactor gives E_f T_m instead
    worst = 0.0
    worst_original = (0.0, Inf)
    for f in fragments, E_f in (f.E_f_L, f.E_f_H)

        norm = quadgk(E -> madland_nix(E, E_f, f.T_m), 0, Inf; rtol = 1e-9)[1]
        worst = max(worst, abs(norm - 1))
        bad = quadgk(E -> madland_nix(E, E_f, f.T_m; correct = false), 0, Inf; rtol = 1e-9)[1]
        worst_original = (max(worst_original[1], bad), min(worst_original[2], bad))
    end
    @printf("∫N(E)dE over [0, ∞): 1 to within %.1e for every fragment with the correct prefactor; ",
        worst)
    @printf("%.2f to %.2f with the original one\n", worst_original[2], worst_original[1])
    worst < NORMALISATION_TOLERANCE ||
        error("the Madland–Nix spectrum is not unit-normalised: $worst")

    good = mass_averaged_spectrum(E_GRID, fragments)
    bad = mass_averaged_spectrum(E_GRID, fragments; correct = false)
    w = trapezoid_weights(E_GRID)
    E_mean = sum(w .* E_GRID .* good) / sum(w .* good)
    E_mean_bad = sum(w .* E_GRID .* bad) / sum(w .* bad)
    @printf("mass-averaged ⟨E⟩ = %.4f MeV (equivalent Maxwellian temperature 2⟨E⟩/3 = %.4f MeV); ",
        E_mean, 2E_mean / 3)
    @printf("with the original prefactor %.4f MeV, %.2f %% higher\n\n", E_mean_bad,
        100 * (E_mean_bad / E_mean - 1))

    # the window normalisation in closed form against quadrature
    T_test, window = 1.3, (0.55, 12.5)
    quad = quadgk(E -> maxwellian(E, T_test), window...)[1]
    assert_close("Maxwellian window integral",
        maxwellian_cdf(window[2], T_test) -
        maxwellian_cdf(window[1], T_test), quad;
        rtol = 1e-10,)

    files = [("Göök, lab", "U5SPGOOK.DAT", PALETTE.blue, :circle),
        ("Vorobyev, lab", "U5SPVORO.DAT", PALETTE.orange, :rect),
        ("Göök, CM light", "U5SPCMLF.DAT", PALETTE.green, :utriangle),
        ("Göök, CM heavy", "U5SPCMHF.DAT", PALETTE.purple, :diamond),]
    results = map(files) do (name, file, colour, marker)
        d = load_measurement(joinpath(DATA, "Date_experimentale", "Spectru_n", file))
        keep = d.σ .> 0
        T, χ²ν = fit_maxwellian(d.x[keep], d.y[keep], d.σ[keep])
        held = maxwellian_cdf(maximum(d.x[keep]), T) - maxwellian_cdf(minimum(d.x[keep]), T)
        @printf("%-15s %3d points, %.2f–%.2f MeV holding %.1f %% of the Maxwellian: T_M = %.4f MeV, χ²/ν = %.2f\n",
            name, count(keep), extrema(d.x[keep])..., 100held, T, χ²ν)
        (; name, E = d.x[keep], N = d.y[keep], σ = d.σ[keep], T, χ²ν, colour, marker)
    end

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1600, 700))
    ax1 = Axis(fig[2, 1], xlabel = L"Neutron energy $E$ [MeV]",
        ylabel = L"$N(E)$ [MeV$^{-1}$]", xscale = log10, yscale = log10,
        xticks = logticks(-1, 1), yticks = logticks(-4, 0),)
    l_g = lines!(ax1, E_GRID, good, color = PALETTE.blue)
    l_b = lines!(
        ax1, E_GRID, bad ./ (sum(w .* bad)), color = PALETTE.red, linestyle = :dash,)
    ylims!(ax1, 1e-4, 1.2)
    text!(ax1, 0.03, 0.05;
        text = rich(
            rich("⟨", it("E"), @sprintf("⟩ = %.3f MeV", E_mean), color = PALETTE.blue),
            "\n",
            rich(@sprintf("%.3f MeV with the original prefactor", E_mean_bad), color = PALETTE.red),),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)

    ax2 = Axis(fig[2, 2], xlabel = L"Neutron energy $E$ [MeV]", ylabel = "Spectrum / Maxwellian fit")
    handles = Any[]
    for r in results
        area = trapezoid_area(r.E, r.N)
        model = windowed_maxwellian.(r.E, r.T, extrema(r.E)...)
        ratio = r.N ./ area ./ model
        errorbars_unstroked!(
            ax2, r.E, ratio, r.σ ./ area ./ model, color = (r.colour, 0.5),
            linewidth = GUIDE_WIDTH, whiskerwidth = 0,)
        push!(handles,
            scatterlines!(ax2, r.E, ratio, color = r.colour, marker = r.marker,
                markersize = MARKERSIZE.dense, linewidth = GUIDE_WIDTH,),)
    end
    hlines!(ax2, [1.0], color = PALETTE.black, linestyle = :dash, linewidth = GUIDE_WIDTH)
    text!(ax2, 13.7, 1.02; text = "Maxwellian fit", align = (:right, :bottom),
        fontsize = ANNOTATION_SIZE,)
    # the centre-of-mass sets are counting-limited above about 6 MeV
    xlims!(ax2, 0, 14)
    ylims!(ax2, 0.6, 1.62)

    # fitted temperatures as a block in the empty lower left of the model panel,
    # above the mean energies; the ratio panel has data wherever a block of this
    # size would go
    fit_line(r) = rich(
        r.name, ": ", it("T"), subscript("M"), @sprintf(" = %.3f MeV, ", r.T),
        it("χ"), superscript("2"), "/", it("ν"), @sprintf(" = %.2f", r.χ²ν), color = r.colour,)
    text!(ax1, 0.03, 0.33;
        text = rich(
            "Maxwellian fits\n", fit_line(results[1]), "\n", fit_line(results[2]), "\n",
            fit_line(results[3]), "\n", fit_line(results[4]),),
        space = :relative, align = (:left, :bottom), fontsize = ANNOTATION_SIZE,)

    Legend(fig[1, 1:2], [[l_g, l_b], handles],
        [["Correct prefactor", "Original prefactor"], [r.name for r in results]],
        ["Madland–Nix", "Maxwellian fits"]; titleposition = :left, nbanks = 2, groupgap = 30,)
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_spectrum"))
end

main()
