# Prompt-fission neutron spectrum of ²³⁵U(n_th,f): the Los Alamos
# (Madland–Nix) model, and Maxwellian fits to four measured spectra.
#
# Madland–Nix with a constant compound-nucleus cross-section and a triangular
# residual-temperature distribution:
#
#   N(E) = 1/(3√(E_f T_m)) [ u₂^{3/2}E₁(u₂) − u₁^{3/2}E₁(u₁) + γ(3/2,u₂) − γ(3/2,u₁) ]
#   u₁ = (√E − √E_f)²/T_m,   u₂ = (√E + √E_f)²/T_m
#
# Ported from Fisiune_4.jl and Fisiune_5.jl.
#
# The prefactor and the χ² convention are both settled by the course's own
# notes (Tudora, *Modele de emisie prompta globala* and *Aplicatie fit cu
# spectru Maxwellian*), which were recovered after this file was first written.
#
# **The prefactor was wrong.** The original wrote
#
#     return (1/3*sqrt(E_F*T_MAX)) * ( ... )
#
# which Julia parses as `(1/3)*sqrt(E_f·T_m)`, because `/` and `*` share
# precedence and associate left to right. The Madland–Nix prefactor is
# **1/(3√(E_f T_m))** — the code multiplied by √(E_f T_m)/3 where it had to
# divide. Since E_f and T_m both vary with fragment mass, the error does *not*
# cancel under the Maxwellian renormalisation applied afterwards; it reweights
# the mass-averaged sum. Both forms are evaluated below so the size of the
# distortion is visible. The course notes print the spectrum with the
# 1/(3√(E_f T_max)) prefactor explicitly, so this is a transcription slip in the
# original rather than a difference of convention.
#
# The level-density constant likewise has a source: ⟨a⟩ = A₀/C with C = 11 MeV
# for an optical-model compound cross-section and **C = 10 MeV for a constant
# one**, which is the case used here.
#
# Two further defects, both removed: E₁(z) and γ(3/2,x) were each re-integrated
# from scratch with `quadgk` at every evaluation — about 34 000 adaptive
# quadratures per run, where `SpecialFunctions` has both in closed form; and
# `Fisiune_5.jl` re-evaluated a normalisation independent of the summation index
# inside the χ² comprehension, running a quadrature and a trapezoid over the
# whole dataset for every data point at every optimiser iteration.
#
# **A retraction.** An earlier version of this header held that `Fisiune_5.jl`
# was wrong to divide χ² by N rather than by ν = N − 1. It was not: the course
# defines χ² = (1/n)Σ(yᵢ − f(xᵢ))²/σᵢ², divided by the number of points, and the
# original was following that definition. The reduced χ² over ν = N − 1 is kept
# here because one parameter is fitted and that is the standard reading, but the
# difference is a convention and the original was not in error. With N in the
# hundreds the two differ by well under a percent; the reference fits quoted
# below are on the course's convention.
#
# For comparison, the fits the course quotes: Hambsch & Kornilov ²³⁵U(n_th,f)
# → T_M = 1.297 MeV at χ² = 2.034; Mannhart ²⁵²Cf → 1.402 and 1.396 MeV at
# χ² = 2.513 and 3.256.
#
# `Fisiune_5.jl` did bound its optimiser at T_M = 0, where T_M^{-3/2} is
# infinite; that is fixed.

using Printf, SpecialFunctions, Optim, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_data.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Exponential integral E₁."
E₁(z) = expint(z)
"Lower incomplete gamma γ(a, x)."
γ_lower(a, x) = gamma(a) * gamma_inc(a, x)[1]

"Maximum residual temperature, T_m = √(C·TXE/A) with C = 10 MeV as in the original."
T_max(A, TXE; C = 10.0) = sqrt(C * TXE / A)

"""
    E_f_pair(A, A_H, TKE)

Average fragment kinetic energy per nucleon for the light and heavy fragment,
from momentum conservation: E_L = TKE·A_H/A, so E_f_L = E_L/A_L.
"""
function E_f_pair(A, A_H, TKE)
    A_L = A - A_H
    return ((A_H / A_L) * TKE / A, (A_L / A_H) * TKE / A)
end

"""
    madland_nix(E, E_f, T_m; correct)

Madland–Nix spectrum at neutron energy `E` for average fragment kinetic energy
per nucleon `E_f` and maximum residual temperature `T_m`. `correct = false`
reproduces the original prefactor.
"""
function madland_nix(E, E_f, T_m; correct = true)
    u₁ = (sqrt(E) - sqrt(E_f))^2 / T_m
    u₂ = (sqrt(E) + sqrt(E_f))^2 / T_m
    bracket = u₂^1.5 * E₁(u₂) - u₁^1.5 * E₁(u₁) + γ_lower(1.5, u₂) - γ_lower(1.5, u₁)
    prefactor = correct ? 1 / (3 * sqrt(E_f * T_m)) : (1 / 3) * sqrt(E_f * T_m)
    return prefactor * bracket
end

"Maxwellian spectrum of temperature `T`, normalised to unit integral."
maxwellian(E, T) = 2 / sqrt(π) * T^(-1.5) * sqrt(E) * exp(-E / T)

"""
    fit_maxwellian(E, N, σN)

Least-squares Maxwellian temperature, with the data renormalised to unit
integral by the trapezoid rule so that shape rather than scale is fitted.
Returns T_M and the reduced χ².
"""
function fit_maxwellian(E, N, σN)
    area = sum((N[1:(end - 1)] .+ N[2:end]) ./ 2 .* diff(E))
    n = N ./ area
    s = σN ./ area
    χ²(T) = sum(((n .- maxwellian.(E, T)) ./ s) .^ 2)
    # bounded: the original allowed T_M = 0, where T^(-3/2) is infinite
    res = optimize(χ², 0.3, 3.0, Brent())
    T = Optim.minimizer(res)
    return T, χ²(T) / (length(E) - 1)
end

"""
    mass_averaged_spectrum(E_grid, fragments; correct)

N(E) averaged over the mass yield, as the original did: for each fragment mass
the light- and heavy-fragment spectra are averaged and weighted by Y(A). This is
the calculation the file exists to perform, and it is where the prefactor error
matters — E_f and T_m both vary with A_H, so a prefactor carrying them cannot be
absorbed into an overall normalisation.
"""
function mass_averaged_spectrum(E_grid, fragments; correct = true)
    N = zeros(length(E_grid))
    for (i, E) in enumerate(E_grid)
        num = 0.0
        den = 0.0
        for f in fragments
            N_L = madland_nix(E, f.E_f_L, f.T_m; correct = correct)
            N_H = madland_nix(E, f.E_f_H, f.T_m; correct = correct)
            num += f.Y * 0.5 * (N_L + N_H)
            den += f.Y
        end
        N[i] = num / den
    end
    return N
end

"Trapezoid weights for a non-uniform grid."
trapz_weights(x) = [i == 1 ? (x[2]-x[1])/2 :
                    i == length(x) ? (x[end]-x[end - 1])/2 :
                    (x[i + 1]-x[i - 1])/2 for i in eachindex(x)]

function main()
    # per-mass TXE, TKE and Y(A) from the yield matrix, exactly as the original
    # file built them, so the mass average below is over the same quantities
    y = load_yields(joinpath(DATA, "Yield", "U5YAZTKE.STR"))
    masses = load_masses(joinpath(DATA, "Defecte_masa", "AUDI2021.csv"))
    A₀, Z₀ = 236, 92
    Δ₀ = Δ(masses, Z₀, A₀)
    S_n = (Δ(masses, 92, 235) + Δ(masses, 0, 1) - Δ₀) / 1000

    fragments = NamedTuple[]
    for a in sort(unique(y.A_H))
        sub = y[y.A_H .== a, :]
        Y = sum(sub.Y)
        Y > 0 || continue
        TKE = sum(sub.TKE .* sub.Y) / Y
        Zp = round(Int, Z₀ * a / A₀ - 0.5)
        δH = Δ(masses, Zp, a)
        δL = Δ(masses, Z₀ - Zp, A₀ - a)
        (δH === nothing || δL === nothing) && continue
        TXE = (Δ₀ - δH - δL) / 1000 + S_n - TKE
        TXE > 0 || continue
        efl, efh = E_f_pair(A₀, a, TKE)
        push!(fragments, (A_H = a, Y = Y, TKE = TKE, TXE = TXE,
            T_m = T_max(A₀, TXE), E_f_L = efl, E_f_H = efh,))
    end
    @printf("mass average over %d fragment masses, A_H %d–%d\n",
        length(fragments), fragments[1].A_H, fragments[end].A_H)
    @printf("  T_m spans %.3f–%.3f MeV, E_f spans %.3f–%.3f MeV\n\n",
        minimum(f.T_m for f in fragments), maximum(f.T_m for f in fragments),
        minimum(min(f.E_f_L, f.E_f_H) for f in fragments),
        maximum(max(f.E_f_L, f.E_f_H) for f in fragments))

    E = 10 .^ range(-1.3, log10(20), length = 260)
    good = mass_averaged_spectrum(E, fragments)
    bad = mass_averaged_spectrum(E, fragments; correct = false)
    w = trapz_weights(E)

    @printf("Madland–Nix prefactor: correct 1/(3√(E_f·T_m)), original (1/3)√(E_f·T_m)\n")
    E_mean_good = sum(w .* E .* good) / sum(w .* good)
    E_mean_bad = sum(w .* E .* bad) / sum(w .* bad)
    @printf("mass-averaged <E>:  correct %.4f MeV,  original prefactor %.4f MeV\n",
        E_mean_good, E_mean_bad)
    @printf("equivalent Maxwellian (2/3)<E>: %.4f vs %.4f MeV\n",
        2E_mean_good/3, 2E_mean_bad/3)
    @printf("evaluated value for ²³⁵U(n_th,f): 1.32 MeV\n")
    @printf("the two differ by %.2f %% in <E> — the prefactor carries E_f and T_m,\n",
        100 * (E_mean_bad / E_mean_good - 1))
    @printf("both of which vary with A_H, so it reweights the mass average and\n")
    @printf("cannot be absorbed into an overall normalisation.\n\n")

    files = [("Göök, lab", "U5SPGOOK.DAT", PALETTE.blue),
        ("Vorobyev, lab", "U5SPVORO.DAT", PALETTE.orange),
        ("Göök, CM light", "U5SPCMLF.DAT", PALETTE.green),
        ("Göök, CM heavy", "U5SPCMHF.DAT", PALETTE.purple),]
    results = map(files) do (name, file, colour)
        d = load_measurement(joinpath(DATA, "Date_experimentale", "Spectru_n", file))
        keep = d.σ .> 0
        T, χ²ν = fit_maxwellian(d.x[keep], d.y[keep], d.σ[keep])
        @printf("%-16s %3d points, T_M = %.4f MeV, χ²/ν = %.2f\n",
            name, count(keep), T, χ²ν)
        (name = name, d = d, keep = keep, T = T, χ²ν = χ²ν, colour = colour)
    end

    fig = Figure(size = (1000, 460))
    ax1 = Axis(fig[2, 1], xlabel = L"Neutron energy $E$ [MeV]",
        ylabel = L"$N(E)$ [arb.]", xscale = log10, yscale = log10,
        xticks = logticks(-1, 1), yticks = logticks(-4, 0),)
    l_g = lines!(ax1, E, good ./ maximum(good), color = PALETTE.blue, linewidth = 1.8)
    l_b = lines!(ax1, E, bad ./ maximum(bad), color = PALETTE.red,
        linewidth = 1.6, linestyle = :dash,)
    ylims!(ax1, 1e-4, 2)
    text!(ax1, 0.03, 0.06;
        text = "shapes coincide after normalising;\nthe error is in the weight each\nfragment mass carries",
        space = :relative, align = (:left, :bottom), fontsize = 14,)

    ax2 = Axis(fig[2, 2], xlabel = L"Neutron energy $E$ [MeV]",
        ylabel = "Spectrum / Maxwellian fit",)
    handles = []
    for r in results
        area = sum((r.d.y[r.keep][1:(end - 1)] .+ r.d.y[r.keep][2:end]) ./ 2 .*
                   diff(r.d.x[r.keep]))
        ratio_data = (r.d.y[r.keep] ./ area) ./ maxwellian.(r.d.x[r.keep], r.T)
        push!(handles,
            scatterlines!(ax2, r.d.x[r.keep], ratio_data,
                color = r.colour, markersize = MARKERSIZE.cloud + 1, linewidth = 1.2,),)
    end
    hlines!(ax2, [1.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    text!(ax2, 13.7, 1.02; text = "Maxwellian fit", space = :data,
        align = (:right, :bottom), fontsize = 13, color = PALETTE.black,)
    ylims!(ax2, 0.6, 1.62)
    # The two centre-of-mass sets are counting-limited above about 6 MeV and
    # their excursions ran off the top of the panel; the ratio is only
    # informative where the measurement has counts.
    xlims!(ax2, 0, 14)

    Legend(fig[1, 1:2], [[l_g, l_b]; handles],
        [["Madland–Nix, correct", "Madland–Nix, original prefactor"];
         [rich(r.name, ", ", it("T"), subscript("M"),
              @sprintf(" = %.2f MeV, ", r.T), it("χ"), superscript("2"), "/",
              it("ν"), @sprintf(" = %.2f", r.χ²ν)) for r in results]],
        orientation = :horizontal, framevisible = false, labelsize = 14,
        nbanks = 2, colgap = 14,)
    rowsize!(fig.layout, 2, Relative(0.80))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_spectrum"))
end

main()
