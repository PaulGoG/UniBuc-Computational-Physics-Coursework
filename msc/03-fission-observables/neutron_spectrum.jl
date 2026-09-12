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
# distortion is visible.
#
# Two further defects, both removed: E₁(z) and γ(3/2,x) were each re-integrated
# from scratch with `quadgk` at every evaluation — about 34 000 adaptive
# quadratures per run, where `SpecialFunctions` has both in closed form; and
# `Fisiune_5.jl` re-evaluated a normalisation independent of the summation index
# inside the χ² comprehension, running a quadrature and a trapezoid over the
# whole dataset for every data point at every optimiser iteration.
#
# `Fisiune_5.jl` also divided χ² by N rather than by ν = N − 1 while labelling
# every figure χ², and bounded its optimiser at T_M = 0 where T_M^{-3/2} is
# infinite.

using Printf, SpecialFunctions, Optim, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))
include(joinpath(@__DIR__, "fission_data.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Exponential integral E₁."
E₁(z) = expint(z)
"Lower incomplete gamma γ(a, x)."
γ_lower(a, x) = gamma(a) * gamma_inc(a, x)[1]

"""
    madland_nix(E, E_f, T_m; correct)

Madland–Nix spectrum at neutron energy `E` for average fragment kinetic energy
per nucleon `E_f` and maximum residual temperature `T_m`. `correct = false`
reproduces the 2018 prefactor.
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
    area = sum((N[1:end-1] .+ N[2:end]) ./ 2 .* diff(E))
    n = N ./ area; s = σN ./ area
    χ²(T) = sum(((n .- maxwellian.(E, T)) ./ s) .^ 2)
    # bounded: the 2018 version allowed T_M = 0, where T^(-3/2) is infinite
    res = optimize(χ², 0.3, 3.0, Brent())
    T = Optim.minimizer(res)
    return T, χ²(T) / (length(E) - 1)
end

function main()
    # the prefactor, at representative fragment parameters
    E_f, T_m = 0.95, 1.05
    E = 10 .^ range(-1.3, log10(20), length = 400)
    good = madland_nix.(E, E_f, T_m)
    bad = madland_nix.(E, E_f, T_m; correct = false)
    ratio = (1 / (3 * sqrt(E_f * T_m))) / ((1 / 3) * sqrt(E_f * T_m))
    @printf("Madland–Nix prefactor: correct 1/(3√(E_f·T_m)), 2018 (1/3)√(E_f·T_m)\n")
    @printf("at E_f = %.2f MeV, T_m = %.2f MeV the two differ by a factor %.4f = 1/(E_f·T_m)\n",
            E_f, T_m, ratio)
    @printf("and since E_f and T_m vary with fragment mass, this does not cancel\n")
    @printf("under renormalisation — it reweights the mass average.\n\n")

    # trapezoid on a logarithmic grid: sum(E.*N)/sum(N) would weight by grid
    # spacing, which is not uniform here
    w = [i == 1 ? (E[2]-E[1])/2 : i == length(E) ? (E[end]-E[end-1])/2 :
         (E[i+1]-E[i-1])/2 for i in eachindex(E)]
    E_mean = sum(w .* E .* good) / sum(w .* good)
    @printf("mean energy of the correct spectrum <E> = %.3f MeV\n", E_mean)
    @printf("equivalent Maxwellian temperature (2/3)<E> = %.3f MeV\n", 2E_mean/3)
    @printf("measured ²³⁵U(n_th,f) values cluster near 1.32 MeV\n\n")

    files = [("Gook, lab", "U5SPGOOK.DAT", PALETTE.blue),
             ("Vorobyev, lab", "U5SPVORO.DAT", PALETTE.orange),
             ("Gook, CM light", "U5SPCMLF.DAT", PALETTE.green),
             ("Gook, CM heavy", "U5SPCMHF.DAT", PALETTE.purple)]
    results = map(files) do (name, file, colour)
        d = load_measurement(joinpath(DATA, "Date_experimentale", "Spectru_n", file))
        keep = d.σ .> 0
        T, χ²ν = fit_maxwellian(d.x[keep], d.y[keep], d.σ[keep])
        @printf("%-16s %3d points, T_M = %.4f MeV, χ²/ν = %.2f\n",
                name, count(keep), T, χ²ν)
        (name = name, d = d, keep = keep, T = T, colour = colour)
    end

    fig = Figure(size = (1000, 460))
    ax1 = Axis(fig[2, 1], xlabel = L"Neutron energy $E$ [MeV]",
        ylabel = L"$N(E)$ [arb.]", xscale = log10, yscale = log10,
        xticks = ([0.1, 1, 10], ["0.1", "1", "10"]))
    l_g = lines!(ax1, E, good ./ maximum(good), color = PALETTE.blue, linewidth = 1.8)
    l_b = lines!(ax1, E, bad ./ maximum(bad), color = PALETTE.red,
        linewidth = 1.6, linestyle = :dash)
    ylims!(ax1, 1e-4, 2)
    text!(ax1, 0.03, 0.06;
        text = "shapes coincide after normalising;\nthe error is in the weight each\nfragment mass carries",
        space = :relative, align = (:left, :bottom), fontsize = 14)

    ax2 = Axis(fig[2, 2], xlabel = L"Neutron energy $E$ [MeV]",
        ylabel = "Spectrum / Maxwellian fit")
    handles = []
    for r in results
        area = sum((r.d.y[r.keep][1:end-1] .+ r.d.y[r.keep][2:end]) ./ 2 .* diff(r.d.x[r.keep]))
        ratio_data = (r.d.y[r.keep] ./ area) ./ maxwellian.(r.d.x[r.keep], r.T)
        push!(handles, scatterlines!(ax2, r.d.x[r.keep], ratio_data,
              color = r.colour, markersize = 5, linewidth = 1.2))
    end
    hlines!(ax2, [1.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    ylims!(ax2, 0.6, 1.5)

    Legend(fig[1, 1:2], [[l_g, l_b]; handles],
        [["Madland–Nix, correct", "Madland–Nix, 2018 prefactor"];
         [@sprintf("%s, T=%.2f", r.name, r.T) for r in results]],
        orientation = :horizontal, framevisible = false, labelsize = 14,
        nbanks = 2, colgap = 14)
    rowsize!(fig.layout, 2, Relative(0.80))
    println("\nwrote ", savefigure(fig, FIGURES, "neutron_spectrum"))
end

main()
