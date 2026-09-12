# Expansion history of a homogeneous isotropic universe from the Friedmann
# equations, and the Hubble parameter fitted to measured H(z).
#
#   H²(z) = H₀² [Ω_r(1+z)⁴ + Ω_m(1+z)³ + Ω_k(1+z)² + Ω_Λ]
#
# Ported from Friedmann_Eq.jl (Dubna 2019, joint work with Alexandru Crăciun).
#
# The original did not run. It referenced a function `Friedmann` that was never
# defined — the definition present is `FriedmannTimeDependentNeutralFluid` — and
# used `fit.param` roughly forty lines before `fit` was assigned. It also
# declared `const G = 6.67408e-11  # Cosmological constant`, which is the
# gravitational constant, not Λ. The companion file additionally read its data
# through a hardcoded `C:\Users\GoG\Desktop\...` path and, inside an ODE
# right-hand side, assigned `u[3] = ...` — mutating the state vector rather than
# returning a derivative.
#
# `data/hubble_parameter.csv` holds 11 measurements of (z, H, σ_H) in
# km s⁻¹ Mpc⁻¹.

using Printf, CSV, DataFrames, LsqFit, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Dimensionless expansion rate E(z) = H(z)/H₀ for a flat universe."
E(z, Ωm, ΩΛ) = sqrt(Ωm * (1 + z)^3 + ΩΛ)

"""
    scale_factor(Ωm, ΩΛ, H₀; a_end, n)

Integrate da/dt = a H(a) backwards and forwards from a = 1 with RK4, in units
where time is measured in 1/H₀. Returns time (Gyr) and scale factor.
"""
function scale_factor(Ωm, ΩΛ; n = 4000, t_span = (-0.95, 1.5))
    # H₀ = 70 km/s/Mpc  ->  1/H₀ = 13.97 Gyr
    inv_H₀_Gyr = 13.968
    f(a) = a <= 0 ? 0.0 : a * sqrt(Ωm / a^3 + ΩΛ)
    ts = range(t_span[1], t_span[2], length = n)
    Δt = step(ts)
    a = Vector{Float64}(undef, n)
    i₀ = argmin(abs.(ts))
    a[i₀] = 1.0
    for i in i₀:(n - 1)                      # forward
        k₁ = f(a[i]);            k₂ = f(a[i] + Δt*k₁/2)
        k₃ = f(a[i] + Δt*k₂/2);  k₄ = f(a[i] + Δt*k₃)
        a[i+1] = a[i] + Δt/6*(k₁ + 2k₂ + 2k₃ + k₄)
    end
    for i in i₀:-1:2                          # backward
        k₁ = f(a[i]);            k₂ = f(a[i] - Δt*k₁/2)
        k₃ = f(a[i] - Δt*k₂/2);  k₄ = f(a[i] - Δt*k₃)
        a[i-1] = max(a[i] - Δt/6*(k₁ + 2k₂ + 2k₃ + k₄), 0.0)
    end
    return collect(ts) .* inv_H₀_Gyr, a
end

function main()
    data = CSV.read(joinpath(@__DIR__, "data", "hubble_parameter.csv"), DataFrame;
                    header = ["z", "H", "σH"])
    @printf("%d H(z) measurements, z from %.2f to %.2f\n",
            nrow(data), minimum(data.z), maximum(data.z))

    # flat LCDM: two free parameters, weighted by the quoted uncertainties
    model(z, p) = p[1] .* sqrt.(p[2] .* (1 .+ z).^3 .+ (1 - p[2]))
    fit = curve_fit(model, data.z, data.H, 1 ./ data.σH .^ 2, [70.0, 0.3])
    H₀, Ωm = fit.param
    σ = stderror(fit)
    χ² = sum(abs2, fit.resid)
    @printf("flat ΛCDM fit:  H₀ = %.2f ± %.2f km/s/Mpc,  Ωm = %.3f ± %.3f\n",
            H₀, σ[1], Ωm, σ[2])
    @printf("χ²/dof = %.2f  (dof = %d)\n", χ² / (nrow(data) - 2), nrow(data) - 2)
    @printf("Planck 2018 for comparison: H₀ = 67.4 ± 0.5, Ωm = 0.315 ± 0.007\n")

    fig = Figure(size = (1000, 440))

    ax1 = Axis(fig[2, 1], xlabel = L"Redshift $z$",
        ylabel = L"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
    zf = range(0, maximum(data.z) * 1.05, length = 200)
    l_fit = lines!(ax1, zf, model(zf, fit.param), color = PALETTE.blue, linewidth = 1.6)
    l_dat = errorbars!(ax1, data.z, data.H, data.σH, color = PALETTE.orange, whiskerwidth = 8)
    scatter!(ax1, data.z, data.H, color = PALETTE.orange, markersize = 9)
    text!(ax1, 0.04, 0.93;
        text = @sprintf("H₀ = %.1f ± %.1f,  Ωm = %.2f ± %.2f", H₀, σ[1], Ωm, σ[2]),
        space = :relative, align = (:left, :top), fontsize = 15, color = PALETTE.blue)

    ax2 = Axis(fig[2, 2], xlabel = L"Time from now $t$ [Gyr]",
        ylabel = L"Scale factor $a(t)$")
    handles = []
    for (Ωm_c, ΩΛ_c, name, colour) in ((0.315, 0.685, L"$\Omega_m = 0.315$, $\Omega_\Lambda = 0.685$", PALETTE.blue),
                                       (1.0, 0.0, L"$\Omega_m = 1$ (Einstein–de Sitter)", PALETTE.orange),
                                       (0.05, 0.95, L"$\Omega_\Lambda$-dominated", PALETTE.green))
        t, a = scale_factor(Ωm_c, ΩΛ_c)
        push!(handles, lines!(ax2, t, a, color = colour, linewidth = 1.6))
    end
    vlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    scatter!(ax2, [0.0], [1.0], color = PALETTE.black, markersize = 9)
    text!(ax2, 0.5, 0.06; text = "Now", space = :relative,
        align = (:center, :bottom), fontsize = 15)
    ylims!(ax2, 0, 3)

    Legend(fig[1, 1:2], [[l_dat]; handles],
        [["H(z) measurements"];
         [L"$\Omega_m = 0.315$", L"$\Omega_m = 1$", L"$\Omega_\Lambda$-dominated"]],
        orientation = :horizontal, framevisible = false, labelsize = 16, colgap = 20)

    rowsize!(fig.layout, 2, Relative(0.85))
    path = savefigure(fig, FIGURES, "friedmann_expansion")
    println("wrote ", path)
end

main()
