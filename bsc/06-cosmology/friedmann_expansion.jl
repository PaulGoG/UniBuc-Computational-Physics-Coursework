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
# km s⁻¹ Mpc⁻¹. `data/supernova_magnitudes.csv` holds 6571 usable host-galaxy
# magnitudes with redshifts, which is what `Friedmann_Magnitude.jl` fitted.
#
# Note on what changed and what did not. The original fitted a *straight line*
# to H(z) (`LineFit(t, p) = p[1]*t .+ p[2]`); flat ΛCDM is fitted here instead,
# which is the model the data are usually read against. The magnitude–redshift
# analysis is kept — it was the entire subject of the companion file — but the
# original's route to it, a `MonteCarloProblem` with randomised initial
# conditions driven through `build_loss_objective`, is replaced by the closed
# form for the luminosity distance. The equation-of-state survey sketched in the
# original's docstring (quintessence, Chaplygin and modified Chaplygin gas) was
# never implemented there and is not implemented here.

using Printf, CSV, DataFrames, LsqFit, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")

"Dimensionless expansion rate E(z) = H(z)/H₀ for a flat universe."
E(z, Ωm, ΩΛ) = sqrt(Ωm * (1 + z)^3 + ΩΛ)

const C_KM_S = 2.99792458e5

"""
    luminosity_distance(z, H₀, Ωm)

d_L = (1+z)·(c/H₀)·∫₀^z dz'/E(z'), in Mpc, for a flat universe. The comoving
integral is done by Simpson's rule on a fixed grid, which is ample here.
"""
function luminosity_distance(z, H₀, Ωm)
    n = 200
    zs = range(0, z, length = n + 1)
    h = step(zs)
    f = [1 / E(zi, Ωm, 1 - Ωm) for zi in zs]
    I = h / 3 * (f[1] + f[end] + 4sum(f[2:2:end-1]) + 2sum(f[3:2:end-2]))
    return (1 + z) * C_KM_S / H₀ * I
end

"Distance modulus plus an effective absolute magnitude."
apparent_magnitude(z, H₀, Ωm, M) = 5 * log10(luminosity_distance(z, H₀, Ωm)) + 25 + M

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

    # --- magnitude-redshift, the subject of Friedmann_Magnitude.jl ------------
    sn = CSV.read(joinpath(@__DIR__, "data", "supernova_magnitudes.csv"), DataFrame;
                  header = ["name", "mag", "z"], skipto = 3, silencewarnings = true)
    zs = Float64[]; ms = Float64[]
    for r in eachrow(sn)
        (ismissing(r.mag) || ismissing(r.z)) && continue
        zv = tryparse(Float64, strip(string(r.z)))
        mv = tryparse(Float64, strip(string(r.mag)))
        (zv === nothing || mv === nothing) && continue
        (0.001 < zv < 2.0 && mv > 0) || continue
        push!(zs, zv); push!(ms, mv)
    end
    @printf("\nmagnitude–redshift: %d usable host galaxies, z from %.4f to %.3f\n",
            length(zs), minimum(zs), maximum(zs))

    mag_model(z, p) = [apparent_magnitude(zi, H₀, clamp(p[2], 0.05, 1.0), p[1]) for zi in z]
    mfit = curve_fit(mag_model, zs, ms, [-20.0, 0.3])
    M_eff, Ωm_sn = mfit.param[1], clamp(mfit.param[2], 0.05, 1.0)
    scatter_rms = sqrt(mean(abs2, mfit.resid))
    railed = Ωm_sn >= 0.999 || Ωm_sn <= 0.051
    @printf("effective absolute magnitude M = %.2f, residual RMS = %.2f mag\n",
            M_eff, scatter_rms)
    @printf("Ωm = %.3f%s\n", Ωm_sn, railed ? "  — AT THE BOUND, i.e. unconstrained" : "")
    @printf("these are host-galaxy magnitudes, not standardised SN Ia peak\n")
    @printf("magnitudes, so with %.2f mag of scatter the data fix the distance\n", scatter_rms)
    @printf("scale but carry essentially no information about Ωm. The trend is\n")
    @printf("real; the cosmological constraint from it is not.\n")
    @printf("(M lands near -19.6, close to the SN Ia absolute magnitude of about\n")
    @printf("-19.3, which is a coincidence of this sample rather than a result.)\n")

    # binned medians, to show the trend through the scatter
    edges = 10 .^ range(log10(0.002), log10(maximum(zs)), length = 18)
    bz = Float64[]; bm = Float64[]
    for i in 1:(length(edges) - 1)
        sel = (zs .>= edges[i]) .& (zs .< edges[i+1])
        count(sel) < 20 && continue
        push!(bz, median(zs[sel])); push!(bm, median(ms[sel]))
    end

    fig = Figure(size = (1200, 440))

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

    ax3 = Axis(fig[2, 3], xlabel = L"Redshift $z$", ylabel = "Host galaxy magnitude",
        xscale = log10, yreversed = true,
        xticks = ([0.002, 0.01, 0.1, 1.0], ["0.002", "0.01", "0.1", "1"]))
    scatter!(ax3, zs, ms, color = (PALETTE.sky, 0.12), markersize = 2)
    zf = 10 .^ range(log10(0.002), log10(maximum(zs)), length = 200)
    l_mag = lines!(ax3, zf, [apparent_magnitude(z, H₀, Ωm_sn, M_eff) for z in zf],
        color = PALETTE.red, linewidth = 2)
    l_bin = scatter!(ax3, bz, bm, color = PALETTE.black, markersize = 9)
    text!(ax3, 0.03, 0.06;
        text = @sprintf("M = %.1f, RMS %.2f mag\nΩm unconstrained", M_eff, scatter_rms),
        space = :relative, align = (:left, :bottom), fontsize = 14)

    Legend(fig[1, 1:3], [[l_dat]; handles; [l_bin, l_mag]],
        [["H(z) measurements"];
         [L"$\Omega_m = 0.315$", L"$\Omega_m = 1$", L"$\Omega_\Lambda$-dominated"];
         ["Binned medians", "Flat ΛCDM"]],
        orientation = :horizontal, framevisible = false, labelsize = 15,
        nbanks = 2, colgap = 16)

    rowsize!(fig.layout, 2, Relative(0.82))
    path = savefigure(fig, FIGURES, "friedmann_expansion")
    println("wrote ", path)
end

main()
