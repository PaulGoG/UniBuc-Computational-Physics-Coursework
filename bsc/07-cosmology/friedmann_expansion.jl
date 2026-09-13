# Expansion history of a homogeneous isotropic universe from the Friedmann
# equations, and flat ΛCDM fitted to two independent datasets.
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
# Data.
#
#   `data/hubble_parameter.csv` holds measurements of (z, H, σ_H) in
#   km s⁻¹ Mpc⁻¹.
#
#   `data/sn_ia_distance_moduli.csv` holds 2376 Type Ia supernova distance
#   moduli against redshift, out to z = 1.91. It is Crăciun's reduction of the
#   raw catalogue the project worked from, selected on `Method == "SNIa"`,
#   deduplicated in redshift and sorted. The raw file was not preserved, and
#   the compilation is not named in anything that survives; the column set —
#   distance modulus keyed by measurement method — is characteristic of a
#   redshift-independent distance compilation of the NED-D kind, but that is an
#   inference and not a citation. Treat the provenance as unverified.
#
#   `data/supernova_magnitudes.csv` holds 6604 host-galaxy apparent magnitudes
#   with redshifts. This is what `Friedmann_Magnitude.jl` actually fitted, and
#   it is kept because that was the whole subject of the companion file, but it
#   is the weaker of the two: host-galaxy magnitudes are not standardised SN Ia
#   peak magnitudes, and the fit to them is reported below only to show what it
#   can and cannot do.
#
# Note on what changed and what did not. The original fitted a *straight line*
# to H(z) (`LineFit(t, p) = p[1]*t .+ p[2]`); flat ΛCDM is fitted here instead,
# which is the model the data are usually read against. The original's route to
# the magnitude analysis, a `MonteCarloProblem` with randomised initial
# conditions driven through `build_loss_objective`, is replaced by the closed
# form for the luminosity distance. The equation-of-state survey sketched in the
# original's docstring (quintessence, Chaplygin and modified Chaplygin gas) was
# never implemented there and is not implemented here.

using Printf, CSV, DataFrames, LsqFit, Statistics
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const C_KM_S = 2.99792458e5

"Dimensionless expansion rate E(z) = H(z)/H₀ for a flat universe."
E(z, Ωm) = sqrt(Ωm * (1 + z)^3 + (1 - Ωm))

"""
    comoving_interpolator(Ωm; zmax, n = 4000)

Tabulate ∫₀^z dz'/E(z') on a fixed grid and return a closure that linearly
interpolates it.

The integral is the whole cost of a distance modulus and depends only on Ωm,
not on H₀, so tabulating it once per trial Ωm turns a fit over N supernovae
into O(N + n) work instead of integrating separately for every point at every
iteration.
"""
function comoving_interpolator(Ωm; zmax, n = 4000)
    zs = range(0, zmax, length = n)
    h = step(zs)
    f = [1 / E(z, Ωm) for z in zs]
    I = similar(f)
    I[1] = 0.0
    for i in 2:n
        I[i] = I[i-1] + h / 2 * (f[i-1] + f[i])
    end
    return function (z)
        z <= 0 && return 0.0
        z >= zmax && return I[end]
        k = min(n - 1, 1 + floor(Int, z / h))
        t = (z - zs[k]) / h
        return I[k] + t * (I[k+1] - I[k])
    end
end

"Luminosity distance in Mpc for a flat universe, given a tabulated comoving integral."
luminosity_distance(z, H₀, Ic) = (1 + z) * C_KM_S / H₀ * Ic(z)

"Distance modulus μ = 5 log₁₀(d_L/Mpc) + 25."
function distance_moduli(zs, H₀, Ωm; zmax)
    Ic = comoving_interpolator(Ωm; zmax)
    return [5 * log10(max(luminosity_distance(z, H₀, Ic), 1e-12)) + 25 for z in zs]
end

const Ωm_BOUNDS = (0.02, 1.0)

"""
    fit_lcdm(z, μ)

Fit flat ΛCDM to distance moduli, returning H₀, Ωm, their standard errors, the
residual RMS, and whether Ωm ended at a bound.

Ωm is bounded: unbounded, a sample that carries no shape information walks off
to unphysical values and reports a confident-looking number for one. Ending at
a bound is the honest signal that the data do not constrain it, so it is
returned rather than hidden.
"""
function fit_lcdm(z, μ)
    zmax = maximum(z) * 1.02
    model(zz, p) = distance_moduli(zz, p[1], clamp(p[2], Ωm_BOUNDS...); zmax)
    f = curve_fit(model, z, μ, [70.0, 0.3];
                  lower = [40.0, Ωm_BOUNDS[1]], upper = [120.0, Ωm_BOUNDS[2]])
    H₀, Ωm = f.param
    railed = Ωm <= Ωm_BOUNDS[1] + 1e-6 || Ωm >= Ωm_BOUNDS[2] - 1e-6
    σ = try
        stderror(f)
    catch
        [NaN, NaN]
    end
    return (; H₀, Ωm, σH₀ = σ[1], σΩm = σ[2],
            rms = sqrt(mean(abs2, f.resid)), railed, n = length(z))
end

"""
    scale_factor(Ωm, ΩΛ; n, t_span)

Integrate da/dt = a H(a) backwards and forwards from a = 1 with RK4, in units
where time is measured in 1/H₀. Returns time (Gyr) and scale factor.
"""
function scale_factor(Ωm, ΩΛ; n = 4000, t_span = (-0.95, 1.5))
    inv_H₀_Gyr = 13.968                       # 1/H₀ for H₀ = 70 km/s/Mpc
    f(a) = a <= 0 ? 0.0 : a * sqrt(Ωm / a^3 + ΩΛ)
    ts = range(t_span[1], t_span[2], length = n)
    Δt = step(ts)
    a = Vector{Float64}(undef, n)
    i₀ = argmin(abs.(ts))
    a[i₀] = 1.0
    for i in i₀:(n - 1)                       # forward
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

"Median μ in logarithmically spaced redshift bins holding at least `minimum_count` points."
function binned_medians(z, μ; nbins = 20, minimum_count = 15)
    edges = 10 .^ range(log10(minimum(z)), log10(maximum(z)), length = nbins + 1)
    bz = Float64[]; bμ = Float64[]
    for i in 1:nbins
        sel = (z .>= edges[i]) .& (z .< edges[i+1])
        count(sel) < minimum_count && continue
        push!(bz, median(z[sel])); push!(bμ, median(μ[sel]))
    end
    return bz, bμ
end

function load_host_galaxy_magnitudes(path)
    raw = CSV.read(path, DataFrame; header = ["name", "mag", "z"],
                   skipto = 3, silencewarnings = true)
    zs = Float64[]; ms = Float64[]
    for r in eachrow(raw)
        (ismissing(r.mag) || ismissing(r.z)) && continue
        zv = tryparse(Float64, strip(string(r.z)))
        mv = tryparse(Float64, strip(string(r.mag)))
        (zv === nothing || mv === nothing) && continue
        (0.001 < zv < 2.0 && mv > 0) || continue
        push!(zs, zv); push!(ms, mv)
    end
    return zs, ms
end

function main()
    # --- H(z) -----------------------------------------------------------------
    hz = CSV.read(joinpath(@__DIR__, "data", "hubble_parameter.csv"), DataFrame;
                  header = ["z", "H", "σH"])
    @printf("%d H(z) measurements, z from %.2f to %.2f\n",
            nrow(hz), minimum(hz.z), maximum(hz.z))

    hz_model(z, p) = p[1] .* sqrt.(p[2] .* (1 .+ z).^3 .+ (1 - p[2]))
    hfit = curve_fit(hz_model, hz.z, hz.H, 1 ./ hz.σH .^ 2, [70.0, 0.3])
    H₀_hz, Ωm_hz = hfit.param
    σ_hz = stderror(hfit)
    χ² = sum(abs2, hfit.resid)
    @printf("flat ΛCDM to H(z):  H₀ = %.2f ± %.2f km/s/Mpc,  Ωm = %.3f ± %.3f\n",
            H₀_hz, σ_hz[1], Ωm_hz, σ_hz[2])
    @printf("χ²/dof = %.2f  (dof = %d)\n", χ² / (nrow(hz) - 2), nrow(hz) - 2)

    # --- SN Ia Hubble diagram -------------------------------------------------
    sn = CSV.read(joinpath(@__DIR__, "data", "sn_ia_distance_moduli.csv"), DataFrame;
                  header = ["mu", "z"], skipto = 2)
    z, μ = sn.z, sn.mu
    @printf("\n%d SN Ia distance moduli, z from %.4f to %.3f\n",
            length(z), minimum(z), maximum(z))

    full = fit_lcdm(z, μ)
    @printf("flat ΛCDM to μ(z):  H₀ = %.2f ± %.2f km/s/Mpc,  Ωm = %.3f ± %.3f\n",
            full.H₀, full.σH₀, full.Ωm, full.σΩm)
    @printf("residual RMS = %.3f mag\n", full.rms)

    # Stability of Ωm under a redshift floor. The whole point: a heterogeneous
    # compilation can return a tight-looking Ωm that moves when the sample does.
    cuts = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12]
    sweep = [(c, fit_lcdm(z[z .> c], μ[z .> c])) for c in cuts]
    println("\nΩm against a redshift floor:")
    for (c, f) in sweep
        @printf("  z > %-5.3f  N = %4d   H₀ = %6.2f   Ωm = %.3f%s\n",
                c, f.n, f.H₀, f.Ωm, f.railed ? "   AT THE BOUND" : "")
    end
    n_railed = count(f -> f[2].railed, sweep)
    spread = maximum(f[2].Ωm for f in sweep) - minimum(f[2].Ωm for f in sweep)
    @printf("Ωm moves by %.2f across these cuts and rails in %d of %d.\n",
            spread, n_railed, length(sweep))
    @printf("Planck 2018 for comparison: H₀ = 67.4 ± 0.5, Ωm = 0.315 ± 0.007\n")
    @printf("H₀ is stable and well determined; Ωm is not. These are compiled\n")
    @printf("distance moduli from mixed sources, not a standardised sample, so\n")
    @printf("the distance scale survives the heterogeneity and the shape does not.\n")

    # --- host-galaxy magnitudes, for fidelity to the companion file -----------
    zg, mg = load_host_galaxy_magnitudes(
        joinpath(@__DIR__, "data", "supernova_magnitudes.csv"))
    zmax_g = maximum(zg) * 1.02
    gmodel(zz, p) = distance_moduli(zz, H₀_hz, clamp(p[2], Ωm_BOUNDS...); zmax = zmax_g) .+ p[1]
    gfit = curve_fit(gmodel, zg, mg, [-20.0, 0.3];
                     lower = [-40.0, Ωm_BOUNDS[1]], upper = [0.0, Ωm_BOUNDS[2]])
    g_rms = sqrt(mean(abs2, gfit.resid))
    g_railed = gfit.param[2] <= Ωm_BOUNDS[1] + 1e-6 || gfit.param[2] >= Ωm_BOUNDS[2] - 1e-6
    @printf("\n%d host-galaxy magnitudes (what Friedmann_Magnitude.jl fitted):\n", length(zg))
    @printf("offset M = %.2f, residual RMS = %.2f mag, Ωm = %.3f%s\n",
            gfit.param[1], g_rms, gfit.param[2], g_railed ? "  — at the bound" : "")
    @printf("%.2f mag of scatter against %.3f for the distance moduli: these fix\n",
            g_rms, full.rms)
    @printf("the distance scale and carry no cosmological information.\n")

    # --- figure ---------------------------------------------------------------
    # Four panels of disjoint content, so each carries its own legend rather
    # than one shared key: a single legend would put the orange of the H(z)
    # points beside the orange of the Einstein–de Sitter curve, and two
    # different "flat ΛCDM fit" entries, as if they meant the same thing.
    fig = Figure(size = (1150, 740))

    Ho(v, e) = rich(it("H"), subscript("0"), @sprintf(" = %.1f ± %.1f", v, e))
    Om(v, e) = rich("Ω", subscript("m"), @sprintf(" = %.2f ± %.2f", v, e))

    ax1 = Axis(fig[1, 1], xlabel = L"Redshift $z$",
        ylabel = L"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
    zf = range(0, maximum(hz.z) * 1.05, length = 200)
    lines!(ax1, zf, hz_model(zf, hfit.param), color = PALETTE.blue, linewidth = 1.8,
        label = "Flat ΛCDM fit")
    errorbars!(ax1, hz.z, hz.H, hz.σH, color = PALETTE.orange, whiskerwidth = 8)
    scatter!(ax1, hz.z, hz.H, color = PALETTE.orange, markersize = MARKERSIZE.dense,
        label = L"$H(z)$ measurements")
    text!(ax1, 0.04, 0.96; text = Ho(H₀_hz, σ_hz[1]), space = :relative,
        align = (:left, :top), fontsize = 15, color = PALETTE.blue)
    text!(ax1, 0.04, 0.87; text = Om(Ωm_hz, σ_hz[2]), space = :relative,
        align = (:left, :top), fontsize = 15, color = PALETTE.blue)
    axislegend(ax1; position = :rb, framevisible = false, labelsize = 14, padding = 4)

    ax2 = Axis(fig[1, 2], xlabel = L"Time from now $t$ [Gyr]",
        ylabel = L"Scale factor $a(t)$")
    Ωlab(v) = rich("Ω", subscript("m"), " = $v")
    for (Ωc, ΩΛc, col, name) in ((0.315, 0.685, PALETTE.blue, Ωlab("0.315")),
                                 (1.0, 0.0, PALETTE.orange, Ωlab("1, Einstein–de Sitter")),
                                 (0.05, 0.95, PALETTE.green, Ωlab("0.05, Λ-dominated")))
        lines!(ax2, scale_factor(Ωc, ΩΛc)..., color = col, linewidth = 1.8, label = name)
    end
    vlines!(ax2, [0.0], color = PALETTE.black, linestyle = :dash, linewidth = 1.0)
    scatter!(ax2, [0.0], [1.0], color = PALETTE.black, markersize = MARKERSIZE.dense)
    # in data coordinates: placed by relative position this sat six Gyr to the
    # right of the line it labels
    text!(ax2, 0.6, 0.06; text = "Now", align = (:left, :bottom), fontsize = 15)
    ylims!(ax2, 0, 3)
    axislegend(ax2; position = :lt, framevisible = false, labelsize = 14, padding = 4)

    ax3 = Axis(fig[2, 1], xlabel = L"Redshift $z$",
        ylabel = L"Distance modulus $m - M$ [mag]", xscale = log10,
        xticks = logticks(-3, 0))
    # blue is the flat LambdaCDM fit and orange the measurements in the H(z)
    # panel, and the two were swapped here, so one figure gave each colour two
    # meanings.
    scatter!(ax3, z, μ, color = (PALETTE.orange, 0.18), markersize = MARKERSIZE.cloud)
    zfit = 10 .^ range(log10(minimum(z)), log10(maximum(z)), length = 250)
    lines!(ax3, zfit, distance_moduli(zfit, full.H₀, full.Ωm; zmax = maximum(z) * 1.02),
        color = PALETTE.blue, linewidth = 2)
    bz, bμ = binned_medians(z, μ)
    scatter!(ax3, bz, bμ, color = PALETTE.black, markersize = MARKERSIZE.dense)
    text!(ax3, 0.04, 0.96; text = Ho(full.H₀, full.σH₀), space = :relative,
        align = (:left, :top), fontsize = 15, color = PALETTE.blue)
    text!(ax3, 0.04, 0.87; text = @sprintf("RMS %.2f mag", full.rms), space = :relative,
        align = (:left, :top), fontsize = 15, color = PALETTE.blue)
    # Built by hand: the data are drawn at markersize 3 and 18 % opacity, which
    # is right for 2376 overlapping points and invisible in a legend swatch.
    axislegend(ax3,
        [MarkerElement(color = PALETTE.orange, marker = :circle,
                       markersize = MARKERSIZE.key),
         MarkerElement(color = PALETTE.black, marker = :circle,
                       markersize = MARKERSIZE.key),
         LineElement(color = PALETTE.blue, linewidth = 2)],
        ["SN Ia", "Binned medians", "Flat ΛCDM fit"];
        position = :rb, framevisible = false, labelsize = 14, padding = 4)

    ax4 = Axis(fig[2, 2], xlabel = L"Redshift floor $z_{\mathrm{min}}$",
        ylabel = rich("Fitted Ω", subscript("m")))
    Ωs = [f[2].Ωm for f in sweep]
    rail = [f[2].railed for f in sweep]
    hlines!(ax4, [Ωm_BOUNDS[2]], color = PALETTE.black, linestyle = :dot, linewidth = 1.2,
        label = rich("Fit bound Ω", subscript("m"), " = 1"))
    lines!(ax4, cuts, Ωs, color = PALETTE.purple, linewidth = 1.8)
    scatter!(ax4, cuts[.!rail], Ωs[.!rail], color = PALETTE.purple,
        markersize = MARKERSIZE.data, label = rich("Fitted Ω", subscript("m")))
    scatter!(ax4, cuts[rail], Ωs[rail], color = PALETTE.red, marker = :xcross,
        markersize = MARKERSIZE.emphasis, label = "At the bound")
    hlines!(ax4, [0.315], color = PALETTE.green, linestyle = :dash, linewidth = 1.6,
        label = "Planck 2018")
    text!(ax4, 0.96, 0.52;
        text = @sprintf("moves by %.2f\nrails in %d of %d cuts", spread, n_railed, length(cuts)),
        space = :relative, align = (:right, :bottom), fontsize = 15, color = PALETTE.purple)
    ylims!(ax4, 0, 1.12)
    axislegend(ax4; position = :lt, framevisible = false, labelsize = 14, padding = 4)

    rowgap!(fig.layout, 1, 26)
    colgap!(fig.layout, 1, 30)

    path = savefigure(fig, FIGURES, "friedmann_expansion")
    println("\nwrote ", path)
end

main()
