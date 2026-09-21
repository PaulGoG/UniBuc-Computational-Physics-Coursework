# Intercomparison of NO₂ differential slant column densities retrieved by DOAS
# from two instruments during a single campaign day, 18 March 2019 at 44° N,
# 28° E in south-east Romania:
#
#   * a ground-based MAX-DOAS running an elevation-scanning sequence;
#   * SWING, the BIRA-IASB UAV-borne whiskbroom imaging spectrometer.
#
# Both retrievals are QDOAS output, fitting NO₂ jointly with O₃, O₄, H₂O and the
# Ring pseudo-absorber.
#
# Ported from SWING_DOAS.jl. The corrections that matter:
#
#   1. **The plotted quantity was mislabelled.** `NO2.SlCol(no2)` is a
#      differential *slant column density* in molec cm⁻², not a concentration
#      and not a mixing ratio. The original axis read "Concentratie NO₂" with no
#      units at all.
#   2. **The retrieval uncertainties were discarded.** `NO2.SlErr` is present in
#      both files and never read, so the comparison carried no error bars.
#   3. **No quality filtering.** `NO2.RMS` spans two orders of magnitude —
#      1.5×10⁻⁴ to 1.8×10⁻² for MAX-DOAS, 3.0×10⁻⁴ to 3.5×10⁻² for SWING — and
#      high-residual fits were plotted on equal footing with good ones. Standard
#      DOAS practice rejects above an RMS threshold; one is applied here.
#
#      This one is not hypothetical. The lab report the original was written for
#      names an unusual NO2 peak at 48 deg and 60 deg around 11:20 as its
#      headline observation on the SWING series. Those two points are the two
#      largest dSCD in that window and carry the worst residuals in the whole
#      file — RMS 0.035 and 0.025 against 5e-4 for the good fits beside them —
#      so the filter deletes them. The elevated NO2 there is real and shows in
#      the 18 deg and 24 deg points, which survive; the named peak is a
#      retrieval failure that was read as physics.
#   4. `88 .- servo_byte` was an undocumented magic constant, and the column kept
#      the name `…position_byte` after being converted to degrees. It is
#      empirically exact — the bytes {28, 34, …, 82} map onto a clean 6° grid
#      matching the MAX-DOAS sequence — so it is kept, named and explained.
#   5. Seconds were rounded rather than carried, producing tick labels of the
#      form `hh:mm:60` on seven rows; `savefig` used a Windows separator into a
#      directory that does not exist, so on Linux it silently wrote files
#      literally named `Grafice\MaxDoas.png`.
#
# A caveat the original did not state: MAX-DOAS looks up from the ground and
# SWING looks down from a UAV, so at equal nominal elevation the two are not
# measuring the same air mass. Converting dSCD to a vertical column would need
# differential air-mass factors, which are not attempted here.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames, Statistics, Dates
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Reject DOAS fits whose residual RMS exceeds this."
const RMS_MAX = 5e-3
"Servo byte b maps to elevation angle 88 - b degrees."
servo_to_elevation(b) = 88 - b

"Fractional UTC hours to a Time, carrying seconds properly."
function fractional_hours(h)
    total = round(Int, h * 3600)
    return Time(0, 0, 0) + Second(total)
end

function main()
    md = CSV.read(joinpath(DATA, "MAXDOAS.csv"), DataFrame; delim = '\t', normalizenames = true)
    sw = CSV.read(joinpath(DATA, "SWING.csv"), DataFrame; delim = '\t', normalizenames = true)

    md_t = fractional_hours.(md[!, "Fractional_time"])
    sw_t = fractional_hours.(sw[!, "Fractional_time"])
    md_col = md[!, "NO2_SlCol_no2_"]
    md_err = md[!, "NO2_SlErr_no2_"]
    sw_col = sw[!, "NO2_SlCol_NO2_"]
    sw_err = sw[!, "NO2_SlErr_NO2_"]
    md_rms = md[!, "NO2_RMS"]
    sw_rms = sw[!, "NO2_RMS"]
    md_elev = md[!, "Elev_viewing_angle"]
    sw_elev = servo_to_elevation.(sw[!, "UAV_servo_sent_position_byte"])

    @printf("MAX-DOAS %d spectra, %s to %s\n", nrow(md), md_t[1], md_t[end])
    @printf("SWING    %d spectra, %s to %s\n", nrow(sw), sw_t[1], sw_t[end])
    @printf("RMS range: MAX-DOAS %.1e–%.1e, SWING %.1e–%.1e\n",
        minimum(md_rms), maximum(md_rms), minimum(sw_rms), maximum(sw_rms))

    md_ok = md_rms .<= RMS_MAX
    sw_ok = sw_rms .<= RMS_MAX
    @printf("RMS ≤ %.0e keeps %d/%d MAX-DOAS and %d/%d SWING spectra\n",
        RMS_MAX, count(md_ok), nrow(md), count(sw_ok), nrow(sw))

    shared = sort(collect(intersect(Set(md_elev[md_ok]), Set(sw_elev[sw_ok]))))
    @printf("elevation angles in common: %s\n", join(shared, ", "))

    for α in shared
        m = md_ok .& (md_elev .== α)
        s = sw_ok .& (sw_elev .== α)
        (count(m) == 0 || count(s) == 0) && continue
        @printf("  α = %2d°: MAX-DOAS ⟨dSCD⟩ = %+.3e ± %.1e (n=%3d),  SWING %+.3e ± %.1e (n=%3d)\n",
            α, mean(md_col[m]), std(md_col[m]), count(m),
            mean(sw_col[s]), std(sw_col[s]), count(s))
    end

    fig = Figure(size = (1300, 900))
    secs(t) = Dates.value(Second(t - Time(0, 0, 0)))
    ticks = [Time(h, 0, 0) for h in 9:15]
    tickspec = (secs.(ticks), [Dates.format(t, "HH:MM") for t in ticks])

    # dSCD in units of 1e16: every tick of this axis otherwise repeated the
    # factor, and Makie wrote the first of them as 1x10^16 rather than 10^16.
    const_scale = 1e16
    dscd_label = L"NO$_2$ dSCD [$10^{16}$ molec cm$^{-2}$]"

    # One colour per elevation angle, shared by the two time series. The figure
    # previously carried a single instrument legend over a left panel that was
    # colour-coded by elevation and held MAX-DOAS alone, so a reader following
    # the legend read its orange points -- MAX-DOAS at 18 degrees -- as SWING.
    elev_colours = Dict(zip(shared,
        (PALETTE.blue, PALETTE.orange, PALETTE.green, PALETTE.purple, PALETTE.red,
            PALETTE.sky,),))

    function series_panel(ax, t, col, elev, ok, name)
        for α in shared
            m = ok .& (elev .== α)
            count(m) == 0 && continue
            scatter!(ax, secs.(t[m]), col[m] ./ const_scale,
                color = elev_colours[α], markersize = MARKERSIZE.dense,)
        end
        text!(ax, 0.02, 0.97; text = name, space = :relative,
            align = (:left, :top), fontsize = 16,)
        ax.xticklabelrotation = π/4
    end

    ax1 = Axis(fig[2, 1], xlabel = "Time [UTC]", ylabel = dscd_label, xticks = tickspec)
    series_panel(ax1, md_t, md_col, md_elev, md_ok, "MAX-DOAS (ground)")

    ax2 = Axis(fig[2, 2], xlabel = "Time [UTC]", ylabel = "", xticks = tickspec)
    series_panel(ax2, sw_t, sw_col, sw_elev, sw_ok, "SWING (UAV)")
    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, grid = false)

    Legend(fig[1, 1:2],
        [MarkerElement(color = elev_colours[α], marker = :circle,
             markersize = MARKERSIZE.key,) for α in shared],
        [latexstring("\\alpha = $(α)^\\circ") for α in shared],
        "Elevation angle";
        orientation = :horizontal, framevisible = false, labelsize = 16,
        titlesize = 16, titleposition = :left, colgap = 22,)

    # The quality cut is the module's result, and it was nowhere in the figure.
    ax3 = Axis(fig[3, 1], xlabel = "Time [UTC]",
        ylabel = "Fit residual RMS", yscale = log10,
        xticks = tickspec, yticks = logticks(-4, -1),)
    ax3.xticklabelrotation = π/4
    # headroom above the worst fits, so that the two notes sit in a clear band
    # rather than on the data or on each other
    ylims!(ax3, 2e-4, 0.9)
    # Colour means elevation angle in the panels above, so the instruments are
    # separated here by shape and tone instead of by hue.
    scatter!(ax3, secs.(md_t), md_rms, color = (PALETTE.black, 0.45),
        markersize = MARKERSIZE.dense, label = "MAX-DOAS (ground)",)
    scatter!(ax3, secs.(sw_t), sw_rms, color = (:grey45, 0.75), marker = :rect,
        markersize = MARKERSIZE.dense, label = "SWING (UAV)",)
    hlines!(ax3, [RMS_MAX], color = PALETTE.red, linestyle = :dash, linewidth = 1.6)
    axislegend(ax3, position = :lt, framevisible = false, labelsize = 14, padding = 2)
    text!(ax3, 0.98, 0.98;
        text = rich("Rejected above RMS ", rsci(RMS_MAX; digits = 0), "\n",
            @sprintf("keeps %d/%d MAX-DOAS and %d/%d SWING",
                count(md_ok), nrow(md), count(sw_ok), nrow(sw))),
        space = :relative, align = (:right, :top), fontsize = 14,
        justification = :right, color = PALETTE.red,)

    # the two spectra the lab report named as a physical peak
    named = findall(i -> sw_rms[i] > RMS_MAX && Time(11, 10) <= sw_t[i] <= Time(11, 25),
        eachindex(sw_rms),)
    scatter!(ax3, secs.(sw_t[named]), sw_rms[named], color = PALETTE.red,
        marker = :xcross, markersize = MARKERSIZE.emphasis,)
    text!(ax3, secs(Time(11, 20)), 0.075;
        text = "The 11:20 peak named by the report:\nthe worst fits in the file",
        align = (:center, :bottom), fontsize = 14, color = PALETTE.red,
        justification = :center,)

    ax4 = Axis(fig[3, 2], xlabel = L"Elevation angle $\alpha$ [deg]",
        ylabel = "", xticks = (shared, [latexstring(string(α)) for α in shared]),)
    md_means = [mean(md_col[md_ok .& (md_elev .== α)]) for α in shared] ./ const_scale
    sw_means = [mean(sw_col[sw_ok .& (sw_elev .== α)]) for α in shared] ./ const_scale
    md_sds = [std(md_col[md_ok .& (md_elev .== α)]) for α in shared] ./ const_scale
    sw_sds = [std(sw_col[sw_ok .& (sw_elev .== α)]) for α in shared] ./ const_scale
    scatterlines!(ax4, shared, md_means, color = PALETTE.black,
        markersize = MARKERSIZE.data, label = "MAX-DOAS (ground)",)
    errorbars!(ax4, shared, md_means, md_sds, color = PALETTE.black, whiskerwidth = 10)
    scatterlines!(ax4, shared, sw_means, color = :grey45, linestyle = :dash,
        markersize = MARKERSIZE.data, marker = :rect, label = "SWING (UAV)",)
    errorbars!(ax4, shared, sw_means, sw_sds, color = :grey45, whiskerwidth = 10)
    ax4.ylabel = L"Mean NO$_2$ dSCD [$10^{16}$ molec cm$^{-2}$]"
    axislegend(ax4, position = :rt, framevisible = false, labelsize = 15, padding = 2)

    rowsize!(fig.layout, 2, Relative(0.44))
    rowgap!(fig.layout, 1, 6)
    rowgap!(fig.layout, 2, 22)
    colgap!(fig.layout, 26)
    println("wrote ", savefigure(fig, FIGURES, "doas_no2_intercomparison"))
end

main()
