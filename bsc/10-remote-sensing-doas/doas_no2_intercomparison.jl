# NO₂ differential slant column densities (dSCD) retrieved by DOAS from two
# instruments operated side by side on the ground on 18 March 2019, at the
# position the files record as 44° N, 28° E, viewing azimuth 185°:
#
#   * a MAX-DOAS running the elevation sequence 3, 6, 12, 18, 30, 48, 60, 89°;
#   * SWING, the BIRA-IASB whiskbroom imaging spectrometer built for UAV
#     operation, here on the ground beside the MAX-DOAS, stepping its servo from
#     60° down to 6° in steps of 6°.
#
# Both files are QDOAS output: NO₂ fitted jointly with O₃, O₄, H₂O and the Ring
# pseudo-absorber, with the fit residual RMS and the retrieval error of every
# column alongside.
#
# Ported from SWING_DOAS.jl in Julia-Workflow-FFUB/Teledetectie_L_4/ on the
# `legacy` branch, which plotted `NO2.SlCol` against time per elevation angle.
# Changed here:
#
#   1. The quantity is a differential slant column density in molec cm⁻², not
#      a concentration; the axes say so.
#   2. The retrieval errors `NO2.SlErr` were never read; they are drawn.
#   3. No quality filter was applied. `NO2.RMS` spans 1.5 × 10⁻⁴ to 3.5 × 10⁻²
#      and fits above RMS_MAX are excluded. The cut is not neutral between the
#      instruments: every SWING spectrum at 6° has RMS between 7 and 13 × 10⁻³,
#      so the whole SWING 6° series goes, with the 12° spectra after 14:00 and
#      six spectra at 48° and 60°; one MAX-DOAS spectrum goes. `main` prints
#      the count per elevation angle.
#   4. The servo byte b is turned into an elevation angle as 88 − b, which puts
#      the bytes 28, 34, …, 82 on the 6° grid of the MAX-DOAS sequence; the
#      constant is the original's and is not otherwise documented.
#   5. Time is read from the files' own hh:mm:ss column, which agrees with the
#      fractional-hour column to the second (asserted), instead of being
#      rebuilt from fractional hours with a rounding that produced hh:mm:60.

include(joinpath(@__DIR__, "..", "..", "activate.jl"))

using Printf, CSV, DataFrames, Statistics, Dates
include(joinpath(@__DIR__, "..", "..", "theme.jl"))

const FIGURES = joinpath(@__DIR__, "figures")
const DATA = joinpath(@__DIR__, "data")

"Fits with residual RMS above this are excluded."
const RMS_MAX = 5e-3
"Fewest accepted spectra for an elevation angle to enter the comparison."
const MIN_SPECTRA = 5
"dSCD are plotted in units of this many molec cm⁻²."
const DSCD_UNIT = 1e16
"Servo byte b maps to elevation angle 88 − b degrees."
servo_to_elevation(b) = 88 - b

"Seconds since midnight of a `Time`."
seconds(t) = Dates.value(t - Time(0)) / 1e9

"""
    read_qdoas(file, column, error, elevation)

One QDOAS output file as a table of time, elevation angle, NO₂ dSCD, its
retrieval error and the fit RMS; `elevation` maps the raw table to angles in
degrees. The hh:mm:ss column is asserted against the fractional-hour column.
"""
function read_qdoas(file, column, error, elevation)
    raw = CSV.read(joinpath(DATA, file), DataFrame; delim = '\t', normalizenames = true)
    t = raw[!, "Time_hh_mm_ss_"]
    eltype(t) <: Time || (t = Time.(String.(t), "HH:MM:SS"))
    s = seconds.(t)
    maximum(abs.(s .- 3600 .* raw.Fractional_time)) <= 1.5 ||
        Base.error("$file: hh:mm:ss and fractional hours disagree by more than a second")
    return DataFrame(
        time = t, seconds = s, elevation = elevation(raw), dscd = raw[!, column],
        error = raw[!, error], rms = raw.NO2_RMS,)
end

"Mean, day spread (1 s.d.), median retrieval error and count of the accepted spectra of `d` at `α`."
function elevation_stats(d, α)
    m = d.ok .& (d.elevation .== α)
    return (
        mean = mean(d.dscd[m]), sd = std(d.dscd[m]), err = median(d.error[m]), n = count(m),)
end

function main()
    md = read_qdoas("MAXDOAS.csv", "NO2_SlCol_no2_", "NO2_SlErr_no2_",
        r -> r[!, "Elev_viewing_angle"],)
    sw = read_qdoas("SWING.csv", "NO2_SlCol_NO2_", "NO2_SlErr_NO2_",
        r -> servo_to_elevation.(r[!, "UAV_servo_sent_position_byte"]),)
    for (name, d) in (("MAX-DOAS", md), ("SWING", sw))
        @printf("%-8s %3d spectra, %s to %s, RMS %.1e–%.1e, median retrieval error %.1e molec cm⁻² (%.0f %% of the median column)\n",
            name, nrow(d), d.time[1], d.time[end], extrema(d.rms)...,
            median(d.error),
            100 * median(d.error ./ abs.(d.dscd)))
    end

    md.ok = md.rms .<= RMS_MAX
    sw.ok = sw.rms .<= RMS_MAX
    @printf("RMS ≤ %.0e keeps %d/%d MAX-DOAS and %d/%d SWING spectra; excluded per elevation angle:\n",
        RMS_MAX, count(md.ok), nrow(md), count(sw.ok), nrow(sw))
    for (name, d) in (("MAX-DOAS", md), ("SWING", sw)), α in sort(unique(d.elevation))

        m = d.elevation .== α
        r = count(m .& .!d.ok)
        r == 0 && continue
        @printf("  %-8s α = %2d°: %2d of %2d, RMS %.1e–%.1e\n", name, α, r, count(m),
            extrema(d.rms[m .& .!d.ok])...)
    end
    all(.!sw.ok[sw.elevation .== 6]) ||
        error("the whole SWING 6° series was expected above RMS_MAX")

    # elevation angles both instruments hold with enough accepted spectra
    accepted(d, α) = count(d.ok .& (d.elevation .== α))
    shared = sort([α
                   for α in intersect(md.elevation, sw.elevation)
                   if accepted(md, α) >= MIN_SPECTRA && accepted(sw, α) >= MIN_SPECTRA])
    println("elevation angles compared: ", join(shared, ", "), "°")
    stats = [(md = elevation_stats(md, α), sw = elevation_stats(sw, α)) for α in shared]
    for (α, s) in zip(shared, stats)
        @printf("  α = %2d°: MAX-DOAS ⟨dSCD⟩ = %.2e, spread %.1e, retrieval error %.1e (n = %2d);  SWING %.2e, spread %.1e, retrieval error %.1e (n = %2d)\n",
            α, s.md.mean, s.md.sd, s.md.err, s.md.n, s.sw.mean, s.sw.sd, s.sw.err, s.sw.n)
    end
    md_means = [s.md.mean for s in stats]
    sw_means = [s.sw.mean for s in stats]
    issorted(md_means, rev = true) ||
        error("the MAX-DOAS mean dSCD does not fall with elevation")
    issorted(sw_means[2:end], rev = true) ||
        error("the SWING mean dSCD does not fall with elevation above $(shared[2])°")
    err_md = median(md.error[md.ok])
    err_sw = median(sw.error[sw.ok])

    # the failed fits at 48° and 60° late in the morning, the largest SWING
    # column of the day among them
    failed_high = findall(i -> !sw.ok[i] && sw.elevation[i] >= 48, eachindex(sw.ok))
    top = argmax(sw.dscd)
    top in failed_high ||
        error("the largest SWING column was expected among the failed high-elevation fits")
    sw.error[top] / sw.dscd[top] > 0.25 ||
        error("the largest SWING column was expected to carry a retrieval error above 25 %")
    median_at(α) = median(sw.dscd[sw.ok .& (sw.elevation .== α)])
    excess = [sw.dscd[i] / median_at(sw.elevation[i]) for i in failed_high]
    @printf("%d failed SWING fits at 48° and 60° between %s and %s, RMS %.3f–%.3f, columns %.0f to %.0f times the medians of their elevation angles; the largest column of the day, %.2e at %d° at %s, is one of them, with a retrieval error of %.0f %%\n",
        length(failed_high), sw.time[failed_high[1]], sw.time[failed_high[end]],
        extrema(sw.rms[failed_high])..., extrema(excess)..., sw.dscd[top], sw.elevation[top],
        sw.time[top], 100 * sw.error[top] / sw.dscd[top])

    # --- figure ---------------------------------------------------------------
    fig = Figure(size = (1600, 1100))
    ticks = [Time(h, 0, 0) for h in 9:15]
    tickspec = (seconds.(ticks), [Dates.format(t, "HH:MM") for t in ticks])
    dscd_label = L"NO$_2$ dSCD [$10^{16}$ molec cm$^{-2}$]"
    colours = Dict(zip(shared,
        (PALETTE.blue, PALETTE.orange, PALETTE.green, PALETTE.purple,
            PALETTE.red, PALETTE.sky,),))

    # one colour per elevation angle, shared by the two time series
    function series_panel!(ax, d, name)
        for α in shared
            m = d.ok .& (d.elevation .== α)
            errorbars_unstroked!(
                ax, d.seconds[m], d.dscd[m] ./ DSCD_UNIT, d.error[m] ./ DSCD_UNIT,
                color = (colours[α], 0.6), linewidth = GUIDE_WIDTH, whiskerwidth = 6,)
            scatter!(ax, d.seconds[m], d.dscd[m] ./ DSCD_UNIT, color = colours[α],
                markersize = MARKERSIZE.dense,)
        end
        text!(ax, 0.02, 0.97; text = name, space = :relative, align = (:left, :top),
            fontsize = ANNOTATION_SIZE,)
        ax.xticklabelrotation = π / 4
    end
    ax1 = Axis(fig[2, 1], xlabel = "Time [UTC]", ylabel = dscd_label, xticks = tickspec)
    series_panel!(ax1, md, "MAX-DOAS")
    ax2 = Axis(fig[2, 2], xlabel = "Time [UTC]", xticks = tickspec)
    series_panel!(ax2, sw, "SWING")
    linkyaxes!(ax1, ax2)
    ylims!(ax1, -0.6, 5.3)
    hideydecorations!(ax2, grid = false)

    # the quality cut: instruments told apart by marker and tone, colour being
    # taken by the elevation angle above
    ax3 = Axis(
        fig[3, 1], xlabel = "Time [UTC]", ylabel = "Fit residual RMS", yscale = log10,
        xticks = tickspec, yticks = logticks(-4, -1),)
    ax3.xticklabelrotation = π / 4
    ylims!(ax3, 1e-4, 1.0)
    scatter!(ax3, md.seconds, md.rms, color = (PALETTE.black, 0.45),
        markersize = MARKERSIZE.dense,)
    scatter!(ax3, sw.seconds, sw.rms, color = (:grey45, 0.75), marker = :rect,
        markersize = MARKERSIZE.dense,)
    hlines!(ax3, [RMS_MAX], color = PALETTE.red, linestyle = :dash, linewidth = GUIDE_WIDTH)
    text!(ax3, 0.98, 0.98;
        text = rich("Excluded above RMS ", rsci(RMS_MAX; digits = 0), ":\n",
            @sprintf("%d of %d MAX-DOAS, %d of %d SWING,", count(.!md.ok), nrow(md),
                count(.!sw.ok), nrow(sw)), "\nthe whole SWING 6° series among them",),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        justification = :right, color = PALETTE.red,)
    scatter!(ax3, sw.seconds[failed_high], sw.rms[failed_high], color = PALETTE.red,
        marker = :xcross, markersize = MARKERSIZE.emphasis,)
    text!(ax3, seconds(Time(11, 45)), 0.03;
        text = @sprintf("%d failed SWING fits at 48° and 60°,\nthe largest column of the day among them",
            length(failed_high)),
        align = (:left, :center), fontsize = ANNOTATION_SIZE, color = PALETTE.red,)

    ax4 = Axis(fig[3, 2], xlabel = L"Elevation angle $\alpha$ [deg]",
        ylabel = L"Mean NO$_2$ dSCD [$10^{16}$ molec cm$^{-2}$]",
        xticks = (shared, string.(shared)),)
    for (d, offset, colour, marker, style) in ((md, -0.5, PALETTE.black, :circle, :solid),
        (sw, +0.5, :grey45, :rect, :dash))
        s = [elevation_stats(d, α) for α in shared]
        errorbars_unstroked!(ax4, shared .+ offset, [x.mean for x in s] ./ DSCD_UNIT,
            [x.sd for x in s] ./ DSCD_UNIT, color = colour, whiskerwidth = 10,
            linewidth = GUIDE_WIDTH,)
        scatterlines!(
            ax4, shared .+ offset, [x.mean for x in s] ./ DSCD_UNIT, color = colour,
            marker = marker, linestyle = style,)
    end
    text!(ax4, 0.97, 0.95;
        text = rich("Bars: 1 s.d. over the day\nMedian retrieval errors\n",
            rsci(err_md; digits = 0), " (MAX-DOAS), ", rsci(err_sw; digits = 0), " (SWING)",),
        space = :relative, align = (:right, :top), fontsize = ANNOTATION_SIZE,
        justification = :right,)

    Legend(fig[1, 1:2],
        [
            [MarkerElement(color = colours[α], marker = :circle, markersize = MARKERSIZE.key)
             for α in shared],
            [
                MarkerElement(color = (PALETTE.black, 0.6), marker = :circle, markersize = MARKERSIZE.key),
                MarkerElement(color = :grey45, marker = :rect, markersize = MARKERSIZE.key),],],
        [[latexstring("\\alpha = $(α)^\\circ") for α in shared], ["MAX-DOAS", "SWING"]],
        ["Elevation angle", "Instrument"]; titleposition = :left,)
    rowgap!(fig.layout, 2, 20)
    println("wrote ", savefigure(fig, FIGURES, "doas_no2_intercomparison"))
end

main()
