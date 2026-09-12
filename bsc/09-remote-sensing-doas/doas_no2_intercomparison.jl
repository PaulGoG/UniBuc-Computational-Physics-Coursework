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
    md_col = md[!, "NO2_SlCol_no2_"];  md_err = md[!, "NO2_SlErr_no2_"]
    sw_col = sw[!, "NO2_SlCol_NO2_"];  sw_err = sw[!, "NO2_SlErr_NO2_"]
    md_rms = md[!, "NO2_RMS"];         sw_rms = sw[!, "NO2_RMS"]
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

    fig = Figure(size = (1040, 470))
    secs(t) = Dates.value(Second(t - Time(0, 0, 0)))
    ticks = [Time(h, 0, 0) for h in 9:15]

    ax1 = Axis(fig[2, 1], xlabel = "Time [UTC]",
        ylabel = L"NO$_2$ dSCD [molec cm$^{-2}$]",
        xticks = (secs.(ticks), [Dates.format(t, "HH:MM") for t in ticks]))
    ax1.xticklabelrotation = π/4

    handles = []
    for (α, colour) in zip(shared, (PALETTE.blue, PALETTE.orange, PALETTE.green,
                                    PALETTE.purple, PALETTE.red, PALETTE.sky))
        m = md_ok .& (md_elev .== α)
        count(m) == 0 && continue
        push!(handles, scatter!(ax1, secs.(md_t[m]), md_col[m],
              color = colour, markersize = 6))
    end

    ax2 = Axis(fig[2, 2], xlabel = L"Elevation angle $\alpha$ [deg]",
        ylabel = L"Mean NO$_2$ dSCD [molec cm$^{-2}$]")
    md_means = [mean(md_col[md_ok .& (md_elev .== α)]) for α in shared]
    sw_means = [mean(sw_col[sw_ok .& (sw_elev .== α)]) for α in shared]
    md_sds = [std(md_col[md_ok .& (md_elev .== α)]) for α in shared]
    sw_sds = [std(sw_col[sw_ok .& (sw_elev .== α)]) for α in shared]
    l_m = scatterlines!(ax2, shared, md_means, color = PALETTE.blue, markersize = 10)
    errorbars!(ax2, shared, md_means, md_sds, color = PALETTE.blue, whiskerwidth = 8)
    l_s = scatterlines!(ax2, shared, sw_means, color = PALETTE.orange,
        markersize = 10, marker = :rect)
    errorbars!(ax2, shared, sw_means, sw_sds, color = PALETTE.orange, whiskerwidth = 8)

    Legend(fig[1, 1:2], [l_m, l_s], ["MAX-DOAS (ground)", "SWING (UAV)"],
        orientation = :horizontal, framevisible = false, labelsize = 17, colgap = 26)
    rowsize!(fig.layout, 2, Relative(0.84))
    println("wrote ", savefigure(fig, FIGURES, "doas_no2_intercomparison"))
end

main()
