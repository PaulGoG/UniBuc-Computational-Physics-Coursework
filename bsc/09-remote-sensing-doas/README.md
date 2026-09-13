# Remote sensing — DOAS, BSc year 4 (2020–2021)

![NO2 intercomparison](figures/doas_no2_intercomparison.png)

NO₂ differential slant column densities retrieved by DOAS from two instruments
on 18 March 2019 at 44° N, 28° E: a ground-based elevation-scanning MAX-DOAS and
SWING, the BIRA-IASB UAV-borne whiskbroom imaging spectrometer. Both are QDOAS
output fitting NO₂ jointly with O₃, O₄, H₂O and the Ring pseudo-absorber.

| Elevation | MAX-DOAS ⟨dSCD⟩ | SWING ⟨dSCD⟩ |
|---|---|---|
| 12° | 2.32 × 10¹⁶ | 7.38 × 10¹⁵ |
| 18° | 1.31 × 10¹⁶ | 2.67 × 10¹⁶ |
| 30° | 7.40 × 10¹⁵ | 1.03 × 10¹⁶ |
| 48° | 4.90 × 10¹⁵ | 3.68 × 10¹⁵ |
| 60° | 4.55 × 10¹⁵ | 1.28 × 10¹⁵ |

molec cm⁻². The MAX-DOAS series falls monotonically with elevation, as the
shortening slant path requires.

Corrections to `SWING_DOAS.jl`:

- **the plotted quantity was mislabelled.** `NO2.SlCol` is a differential *slant
  column density* in molec cm⁻², not a concentration and not a mixing ratio. The
  axis read "Concentratie NO₂" with no units
- **the retrieval uncertainties were discarded.** `NO2.SlErr` is present in both
  files and was never read, so the comparison had no error bars
- **no quality filtering.** `NO2.RMS` spans two orders of magnitude in both
  instruments and high-residual fits were plotted on equal footing with good
  ones. An RMS ≤ 5 × 10⁻³ cut keeps 376/377 MAX-DOAS and 445/510 SWING spectra
- `88 - servo_byte` was an undocumented magic constant, and the column kept the
  name `…position_byte` after becoming degrees. It is empirically exact — the
  bytes map onto a clean 6° grid matching the MAX-DOAS sequence — so it is kept,
  named and explained
- seconds were rounded rather than carried, producing `hh:mm:60` labels on seven
  rows; `savefig` used a Windows separator into a directory that does not exist,
  so on Linux it silently wrote files literally named `Grafice\MaxDoas.png`

## What the quality cut removes

The lab report this code was written for survives, and its headline observation
is worth recording because the filter above deletes it.

Reporting on the SWING series, it singles out *"excepție făcând unghiurile de 48
și 60 de grade în jurul orei 11:20 unde se observă un pic neobișnuit al
concentrației de NO₂"* — an unusual NO₂ peak at 48° and 60° around 11:20. The
feature is real in the data, and it is an artefact:

| Time | Elevation | dSCD [molec cm⁻²] | RMS | |
|---|---|---|---|---|
| 11:20:35 | 60° | 6.16 × 10¹⁶ | 0.0350 | **rejected** |
| 11:16:14 | 48° | 4.17 × 10¹⁶ | 0.0254 | **rejected** |
| 11:13:31 | 18° | 1.96 × 10¹⁶ | 0.000792 | kept |
| 11:18:57 | 18° | 1.71 × 10¹⁶ | 0.000758 | kept |
| 11:12:59 | 24° | 1.05 × 10¹⁶ | 0.000506 | kept |

The two points singled out are the two largest dSCD in that window and they
carry the worst spectral residuals in the entire SWING file, five to seven times
the 5 × 10⁻³ threshold and thirty to seventy times the residuals of the good
fits beside them. Four of the twenty-two spectra in the window are rejected;
these are two of them.

There *is* elevated NO₂ around 11:10–11:25 — the 18° and 24° points at
1–2 × 10¹⁶ are three to four times their own medians, they survive the cut, and
low elevations look through more polluted air, which is what one expects. The
peak the report named is not it. That is the whole case for quality filtering,
made concrete: without it a retrieval failure reads as a physical signal, and
the report drew exactly that conclusion.

## Caveat

MAX-DOAS looks up from the ground and SWING looks down from a UAV, so at equal
nominal elevation they are not sampling the same air mass. Turning dSCD into a
vertical column needs differential air-mass factors, which is not attempted
here.
