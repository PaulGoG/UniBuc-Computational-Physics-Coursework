# Geometrical optics

Paraxial lens design and prism dispersion, ported from Octave and from measured
laboratory data.

## `paraxial_systems.jl`

![Cardinal elements](figures/paraxial_systems.png)

Ray-transfer analysis of two classical designs. The ray state is `[n·u ; y]`
and the system matrix is accumulated from the last surface backwards.

| | Erfle eyepiece | Tessar objective |
|---|---|---|
| Surfaces | 9 | 7 |
| Σd | 45.20 mm | 15.65 mm |
| Focal length | **25.39 mm** | **50.79 mm** |
| Front focal, from V₁ | 9.97 mm | 43.11 mm |
| Back focal, from V_N | 15.33 mm | 44.06 mm |
| H₁H₂ separation | 19.71 mm | 1.25 mm |

Two independent checks that the engine and the prescriptions are right: a thin
lens of known power comes back at f = 100.000000 mm against the lensmaker value
with both principal planes at the vertex; and the Tessar lands at **50.79 mm**,
a classic normal lens for the format, which is what a design of this scale
should give.

`"Erlfe"` is a misspelling of **Erfle**, the 1921 wide-field eyepiece; the
nine-surface prescription resolves to the six-element Type II form. `Tessar.m`
is the genuine 1902 Rudolph four-element design.

Corrections to the originals:

- `f = -1/S(1,2)` was labelled "Convergenta sistemului". That is the focal
  length; the convergence, or optical power, is `P = 1/f = -S(1,2)`.
- The principal-plane separation was written `sum(d) - abs(zH1) - abs(zH2)`,
  which is correct only when both principal planes fall on the same side. The
  files' own diagrams place H₁ at `-zH1` and H₂ at `sum(d)+zH2`, so the
  separation is `sum(d) + zH1 + zH2`. **For these two prescriptions the defect
  is latent**: both zH₁ and zH₂ come out negative, so the two formulae agree to
  the digit. It would bite on any design with a principal plane on the other
  side.
- `Tessar.m` acknowledged in a comment that plotting fails for `z1 = zf1` — an
  object at the front focal plane sends the image to infinity — but never
  guarded it. `conjugate` returns `Inf` explicitly here.

Both remain strictly paraxial, as in 2018: no real ray trace, no aberrations,
and a single refractive index per glass. That last point is the real limitation,
because the cemented doublet in the Tessar and the cemented doublet and triplet
in the Erfle exist precisely to achromatise, and without dispersion data that
cannot be evaluated at all. Neither is an aperture stop modelled, so the Tessar
has no f-number.

## `prism_dispersion.jl`

![Prism dispersion](figures/prism_dispersion.png)

Refractive index against wavelength for a prism, measured by minimum deviation
on a goniometer: six lines from 4050 to 6700 Å, with both goniometer readings,
the minimum deviation and the index the student derived. The data come from the
first- and second-year Origin lab project `PrismaOptica.opj`; **there was no
companion code — this analysis did not exist in 2018.**

It is written now because `msc/02-experimental-methods/hg_spectroscope_calibration.jl`
has to *assume* that a spectroscope's drum reading is linear in refractive
index, and fits a Cauchy relation to the drum on that assumption. Here the index
is measured directly, so the relation can be tested against `n` itself.

Three internal checks, all of which the data pass:

| Check | Result |
|---|---|
| δ_min = (α₂ − α₁)/2 | holds to 1.4 × 10⁻¹⁴ ° |
| apex angle from each line | **59.9000 ± 0.0001°**, spread 3.6 × 10⁻⁴ ° |
| n falls with λ (normal dispersion) | yes |

The second is the strong one: the apex angle is a property of the prism, not of
the light, so six lines give six independent determinations of one number — and
they agree to the fourth decimal.

The Cauchy fit gives **A = 1.5071, B = 0.00354 µm², R² = 0.909**. A lands within
0.5 % of BK7 borosilicate crown (1.5046), so the glass is identified. But R² =
0.909 on six points is poor, and the residuals do not scatter: they run
+ + − − + + across the series, a smooth arc of amplitude 2 × 10⁻³. A
two-parameter Cauchy relation is too few across this range.

That is the answer to the question the spectroscope calibration could not
settle. Its residuals show the same systematic structure, and this says the
cause is the dispersion model rather than the linearity assumption or the
instrument.

## Moved out

`hg_spectroscope_calibration.jl` lived here until its source file turned up in
the MSc *Metode Experimentale in Fizica* directory, beside the submitted reports
for that course's interferometry and polarimetry labs. It is
`msc/02-experimental-methods/` now.
