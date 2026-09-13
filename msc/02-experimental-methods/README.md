# Experimental methods in physics — MSc year 1 (2021–2022)

Optical-metrology coursework from *Metode Experimentale in Fizica*.

## `hg_spectroscope_calibration.jl`

![Hg calibration](figures/hg_spectroscope_calibration.png)

Prism spectroscope calibrated against eleven Hg I lines under the Cauchy
relation `x = A + B/λ²`, where `x` is the drum reading, on the working
assumption that the drum is linear in the refractive index — which holds near
minimum deviation.

**A transcription error in the reference data.** The original listed 5789.66 Å
beside 5790.65 Å at adjacent drum divisions 35 and 36. No Hg I line exists at
5789.66; the second member of the yellow doublet is **5769.60 Å**, and a 1 Å
splitting across one division is impossible at the ~13 Å/division local
dispersion. Corrected here — though in fairness it barely moves the fit
(R² 0.97174 against 0.97205), since it is one point of eleven and the two drum
readings differ by one division. It is worth fixing because it is a wrong
physical constant, not because it changed the answer.

What does matter is that the fit is mediocre: R² = 0.972 with an RMS residual of
11.2 divisions on a 2–200 division range, and the residual panel shows clear
systematic structure rather than scatter. The linearity assumption holds only
near minimum deviation, and a two-parameter Cauchy relation is evidently not
enough across 4000–7100 Å.

`bsc/06-geometrical-optics/prism_dispersion.jl` fits the same Cauchy relation to
*measured* refractive indices from a prism goniometer, which is the direct test
of the assumption this calibration has to make.

## What changed from the original

The original computed A and B, never printed them, never displayed the figure
and had its `savefig` commented out, so running it as a script produced nothing
whatsoever. It also used Levenberg–Marquardt with a finite-difference Jacobian
from `p0 = [1, 1]`, nine orders of magnitude away from B ≈ 4 × 10⁹, for a model
that is linear in its parameters and solves in one least-squares step.

## Provenance

Filed under `bsc/06-geometrical-optics/` until the source turned up as
`CALIBRARE.jl` in this course's directory, beside the submitted reports for its
interferometry and polarimetry labs. The repository orders modules by degree and
then chronologically, so it belongs here, in MSc year 1 semester 1.
