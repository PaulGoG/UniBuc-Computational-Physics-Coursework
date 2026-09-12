# Nuclear and particle physics — MSc year 2

## `breit_wigner_interference.jl`

![Breit-Wigner interference](figures/breit_wigner_interference.png)

Two Breit–Wigner resonances with the same quantum numbers, so their amplitudes
add before squaring.

| | Peak [MeV] | FWHM [MeV] |
|---|---|---|
| Resonance 1 alone | 2300.00 | 150.00 |
| Resonance 2 alone | 2340.00 | 320.00 |
| Interference, φ = 30° | 2308.40 | 193.59 |
| Interference, φ = 45° | 2311.90 | 191.00 |

The isolated resonances return their input parameters exactly, which validates
the width measurement. The interference line sits at neither 2300 nor 2340 and
is narrower than the broad resonance — its apparent position and width depend on
the relative phase, which is the point of the exercise.

The physics in the original was right; the width measurement was not. FWHM was
located by scanning for samples satisfying
`minimum(y)*0.005 >= abs(maximum(y)/2 - i)`, an absolute tolerance keyed to the
array minimum, which finds nothing at all if the sampling straddles the half
maximum and then silently returns whichever near-misses came first and last.
Replaced by linear interpolation between the bracketing samples.

## `hubble_parameter.jl`

![Hubble parameter](figures/hubble_parameter.png)

H₀ from three galaxies, with redshifts from five reference lines and distances
from apparent size relative to a standard.

**H₀ = 73.1 ± 11.4 km s⁻¹ Mpc⁻¹**, a Hubble time of 13.4 Gyr — consistent with
both Planck (67.4 ± 0.5) and SH0ES (73.0 ± 1.0). The measured recession
velocities are good: 745, 936 and 2827 km s⁻¹ against catalogue values near 727,
897 and 2820.

Two corrections:

- **A one-character typo.** `Suma_σ² =+ σ_z[i]^2` parses as
  `Suma_σ² = +σ_z[i]^2` — an assignment, not `+=`. The running sum of variances
  was overwritten every iteration and kept only the last term, so every averaged
  redshift uncertainty was too small by roughly √n, and the fit weights derived
  from them were wrong.
- **The distance errors were left out of the fit.** They dominate here: σ_d/d
  runs to 30 %. Weighting by velocity errors alone gives H₀ = 89.1 ± 2.0, because
  the distant, poorly-measured galaxy is over-weighted. Using the effective
  variance σ_v² + H₀²σ_d², iterated to convergence, gives 73.1 ± 11.4 — a
  different central value *and* an honest uncertainty.
