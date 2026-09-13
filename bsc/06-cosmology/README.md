# Cosmology — Dubna 2019

Joint work with [Alexandru Crăciun](https://github.com/Craciun-Alexandru),
presented together at JINR Dubna in 2019.

## `friedmann_expansion.jl`

![Friedmann expansion](figures/friedmann_expansion.png)

Flat ΛCDM fitted to two independent datasets, and the expansion history
integrated for three cosmologies.

### H(z)

| | Fitted | Planck 2018 |
|---|---|---|
| H₀ | **73.3 ± 5.1** km s⁻¹ Mpc⁻¹ | 67.4 ± 0.5 |
| Ω_m | **0.267 ± 0.064** | 0.315 ± 0.007 |

χ²/dof = 0.73 on 10 degrees of freedom, from **twelve** measurements of the
Hubble parameter. Both parameters agree with Planck inside the (large)
uncertainties of this small dataset; the central H₀ sits on the local-measurement
side of the Hubble tension, which is what H(z) compilations of this kind tend to
give.

### The SN Ia Hubble diagram

2376 Type Ia supernova distance moduli out to z = 1.91:

| | Fitted | Planck 2018 |
|---|---|---|
| H₀ | **72.33 ± 0.68** km s⁻¹ Mpc⁻¹ | 67.4 ± 0.5 |
| Ω_m | 0.535 ± 0.049 | 0.315 ± 0.007 |

Residual RMS 0.578 mag.

**H₀ is well determined and Ω_m is not, and the fit does not say so on its own.**
Ω_m comes back with a 9 % formal error, four sigma from Planck, and the only way
to see that the number is worthless is to move the sample. Raising a floor on
redshift walks it from 0.535 to the bound:

| z floor | N | H₀ | Ω_m |
|---|---|---|---|
| 0 | 2376 | 72.33 | 0.535 |
| 0.01 | 2308 | 72.11 | 0.547 |
| 0.02 | 2146 | 69.97 | 0.671 |
| 0.03 | 1993 | 66.78 | 0.890 |
| 0.05 | 1828 | 64.71 | **1.000 — at the bound** |
| 0.12 | 1650 | 64.18 | **1.000 — at the bound** |

Ω_m moves by 0.48 and rails in three of eight cuts, while H₀ drifts by 8 km s⁻¹
Mpc⁻¹ in the same direction — the two are degenerate, and the low-redshift
points are anchoring both. These are compiled distance moduli from mixed
sources, not a standardised sample: the distance scale survives the
heterogeneity and the shape does not. That is why standardised SN Ia compilations
exist, and the fourth panel is the honest way to report it.

The fit is bounded in Ω_m deliberately. Unbounded, a sample carrying no shape
information wanders off to unphysical values and still reports a confident error
bar; ending at a bound is the signal that the data do not constrain the
parameter, so the code returns that flag rather than hiding it.

### Host-galaxy magnitudes, for comparison

`Friedmann_Magnitude.jl` fitted 6571 host-galaxy apparent magnitudes, not
distance moduli. Refitting them gives an effective offset M = −19.57 with a
residual RMS of **1.35 mag** and Ω_m at its bound — against **0.578 mag** for the
distance moduli. Host-galaxy magnitudes are not standardised SN Ia peak
magnitudes, so they fix the distance scale and carry no cosmological information
at all. Kept because it is what the companion file actually did.

## Data

| File | Contents |
|---|---|
| `hubble_parameter.csv` | 12 measurements of (z, H, σ_H) |
| `sn_ia_distance_moduli.csv` | 2376 SN Ia distance moduli, z = 0.0008–1.91 |
| `supernova_magnitudes.csv` | 6604 host-galaxy apparent magnitudes with redshifts |

The distance moduli are Crăciun's reduction of the raw catalogue the project
worked from — selected on `Method == "SNIa"`, deduplicated in redshift, sorted —
recovered from the project directory alongside his `deqSolver.jl`. **The raw file
was not preserved and the compilation is not named in anything that survives.**
The column set, a distance modulus keyed by measurement method, is characteristic
of a redshift-independent distance compilation of the NED-D kind, but that is an
inference and not a citation. Treat the provenance as unverified; it is the
reason the Ω_m instability above should be read as a property of an unidentified
heterogeneous compilation rather than as a result about cosmology.

## What changed from the original

**The original did not run.** `Friedmann_Eq.jl` called a function `Friedmann`
that was never defined — the definition present is
`FriedmannTimeDependentNeutralFluid` — and used `fit.param` some forty lines
before `fit` was assigned. It also declared

```julia
const G = 6.67408*1e-11        #Cosmological constant
```

which is the gravitational constant, not Λ. `Friedmann_Magnitude.jl` read its
data through a hardcoded `C:\Users\GoG\Desktop\...` path and, inside an ODE
right-hand side, wrote `u[3] = ...` — mutating the state vector instead of
returning a derivative.

The original fitted a **straight line** to H(z) (`LineFit(t, p) = p[1]*t .+ p[2]`);
flat ΛCDM is fitted here instead. The magnitude analysis is kept, but its route
is not: the original drove a `MonteCarloProblem` with randomised initial
conditions through `build_loss_objective`, which is replaced by the closed form
for the luminosity distance. The equation-of-state survey sketched in the
original's docstring — quintessence, Chaplygin and modified Chaplygin gas — was
never implemented there and is not implemented here.

The comoving integral is tabulated once per trial Ω_m and interpolated rather
than integrated per supernova per iteration, which is what makes a 2376-point
fit and an eight-point stability sweep run in seconds.
