# UniBuc Computational Physics Coursework — `legacy`

Verbatim snapshot of the code written during my studies at the Faculty of
Physics, University of Bucharest. It is kept as a reference point: `main`
reorganises, translates and repairs this material, and every change made there
reads as a diff against this branch.

Nothing here has been corrected. Comments, identifiers and printed output are in
Romanian, the numerical methods stand as they were understood at the time, and
the programming is script-style throughout. Several programs contain mistakes
that were never caught; they are left in place deliberately.

## Layout

The two directories mirror the archives this material was recovered from.

| Directory | Languages | Period |
|---|---|---|
| `Code_Archive/` | C, C++, MATLAB/Octave | 2018–2020 |
| `Julia-Workflow-FFUB/` | Julia | 2019–2023 |

### `Code_Archive/`

| Path | Contents |
|---|---|
| `Numerical_Methods/` | Adams predictor–corrector, Runge–Kutta 4, ODE systems, pendulum family, <sup>238</sup>U decay chain |
| `Old_2018/C_C++/` | Harmonic, damped and driven oscillators, Euler and RK integrators, phase-space portraits, Compton scattering, Ising model, strange attractors |
| `Old_2018/Octave/` | Attractor distributions, histogramming, pseudorandom sampling |
| `Optics/` | Erfle eyepiece and Tessar objective ray traces |

### `Julia-Workflow-FFUB/`

Directory suffixes record the course stage: `_L_n` is year *n* of the BSc,
`_M_n` is year *n* of the MSc.

| Path | Contents |
|---|---|
| `Dubna_2019_Cosmology/` | Friedmann equations, distance modulus against a supernova catalogue |
| `Examen_PDF_MN_II_L_4/` | RK4 for scalar and vector ODEs, method of lines for a solitonic PDE |
| `FPECA_M_2/` | Breit–Wigner resonance, Hubble diagram |
| `Fisiune_M_2/` | Neutron-induced fission of <sup>235</sup>U: yields, *Q*-values, excitation energies, neutron multiplicities and spectra |
| `IRM_M_1/` | Bethe–Bloch stopping power, α attenuation, NaI(Tl) activation fit, neutron reactions, natural-Cd cross-section |
| `Radionuclizi_M_1/` | Mass defects and decay chains from the AUDI and Möller tables |
| `Single_Files/` | Hg lamp calibration, aircraft-accident frequency, linear-system maximisation, Monte Carlo sampling |
| `Statistics_symul_L_4/` | Bayesian updating of a Beta posterior over repeated coin tosses |
| `Teledetectie_L_4/` | MAX-DOAS and SWING differential optical absorption spectroscopy |

## Provenance and attribution

The material was recovered from two personal repositories, now private and
archived. Their commit histories are not carried over: this branch is a single
import of the final state of both, so that the public repository starts without
the generated artefacts and bulk data files that had accumulated in them.

`Dubna_2019_Cosmology/` is joint work with
[Alexandru Crăciun](https://github.com/Craciun-Alexandru), presented together at
JINR Dubna in 2019.

`Statistics_symul_L_4/CoinTossAveraging.gif` is the only file altered in this
import: the original 1001-frame animation was 14.6 MB, and it was re-encoded to
640×450 at every fourth frame to keep the repository a reasonable size. The
script that produces it is unchanged and will regenerate the original.
