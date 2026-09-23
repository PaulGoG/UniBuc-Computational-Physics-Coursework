# UniBuc Computational Physics Coursework

[![CI](https://github.com/PaulGoG/UniBuc-Computational-Physics-Coursework/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PaulGoG/UniBuc-Computational-Physics-Coursework/actions/workflows/CI.yml?query=branch%3Amain)

Computational physics written during my studies at the Faculty of Physics,
University of Bucharest — BSc 2017–2021, MSc 2021–2023. Originally C, C++,
Octave and Julia; rewritten here in Julia, organised by degree and in
chronological order.

Everything is a script. There is no package scaffolding, and that is deliberate:
these were scripts when they were written and they remain scripts, just ones
that run, say what they compute, and check themselves against something.

## Layout

```
.
├── Project.toml                    single root environment, bounded by [compat]
├── activate.jl                     activates and instantiates it
├── theme.jl                        shared figure style and colourblind-safe palette
├── .JuliaFormatter.toml            formatting rules, SciMLStyle
├── .github/workflows/CI.yml        formatting check, and one figure rendered from scratch
├── CITATION.cff                    citation metadata, rendered by GitHub
├── LICENSE                         MIT
├── bsc/
│   ├── 01-oscillators-and-integrators/   year 1, 2017–2018, from C++
│   ├── 02-chaos-and-attractors/          year 1, 2017–2018, from C++ and Octave
│   ├── 03-compton-scattering/            year 1 semester 2, 2018, from C++
│   ├── 04-ising-model/                   summer internship, 2018, from C++
│   ├── 05-numerical-methods/             year not recorded, from Octave
│   ├── 06-geometrical-optics/            year not recorded, Octave and my goniometer readings
│   ├── 07-cosmology/                     July 2019, JINR Dubna (LIT)
│   ├── 08-numerical-methods-ii/          year 4, 2020–2021
│   ├── 09-bayesian-statistics/           year 4, 2020–2021
│   └── 10-remote-sensing-doas/           year 4, 2020–2021
└── msc/
    ├── 01-radiation-matter-interaction/  year 1, 2021–2022
    ├── 02-experimental-methods/          year 1, 2021–2022
    ├── 03-radionuclides/                 year 1, 2021–2022
    ├── 04-fission-observables/           year 2, 2022–2023
    ├── 05-nuclear-particle-physics/      year 2, 2022–2023
    └── 06-aircraft-hazard-screening/     late 2023, not a course assignment
```

The BSc is Engineering Physics, four years and 240 ECTS, so BSc year 4 is the
2020–2021 academic year and the MSc runs 2021–2023. Years come from the archive
folder names, which carry the year of study (`Old_2018`, `…_L_4`, `…_M_1`,
`…_M_2`); the two modules whose folders carry no marker say so rather than
guess.

Each directory has its own README with the figures, the results, and what was
wrong with the original.

## By subject

The layout above is by degree and then chronological, which is how the work
happened. Several subjects span that ordering, so they are indexed here too.

| Subject | Where |
|---|---|
| **Numerical integration of ODEs** | `bsc/01` (integrator order, symplectic vs explicit), `bsc/02` (applied to chaotic flows), `bsc/05` (Adams, RK4, pendulum), `bsc/08` (non-unique IVPs, method of lines) |
| **Cosmology** | `bsc/07` (Friedmann, the Hubble diagram, the acceleration discovery), `msc/05` (`hubble_parameter.jl`, H₀ from galaxy recession) |
| **Optics** | `bsc/06` (paraxial ray transfer, prism dispersion), `msc/02` (prism spectroscope calibration) — the two halves of the same Cauchy relation, one assumed and one measured |
| **Monte Carlo and statistics** | `bsc/03` (Compton chain, Klein–Nishina sampling), `bsc/04` (Ising, Metropolis), `bsc/09` (Bayesian updating, rejection sampling) |
| **Nuclear structure and reactions** | `msc/01` (radiation–matter), `msc/03` (masses, separation, pairing), `msc/04` (fission observables) |
| **Atmospheric remote sensing** | `bsc/10` (DOAS NO₂ columns from two instruments) |
| **Radiological assessment** | `msc/06` (external-hazard screening), and [AtmosphericDispersion.jl](https://github.com/PaulGoG/AtmosphericDispersion.jl) |

The numerical-methods pair, `bsc/05` and `bsc/08`, is one course taught twice:
the first in Octave, the second in Mathematica and Octave in the fourth year.
Only the second is dated — its archive folder is marked `L_4` — so the first is
placed before it by the course numbering alone. They are kept apart because what
they show is different: `05` is about the schemes themselves, `08` about what
happens when a problem is ill-posed.

## Running it

```bash
julia bsc/04-ising-model/ising_annealing.jl
```

Any script runs the same way, from any directory: each activates and
instantiates the root environment as its first statement, so the first run is
slow. Each writes its figures into its own `figures/` directory and prints its
numerical results to the terminal. `julia -i activate.jl` opens a REPL in the
same environment.

Formatting is checked against the committed `.JuliaFormatter.toml`:

```bash
julia -e 'include("activate.jl"); using JuliaFormatter; format(".")'
```

## What the rewrite found

Every script is now checked against something — a closed-form solution, a
conserved quantity, an evaluated nuclear datum, or a literature value. Doing
that turned up a substantial number of defects in the originals. The ones that
changed results:

| Where | Defect |
|---|---|
| `OscilatorArmonique/Amortizat/Fortat.cpp`, `AtractorI/helpers.cpp` | an `ifstream` opened on a file that an unflushed `ofstream` was still writing; a failed extraction zeroes its target since C++11, so the tail of every separation and extremum series was zero, and all of it for a run whose output fits one stream buffer. Compiled with its paths repaired, the attractor program does print non-zero distances, so the design condemned the tails, not necessarily every run |
| `Runge_Kutta_4_pendul_dampat_fortat.m` | the θ-stages used the ω-slopes, making the scheme **first-order, not fourth**. Measured order: 1.005 against 3.996 |
| `IsingFinal.cpp` | rows allocated `n+2` wide, columns only `n`, while indexing `a[i][n+1]` — a heap overflow; ghost cells never refreshed; every bond double-counted, so Metropolis ran at effective temperature **T/2**; cooling per flip rather than per sweep |
| `Fisiune_4.jl` | Madland–Nix prefactor written `(1/3*sqrt(E_F*T_MAX))`, which multiplies by √(E_f T_m) where the model divides by it; each fragment's spectrum then integrates to E_f T_m instead of one, reweighting the mass average by 4.8 % in ⟨E⟩ |
| `Fisiune_2.jl` | Y(N) summed inside the (A,Z) loop over a quantity depending only on N: **805 % total yield instead of 100 %** |
| `Frecventa_accident_aviatic.jl` | two crash-frequency estimators, neither of them an areal density, differing by one power of length — **neither was a frequency**, yet both were screened against per-year thresholds |
| `Hubble.jl` | `Suma_σ² =+ σ_z[i]^2` parses as an assignment, not `+=`, so every averaged uncertainty kept only its last term |
| `EfectulCompton.cpp` | scattering angle sampled uniformly in θ — neither Klein–Nishina nor isotropic, which is uniform in cos θ |
| `Dez_U_28.m` | the ²³⁸U half-life used as the mean lifetime |
| `Spatiul_FazelorEcDiff.cpp` | `EulerImplicit` advances both components from old values — it is explicit Euler, and the misnomer had spread into six output filenames |
| `Radionuclizi_4_1.jl` | Bateman time grid stepping in units of 450 Myr for a transient peaking at 2.4 years, collapsing the whole figure onto x = 0 |
| `Calibrare_Hg.jl` | 5789.66 Å listed where no Hg I line exists; the yellow doublet member is 5769.60 Å |
| `Breit_Wigner.jl` | the full width at half maximum located by scanning for samples within an absolute tolerance keyed to the array minimum, which returns whichever near-misses came first and last |
| `SWING_DOAS.jl` | no quality filter and no error bars: the largest SWING column of the day is a failed fit with RMS 0.035 and a 37 % retrieval error, and every SWING spectrum at 6° fails the residual cut |

And one that is a property of the problem rather than the code:
`ODE_RK4.jl` states the initial-value problem `y' = √|sin y|, y(0) = 0`, which is
**not well posed** — √|sin y| is not Lipschitz at any multiple of π, so
uniqueness fails. A perturbation of 10⁻¹² changes y(100) from 0 to 59.7, the
value of the solution that never waits, obtained by quadrature. The 2021 file
set `y[1] = 1`, contradicting its own comment.

## Some results

| | |
|---|---|
| Largest Lyapunov exponents, Benettin with blocked errors | Lorenz 0.906 ± 0.014 (literature 0.9056), Rössler 0.0719 ± 0.0016 (0.0714) |
| Observed orders of Euler and RK4 | 1.000 ± 0.001 and 3.994 ± 0.010 |
| Maximum binding energy per nucleon | 8.7946 MeV at **⁶²Ni**, over ⁵⁶Fe at 8.7904 |
| ⟨TKE⟩ in ²³⁵U(n_th,f) | 170.55 ± 0.05 MeV |
| S_n(²³⁶U) | 6.546 MeV |
| Prompt-neutron spectrum, Madland–Nix mass average | ⟨E⟩ = 2.103 MeV, Maxwellian-equivalent T = 1.402 MeV |
| Prompt-neutron pair multiplicity, yield-weighted | 2.538 (measured 2.40–2.46 across four datasets) |
| H₀ from three galaxies | 70.6 ± 12.2 km s⁻¹ Mpc⁻¹, the reference galaxy's size error being common to two of them |
| H₀ from H(z), flat ΛCDM | 73.3 ± 5.1 km s⁻¹ Mpc⁻¹, Ω_m = 0.267 ± 0.064 |
| Tessar objective focal length | 50.79 mm |
| Prism glass, two-term Cauchy | A = 1.5071 ± 0.0026, B = 0.00354 ± 0.00056 µm² |
| ¹²⁸I γ line, from a NaI(Tl) calibration | 441.1 ± 7.0 keV (accepted 442.9) |
| Aircraft crash frequency, 5 km from the route, R = 50 m | 1.09 × 10⁻⁸ yr⁻¹ (forward cosine kernel) |

## Branches

The [`legacy`](../../tree/legacy) branch holds the material exactly as it was
submitted — Romanian identifiers
and comments, the defects above intact — as the reference this branch is read
against. It is a single import: the commit histories of the two private
repositories this came from are not carried over.

The fission work in `msc/04` is the coursework ancestor of a deterministic
prompt-emission model that is still private and due for release; the same is
true of the EXFOR retrieval and fragment-temperature work it leans on.

## Attribution

`bsc/07-cosmology/` is joint work with
[Alexandru Crăciun](https://github.com/Craciun-Alexandru), presented together at
JINR Dubna in 2019. Its `sn_ia_distance_moduli.csv` is his reduction of the
NED-D compilation of redshift-independent distances (Steer et al. 2017,
doi:10.3847/1538-3881/153/1/37) and is used here with that attribution.

## How to cite

`CITATION.cff` carries the metadata; GitHub renders it as a citation block.

```bibtex
@software{Gogita_UniBuc_Computational_Physics_Coursework,
  author = {Gogîță, Paul-Adrian},
  title  = {UniBuc Computational Physics Coursework},
  year   = {2026},
  url    = {https://github.com/PaulGoG/UniBuc-Computational-Physics-Coursework}
}
```

## Licence

MIT, see `LICENSE`.
