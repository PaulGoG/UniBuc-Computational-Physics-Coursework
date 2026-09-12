# UniBuc Computational Physics Coursework

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
├── Project.toml, Manifest.toml     single root environment, resolved versions pinned
├── activate.jl                     activates and instantiates it
├── theme.jl                        shared figure style and colourblind-safe palette
├── bsc/
│   ├── 01-oscillators-and-integrators/   2018, from C++
│   ├── 02-chaos-and-attractors/          2018, from C++ and Octave
│   ├── 03-ising-and-compton/             2018, from C++
│   ├── 04-numerical-methods/             from Octave
│   ├── 05-geometrical-optics/            from Octave
│   ├── 06-cosmology/                     2019, JINR Dubna
│   ├── 07-numerical-methods-ii/          2021
│   ├── 08-bayesian-statistics/           2021
│   └── 09-remote-sensing-doas/           2021
└── msc/
    ├── 01-radiation-matter-interaction/
    ├── 02-radionuclides/
    ├── 03-fission-observables/
    ├── 04-nuclear-particle-physics/
    └── 05-aircraft-hazard-screening/
```

Each directory has its own README with the figures, the results, and what was
wrong with the original.

## Running it

```bash
julia activate.jl                                    # instantiate the environment
julia --project=. bsc/03-ising-and-compton/ising_annealing.jl
```

Any script runs the same way. Each writes its figures into its own `figures/`
directory and prints its numerical results to the terminal.

## What the rewrite found

Every script is now checked against something — a closed-form solution, a
conserved quantity, an evaluated nuclear datum, or a literature value. Doing
that turned up a substantial number of defects in the originals. The ones that
changed results:

| Where | Defect |
|---|---|
| `OscilatorArmonique/Amortizat/Fortat.cpp`, `AtractorI/helpers.cpp` | an `ifstream` opened on the same file as an unflushed `ofstream`, so every extraction failed and — since C++11 — set its target to zero. **Every separation and extremum result those four programs printed was identically zero** |
| `Runge_Kutta_4_pendul_dampat_fortat.m` | the θ-stages used the ω-slopes, making the scheme **first-order, not fourth**. Measured: 1.01 against 3.83 |
| `IsingFinal.cpp` | rows allocated `n+2` wide, columns only `n`, while indexing `a[i][n+1]` — a heap overflow; ghost cells never refreshed; every bond double-counted, so Metropolis ran at effective temperature **T/2**; cooling per flip rather than per sweep |
| `Fisiune_4.jl` | Madland–Nix prefactor written `(1/3)*sqrt(E_f*T_m)` where the model needs `1/(3*sqrt(E_f*T_m))` — multiplying where it must divide |
| `Fisiune_2.jl` | Y(N) summed inside the (A,Z) loop over a quantity depending only on N: **805 % total yield instead of 100 %** |
| `Frecventa_accident_aviatic.jl` | two crash-frequency estimators, neither of them an areal density, differing by one power of length — **neither was a frequency**, yet both were screened against per-year thresholds |
| `Hubble.jl` | `Suma_σ² =+ σ_z[i]^2` parses as an assignment, not `+=`, so every averaged uncertainty kept only its last term |
| `EfectulCompton.cpp` | scattering angle sampled uniformly in θ — neither Klein–Nishina nor isotropic, which is uniform in cos θ |
| `Dez_U_28.m` | the ²³⁸U half-life used as the mean lifetime |
| `Spatiul_FazelorEcDiff.cpp` | `EulerImplicit` advances both components from old values — it is explicit Euler, and the misnomer had spread into six output filenames |
| `Radionuclizi_4_1.jl` | Bateman time grid stepping in units of 450 Myr for a transient peaking at 2.4 years, collapsing the whole figure onto x = 0 |
| `Calibrare_Hg.jl` | 5789.66 Å listed where no Hg I line exists; the yellow doublet member is 5769.60 Å |

And one that is a property of the problem rather than the code:
`ODE_RK4.jl` states the initial-value problem `y' = √|sin y|, y(0) = 0`, which is
**not well posed** — √|sin y| is not Lipschitz at the origin, so uniqueness
fails. A perturbation of 10⁻¹² changes y(100) from 0 to 59.7. The 2018 file set
`y[1] = 1`, contradicting its own comment.

## Some results

| | |
|---|---|
| Lorenz largest Lyapunov exponent | 0.9018 (literature 0.9056) |
| Maximum binding energy per nucleon | 8.7946 MeV at **⁶²Ni**, over ⁵⁶Fe at 8.7904 |
| ⟨TKE⟩ in ²³⁵U(n_th,f) | 170.55 MeV |
| S_n(²³⁶U) | 6.546 MeV |
| Prompt-neutron spectrum temperature | 1.322 and 1.336 MeV (lab), 0.796 and 0.839 (CM) |
| H₀ from three galaxies | 73.1 ± 11.4 km s⁻¹ Mpc⁻¹ |
| H₀ from H(z), flat ΛCDM | 73.3 ± 5.1 km s⁻¹ Mpc⁻¹, Ω_m = 0.267 ± 0.064 |
| Tessar objective focal length | 50.79 mm |
| ¹²⁸I γ line, from a NaI(Tl) calibration | 441.1 ± 12.1 keV (accepted 442.9) |

## Branches

`legacy` holds the material exactly as it was submitted — Romanian identifiers
and comments, the defects above intact — as the reference this branch is read
against. It is a single import: the commit histories of the two archives this
came from are not carried over.

## Attribution

`bsc/06-cosmology/` is joint work with
[Alexandru Crăciun](https://github.com/Craciun-Alexandru), presented together at
JINR Dubna in 2019.

The BSc thesis and its code live in their own repository. The MSc thesis is not
here.

## Licence

MIT, see `LICENSE`.
