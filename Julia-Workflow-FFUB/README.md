# Julia-Workflow-FFUB

Julia coursework of the BSc fourth year and the MSc, 2019–2023, as submitted at
the Faculty of Physics, University of Bucharest (FFUB). Comments and identifiers
are in Romanian, and the scripts read their data files by relative path from
their own directory. The names are English since 2026; the mapping to the
submitted names is in the [root README](../README.md#file-names). Directory
suffixes record the year of study: `_L_n` is BSc year *n*, `_M_n` MSc year *n*.

## `Dubna_2019_Cosmology/`

Joint work with Alexandru Crăciun at JINR Dubna, July 2019.

| File | Contents |
|---|---|
| `Friedmann_Eq.jl` | The Friedmann equation for `H(z)`, analytic and numerical, fitted to the `H(z)` data |
| `Friedmann_Magnitude.jl` | Distance modulus against redshift, the cosmological parameters fitted with Optim; the data path is an absolute Windows one |
| `Data.csv`, `goodData.csv` | The supernova distance moduli, full and cleaned |

## `Numerical_Methods_II_Exam_L_4/`

| File | Contents |
|---|---|
| `ODE_RK4.jl` | `y' = √|sin y|`, `y(0) = 0`, by RK4: an initial-value problem that is not well posed |
| `ODE_system_RK4.jl` | A second-order system, `y' = z`, `z' = −y sin y`, by RK4 |
| `SolitonicEq_MOL.jl` | The nonlinear Schrödinger equation by the method of lines |

## `Nuclear_Particle_Physics_M_2/`

| File | Contents |
|---|---|
| `Breit_Wigner.jl` | Two interfering Breit–Wigner resonances and their widths |
| `Hubble.jl` | `H₀` from the recession of three galaxies with propagated uncertainties |

## `Fission_M_2/`

Neutron-induced fission of ²³⁵U, five assignments on one yield matrix.

| File | Contents |
|---|---|
| `Fission_1.jl` | Energy released per fragmentation, `Q(A, Z)`, from the mass table over three charges about `Z_p(A)` |
| `Fission_2.jl` | The usual distributions of prompt fission: `Y(A)`, `Y(Z)`, `Y(N)`, `Y(TKE)` |
| `Fission_3.jl` | Partition of the total excitation energy between the light and heavy fragments by a scission-point model, against measured multiplicities |
| `Fission_4.jl` | Prompt-neutron spectrum in the global Los Alamos (Madland–Nix) treatment, one most probable fragmentation |
| `Fission_5.jl` | Maxwellian fits to the measured prompt-neutron spectra |
| `Dump.jl` | Functions shared by the five scripts |
| `Data_files/Yield/U5YAZTKE.STR` | The yield matrix `Y(A, Z, TKE)` of ²³⁵U(n_th,f) |
| `Data_files/Mass_defects/AUDI2021.csv` | The AME2020 mass table |
| `Data_files/Auxiliary_parametrisations/B2MOLLER.ANA`, `SZSN.GC` | Ground-state deformations of Möller and the Gilbert–Cameron shell corrections |
| `Data_files/Experimental/Neutron_multiplicity/` | Measured `ν(A)` of Göök, Maslin, Nishio and Vorobyev |
| `Data_files/Experimental/Neutron_spectrum/` | Measured prompt-neutron spectra of Göök and Vorobyev |

## `Radiation_Matter_Interaction_M_1/`

| File | Contents |
|---|---|
| `Bethe_Bloch.jl` | Bethe–Bloch stopping power of a charged particle in silicon |
| `Alpha_attenuation.jl` | Residual energy of α particles through Mylar absorbers, a linear fit |
| `NaITl_activation_fit.jl` | Energy calibration of a NaI(Tl) detector and the half-lives of activation products |
| `Neutron_reactions.jl` | Neutron activation counts against time, the Schweidler law linearised |
| `Natural_Cd_cross_section.jl` | The neutron cross-section of natural cadmium from its isotopes |

## `Radionuclides_M_1/`

| File | Contents |
|---|---|
| `Radionuclides_1.jl` | Mean binding energy per nucleon from the mass tables |
| `Radionuclides_2_1.jl` | Separation energies of nucleons and clusters |
| `Radionuclides_2_2.jl` | The separation energy of one nucleus from another, read from the console |
| `Radionuclides_3.jl` | Pairing energies |
| `Radionuclides_4_1.jl` | Radioactive decay series by the Bateman equations |
| `Radionuclides_4_2.jl` | Activity of an ²⁷Al sample under a pulsed neutron flux |
| `Radionuclides_5.jl` | Shell corrections against the Möller table |
| `AUDI95.csv`, `AUDI2021.csv`, `MOLLER.csv` | The AME1995 and AME2020 mass tables and Möller's microscopic corrections |

## `Single_Files/`

| File | Contents |
|---|---|
| `Hg_calibration.jl` | Calibration of a prism spectroscope on the lines of a mercury lamp |
| `Aircraft_crash_frequency.jl` | Aircraft crash frequency at a site, two estimators |
| `MaximizeLinearSystemEq.jl` | A linear programme through JuMP and GLPK |
| `MonteCarloDistribution.jl` | Rejection sampling from a bi-Gaussian |

## `Statistics_Simulation_L_4/`

| File | Contents |
|---|---|
| `Statistics.jl` | Bayesian updating of a Beta posterior over repeated coin tosses, animated |
| `CoinTossAveraging.gif` | The animation, re-encoded at every fourth frame for the import |

## `Remote_Sensing_L_4/`

| File | Contents |
|---|---|
| `SWING_DOAS.jl` | NO₂ slant columns from a MAX-DOAS and the SWING instrument, 18 March 2019 |
| `MAXDOAS.csv`, `SWING.csv` | The two instruments' retrievals |
