# UniBuc Computational Physics Coursework — `legacy`

Snapshot of the code written during my studies at the Faculty of Physics,
University of Bucharest. It is kept as a reference point: `main` reorganises,
translates and repairs this material, and every change made there reads as a
diff against this branch.

Nothing here has been corrected. Comments, identifiers and printed output are in
Romanian, the numerical methods stand as they were understood at the time, and
the programming is script-style throughout. Several programs contain mistakes
that were never caught; they are left in place deliberately. The one change
since submission is to the names: files and directories carry English names,
and the `include` lines and data-path strings that name them were changed to
match. The mapping is [at the end](#file-names). Each of the two archives has a
README listing its files.

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
| `Old_2018/C_C++/` | Harmonic, damped and driven oscillators, Euler and RK integrators, phase-space portraits, Compton scattering, Ising model, strange attractors (`Attractors/`) |
| `Old_2018/Octave/` | Attractor distributions, histogramming, pseudorandom sampling |
| `Optics/` | Erfle eyepiece and Tessar objective ray traces |

### `Julia-Workflow-FFUB/`

Directory suffixes record the course stage: `_L_n` is year *n* of the BSc,
`_M_n` is year *n* of the MSc.

| Path | Contents |
|---|---|
| `Dubna_2019_Cosmology/` | Friedmann equations, distance modulus against a supernova catalogue |
| `Numerical_Methods_II_Exam_L_4/` | RK4 for scalar and vector ODEs, method of lines for a solitonic PDE |
| `Nuclear_Particle_Physics_M_2/` | Breit–Wigner resonance, Hubble diagram |
| `Fission_M_2/` | Neutron-induced fission of <sup>235</sup>U: yields, *Q*-values, excitation energies, neutron multiplicities and spectra |
| `Radiation_Matter_Interaction_M_1/` | Bethe–Bloch stopping power, α attenuation, NaI(Tl) activation fit, neutron reactions, natural-Cd cross-section |
| `Radionuclides_M_1/` | Binding, separation and pairing energies, decay chains and shell corrections from the AME and Möller tables |
| `Single_Files/` | Hg lamp calibration, aircraft-accident frequency, linear-system maximisation, Monte Carlo sampling |
| `Statistics_Simulation_L_4/` | Bayesian updating of a Beta posterior over repeated coin tosses |
| `Remote_Sensing_L_4/` | MAX-DOAS and SWING differential optical absorption spectroscopy |

## Provenance and attribution

The material was recovered from two personal repositories, now private and
archived. Their commit histories are not carried over: this branch is a single
import of the final state of both, so that the public repository starts without
the generated artefacts and bulk data files that had accumulated in them.

`Dubna_2019_Cosmology/` is joint work with
[Alexandru Crăciun](https://github.com/Craciun-Alexandru), presented together at
JINR Dubna in 2019.

`Statistics_Simulation_L_4/CoinTossAveraging.gif` is the only file altered in
this import: the original 1001-frame animation was 14.6 MB, and it was
re-encoded to 640×450 at every fourth frame to keep the repository a reasonable
size. The script that produces it is unchanged and will regenerate the original.

## File names

The names as submitted, so that `main`'s account of what was wrong reads
against the files it names, and so that the archive folders can be found. The
year-of-study suffixes (`_L_n`, `_M_n`) and `Old_2018` are kept as they were.
Beyond the names, the only edits are the two `#include` lines of the attractor
program and the data paths inside `Fission_M_2/`, which point at the renamed
directories.

| Now | As submitted |
|---|---|
| `Numerical_Methods/Adams_Predictor_Corrector.m` | `Adams_Predictor_Corector.m` |
| `Numerical_Methods/U238_decay.m` | `Dez_U_28.m` |
| `Numerical_Methods/Damped_pendulum.m` | `Pendul_dampat.m` |
| `Numerical_Methods/Damped_forced_pendulum.m` | `Pendul_dampat_fortat.m` |
| `Numerical_Methods/Simple_pendulum.m` | `Pendul_simplu.m` |
| `Numerical_Methods/RK4_damped_forced_pendulum.m` | `Runge_Kutta_4_pendul_dampat_fortat.m` |
| `Numerical_Methods/ODE_system_RK4.m` | `Sys_ODE_RK4.m` |
| `Old_2018/C_C++/Attractors/` | `AtractorI/` |
| `Attractors/Attractors.cpp`, `Integrators.cpp`, `Integrators.h` | `Atractori.cpp`, `Integratori.cpp`, `Integratori.h` |
| `Attractors/Thread_demo.c` | `badica.c` |
| `Old_2018/C_C++/Euler_ODE.cpp` | `DiffEqEuler.cpp` |
| `Old_2018/C_C++/Compton_effect.cpp` | `EfectulCompton.cpp` |
| `Old_2018/C_C++/Ising_model.cpp` | `IsingFinal.cpp` |
| `Old_2018/C_C++/Damped_oscillator.cpp` | `OscilatorAmortizat.cpp` |
| `Old_2018/C_C++/Harmonic_oscillator.cpp` | `OscilatorArmonique.cpp` |
| `Old_2018/C_C++/Forced_oscillator.cpp` | `OscilatorFortat.cpp` |
| `Old_2018/C_C++/ODE_exam_problem.cpp` | `PbEcDiffExamen.cpp` |
| `Old_2018/C_C++/Runge_Kutta_trial.cpp` | `RKtrial.cpp` |
| `Old_2018/C_C++/Phase_space_ODE.cpp` | `Spatiul_FazelorEcDiff.cpp` |
| `Old_2018/Octave/Attractor_distribution.m` | `AtractorDistributie.m` |
| `Old_2018/Octave/Attractor_full.m` | `AtractorFull.m` |
| `Old_2018/Octave/Interval_histograms.m` | `HistogrameAmbalate.m` |
| `Optics/Erfle.m` | `Erlfe.m` |
| `Numerical_Methods_II_Exam_L_4/` | `Examen_PDF_MN_II_L_4/` |
| `Nuclear_Particle_Physics_M_2/` | `FPECA_M_2/` |
| `Fission_M_2/`, with `Fission_1.jl` … `Fission_5.jl` | `Fisiune_M_2/`, `Fisiune_1.jl` … `Fisiune_5.jl` |
| `Fission_M_2/Data_files/Experimental/Neutron_multiplicity/`, `Neutron_spectrum/` | `Date_experimentale/Multiplicitate_n/`, `Spectru_n/` |
| `Fission_M_2/Data_files/Mass_defects/`, `Auxiliary_parametrisations/` | `Defecte_masa/`, `Parametrizari_auxiliare/` |
| `Radiation_Matter_Interaction_M_1/` | `IRM_M_1/` |
| `Radiation_Matter_Interaction_M_1/Alpha_attenuation.jl` | `Atenuare_alpha.jl` |
| `Radiation_Matter_Interaction_M_1/Bethe_Bloch.jl` | `Calcul_Bethe_Bloch.jl` |
| `Radiation_Matter_Interaction_M_1/NaITl_activation_fit.jl` | `FitActivareNaITl.jl` |
| `Radiation_Matter_Interaction_M_1/Neutron_reactions.jl` | `ReactiiNeutronice.jl` |
| `Radiation_Matter_Interaction_M_1/Natural_Cd_cross_section.jl` | `SigmaCdNatural.jl` |
| `Radionuclides_M_1/`, with `Radionuclides_1.jl` … `Radionuclides_5.jl` | `Radionuclizi_M_1/`, `Radionuclizi_1.jl` … `Radionuclizi_5.jl` |
| `Single_Files/Hg_calibration.jl` | `Calibrare_Hg.jl` |
| `Single_Files/Aircraft_crash_frequency.jl` | `Frecventa_accident_aviatic.jl` |
| `Statistics_Simulation_L_4/Statistics.jl` | `Statistics_symul_L_4/Statis.jl` |
| `Remote_Sensing_L_4/` | `Teledetectie_L_4/` |
