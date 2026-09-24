# Code_Archive

C, C++ and Octave coursework of the BSc years, 2018–2021, as submitted. Comments
and identifiers are in Romanian; several programs read files they have not
finished writing, or hold absolute Windows paths, and none of that is repaired
here. The names are English since 2026, and the mapping to the submitted names
is in the [root README](../README.md#file-names).

## `Numerical_Methods/`

Octave, the first numerical-methods course.

| File | Contents |
|---|---|
| `Adams_Predictor_Corrector.m` | Four-step Adams predictor–corrector for a first-order ODE, started with Runge–Kutta |
| `Simple_pendulum.m` | The simple pendulum by the Euler method |
| `Damped_pendulum.m` | The damped pendulum by Euler–Cromer |
| `Damped_forced_pendulum.m` | The damped, driven pendulum by Euler–Cromer |
| `RK4_damped_forced_pendulum.m` | The damped, driven pendulum by a fourth-order Runge–Kutta scheme whose angle stages take the wrong slopes, so that it is first order in practice |
| `ODE_system_RK4.m` | Runge–Kutta 4 for a system of ODEs |
| `U238_decay.m` | Decay of ²³⁸U, `dN/dt = −N/τ`, by Euler, with the half-life used where the mean lifetime belongs |

## `Old_2018/C_C++/`

C++ of the first year, written for Code::Blocks on Windows. The oscillator and
attractor programs write their trajectories to files and read them back within
the same run, before the streams are flushed.

| File | Contents |
|---|---|
| `Harmonic_oscillator.cpp`, `Damped_oscillator.cpp`, `Forced_oscillator.cpp` | Two nearby trajectories of the oscillator, their separation in time and the intervals between extrema |
| `Euler_ODE.cpp` | Euler integration of a first-order ODE |
| `Runge_Kutta_trial.cpp` | A first fourth-order Runge–Kutta integrator |
| `ODE_exam_problem.cpp` | An examination problem: an ODE integrated with its energies written out |
| `Phase_space_ODE.cpp` | Phase-space portraits under explicit Euler and a scheme labelled implicit that is explicit as well |
| `Compton_effect.cpp` | Monte Carlo chain of Compton scatterings down to an absorption threshold, the angle sampled uniformly in θ |
| `Ising_model.cpp` | Two-dimensional Ising model by Metropolis with cooling |
| `Attractors/Attractors.cpp` | Driver for the Lorenz and Rössler systems: trajectories, separation of nearby trajectories, intervals between extrema |
| `Attractors/Integrators.h`, `Integrators.cpp` | The two flows and their Runge–Kutta integration, written to file |
| `Attractors/helpers.h`, `helpers.cpp` | Separation and extremum intervals, read back from the trajectory files |
| `Attractors/Thread_demo.c` | A POSIX threads exercise unrelated to the attractors, kept where it was filed |

## `Old_2018/Octave/`

Post-processing of the C++ output above.

| File | Contents |
|---|---|
| `Attractor_full.m` | The attractor trajectories from the C++ output files |
| `Attractor_distribution.m` | Histogram of the intervals between extrema |
| `Interval_histograms.m` | Six interval histograms on one figure |
| `Random.m`, `Test.m` | Scratch files: a matrix syntax check and a reading of the separation files |

## `Optics/`

Ray-transfer matrices for two lens systems, distances in millimetres.

| File | Contents |
|---|---|
| `Erfle.m` | Erfle eyepiece: cardinal points and focal length |
| `Tessar.m` | Tessar objective: the same |
