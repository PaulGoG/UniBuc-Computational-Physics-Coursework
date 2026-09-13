# Oscillators and integrators — BSc year 1 (2017–2018)

First-year and second-year work, originally written in C++ and ported to Julia.
Three programs on the harmonic oscillator and on the accuracy of the elementary
integration schemes used to solve it.

## `integrator_order_comparison.jl`

Forward Euler against classical RK4 on `y' = 3e^{-x} - 0.4y`, whose exact
solution `y = (y₀+5)e^{-0.4x} - 5e^{-x}` is available in closed form. Sweeping
the step size over sixteen powers of two recovers the observed orders of
convergence: **1.002** for Euler and **4.121** for RK4.

`DiffEqEuler.cpp`'s autonomous nonlinear problem `y' = y + sin(0.4y)` is
integrated alongside it against a finely resolved RK4 reference, since it has no
closed form: orders **0.979** and **3.955**.

![Order of convergence](figures/integrator_order_comparison.png)

Ported from `RKtrial.cpp` and `DiffEqEuler.cpp`. Both Butcher tableaux in the
originals were correct, but neither program compared against the exact solution
or swept the step size, so the O(h) against O(h⁴) contrast they were written to
show was never actually demonstrated. `DiffEqEuler.cpp` additionally used one
variable as both the initial abscissa and the length of the integration
interval, and never advanced the abscissa at all.

## `symplectic_vs_explicit_euler.jl`

Explicit Euler against the two symplectic (Euler–Cromer) variants on
`ẍ = -Ω²x`. The three schemes agree to first order in Δt yet behave completely
differently over long times: explicit Euler spirals outward without bound while
both symplectic variants keep their orbits closed and their energy oscillating
about the initial value.

![Phase portraits and energy drift](figures/symplectic_vs_explicit_euler.png)

The amplification factor of explicit Euler on this system is √(1 + Ω²Δt²) per
step, so the energy grows by (1 + Ω²Δt²)^N. Over 800 steps at Δt = 0.05 that
predicts a factor 7.371 — the measured growth is 7.371.

Ported from `Spatiul_FazelorEcDiff.cpp` and `PbEcDiffExamen.cpp`. The original
called its first scheme `EulerImplicit`, but it advances both components from
the old values, which is explicit Euler; the misnomer propagated into six output
filenames, so everything labelled "implicit" in that project was in fact about
the explicit method. It also worked in single precision on a study whose entire
subject is accumulated integration error, and drew its initial conditions as
integers on a 101×101 lattice with modulo bias.

## `oscillator_trajectory_separation.jl`

Two trajectories starting 10⁻⁷ apart in phase space. The three parameter sets of
the 2018 programs are integrated as written, alongside two non-degenerate
replacements.

| Case | Fitted decay rate | Input δ |
|---|---|---|
| Undamped (2018) | −0.00000 | 0 |
| Critically damped, δ = ω₀ = 1 (2018) | 0.83410 | 1 |
| Driven, δ = 0 (2018) | −0.00000 | 0 |
| Under-damped, δ = 0.15 | 0.14996 | 0.15 |
| Driven, δ = 0.3, ω = 1.6 | 0.29999 | 0.3 |

**The 2018 undamped and 2018 driven cases give identical separation**, to every
digit. The driving is the same along both trajectories and cancels in their
difference, so `OscilatorFortat.cpp` could not have learned anything its
undamped sibling did not already show — even if its file handles had worked.
The critically damped case fits a rate of 0.834 rather than 1, because at
δ = ω₀ the solution carries a `(A + Bt)e^{−δt}` factor and a pure exponential fit
splits the difference.

![Trajectory separation](figures/oscillator_trajectory_separation.png)

Because the system is linear, the difference of two solutions satisfies the
*homogeneous* equation — the driving is identical along both trajectories and
cancels — so the separation decays with envelope `e^{-δt}` and never grows. For
δ = 0 the difference vector rotates rigidly and |Δ| is exactly constant; for
δ > 0 it oscillates inside the envelope by a relative amount of order δ/ω₀.
Fitting the decay rate recovers the input damping to five digits: 0.00000,
0.14996 and 0.29999 against 0, 0.15 and 0.30.

That is the point of the exercise, and it is what makes the return-time
statistics the 2018 programs were built around meaningless here. They only carry
information once the dynamics are chaotic, which is `02-chaos-and-attractors`.

Ported from `OscilatorArmonique.cpp`, `OscilatorAmortizat.cpp` and
`OscilatorFortat.cpp`. All three opened an `ifstream` on the same file as an
unflushed `ofstream`, so every extraction failed and — since C++11 — silently
set its target to zero. **The separation and extremum-interval results those
three programs printed were identically zero; no real number was ever
produced.** The analysis is done in memory here, so the failure cannot recur.
The parameter choices were also degenerate: the damped case used δ = ω₀ = 1,
exactly critical damping, which does not oscillate at all, and the driven case
used δ = 0, so it never reaches a steady state.
