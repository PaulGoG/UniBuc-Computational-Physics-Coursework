# Numerical methods — BSc, year not recorded

Octave coursework on multistep methods, Runge–Kutta, and the pendulum. Ported to
Julia. Every script here is checked against a closed-form solution, which none of
the originals did.

## `adams_predictor_corrector.jl`

![Adams PECE](figures/adams_predictor_corrector.png)

Four-step Adams–Bashforth predictor with three-step Adams–Moulton corrector, in
PECE form, on `y' = -y + 2 sin x` with exact solution `y = sin x - cos x + e^{-x}`.
Starting values from RK4. Observed order **4.44**.

The Adams coefficients in the original were correct. What it threw away was the
free local-error estimate: the difference between predictor and corrector is
proportional to the local truncation error, and is the entire reason to run a
predictor–corrector pair rather than a corrector alone. Restored here, it tracks
the true global error across the whole interval.

## `pendulum_integrators.jl`

![Pendulum integrators](figures/pendulum_integrators.png)

![Pendulum integrators compared](figures/pendulum_integrators.gif)

The order study measures the defect; this shows it. Both pendulums are given the
same coarse step of 0.05 s, and the one whose θ-stages were fed the ω-slopes
drifts visibly in phase within a few swings. At the 0.002 s step used elsewhere
the two are indistinguishable, which is exactly why the error went unnoticed.

**The original "RK4" was first-order.** For the coupled system `θ̇ = ω`,
`ω̇ = F(θ, ω, t)`, the θ-stage increments must use the θ-slopes. The original
wrote

```octave
k2 = frhs(omega(step) + dt*k1/2, theta(step) + dt*k1/2, time(step) + dt/2);
pk1 = omega(step);  pk2 = omega(step) + dt*pk1/2;  ...
```

so `k1`, which is dω/dt, advanced *both* components, and the θ-stages propagated
ω as though dω/dt = ω. The θ update expands to `dt·ω + dt²·ω/2` where RK4 needs
`dt·ω + dt²·F/2`. Both versions are integrated here and their orders measured
against a fine reference:

| Stage coupling | Observed order |
|---|---|
| Correct RK4 | **3.83** |
| As written originally | **1.01** |

The middle panel contrasts the full pendulum `θ̈ = -(g/L) sin θ` with the
small-angle approximation from θ₀ = 2.5 rad — they diverge by up to 4.9 rad,
since the true period lengthens with amplitude. It uses g/L = 9.8, the value of
`Pendul_simplu.m` (L = 1 m, g = 9.8); the driven files set L = 9.8 m, hence the
g/L = 1 of the Poincaré panel. `Pendul_simplu.m` is headed
"pendulul matematic" but integrates the linearised equation, not the pendulum.

The right panel is the Poincaré section of the driven damped pendulum at
Giordano's chaotic parameters (q = 0.5, F_D = 1.2, Ω_D = 2/3), sampled once per
drive period after discarding the transient. The original files plotted the wrapped
angle against time and never took a section, which is what this parameter set
exists for.

## `uranium238_decay.jl`

![U-238 decay](figures/uranium238_decay.png)

`dN/dt = -N/τ` by forward Euler against `N₀e^{-t/τ}`.

The original wrote `tau = 4.4e9; % timpul mediu de viata U238`. That is the
**half-life** of ²³⁸U, not the mean lifetime: τ = T½/ln 2 = 6.446 Gyr. Using one
for the other decays the sample far too fast — at t = τ the surviving fraction
is 0.3674 rather than 0.2357, a 36 % error. Both curves are shown. On the
original's grid (Δt = 10⁷ yr, 1000 points) the Euler solution recovers a
half-life of 4.470 Gyr against the true 4.468.

It is also not a decay chain despite the filename. ²³⁸U does decay through
fourteen members to ²⁰⁶Pb, but none of that is modelled; the single-isotope law
is what the file solves.

## `ode_system_rk4.jl`

![ODE system](figures/ode_system_rk4.png)

RK4 on a three-component linear system with a closed-form solution. Components 1
and 2 form a driven oscillator `y₁'' + y₁ = 1 - 2eᵗ`; component 3 is a
quadrature of it. Observed order **3.98**, maximum error 6.4 × 10⁻⁴ over
t ∈ [1, 8], where the solution itself reaches −3000.

The tableau was correct. The problems were structural: an orientation check
written as `m = size(alpha); if m==1 ...` where `size` returns `[1 3]`, so the
condition is `[1 0]` and — since `if` requires every element non-zero — the
transpose never fired; the code worked only because assigning into a column
silently reoriented. The routine also called its right-hand side by name instead
of taking it as an argument, dumped a 101×4 matrix to the console twice through
missing semicolons, and demultiplexed its state by column-major linear indexing.
