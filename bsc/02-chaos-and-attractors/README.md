# Chaos and attractors — BSc year 1 (2017–2018)

The Lorenz and Rössler systems integrated with RK4, their largest Lyapunov
exponent from the divergence of neighbouring trajectories, and the
extremum-interval statistic of the separation that the 2018 program computed.
Ported from `Integrators.cpp`, `Attractors.cpp` and `helpers.cpp` in
`Code_Archive/Old_2018/C_C++/Attractors/` on the `legacy` branch, with the Octave
post-processing of `Code_Archive/Old_2018/Octave/`. `attractors_core.jl` holds
the vector fields, parameters and the RK4 step both scripts share.

## `lorenz_rossler_attractors.jl`

![Lorenz and Rössler attractors](figures/lorenz_rossler_attractors.png)

![Attractors traced out](figures/lorenz_rossler_attractors.gif)

Projections of the two attractors after a transient, 2.5 × 10⁵ points each at
h = 0.002. The non-trivial Lorenz fixed points (±√(β(ρ−1)), ±√(β(ρ−1)), ρ−1)
lie at ±8.4853 for β = 8/3; the 2018 code wrote β = 2.66, which puts them at
±8.4747. The animation traces the trajectory itself: on the Lorenz attractor the
state circles one lobe an unpredictable number of times before crossing to the
other; on the Rössler attractor it spirals outward in a near-plane and is folded
back.

The RK4 stage coupling of the original is correct. Besides β, its Rössler branch
is unreachable (the selector is assigned and never read), its output paths are
absolute Windows paths, and its parameters are `#define` macros, so `a`, `b`
and `c` textually replace any identifier of those names in the translation unit.

## `lyapunov_and_return_times.jl`

![Lyapunov divergence and return times](figures/lyapunov_and_return_times.png)

**Lyapunov exponents** by the Benettin algorithm: the perturbation is rescaled
to d₀ = 10⁻⁷ every τ = 0.1 time units and λ₁ is the mean of ln(d/d₀)/τ. The
increments are correlated, so the standard error is taken from blocked averages
(Flyvbjerg and Petersen, J. Chem. Phys. **91**, 461 (1989),
doi:10.1063/1.457480), at the longest blocks that still leave 32 of them: 48
blocks of 410 time units for Lorenz and 61 of 1640 for Rössler. The increments
are anti-correlated over a turn of the attractor, so the blocked estimate first
rises, to 0.014 for Lorenz at blocks of 0.4 time units, and then falls towards
its plateau. For Lorenz it is still falling at the longest blocks (0.0014,
0.0014 and 0.0012 at 102, 205 and 410 time units), so the quoted error is an
upper estimate: eight independent runs of the same length scatter by 0.0006.
The script asserts that the estimate is not still rising at the last level,
and agreement within three standard errors with the values of Sprott, *Chaos
and Time-Series Analysis*, Oxford University Press (2003).

| System | Run length | λ₁ | Reference |
|---|---|---|---|
| Lorenz | 2 × 10⁴ | 0.9064 ± 0.0012 | 0.9056 |
| Rössler | 10⁵ | 0.0719 ± 0.0005 | 0.0714 |

**The un-renormalised separation**, which is what the 2018 program followed,
grows from 10⁻⁷ to order one in 16 time units for Lorenz and in 219 for
Rössler, whose exponent is thirteen times smaller, and saturates at the size of
the attractor. The left panel runs each system to just past the time its
separation reaches one: Lorenz has saturated by then, Rössler is at 2.2 and
still growing.

**Extremum intervals.** `helpers.cpp` drops 30 time units at each end of the
separation series, accepts an interior extremum when it differs from the last
accepted one by more than 10⁻⁴ in relative terms, and records the min → max
and max → min intervals in two columns; the Octave scripts histogrammed them.
On the run that program made, 200 000 steps of 0.001 over 200 time units from
(1, 1, 1) for Lorenz only (its Rössler branch is unreachable), the definition
gives 269 rising intervals with mean 0.237 ± 0.007 and 270 falling ones with
mean 0.282 ± 0.008. This statistic describes the oscillation of the separation
on the saturated attractor, not the exponential growth; return-time statistics
of a linear oscillator, by contrast, carry no information
(`01-oscillators-and-integrators`).

The program reads its own output back through an `ifstream` opened on the file
its `ofstream` is still writing, which the C++ standard does not guarantee.
Compiled with local paths and run with the inputs above, it does produce
140 001 distances, none of them zero, and 301 interval rows with mean 0.28, so
the design alone does not produce zeros; the output of the 2018 runs is unknown.
`Test.m` and `Attractor_distribution.m` are not ported. The first reads
`DistantaLorentz.txt`, `TimpiUCLorentz.txt` and `DateExtreme.txt`, which no
program in the archive writes, with histogram edges `0:0.05:0.8` that match the
intervals above; the second seeds a six-parameter fit from `rand` and does not
check its convergence flag.
