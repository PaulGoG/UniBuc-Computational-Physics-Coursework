# Ising model — summer internship, 2018

## `ising_annealing.jl`

![Ising annealing](figures/ising_annealing.png)

Two-dimensional square-lattice Ising model, `H = -J Σ⟨ij⟩ sᵢsⱼ`, 64 × 64 with
periodic boundaries, driven to its ground state by Metropolis single-spin-flip
dynamics under a temperature cooled geometrically from 5 to 0.05 J/k_B over
3000 sweeps, from a seeded stable random stream. The run first reaches the
exact ground state `E/N = -2J` at `|m| = 1` at T = 0.93 J/k_B. A single flip
out of it costs 8J, with a Boltzmann factor of 2 × 10⁻⁴ at that temperature, so
60 later sweeps still end with a spin flipped, the last at T = 0.57 J/k_B; from
there the ground state is held to the end. The script asserts the final state,
and the 60 against the 51 ± 7 that the Boltzmann factor gives over the cooling
schedule. The
snapshots show the disordered, critical and ordered regimes against the exact
Onsager temperature `T_c = 2/ln(1+√2) = 2.269 J/k_B`. Before the anneal the
script checks the Metropolis kernel itself: the flip energy
`ΔE = 2 J sᵢⱼ Σ_nb` is asserted equal to the difference of two full energy
evaluations at 400 sites of random lattices, under periodic and under open
boundaries. Spins are drawn black (+1) and white (−1);
the energy and magnetisation traces keep their own colours.

![Annealing animation](figures/ising_annealing.gif)

The anneal itself, a 64 × 64 lattice sampled every 20 sweeps, with the energy
and magnetisation traced out beside it. Where the snapshots show three states,
the animation shows the domains forming: small and short-lived well above
`T_c`, growing and merging as the temperature passes through it, then freezing
into a single spanning domain. At this cooling rate the magnetisation stays
below 0.3 until T ≈ 1.8 J/k_B, well below the Onsager temperature, which is
marked.

## What changed from the original

Ported from `Ising_model.cpp` in `Code_Archive/Old_2018/C_C++/` on the `legacy`
branch, which has five defects, every one of them affecting the result:

- rows allocated `n+2` wide but columns only `n`, while the code indexed
  `a[i][n+1]`, a heap overflow of two `short`s past every row
- ghost cells bordered once before the Monte Carlo loop and never refreshed,
  so the periodic images were wrong after the first boundary flip
- `H` summed every bond twice and `deltaH` returned `4·J·s·Σ` where the correct
  value is `2·J·s·Σ`, so `exp(-ΔE/T)` ran at an effective temperature of T/2
- `T = T/step` applied per attempted flip rather than per sweep, so the
  temperature underflowed to zero within a few thousand attempts and the
  dynamics degenerated to greedy descent
- `unsigned int` loop counter tested against an `unsigned long long` bound, an
  infinite loop for any bound ≥ 2³², and the prompt asked the user for "un
  numar foarte mare"

Periodic boundaries are handled here by modular indexing, without ghost cells.
The original also misspelled its input and output file names as `Icing`.
