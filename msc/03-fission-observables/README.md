# Fission observables — MSc year 2

Neutron-induced fission of ²³⁵U at thermal energy, from the Straede
Y(A_H, Z_H, TKE) matrix, the AME mass evaluation, Möller–Nix deformations,
Gilbert–Cameron level densities, and four measured ν(A) and neutron-spectrum
datasets.

## `fission_data.jl`

Loaders and indices for all of the above. The originals accessed these tables by
boolean masking inside nested loops; `Fisiune_2.jl:TKE_A` alone made roughly
8000 full passes over the 20 705-row yield matrix.

## `fission_q_value.jl`

![Q value](figures/fission_q_value.png)

Q(A_H) averaged over the isobaric charge distribution. Mean **184.6 MeV**, range
165.6–196.0, peaking at **A_H = 130** — at the ¹³²Sn shell closure, where it
must.

## `fragment_yields.jl`

![Fragment yields](figures/fragment_yields.png)

| Quantity | Value |
|---|---|
| Y(A) peak | A_H = **140**, 6.71 % |
| Y(Z) peak | Z_H = **54** (Xe) |
| Even–odd staggering δ | 0.121 |
| ⟨TKE⟩ | **170.55 MeV** |
| S_n(²³⁶U) | **6.546 MeV** |
| TXE(A) | 20.7–44.0 MeV |
| KE_L(A) | 77.9–101.1 MeV |
| KE_H(A) | 49.0–80.2 MeV |

The single-fragment kinetic energies follow from momentum conservation,
KE_L = TKE·A_H/A₀; they close against TKE to 2.8 × 10⁻¹⁴ MeV.

⟨TKE⟩ and S_n both reproduce the evaluated values to three digits.

**Y(N) was over-counted.** The original summed over all rows with a given
N = A_H − Z_H *inside* the double loop over (A_H, Z_H), so the sum was re-added
once per (A, Z) pair mapping to that N — and with five charge splits per mass,
several do. The measured effect: total yield **805.4 %** instead of 100.2 %, an
eightfold inflation, with the peak displaced from N = 87 to N = 86. The
normalisation applied afterwards hid the absolute error but not the distortion.

Also corrected: uncertainties added linearly in three places while the comments
stated quadrature; and `Sortare_distributie` assumed its abscissa contained
every integer between its extremes with no gaps or duplicates — the author hit
this, and five commented-out `deleteat!` lines in `Fisiune_3.jl` document it.

## `neutron_spectrum.jl`

![Neutron spectrum](figures/neutron_spectrum.png)

Madland–Nix spectrum, **averaged over the mass yield** as the original did —
per-mass TXE, TKE, T_m = √(10·TXE/A) and E_f, with the light- and
heavy-fragment spectra averaged and weighted by Y(A) over 41 masses — and
Maxwellian fits to four measured spectra.

The mass average is what makes the prefactor error measurable. T_m spans
0.936–1.365 MeV and E_f spans 0.310–1.273 MeV across the fragments, so a
prefactor carrying both cannot be absorbed into an overall normalisation:

| | ⟨E⟩ | (2/3)⟨E⟩ |
|---|---|---|
| Correct prefactor | 2.1025 MeV | 1.4017 MeV |
| 2018 prefactor | 2.2025 MeV | 1.4683 MeV |

a **4.76 %** shift in the mean energy. Real, and in the direction that matters,
but modest — worth stating precisely rather than left as an assertion. The
correct equivalent Maxwellian temperature of 1.40 MeV sits above the evaluated
1.32, which is the Los Alamos model with C = 10 rather than a coding issue.

| Dataset | T_M [MeV] | χ²/ν |
|---|---|---|
| Göök, lab | 1.3222 | 54.9 |
| Vorobyev, lab | 1.3362 | 3.35 |
| Göök, CM light | 0.7964 | 4.98 |
| Göök, CM heavy | 0.8388 | 1.83 |

The laboratory temperatures land on the evaluated ~1.32 MeV for ²³⁵U(n_th,f),
and the centre-of-mass values sit well below them, as removing the fragment
motion requires.

**The Madland–Nix prefactor was wrong.** The original wrote

```julia
return (1/3*sqrt(E_F*T_MAX)) * ( ... )
```

which parses as `(1/3)·√(E_f·T_m)`, since `/` and `*` share precedence and
associate left to right. The prefactor is **1/(3√(E_f T_m))** — it multiplied
where it had to divide. Because E_f and T_m both vary with fragment mass, the
error does not cancel under the Maxwellian renormalisation applied afterwards;
it reweights the mass average.

E₁(z) and γ(3/2,x) were also re-integrated with `quadgk` at every evaluation,
some 34 000 adaptive quadratures per run, where `SpecialFunctions` has both in
closed form. And `Fisiune_5.jl` divided χ² by N while labelling it χ², and
bounded its optimiser at T_M = 0 where T_M^{-3/2} is infinite.

## `neutron_multiplicity.jl`

![Neutron multiplicity](figures/neutron_multiplicity.png)

ν(A) from the energy balance `Fisiune_3.jl` used,

```
ν_pair = (TXE − q) / (⟨ε⟩ + ⟨S_n⟩ + p),   ⟨ε⟩ = 4T_m/3,  T_m = √(TXE/a_tot)
p = 6.71 − Z²·0.156/A = 1.116 MeV        q = 0.75 + Z²·0.088/A = 3.905 MeV
```

against Göök, Maslin, Nishio and Vorobyev.

⟨ν_pair⟩ = **2.883** against the evaluated **2.42**, +19 %. The residual is the
model: prompt γ emission competes for the same excitation energy and carries off
6–7 MeV per fission, which this balance does not account for. The measured sets
average 1.16–1.32 neutrons per fragment, consistent with 2.42 for the pair.

The pair multiplicity is split between fragments by the ratio of their
level-density parameters. The original took that ratio from its scission-point
deformation energies; that model is not reproduced, and the level-density ratio
is the statistical-equilibrium limit of the same quantity.

The defect fixed from `Fisiune_3.jl`:

```julia
for Z_H in minimum(dGC.n):maximum(dGC.n)
```

bounded the loop over *fragment proton number* by the index column of the
Gilbert–Cameron table, so Z_H ran 11 to 150. Unphysical pairs such as
(A_H = 130, Z_H = 100) survived the guards and polluted the level-density
distribution and everything averaged over it. The charge range now comes from
the yield matrix, which spans Z_H = 44–63.

The scission-point deformation model of `Fisiune_3.jl` is **not** reproduced. It
hardcodes six undocumented parametrisation constants — `a = 0.58/13`,
`b = -28*a`, `-0.58/6`, `0.6/15`, `-50*a` — with no source given, and those are
the most model-dependent numbers in the archive.

`Dump.jl` is deleted: it referenced `distributie_bidym`, `Y`, `y_A`, `DataFrame`
and `CSV`, none of which it defined or imported, so it could not run and was
never included by anything. Its one useful routine, the (A,Z)-resolved
⟨TKE⟩(A,Z), is subsumed by the marginalisation in `fragment_yields.jl`.
