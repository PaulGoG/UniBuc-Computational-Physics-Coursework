# Fission observables — MSc year 2 (2022–2023)

Neutron-induced fission of ²³⁵U at thermal energy, from a Y(A_H, Z_H, TKE)
matrix, the AME2020 masses, the FRDM ground-state deformations, the
Gilbert–Cameron shell corrections, and four measured ν(A) and four measured
neutron-spectrum sets. Ported from `Fission_1.jl` to `Fission_5.jl` in
`Julia-Workflow-FFUB/Fission_M_2/` on the `legacy` branch. `fission_data.jl`
reads every table once into typed containers; `fission_core.jl` holds the
compound system, the isobaric charge distribution, Q values and separation
energies from the mass table, the reduction of the yield matrix to its
measured cells, and the yield-weighted mean with its uncertainty. Every
script asserts an identity or a closed form.

## `fission_q_value.jl`

![Q value](figures/fission_q_value.png)

Q(A_H, Z_H) = Δ(²³⁶U) − Δ(A_H, Z_H) − Δ(A₀ − A_H, Z₀ − Z_H) for heavy-fragment
masses 118 to 160, averaged at each mass over the three charges nearest
Z_p(A_H) = Z_UCD(A_H) − 0.5 with the weights of a Gaussian isobaric charge
distribution of rms width 0.6, the mass-excess uncertainties propagated in
quadrature. The maximum is 196.03 ± 0.01 MeV at A_H = 130, the ¹³²Sn shell
closure; the yield-weighted mean over the 41 masses with a yield is
⟨Q⟩ = 186.62 MeV, while the arithmetic mean over the 43 mass points, which is
what the original reported, is 183.70 MeV. The polarisation and width are
rounded values of the Wahl Z_p systematics (At. Data Nucl. Data Tables **39**,
1 (1988), doi:10.1016/0092-640X(88)90016-2).

## `fragment_yields.jl`

![Fragment yields](figures/fragment_yields.png)

What the matrix measures is Y(A_H, TKE) over 4141 cells, 118 ≤ A_H ≤ 158 and
100 ≤ TKE ≤ 200 MeV. The five charge rows of each cell split one measured yield
by fractions that are identical across TKE for a given mass, to 3 × 10⁻⁶,
which the script asserts: Y(Z), Y(N) and the even–odd staggering
δ = 0.121 describe the charge split imposed on the matrix, not a measurement.

| Quantity | Value |
|---|---|
| Y(A) peak | A_H = 140, 6.71 ± 0.14 % |
| Y(Z) peak | Z_H = 54 (imposed split) |
| ⟨A_H⟩ | 139.323 ± 0.019 |
| ⟨A_L⟩ | 96.677 ± 0.019 |
| ⟨TKE⟩ | 170.549 ± 0.047 MeV |
| ⟨Q⟩ | 186.624 ± 0.023 MeV |
| ⟨TXE⟩ | 22.621 ± 0.041 MeV |
| ⟨KE_L⟩, ⟨KE_H⟩ | 100.55, 69.99 MeV |
| S_n(²³⁶U) | 6.546 MeV |
| ⟨TKE⟩(A), TXE(A) | 148.3–178.7, 22.0–40.7 MeV |

Every total is a weighted mean over the measured cells of a quantity fixed
per cell — the cell's A_H, A_L or TKE, the mass's Q(A), and
Q(A) + S_n − TKE for TXE — so its variance is Σ_c [(q_c − ⟨q⟩) σY_c/ΣY]², with
the rows of one cell added linearly (they share one relative error) and the
cells in quadrature. Propagated at the mass level instead, after averaging TKE
within each mass, σ⟨TKE⟩ comes out 0.022 MeV: the within-mass spread of TKE is
lost. The identities ⟨A_H⟩ + ⟨A_L⟩ = 236 and ⟨TXE⟩ = ⟨Q⟩ + S_n − ⟨TKE⟩ are
asserted; the single-fragment energies follow from momentum conservation,
KE_L = TKE·A_H/A₀. The mass-excess uncertainties behind Q are not carried into
the totals; they are at most 0.09 MeV per mass.

`Fission_2.jl` summed Y(N) over all rows with a given N inside the double
loop over (A_H, Z_H), once per pair mapping to that N: its total came to
805 % instead of 100 % and its peak moved from N = 87 to 86, reproduced by the
script. Three of its uncertainty accumulations added linearly where the
comments said quadrature.

## `neutron_spectrum.jl`

![Neutron spectrum](figures/neutron_spectrum.png)

The Madland–Nix spectrum (Nucl. Sci. Eng. **81**, 213 (1982),
doi:10.13182/NSE82-5) with a constant compound-nucleus cross-section and a
triangular residual-temperature distribution,

```
N(E) = 1/(3√(E_f T_m)) [u₂^{3/2} E₁(u₂) − u₁^{3/2} E₁(u₁) + γ(3/2, u₂) − γ(3/2, u₁)],
u₁,₂ = (√E ∓ √E_f)² / T_m,
```

averaged over the mass yield as `Fission_4.jl` did: per heavy-fragment mass,
TXE = Q(A, Z_p) + S_n − ⟨TKE⟩(A) at the single most probable charge,
T_m = √(C·TXE/A₀) with C = 10 MeV as in that file, E_f per nucleon from
momentum conservation, the light- and heavy-fragment spectra averaged and
weighted by Y(A) over 41 masses. T_m spans 0.94–1.37 MeV and E_f 0.31–1.27 MeV.

The closed form integrates to one over [0, ∞) for every E_f and T_m, which the
script asserts for each fragment to 10⁻⁹. `Fission_4.jl` wrote the prefactor
as `(1/3*sqrt(E_F*T_MAX))`, which Julia parses as (1/3)√(E_f T_m): its
spectrum integrates to E_f T_m, between 0.31 and 1.33 across the fragments,
so the error reweights the mass average rather than cancelling under the
renormalisation applied afterwards. The mass-averaged mean energy is
⟨E⟩ = 2.103 MeV (equivalent Maxwellian temperature 2⟨E⟩/3 = 1.402 MeV) against
2.203 MeV with the original prefactor, 4.8 % higher.

Maxwellian fits to the four measured spectra, data and model both normalised
to unit area over the measured window and χ² reduced by ν = N − 1:

| Dataset | Window [MeV] | Fraction of the Maxwellian in the window | T_M [MeV] | χ²/ν |
|---|---|---|---|---|
| Göök, laboratory | 0.55–12.5 | 84 % | 1.295 | 3.95 |
| Vorobyev, laboratory | 0.22–16.7 | 96 % | 1.360 | 2.61 |
| Göök, centre of mass, light fragment | 0.05–6.95 | 99 % | 0.800 | 5.04 |
| Göök, centre of mass, heavy fragment | 0.05–6.95 | 99 % | 0.841 | 1.88 |

The window matters: a Maxwellian normalised over all energies compared with
the Göök laboratory spectrum, which holds 84 % of it, offsets the whole ratio
by a fifth and returns T_M = 1.322 MeV at χ²/ν = 55. The centre-of-mass sets
sit well below the laboratory ones, as removing the fragment motion requires,
and are counting-limited above about 6 MeV. `Fission_5.jl` divided χ² by N
rather than ν, a convention, and bounded its optimiser at T_M = 0, where
T_M^{−3/2} is infinite; both files re-integrated E₁, γ(3/2, x) and the window
integral by quadrature at every evaluation.

## `neutron_multiplicity.jl`

![Neutron multiplicity](figures/neutron_multiplicity.png)

ν(A) from the energy balance of `Fission_3.jl`,

```
ν_pair = (TXE − q) / (⟨ε⟩ + ⟨S_n⟩ + p),   ⟨ε⟩ = 4T_m/3,   T_m = √(TXE/(a_L + a_H)),
p = 6.71 − 0.156 Z₀²/A₀ = 1.115 MeV,     q = 0.75 + 0.088 Z₀²/A₀ = 3.906 MeV,
```

with the Gilbert–Cameron level-density parameter a = A[0.00917 (S_Z + S_N) +
0.142] (Can. J. Phys. **43**, 1446 (1965), doi:10.1139/p65-139). The pair
multiplicity is split by the heavy fragment's share R of the excitation
energy. `Fission_3.jl` takes R from the scission-point deformation: each
fragment pays the liquid-drop energy of deforming from its FRDM ground-state
β₂ to a scission β(Z), piecewise linear through (28, 0), (41, 0.58),
(44, 0.58), (50, 0) and (65, 0.6), and the remainder is shared in the ratio of
the level-density parameters. The statistical limit R = a_H/(a_L + a_H) is
computed alongside. Everything is evaluated per (A_H, Z_H) over the charges
the mass table lists, averaged over Z_H with the Gaussian isobaric weights,
and weighted by the mass yield.

| | Yield-weighted |
|---|---|
| ⟨ν_pair⟩ | 2.538 |
| ⟨ν_pair⟩, unweighted over masses | 2.827 |
| ⟨TXE⟩ | 22.56 MeV |
| R, scission deformation | 0.398 |
| R, level densities alone | 0.463 |

Against the measured sets, whose own yield-weighted per-fragment means give
2.40–2.46 neutrons per pair:

| Dataset | Points | 2 ν̄ measured | Deformation partition, χ²/N, rms | Level densities, χ²/N, rms |
|---|---|---|---|---|
| Göök | 79 | 2.460 | 76, 0.36 | 192, 0.48 |
| Maslin | 48 | 2.434 | 31, 0.37 | 24, 0.39 |
| Nishio | 59 | 2.403 | 8.9, 0.17 | 17, 0.23 |
| Vorobyev | 97 | 2.409 | 67, 0.33 | 112, 0.40 |

The deformation partition of the original reproduces the sawtooth better on
three of the four sets, and it is what the figure draws; the level-density
partition is drawn dashed. The model sits 3–6 % above the measured pair
values: prompt γ emission competes for the same excitation energy and this
balance does not account for it. The script asserts that β(Z) vanishes at the
closed shells Z = 28 and 50, that 0 < R < 1 for every split, and that
ν_H + ν_L = ν_pair at every mass.

`Fission_3.jl` bounded its charge loop by the index column of the
Gilbert–Cameron table, 11 to 150; charges far from Z_p enter with a Gaussian
weight below 10⁻⁵ beyond |Z − Z_p| = 3, so the bound cost time, not results.
Its `Energie_separare` returned `[S, σ]` on success and a bare `NaN` on
failure, indexed by the caller regardless. `Dump.jl` referenced names it
neither defined nor imported and was never included; it is not ported.

**Constants.** The energy balance and its constants are those of the
TXE-partition application of the fission course this was written for
(A. Tudora, University of Bucharest, 2022), which the 2023 files carried
without attribution. The liquid-drop deformation energy with 15.4941, 17.9439,
1.7826, 0.7053 and 1.1529 MeV is the Myers–Swiatecki parameter set of 1967
(Ark. Fys. **36**, 343, the Lysekil symposium); the scission deformation is
the application's piecewise-linear β(Z) through (28, 0), (41, 0.58),
(44, 0.58), (50, 0) and (65, 0.6); p = 6.71 − 0.156 Z₀²/A₀ and
q = 0.75 + 0.088 Z₀²/A₀ are its systematics of the prompt γ energy
E_γ = p ν + q. The level-density constant of the spectrum, T_m = √(C·TXE/A₀)
with C = 10 MeV, is the application's ⟨a⟩ = A₀/C for an energy-independent
compound-nucleus cross-section (C = 11 MeV with an optical-model one).

## Data

| File | Content | Source |
|---|---|---|
| `Yield/U5YAZTKE.STR` | Y(A_H, Z_H, TKE) [%] and its error, 20 705 rows, five charges per (A_H, TKE) cell | Header `5Z/A Y(A,TKE) Straede`: the Y(A, TKE) of Straede, Budtz-Jørgensen and Knitter, Nucl. Phys. A **462**, 85 (1987), doi:10.1016/0375-9474(87)90381-2 (EXFOR 23591), with a charge split imposed on it for the course: a Gaussian isobaric distribution about Z_p(A) with the polarisation ΔZ(A) and width of Wahl's Z_p model (At. Data Nucl. Data Tables **39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2), five charges per mass; the digitisation of the measurement is not recorded |
| `Experimental/Neutron_multiplicity/U5NUA*.DAT` | ν(A) and its error per fragment mass | Göök, Hambsch, Oberstedt and Vidali, Phys. Rev. C **98**, 044615 (2018), doi:10.1103/PhysRevC.98.044615; Maslin, Rodgers and Core, Phys. Rev. **164**, 1520 (1967), doi:10.1103/PhysRev.164.1520; Nishio, Nakagome, Yamamoto and Kimura, Nucl. Phys. A **632**, 540 (1998), doi:10.1016/S0375-9474(98)00008-6; Vorobyev et al., EPJ Web Conf. **8**, 03004 (2010), doi:10.1051/epjconf/20100803004. The files carry the first author's name only; the digitisation is not recorded |
| `Experimental/Neutron_spectrum/U5SP*.DAT` | Neutron spectra with errors: Göök laboratory frame, Göök centre of mass for the light and heavy fragment, Vorobyev laboratory frame | The Göök sets from the measurement of Phys. Rev. C **98**, 044615 (2018); the Vorobyev set from the same group's ²³⁵U measurements (EPJ Web Conf. **8**, 03004 (2010) describes them); the exact publication and digitisation are not recorded in the files |
| `Mass_defects/AUDI2021.csv`, `AUDI95.csv` | AME2020 and AME1995 mass excesses | As in `msc/03-radionuclides`, byte-identical copies |
| `Auxiliary_parametrisations/B2MOLLER.ANA` | Z, A, ground-state β₂ | FRDM, Möller, Nix, Myers and Swiatecki, At. Data Nucl. Data Tables **59**, 185 (1995), doi:10.1006/adnd.1995.1002, via `moller.dat` of RIPL-1 as the header says |
| `Auxiliary_parametrisations/SZSN.GC` | Nucleon number, shell corrections S(N) and S(Z) [MeV] | Gilbert and Cameron, Can. J. Phys. **43**, 1446 (1965), doi:10.1139/p65-139; the extraction is not recorded |
