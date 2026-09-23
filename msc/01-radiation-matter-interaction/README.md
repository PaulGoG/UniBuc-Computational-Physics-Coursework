# Interaction of radiation with matter — MSc year 1 (2021–2022)

Four scripts from the laboratory course, ported from the Julia files in
`Julia-Workflow-FFUB/IRM_M_1/` on the `legacy` branch. `radiation_matter_core.jl`
holds the linear least squares with a general data covariance that two of them
share. Every result below is asserted in the script that produces it, against a
tabulated value, a statistical identity or an evaluated datum.

## `bethe_bloch_stopping_power.jl`

![Bethe-Bloch](figures/bethe_bloch_stopping_power.png)

Electronic stopping power of a 5 MeV α in silicon from the Bethe–Bloch formula
with the Sternheimer–Barkas mean excitation potential, I = 172.3 eV:
1473 MeV cm⁻¹, or 632 MeV cm² g⁻¹. ASTAR (Berger, Coursey, Zucker and Chang,
NIST Standard Reference Database 124, doi:10.18434/T4NC7P) gives
617.4 MeV cm² g⁻¹, so the formula sits 2.4 % high, as it must without the shell
and Barkas corrections at β = 0.052. Integrating 1/(dE/dx) from 0.5 to 5 MeV
gives a CSDA range of 21.6 µm against 22.0 µm from the ASTAR ranges at the two
energies. Both comparisons are asserted within 5 %.

The formula was transcribed correctly in `Calcul_Bethe_Bloch.jl`. The target
material is named nowhere in that file, A is written as 28 rather than 28.085,
and the logarithm's argument types the same sub-expression twice instead of
squaring it, which hides that the second factor is W_max.

## `alpha_attenuation_mylar.jl`

![Alpha attenuation](figures/alpha_attenuation_mylar.png)

Residual energy of α particles after each of five absorber settings of stacked
Mylar, fitted with `ε(x) = ε₀ − S x` by weighted least squares on all five
measurements: S = 380 ± 150 keV per setting, ε₀ = 4900 ± 450 keV,
χ² = 0.31 on 3 degrees of freedom.

`Atenuare_alpha.jl` differenced the data first, Δε = ε(0) − ε(x), and fitted
the four differences unweighted with a free intercept, drawing error bars that
did not enter the fit and never printing the coefficients. Differencing
correlates the points through the shared ε(0). With that covariance written out
and no intercept — Δε(0) = 0 by construction — the differenced fit returns the
same slope, error and χ² as the five-point fit, which the script asserts. A
free intercept in the differenced model is a redundant parameter: it comes
back as −100 ± 1200 keV and widens the slope error by 1.66×.

**Abscissa.** The absorber setting is recorded as the original labelled it,
"x (μm)", and the script is the only surviving record of the measurement.
Read as micrometres the slope would be 375 keV µm⁻¹, about three times the
electronic stopping power of Mylar between 4 and 5 MeV (115–135 keV µm⁻¹,
ASTAR), and a 1.85 MeV residual after 8 µm is not compatible with the 28.8 µm
CSDA range of a 5 MeV α in Mylar, so the setting is more likely a foil count
than a thickness. The axis is labelled as the setting, without a unit.

## `natural_cd_cross_section.jl`

![Cd cross-section](figures/natural_cd_cross_section.png)

Abundance-weighted (n,γ) cross-section of natural cadmium at 0.25 eV from the
eight stable isotopes: 2203 b, of which ¹¹³Cd supplies 99.93 % through its
0.178 eV resonance. The abundances and cross-sections are those entered in
`SigmaCdNatural.jl`, which cites the IAEA neutron cross-section atlas for the
latter; the abundances are asserted to sum to 100 % within 0.1.

The energy is part of the result. The conventional number for natural Cd is
2520 b at the 2200 m s⁻¹ thermal point (0.0253 eV, where σ(¹¹³Cd) = 20 600 b),
and a bare "2203 b" would be read as that. The original never printed the
result — a bare top-level expression, so running the file produced no output —
and its arrays were positional, with nothing tying σ = 18 000 b to ¹¹³Cd.

## `neutron_activation_halflives.jl`

![Activation half-lives](figures/neutron_activation_halflives.png)

Half-lives of three activation products from counts accumulated in successive
equal acquisition intervals, by weighted least squares of ln N = a − λt with
Poisson weights on every recorded point, and the NaI(Tl) energy calibration
that identifies the ¹²⁸I photopeak. Reference half-lives from NUBASE2020
(Kondev et al., Chin. Phys. C **45**, 030001 (2021),
doi:10.1088/1674-1137/abddae).

| Product | Points | T½ [min] | χ²/ν | NUBASE2020 [min] |
|---|---|---|---|---|
| ²⁸Al | 5 | 2.290 ± 0.045 | 4.2/3 | 2.245 ± 0.005 |
| ²⁷Mg | 5 | 9.30 ± 0.72 | 46.7/3 | 9.435 ± 0.027 |
| ¹²⁸I | 4 | 22.6 ± 1.4 | 0.4/2 | 24.99 ± 0.02 |

²⁸Al and ¹²⁸I agree with the evaluated values within 1.0σ and 1.7σ. The ²⁷Mg
series is not Poisson-compatible with a single exponential (p = 4 × 10⁻¹⁰): its
second count is 3.5σ below its neighbours; the error scaled by √(χ²/ν) is
2.8 min, and the central value still lands on the evaluated one.

`ReactiiNeutronice.jl` and `FitActivareNaITl.jl` linearised the decay law as
ln(N₀/Nᵢ) = λ(tᵢ − t₀) with the first count as reference; the first fitted
unweighted and never called `stderror`, the second used weights Nᵢ, which
drops the 1/N₀ term. Every ln(N₀/Nᵢ) shares N₀, so those points carry the
covariance 1/Nᵢ δᵢⱼ + 1/N₀; with a free intercept the common term is absorbed
and the slope and its variance are exactly those of a weighted fit to N₂ … Nₙ
alone — the reference count, the most precise of the series, drops out of λ.
The script asserts this identity and prints the corresponding values
(2.354 ± 0.076, 14.2 ± 2.6 and 21.8 ± 2.1 min) beside the all-point ones. The
first ¹²⁸I count has no recorded time in the original, which lists 520, 780 and
1040 s for the other three; it is placed one interval earlier, at 260 s, and
the slope does not depend on that choice. Equal acquisition intervals matter:
the factor (1 − e^{−λΔ}) relating counts to activity is then common to every
point and changes only the intercept.

The calibration is an unweighted straight line through four lines (511, 1173,
1274 and 1332 keV at channels 77, 172, 185 and 195): a = −27.6 ± 11.2 keV,
b = 6.996 ± 0.068 keV per channel, ρ(a, b) = −0.96, residual standard
deviation 6.4 keV. The studied peak at channel 67 — below the lowest
calibration point, an extrapolation — comes out at 441.1 ± 7.0 keV with the
full covariance; dropping the covariance term gives ± 12.1 keV. ¹²⁸I β⁻ decay
feeds the 442.9 keV 2⁺ → 0⁺ transition of ¹²⁸Xe (Elekes and Timar, Nucl. Data
Sheets **129**, 191 (2015), doi:10.1016/j.nds.2015.09.002), 0.25σ away. The
original printed the energy without an uncertainty.
