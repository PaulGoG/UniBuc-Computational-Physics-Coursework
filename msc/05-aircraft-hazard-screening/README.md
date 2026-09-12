# Aircraft crash hazard screening

![Impact frequency](figures/aircraft_crash_frequency.png)

Annual frequency of an aircraft crash onto an installation sited a perpendicular
distance y₀ from an air corridor, screened against the 10⁻⁵ and 10⁻⁷ yr⁻¹
design-basis thresholds.

## Formulation

Following the structure of DOE-STD-3014 and NUREG-0800 §3.5.1.6:

```
F = N · P · ∫ f(s(x)) dx · A_eff            [yr⁻¹]

  N      flights per year on the route                      [yr⁻¹]
  P      probability of loss of control per km of flight    [km⁻¹]
  f(s)   crash-location probability density per unit AREA   [km⁻²]
  A_eff  effective target area of the site                  [km²]
```

`yr⁻¹ · km⁻¹ · km⁻² · km · km² = yr⁻¹`, so the result is a frequency by
construction. The crash-location density is the object that matters: for an
exponential radial density `p(r) = g e^{−gr}` per unit r and isotropic direction,
the density per unit **area** is

```
f(r) = p(r) / (2πr) = g e^{−gr} / (2πr)     [km⁻²]
```

It is the `1/(2πr)` that turns a radial density into an areal one.

## Results

| y₀ | R = 50 m | R = 200 m |
|---|---|---|
| 5 km | 1.37 × 10⁻⁸ yr⁻¹ | 2.20 × 10⁻⁷ yr⁻¹ |
| 10 km | 3.19 × 10⁻⁹ | 5.10 × 10⁻⁸ |
| 25 km | 6.56 × 10⁻¹¹ | 1.05 × 10⁻⁹ |
| 50 km | 1.49 × 10⁻¹³ | 2.39 × 10⁻¹² |

At the nearest approach considered, 5 km, the 10⁻⁷ yr⁻¹ threshold corresponds to
a site radius of 135 m and the 10⁻⁵ threshold to 1349 m.

## What was wrong

The 2018 file computed two estimators,

```julia
integrandRadial          = (x/s) e^{-gs}
integrandRadialUnghiular = (x y₀/s²) e^{-gs}     # scaled by g/2
```

**Neither is an areal density.** The first integrates to a length, the second to
a dimensionless number, so after both were multiplied by the same prefactor
`P·N·π(R·10⁻³)²` the two results differed by one power of length and *neither*
was a frequency — yet they were plotted on the same axis and compared against
the same yr⁻¹ thresholds. The first also omitted the normalisation `g` of the
exponential density which the second included, a further factor 1/g = 4.35
between them.

Two further corrections:

- **The route was integrated on one side only**, 0 to x₀, while the comment
  called x₀ the route *length*. A site beside the middle of a corridor receives
  crashes from both directions, so the one-sided result is low by a factor
  approaching two.
- **The effective target area was the bare geometric disc** πR². DOE-STD-3014
  adds skid and shadow contributions, so the geometric footprint is a lower
  bound on the target the aircraft actually presents. It is used here and
  labelled as such rather than silently.
