# Cosmology — Dubna 2019

Joint work with [Alexandru Crăciun](https://github.com/Craciun-Alexandru),
presented together at JINR Dubna in 2019.

## `friedmann_expansion.jl`

![Friedmann expansion](figures/friedmann_expansion.png)

Flat ΛCDM fitted to eleven measurements of the Hubble parameter, and the
expansion history integrated for three cosmologies.

| | Fitted | Planck 2018 |
|---|---|---|
| H₀ | **73.3 ± 5.1** km s⁻¹ Mpc⁻¹ | 67.4 ± 0.5 |
| Ω_m | **0.267 ± 0.064** | 0.315 ± 0.007 |

χ²/dof = 0.73 on 10 degrees of freedom. Both parameters agree with Planck inside
the (large) uncertainties of this small dataset; the central H₀ sits on the
local-measurement side of the Hubble tension, which is what H(z) compilations of
this kind tend to give.

**The original did not run.** `Friedmann_Eq.jl` called a function `Friedmann`
that was never defined — the definition present is
`FriedmannTimeDependentNeutralFluid` — and used `fit.param` some forty lines
before `fit` was assigned. It also declared

```julia
const G = 6.67408*1e-11        #Cosmological constant
```

which is the gravitational constant, not Λ. `Friedmann_Magnitude.jl` read its
data through a hardcoded `C:\Users\GoG\Desktop\...` path and, inside an ODE
right-hand side, wrote `u[3] = ...` — mutating the state vector instead of
returning a derivative.

`data/supernova_magnitudes.csv` (20 588 rows) is retained but not yet used; the
magnitude–redshift diagram is future work.
