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

## Magnitude–redshift

The second analysis, and the whole subject of `Friedmann_Magnitude.jl`: 6571
usable host-galaxy magnitudes against redshift, fitted with the flat-ΛCDM
luminosity distance.

The distance scale is well determined — the effective absolute magnitude comes
out **M = −19.57** with a residual RMS of 1.35 mag — but **Ωm rails at its
bound and is unconstrained**. That is the honest result: these are host-galaxy
magnitudes, not standardised SN Ia peak magnitudes, so with 1.35 mag of scatter
the data fix the distance scale and carry essentially no cosmological
information. The binned medians show the trend clearly; the constraint on Ωm
does not exist.

M landing near −19.6, close to the SN Ia absolute magnitude of about −19.3, is a
coincidence of this sample rather than a result.

## What changed from the original

The original fitted a **straight line** to H(z) (`LineFit(t, p) = p[1]*t .+ p[2]`);
flat ΛCDM is fitted here instead. The magnitude–redshift analysis is kept, but
its route is not: the original drove a `MonteCarloProblem` with randomised
initial conditions through `build_loss_objective`, which is replaced by the
closed form for the luminosity distance. The equation-of-state survey sketched
in the original's docstring — quintessence, Chaplygin and modified Chaplygin gas
— was never implemented there and is not implemented here.
