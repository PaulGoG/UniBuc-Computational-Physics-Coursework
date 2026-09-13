# Compton scattering — BSc year 1, semester 2 (2018)

## `compton_scattering_chain.jl`

![Compton scattering chain](figures/compton_scattering_chain.png)

A 1 MeV photon followed through successive Compton scatterings down to a 20 keV
absorption threshold, over 20 000 histories. Takes a mean of 28.4 scatterings.

![Energy degradation](figures/compton_degradation.gif)

The ensemble degrading, and the reason the sampling matters. The left panel is
the photon-energy distribution after each scattering; the right is the
Klein–Nishina cross-section **at the mean energy of the survivors at that
moment**. At 1 MeV it is strongly forward-peaked; by 30 keV it has relaxed
towards the symmetric Thomson form. The distribution the sampler must draw from
therefore changes at every step, which is precisely what a single fixed
distribution — let alone a uniform one — cannot represent.

## What changed from the original

Ported from `EfectulCompton.cpp`. The Compton relation, the recoil-angle formula
and the momentum-triangle sine rule were all individually correct. The sampling
was not: **the scattering angle was drawn uniformly on [0°, 180°]**, which is
neither Klein–Nishina nor even isotropic scattering — isotropy is uniform in
cos θ, not in θ. The right-hand panel of the static figure shows the difference;
the mean scattering angle is 64.7° under Klein–Nishina against 89.9° for the
original sampling, and every energy distribution the original produced was
unphysical.

Two further defects are simply gone rather than fixed. The original computed the
scattered energy a second time from the recoil angle by the sine rule and
averaged the two, but that second value is an algebraic identity of the first
and carries no independent information — so the recoil angle is no longer
computed at all. And its "abaterea standard" was a relative deviation of a
single sample from a running mean, cast to `int` and taken mod 100, with no sum
of squares anywhere.

The angle is now sampled from the Klein–Nishina differential cross-section by
rejection, recomputing the envelope at each scattering because the distribution
depends on the current photon energy.
