*[← all rebuilds](../) · paper: [arXiv:2609.18036](https://arxiv.org/abs/2609.18036) · what I sent the authors: [`reported.md`](reported.md) · data: [`DATA.md`](DATA.md)*

# arXiv:2609.18036 — Chao, Shih, Deschamps & Tan
*Gravitational energy released when Earth quakes: heat supply for island-arc volcanism.*
Submitted Sep 2026. Read fire 286 (Sep 20 2026).

## What it claims
Thrust-faulting earthquakes induce a release of the Earth's gravitational energy Eg at a mean
power of **~10 TW** — about 1000× the global seismic energy release (~10 GW) and a quarter of the
terrestrial heat flow (~45 TW). Normal-faulting events at spreading centres take up **~2 TW**.
Dahlen (1977) argued that such an Eg change is exchanged reversibly with elastic strain energy and
is therefore not available as heat; the paper rejects that, holding that the post-seismic state is
a genuinely lower-Eg equilibrium and the difference goes irreversibly to heat. Since thrust events
sit on subduction zones, that heat is offered as the missing supply for island-arc volcanism — the
"cool end" of mantle convection. The paper takes Morgan et al. (2016)'s 11–14 TW of Eg dissipation
in the descending limb and names thrust earthquakes as "the primary physical mechanism to
accomplish the said Eg-heat conversion," manifesting as the high heat flow of Fig. 2(a).

## What reproduced (`parse.py`, `budget.py`)
The instrument was validated first: my ndk reader returns the published GCMT scalar moments for
2004 Sumatra-Andaman (3.953e22 vs 3.95e22 N m) and 2011 Tohoku (5.312e22 vs 5.31e22) to 0.1%.

From GCMT 1976–2020 (56,832 events, 45 yr), with the paper's own criterion that thrust means
`Mrr > 0`:

| quantity | this rebuild | the paper |
|---|---|---|
| global seismic power, Kanamori $E_s = M_0/2\times10^4$ | **11.8 GW** | ~10 GW |
| thrust Eg release, using Dahlen's 1000× | **9.4 TW** | ~10 TW |
| normal Eg gain, same ratio | **2.3 TW** | ~2 TW |

Total moment rate 7.42e21 N m/yr; thrust carries 80% of it. The single ratio $|\Delta E_g|/M_0$
needed is 0.050 from Dahlen's factor of 1000, and 0.076 from the paper's own Sumatra figure of
3e21 J. **Their three headline numbers are arithmetically sound and I could not shake them.**

## The finding: the conversion is placed where the energy is not (`depth.py`, `fig.py`)
The paper's mechanism requires the thrust earthquakes to convert 10–14 TW. Those earthquakes are
not distributed through the mantle. Moment-weighted, from the same catalogue:

| centroid depth | % of thrust moment |
|---|---|
| 0–30 km | 76.5 |
| 30–50 km | 14.0 |
| 50–70 km | 4.5 |
| below 70 km | 5.0 |

**95.0% of the thrust seismic moment is shallower than 70 km; the moment-weighted mean depth is
38.4 km.** So the proposed converter operates in the top ~2% of the mantle's depth.

Now the supply side, with no free parameters: 3.0 km²/yr of plate subducted (equal to the seafloor
creation rate), 100 km of thermal lithosphere, Δρ = 50 kg/m³ (≈500 K colder at α = 3e-5).
The gravitational power released in descending the top $z$ km is $\Delta\rho\, g\, z\, \dot V$:

| descent through | power |
|---|---|
| top 70 km | **0.33 TW** |
| top 100 km | 0.47 TW |
| top 660 km | 3.1 TW |
| whole mantle, 2890 km | **13.5 TW** |

The last line is the check that this arithmetic is right: it lands on Morgan et al.'s 11–14 TW,
which is a *whole-descending-limb* figure. **To have 10–14 TW you have to go to the base of the
mantle. To have thrust earthquakes you have to stay in the top 70 km, where 0.33 TW has been
released.** The paper takes the first number and spends it at the second depth, and offers no
mechanism for carrying the energy up.

And the number that is actually available there is the interesting one. Island-arc magma
production is 1.6–8 km³/yr (40 km³/km/Myr over ~4×10⁴ km of arc; or 0.6 ± 0.2 km³/yr extrusive at
6:1–13:1 intrusive:extrusive), and at 2800 kg/m³ × (1200 J/kg/K × 1200 K + 4×10⁵ J/kg) that is
**0.3–1.3 TW**. The gravitational power passing through the seismogenic depth range, 0.33 TW, is
already in that band. *The paradox may not need the extra factor of thirty.* That is a gift as
much as a correction: the mechanism can be right and the budget still be shallow.

Figure: `where_the_moment_is.png`.

## What I tried and could NOT establish — the heat-flow test failed its own control
The paper says the 10 TW "manifest[s] as the high heat-flow anomalies at subduction zones." That
is testable, so I tested it with the Global Heat Flow Database (R2024 v.2026.03, 89,583 usable
points) and located subduction by the paper's own thrust criterion. **Two instruments, both dead:**

1. `budget3.py` — annulus-by-annulus integral of the heat-flow excess above a far-field baseline.
   Thrust gave **6.02 TW**, normal gave **6.15 TW**. Indistinguishable, where the paper says +10
   and −2. It was measuring the continent/ocean contrast: the "excess" persists to 1000 km, and
   the far-field baseline (52 mW/m²) is mostly old continental shield.
2. `budget4.py` — oceanic points only, each normalised against the heat flow its own seafloor age
   predicts, age taken from the recorded water depth via GDH1. Median residual **−56 mW/m²**
   everywhere, i.e. the model over-predicts globally (hydrothermal advection, sediment loading,
   flexure — a depth-to-age inversion is not sound near a trench). Thrust and normal residuals
   again track each other at every distance (−77/−70, −69/−64, −50/−52 mW/m²).

So: **I do not know whether the surface heat flow refutes the 10 TW, and this database plus these
methods cannot tell me.** The one thing that survives is a like-for-like contrast within a single
broken instrument: near-field heat flow around normal-faulting (spreading) seismicity is as high
as or higher than around thrust seismicity, which is the opposite of what a +10/−2 TW split
predicts — but a contrast measured with a ruler that failed its control is a hint, not a result.

## What I don't know
- Whether Dahlen's reversibility argument is right. That is the physics the paper turns on and it
  is not something arithmetic settles. My finding is orthogonal: even granting the paper its Eg,
  the *location* of the conversion does not have the energy passing through it.
- Whether ΔEg computed by their Eq. (2) — a whole-Earth integral over the coseismic displacement
  field — can be localised at all. The paper needs it whole-Earth in magnitude and local in
  deposition, and those may not be the same quantity.
- Steady state: 10 TW is 21% of the global heat budget. If it is real and shallow, it has to come
  out somewhere, and I could not find where. See above — that is a failure of my instrument, not
  evidence either way.

## Files
`parse.py` (validated ndk reader) · `budget.py` (their three numbers) · `depth.py` (the finding)
· `fig.py` → `where_the_moment_is.png` · `heatflow.py`, `budget3.py`, `budget4.py` (the two
failed instruments, kept because a dead check is worth more than a silent one)

## Data
GCMT 1976–2020, globalcmt.org. GHFDB Release 2024 v.2026.03, GFZ Data Services,
doi:10.5880/fidgeo.2024.014 (CC-BY-4.0).
