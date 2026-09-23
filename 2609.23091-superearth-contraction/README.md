# Super-Earth interiors shrink by about 10% as they crystallise — rebuilt

**Paper:** [arXiv:2609.23091](https://arxiv.org/abs/2609.23091) — Lichtenberg, Attia,
Nicholls, Sastre, Bower, Stuitje, Pascal & Soucasse (2026-09-19).

**Verdict: the paper's arithmetic holds and is internally consistent. I did not
reproduce their number — their equation of state is a companion paper not yet
public. What I did instead was translate their headline into the quantity that
actually does the work, which they do not state, and which makes their own open
question much sharper.**

---

## The claim

From a standalone structure sweep at 1 M⊕ with an Earth-like core-mass fraction,
stepping the adiabat from molten to fully crystallised: the interior radius falls
from about **7150 km to about 6340 km**, a contraction of **~11 %**, with the
mantle shell thinning ~20 % and the core radius shrinking ~3 %. Across a coupled
1–10 M⊕ grid the contraction falls from ~11 % to ~9.5 %. They note this is
roughly **twice** the ~5 % Bower et al. (2019) found for an Earth-mass magma
ocean, and suggest the silicate equation of state is the likely origin.

## What is and is not reproducible

The evolution model is five coupled codes, and the mantle equation of state
(PALEOS, Attia et al. 2026) is *submitted, not yet available*. So there is no
honest way for an outsider to reproduce 11 %.

The **structure problem** is another matter: it is fully specified in their
Eq. 1, and the three numbers they quote are enough to answer a question they
never ask.

## Validating the instrument first

The paper's own Appendix A verifies its solver against an n = 1 polytrope and a
two-layer constant-density sphere. I used the same two exact tests, plus Earth.

| test | mine | exact / PREM | rel. |
|---|---|---|---|
| n = 1 polytrope, R | 6371.0531 km | 6371.0000 km | 8.3 × 10⁻⁶ |
| n = 1 polytrope, P_c | 139.6964 GPa | 139.6967 GPa | 2.3 × 10⁻⁶ |
| n = 1 polytrope, ρ(r) | — | — | max 5.1 × 10⁻⁶ |
| 2-layer sphere, R | 6512.7210 km | 6512.7211 km | 2.1 × 10⁻⁸ |
| 2-layer sphere, P_c | 380.4287 GPa | 380.4290 GPa | 5.8 × 10⁻⁷ |
| Earth, **pure ε-Fe core**, P_c | 422.5 GPa | 363.9 GPa | **0.161 — FAIL** |
| Earth, core × 0.90, R_cmb | 3469.8 km | 3480 km | 2.9 × 10⁻³ |
| Earth, core × 0.90, P_cmb | 135.07 GPa | 135.8 GPa | 5.4 × 10⁻³ |
| Earth, core × 0.90, P_c | 365.36 GPa | 363.9 GPa | 4.0 × 10⁻³ |

The polytrope radius lands 8 × 10⁻⁶ from the exact Lane–Emden answer; the paper
reports 1.4 × 10⁻³ for its own full solver on the same test.

**The Earth failure is kept in the suite rather than tuned away, because it is
the diagnosis.** A pure ε-Fe core is too dense to be Earth's: the core comes out
too small and the central pressure 16 % too high. Applying the one published
correction at its accepted value — Earth's core carries ~10 wt % light elements,
so ~10 % less dense than pure iron — moves **R_cmb, P_cmb and P_c each to within
0.5 % of PREM at once**. Three quantities, one parameter, nothing tuned. R stays
1.2 % small because my mantle is pure bridgmanite with none of the lighter
upper-mantle phases, which is the expected direction.

Free bonus check: the solid reference planets give R(10 M⊕)/R(1 M⊕) = 1.859
against the standard rocky mass–radius scaling M^0.27 → 1.862.

Run it: `python3 validate.py`.

## Finding 1 — eleven per cent of radius is fifty per cent of density

No equation of state is needed. Mantle mass is conserved, so its mean density is
its mass over its shell volume, and

    rho_solid / rho_molten = (R_m^3 - Rc_m^3) / (R_s^3 - Rc_s^3)

Put their own three numbers in (`geometry.py`):

| | radius contraction | ⟺ mean mantle densification |
|---|---|---|
| this paper, 1 M⊕ | 11.3 % | **50 %** |
| this paper, 10 M⊕ | 9.5 % | 40 % |
| Bower et al. (2019) | 5.0 % | 18 % |

Varying the solid core radius from 2800 to 4000 km moves the 1 M⊕ figure only
from 1.47 to 1.55, so the conclusion does not depend on a core radius I had to
assume. Their numbers are also mutually consistent: the stated "~20 % shell
thinning" independently picks out a solid core radius of 3538 km, against
Earth's 3480 km.

**So the factor of two they note between themselves and Bower et al. is a factor
of 2.8 in the quantity that produces it.** Radius is a compressed readout of
density: at fixed core, d(contraction)/d(densification) = 0.2752 measured, against
the thin-shell prediction (1 − (R_c/R)³)/3 = 0.2782 — two roads, 1.1 % apart.

## Finding 2 — that 50 % is not a melt–solid contrast

The published **isochemical** melt–solid density contrast is nowhere near 50 %.
It falls from roughly 20 % in the shallow mantle to about **4 %** at core–mantle
boundary pressure (Karki et al. 2018, GRL 45), and Petitgirard et al. (2015, PNAS
112) put MgSiO₃ glass just **1.6 %** below bridgmanite at 133 GPa. Most of a
mantle's mass sits where the contrast is smallest.

These are different quantities, and the gap between them is a feedback: a less
dense mantle inflates the planet, interior pressures drop, and the mantle
decompresses further. `contraction.py` measures that amplification instead of
arguing about it — it takes the validated solid planet, makes the mantle a
uniform fraction δ less dense *at the same pressure and temperature*, which is
exactly what the melt literature reports, and re-solves the structure.

Sweeping δ at four masses (`contraction.py`, `run_contraction.log`):

| δ (same-pressure deficit) | 1 M⊕ | 3 M⊕ | 5 M⊕ | 10 M⊕ |
|---|---|---|---|---|
| 5 % | 1.63 % | 1.83 % | 1.94 % | 2.15 % |
| 10 % | 3.33 % | 3.72 % | 3.95 % | 4.36 % |
| 20 % | 6.94 % | 7.71 % | 8.18 % | 8.98 % |
| 30 % | 10.89 % | 12.03 % | 12.73 % | 13.92 % |

Inverting for their number:

| | target | required δ | R_molten | mean densification | amplification |
|---|---|---|---|---|---|
| 1 M⊕ | 11.3 % | **30.98 %** | 7098 km | **50.8 %** | 1.64× |
| 10 M⊕ | 9.5 % | 21.08 % | 12934 km | 38.9 % | 1.84× |

Two things to notice. The mean densification, 50.8 % and 38.9 %, lands on the
EoS-free arithmetic of Finding 1 (50.0 % and 39.9 %) to about one point, by a
completely different route. And R_molten = 7098 km sits 0.7 % from their
7150 km, from an equation of state that shares nothing with theirs.

**So their 11.3 % requires a mantle ~31 % less dense in the molten state than
in the solid state at the same pressure.** The published contrast is 20 %
falling to 4 %.

### The forward direction

`literature_melt.py` builds δ(P) = δ₀·exp(−P/P\*) from the two published
anchors and solves the structure with it directly, so no weighting argument is
needed:

| δ(0) | δ(136 GPa) | R_molten | contraction |
|---|---|---|---|
| 0.20 | 0.040 | 6568.8 km | **4.16 %** |
| 0.20 | 0.016 | 6507.5 km | 3.25 % |
| 0.25 | 0.040 | 6630.6 km | 5.05 % |
| 0.15 | 0.040 | 6509.7 km | 3.29 % |
| 0.30 | 0.060 | 6747.3 km | 6.69 % |
| 0.35 | 0.080 | 6875.2 km | 8.43 % |

The central case gives **4.16 %** — which is essentially the ~5 % of Bower et
al. (2019) that this paper is twice. Even the most generous row, a contrast
twice the published one at both ends, reaches only 75 % of 11.3 %.

The paper attributes a further 1–2 points of its 11.3 % to thermal contraction
of the superheated melt; that is their accounting, used as theirs, and is not
in the table above.

## Finding 3 — one sentence runs backwards

The paper reports the contraction falling with mass (11.3 % → 9.5 % over
1 → 10 M⊕) and explains it as *"a more strongly compressed massive interior
returns a smaller fractional radius change for the same melt-to-solid density
contrast."*

Holding the contrast fixed is exactly what a uniform δ does, so that is a
controlled experiment I can run. **It goes the other way**: at δ = 30.98 %, the
contraction is 11.30 % at 1 M⊕ and 14.42 % at 10 M⊕, a fitted mass exponent of
**+0.105** against their **−0.075**.

Part is geometry: a heavier rocky planet carries a proportionally smaller core
(R_cmb/R = 0.509 at 10 M⊕ against 0.553 at 1 M⊕), so the thin-shell coefficient
rises from 0.277 to 0.289. That is 4 % of a 31 % effect. The rest is the
compression feedback, which strengthens with mass.

`robustness.py` tries nine ways to flip the sign — a core geotherm three times
steeper, a core based at 6000 K, mantle adiabat exponents 0.25 and 0.45,
core-mass fractions 0.20 and 0.50, and mantle K₀′ of 3.5 and 5.0:

```
  9/9 variants give a RISING trend; range 1.249 to 1.374.
  the paper's 0.841 is outside that range by 49 % at the closest variant.
```

**The trend they report is real; the reason offered for it appears inverted.**
At a fixed contrast the massive interior returns a *larger* fractional radius
change. What must actually fall with mass is the contrast itself, squeezed by
the higher pressures — which is the same physics that makes the contrast small
in the first place, and it has to overcome a ~30 % geometric and compressive
gain to show up at all. That is checkable directly in their own PALEOS tables.

## What `dull.py` caught, and what it missed

It flagged 1152 floating-point underflows in `literature_melt.py:44`. Two
minutes of checking showed they come from the shooting bracket probing absurd
central pressures that are then rejected; converged solutions sit 113–180×
below where the underflow begins. A false alarm, cleared by looking.

It missed the one that mattered. I wrote the mean-density ratio the wrong way
up and printed a mantle getting 50 % **less** dense as it froze. Every check was
green — a well-conditioned wrong formula is precisely the blind spot that
module's own test suite asserts. What caught it was the sign looking strange and
`geometry.py` disagreeing. The second road is still the only road.

## Files

| | |
|---|---|
| `structure.py` | hydrostatic structure in mass coordinate; Vinet + Mie-Grüneisen-Debye |
| `validate.py` | polytrope, constant-density sphere, Earth — the known answers |
| `geometry.py` | the EoS-free translation of their headline |
| `contraction.py` | δ → contraction, 1–10 M⊕, and the inversion |

All four run under [`../tools/dull.py`](../tools/README.md).

## Sources
- ε-Fe Vinet: Dewaele et al. 2006, PRL **97**, 215504.
- MgSiO₃ bridgmanite: Stixrude & Lithgow-Bertelloni 2011, GJI **184**, 1180.
- PREM: Dziewonski & Anderson 1981, PEPI **25**, 297.
- Melt–solid contrast: Karki et al. 2018, GRL **45**; Petitgirard et al. 2015,
  PNAS **112**, 14186.
