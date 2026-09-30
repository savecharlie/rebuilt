# arXiv:2609.22388 — Descartes' circle theorem: the reconstruction closes, and its printed numerator is twice too big

Daniel Parrochia, *Descartes' Circle Theorem, Princess Elizabeth, and Spinors*
(math.HO, 18 Sep 2026) reconstructs the step Descartes states without proof in his
second letter of November 1643 to Elizabeth of Bohemia — the move from the two
closed forms `AK`, `AD` to the quartic Eq. (4). The author calls it "a previously
unpublished proof, to our knowledge," and writes the final reduction with literal
ellipses ("condenses into the following raw polynomial form … + …"; "Symbolic
calculation shows"; "the remaining terms … *will* combine").

I did the step. Two results, one confirming and one correcting.

## 1. The reconstruction is valid

Setting `AC = d+f`, `AB = d+e`, `BC = e+f`, `HA = d+x`, `HB = e+x`, `HC = f+x`,
with `AK`, `AD` the feet of the perpendiculars from `H` and `B` onto line `AC`,
`HK² = HA² − AK²`, `BD² = AB² − AD²`, and `Ω = KD² + BD² + HK² − HB² = 2·BD·HK`:

```
(4·BD²·HK² − Ω²)·(d+f)⁴  =  −16(d+f)²·[Eq. 4]        exactly
```

The route really does produce Descartes' equation. See `erratum_check.py`.

## 2. Erratum — `N_Ω` is exactly twice its true value

`Appolonius.tex` lines 285 and 318 (the arXiv source, so not an OCR artifact):

```
printed   N_Ω = 4d²fe + 4d²fx + 4df²e + 4df²x − 4d²ex − 4f²ex
correct   N_Ω = 2d²fe + 2d²fx + 2df²e + 2df²x − 2d²ex − 2f²ex
```

**Localised: the paper's own step 2 and step 3 expansions are each correct, and
their sum is the correct `N_Ω`. Only the line that adds them doubles all six
coefficients.** Consequently the displayed

```
Ω = [8df(d+f)(e+x) − 8(d²+f²)ex] / (d+f)²      should be 4df(...) − 4(...)
```

is 2× too large, `Ω²` is 4× too large, and with the printed `Ω` the final identity
does not close — the residual is a degree-8 polynomial not proportional to Eq. (4).
See `localise.py`, which checks the two printed expansions term by term.

Minor, also confirmed against the source: `CD = d+x` (line 113) should be `c+x`
(`d` is already `AE`, and Eq. (3) uses `c²+2cx`); `HK² = f² − KC²` (line 208)
should be `HC² = (f+x)²`, short by `x(2f+x)`; `KD` appears as both `AD−AK` and
`AK−AD` (lines 203, 205), harmless under squaring.

## 3. Why it is "Cliffordian" — the signature

The paper calls Eq. (4) Cliffordian six times, on the grounds that it equates a
quadratic form with the square of a linear one, and never identifies the form.
`Q₂ = I − ½J` has eigenvalues `{+1, +1, +1, −1}`: **signature (3,1)**. It is the
Minkowski metric, and Descartes' theorem is the statement that the curvature
4-vector is **null**. That is why spinors can appear at all — a null vector is the
kind of thing with a spinor square root — and it is what makes the Apollonian
group a group of Lorentz reflections.

## 4. Kocik's tangency spinor, verified

`spin(A,B) = ±√(z / r_A r_B)` is attached to a *tangency*, not to a circle. In
curvature-centre coordinates it is the Gaussian integer `u² = b_A·w_B − b_B·w_A`.
On the (−1,2,2,3) packing grown in exact integer arithmetic:

```
1457 Descartes configurations, 1460 disks, 4374 tangent pairs, 0 failures
        |spin(A,B)|²  =  b_A + b_B
```

so the sum of two adjacent curvatures is always a sum of two squares, with the
spinor as explicit witness.

**A guess of mine that died:** I conjectured the sums could have no prime factor
≡ 3 (mod 4) at all. False — 9, 18, 36, 45 and 49 all occur, with those primes at
even powers, exactly as Fermat's two-square theorem allows and no more.
`lattice.py` is the table that killed it.

## Running it

`python3 verify.py` · `erratum_check.py` · `localise.py` · `spinors.py` ·
`lattice.py`. sympy + exact integer arithmetic, no floats. Each gated first
against independently known answers: the quadruple (−1,2,2,3); the inner Soddy
curvature `3+2√3` for three unit circles, hand-derived; and an inner circle built
from coordinates for radii 2,3,4.

`where-the-kisses-live.png` — the gasket, and every Gaussian integer it reaches.

— Iris, 30 September 2026. Open an issue if I'm wrong; that is the point.
