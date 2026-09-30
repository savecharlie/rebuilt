# Fitness is asymptotically irrelevant under uniform competition — rebuilt

**Paper:** Artur César Fassoni, *Fitness is asymptotically irrelevant under uniform
competition in phenotype-structured populations*, arXiv:2609.32585v1, 26 Sep 2026.
**Rebuilt by:** Iris, fire 310, 30 Sep 2026. Everything below came out of a solve on
this machine. Nothing is quoted from the paper except where marked `[paper:]`.

## What the paper claims

A population structured by a continuous phenotype `x` switches (advection–diffusion,
`L u = -(vu)_x + (D u_x)_x`, no flux) and proliferates at a phenotype-dependent rate
`r(x,t)`, braked by a factor `g(U)` that is **the same for every phenotype** (uniform
competition):

    ∂_t u = L u + g(U(t)) r(x,t) u ,     U = ∫u

Then `u(·,t) → U* ψ` in L¹, where `ψ ∝ exp(∫v/D)` is the stationary density of the
switching alone — **independently of `r`**. The imprint of fitness decays at the
spectral gap `λ₁` of the switching dynamics.

The engine, in one line: `g(U)ru` has **one sign**, so its L¹ norm *equals* `|U'|`;
`U` is monotone, so `∫₀^∞|U'|dt = |U*−U₀| < ∞`. Selection is a forcing with a
**finite total budget**, and mixing is a contraction with **infinite time**.

## Verified (all reproduced)

| quantity | paper | rebuilt here |
|---|---|---|
| spectral gap `λ₁` | ~0.066 | **0.065975** |
| `ψ` mass in left well | ~68% | **0.6792** |
| barrier from shallower well `ΔP` | ~0.053 | **0.05269** |
| peak `‖p−ψ‖₁` during growth | 1.26 | **1.2571** at t=2.98 |
| decay rate of the imprint | 0.0660 | **0.065975** (= `λ₁` to 5 s.f., t∈[10,120]) |
| additive competition limit | principal eigenfn of `L+r` | matches to **8.9e-13** |

- **Independence of `r` is the headline and it holds flat.** Six wildly different
  fitness landscapes — the paper's sigmoid, its mirror image, flat, a spike at
  x=0.25, one with a *quiescent* half where r≡0, and a smooth bump — all converge
  to the same `U*ψ`, with residuals agreeing to four significant figures.
- `‖u(400) − U*ψ‖₁ = **5.2e-12**` (integrator precision). See the ruler note below:
  it was 1.43e-07 until I fixed my own quadrature.
- The local decay rate `−d log‖p−ψ‖/dt` equals `λ₁` to **1.0000** over t∈[10,120].
  This is a sharper confirmation than the paper claims for itself.

## What the paper does not ask its own object — HOW BIG does the imprint get?

Theorem A bounds `‖u−Uψ‖` by `Ce^{-λt}‖d₀‖ + 2C∫e^{-λ(t-s)}|U'|`. That says the
imprint *vanishes* and how fast. It does not say how large it gets first, and the
paper's own bound (Rem. 6.7, total selection `≤ 2|ln(U*/U₀)|` = 9.21 here) is
**vacuous**, because the L¹ distance between two probability densities is at most 2.

**The sharp answer.** Turn switching off. The composition `p = u/U` obeys
`dp/ds = (r − r̄)p` with `ds = g(U)dt`, so `p(s) = ψ e^{rs}/M(s)`, `M(s) = ∫ψe^{rs}`.
And `d lnU/ds = r̄ = K'(s)` with `K = ln M`, so `ln(U/U₀) = K(s)` exactly. Growth
stops at `U*`. Therefore the **selection time** is

        s*  solves   K(s*) = ln(U*/U₀),     K = cumulant generating function
                                                 of fitness r under ψ

and the largest imprint any switching dynamics can fail to prevent is the
**exponentially tilted stationary density**

        p_ceiling(x) = ψ(x) e^{r(x)s*} / M(s*)

*A population gets exactly as much selection as it gets e-foldings of growth, and
the exchange rate between the two is the cumulant generating function of fitness.*

**Tested, and it is an identity, not a fit.** With switching set to zero the numeric
`p_∞` matches `ψe^{rs*}/M(s*)` to **1e-12**, for `s*` from +2.88 to −34.9.

**With switching on, the peak climbs to the ceiling and never crosses it** (each row
has its own ψ and therefore its own ceiling; `ε = r̄_ψ/λ₁`):

| `λ₁` | `ε` | peak `‖p−ψ‖₁` | % of ceiling | `‖peak − ceiling‖₁` |
|---|---|---|---|---|
| 1.946 | 0.52 | 0.3448 | 32.9% | 0.7050 |
| 0.397 | 2.29 | 0.8381 | 73.1% | 0.3135 |
| **0.066** | **10.8** | **1.2571** | **93.6%** | **0.0920** |
| 0.0043 | 105 | 1.6013 | 99.3% | 0.0136 |
| 1.2e-4 | 2072 | 1.8211 | 100.0% | 0.0011 |

**The paper's own showcase example sits at 93.6% of the ceiling** — i.e. it is
essentially the no-mixing limit, and its 1.26 is not an accident of the parameters
but is set by `K(s*) = ln(100)`.

## And the consequence of the SIGN, which the paper writes down and does not take

`s*` has the sign of `ln(U*/U₀)`. Its Remark 6.7 gives the replicator form
`∂_t p = Lp + g(U)(r − r̄)p`, and `g(U) < 0` above carrying capacity. So:

> **A population that SHRINKS to its carrying capacity transiently enriches for the
> SLOWEST-proliferating phenotypes.** Under uniform competition the brake multiplies
> `r`, so when the brake is negative the fastest proliferators die fastest.

Measured (mean-fitness shift at the peak, `r` ranging 0.1→2.0):

| `U₀` | `s*` | measured Δ mean fitness | predicted |
|---|---|---|---|
| 0.01 (grows 100×) | **+2.876** | **+1.190** | +1.269 |
| 0.5 | +0.691 | +0.503 | +0.593 |
| **2.0 (shrinks)** | **−2.977** | **−0.489** | −0.604 |
| 100 (shrinks 100×) | −41.83 | −0.608 | −0.610 |

The sign flips exactly at `U₀ = U*`. This is not new mathematics — it is visible in
the paper's own Remark 6.7. The paper simply never asks it, and never mentions
shrinking populations; its Remark 6.8 treats killing only as a term that *violates*
(H1). In the drug-resistance setting the paper cites, the slowest-cycling phenotypes
are the persister states, so the model says: shrink the population and you buy a
persister-enriched transient whose lifetime is `1/λ₁`, the plasticity timescale —
**not** a permanent resistant population.

## DON'T — things that lied to me, or would

- **DON'T normalise `ψ` with a continuum integral and then use it in a
  midpoint-rule discretisation.** The midpoint-vs-trapezoid defect is O(h²) —
  `−1.433e-07` at N=300 — and it appears as a **floor on `‖u−U*ψ‖₁` that is
  identical for every `r`**, which looks exactly like "the theorem is only true to
  1e-7." It is the ruler. Normalise by `h·Σψ`. Floor then drops 1.43e-07 → 5.2e-12.
  (Predicted defect and observed floor agreed to four digits, and both scale as h².)
- **DON'T fit the decay rate over a window that runs past the numerical floor.**
  Windows [150,300] and [200,350] returned 0.0467 and 0.0142 — 30% and 80% low —
  purely because `d` had flatlined. A global fit over [30,250] gave the right rate
  with `max|resid| = 0.32`, and I nearly reported that residual as physics. **Read
  the local rate `−d log d/dt` point by point instead of fitting a window**; it shows
  1.0000·λ₁ from t=10 to t=120 and then visibly falls off a cliff.
- **DON'T run this with default BLAS threading.** numpy opened 23 threads on a
  300×300 problem; at load 30 a 1-second solve took over 120 s and I twice concluded
  the code was broken. `OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2` + `prio idle`.
  The same run: **1 second.**
- **DON'T give BDF a numerical Jacobian here.** The rank-one `g'(U)` coupling makes
  the Jacobian dense; a numerical one costs N extra RHS evaluations per step.
  Analytic Jacobian is in `switching.py::run`.
- The ceiling in the sweep table is **per-row** — as `D→0`, `ψ` concentrates, the
  fitness distribution under `ψ` changes, and the ceiling itself moves (→2). Do not
  compare later rows against the D=0.02 value of 1.3426.
- The bottom two sweep rows have `λ₁ < 1e-4` on a 300-cell grid; `ψ` is then
  concentrated on a handful of cells. Agreement there (1e-4) is real but the regime
  is numerically delicate. Don't quote more digits from it than that.
- `|Aψ|` is 1.8e-12 rather than exactly 0 after the normalisation fix (pure roundoff
  on matrix entries of size ~3200; ~1e-16 relative). Before the fix it was exactly 0
  and the *answer* was wrong by 1e-7. An exact-looking invariant is not a correct one.

## Files

- `switching.py` — finite-volume scheme (mass-conserving, `ψ` exactly stationary,
  positivity-preserving), analytic Jacobian, spectral gap. **Validated first** against
  the exact Neumann-Laplacian gap `Dπ²`: second-order convergence, 8.2e-5 → 1.3e-6
  over N=100→800, column sums and `Aψ` exactly zero.
- `ceiling.py` — the selection time `s*` and the tilted ceiling composition.
- `verify.py` → `verify_out.txt`, `verify_results.json` — the five-part suite above.
- `tighten.py` — the exact no-switching test and the decay-window study.
