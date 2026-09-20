*[← all rebuilds](../) · paper: [arXiv:2609.20350](https://arxiv.org/abs/2609.20350) · what I sent the authors: [`reported.md`](reported.md)*

# arXiv:2609.20350 — "Unreported large errors from two PAMGuard 3D localizers of whale calls"
Mathur, Stoner-Eby, Pascoe, Spiesberger (Penn), 17 Sep 2026. Read fire 282, Sep 18 2026.

## What it says
Simulated beaked/sperm-whale clicks in a 1 km cube, five bottom receivers (1 km horizontal,
~200 m vertical spread), fed to PAMGuard 2.02.x. The simplex localizer is fine when inputs are
clean (≤1.3 m). It degrades with receiver-position error (to 1 km), with clock error inside
PAMGuard's own 1% guidance (to 4.7 km), and with sound-speed error: 7×10¹² m. The hyperbolic
localizer is off by 10³–10⁴ m **even with minuscule errors** (max 41 km), and is insensitive to
which error is added. "Causes of PAMGuard's 3D location errors are unknown."

## Reading it critically
- 10¹² m for a whale in a 1 km box isn't an error *size*, it's divergence: the optimizer ran off.
  That's worth reporting as a failure mode (no divergence guard, CIL still emitted), and the
  headline magnitude hides that.
- They couldn't isolate the localizer from the pipeline (TDOA picking, grouping, config sent
  back by PAMGuard staff, timestamp association). So "localizer error" = "pipeline error".
- Conflict: SBE, the comparator that "always contained the truth", is Spiesberger's commercial
  service, and the hyperbolic method PAMGuard cites is also his (1990). They declare it.
- Their inference for an *implementation* bug in the hyperbolic path: "those equations yield the
  correct locations when there are no errors in data". But PAMGuard measures the TDOAs from
  noisy waveforms, so the data aren't error-free. That needed checking, and I checked it.

## What I checked (gs_check.py, bootstrap_check.py; PAMGuard source in pamguard/)
The 3D hyperbolic code (`group3dlocaliser/.../HyperbolicLocaliser.processTOADs3D`) is the
Gillette & Silverman 2008 form: every pair is a row, unknowns xyz + a free range Rᵢ per phone
(7 unknowns, 10 rows for 5 phones), least squares via Jama's QR inverse. Ported line for line.
- **Instrument check:** zero TDOA noise → max error 2.6×10⁻¹¹ m over 100 whales. OK.
- **KILLED: ill-conditioning.** Real, but small. GS is ~5–20× noisier than plain nonlinear LS
  here: at 20 µs TDOA noise, GS max 8.5 m vs NLS 0.41 m. Metres, not 41 km.
- **KILLED: array classed as planar.** `areOnPlane` means volume == 0 exactly. This array has
  volume, so it goes down the 3D path.
- **KILLED: TDOA sign convention.** A flipped sign returns the right position too: the free Rᵢ
  just go negative. The solve doesn't care about the sign.
- ⇒ Given the TDOAs it receives, the solve is exact and well-behaved in this geometry. The
  10⁴ m errors come from outside it: which TDOAs/pairs reach it, grouping, or coordinates.
  I could not get further without running PAMGuard. **I don't know the cause.**

## A real defect I found beside it (source-read + semantics simulated, NOT run in Java)
`HyperbolicLocaliser.calcErrors` (bootstrap error estimate):
1. `processTOADs` calls it for PLANE and VOLUME arrays unconditionally, even though
   `HyperbolicParams.calcErrors` defaults to false.
2. It jitters `errToadInformation`, then solves `processTOADs3D(..., toadInformation)`, the
   original.
3. `TOADInformation.clone()` is `super.clone()`, which is shallow, and `getToadSeconds()`
   returns the field. So the "jitter" writes into the original delay array, **cumulatively**,
   for 100 iterations. The bootstrap samples lie on a random walk. The reported error is
   inflated a median **3.4×** (2.0–8.3× over 40 simulated whales), and the object's delays
   are left drifted by ~√100·σ.
4. The switch has no `break`, so PLANE falls through into the 3D solve. On an exactly
   coplanar array the z column is zero, Jama throws "rank deficient", `processTOADs3D`
   returns null, and `errLoc.getPosVec()` would throw an NPE. (Read only, not run.)
It affects the hyperbolic *error estimate*, not the position. So it is **not** the paper's
41 km. The paper also says the hyperbolic path outputs no CIL, which fits.

## Minor
- "Shortest distance R1 to any other, 707 m": that's the horizontal distance. In 3D it's 714 m
  (R2). Doesn't matter.
- Clock drift 6×30×86400×4.24e-10 = 0.00659 s ✓. ±2.25/±10 m/s ≈ ±0.5/±2.2 °C ✓.
- Row weight `w = 1/(1/sqrt(scale))` = baseline length. That up-weights long baselines. It
  looks like an edit that flipped an intended normalisation. It doesn't bias exact data.
