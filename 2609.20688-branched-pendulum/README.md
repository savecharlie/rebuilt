# arXiv:2609.20688 — Toda & Ooshida, *Experimental detection of energy transfer into the antiphase mode in a branched double pendulum*

Read Sep 22 2026. Rebuild in this directory. Essay: `../essays/2609.20688-branched-pendulum.md`.

## What it is
A parent link hangs from a pivot, and two identical child links hang from the parent. Swapping the children
is a symmetry, so the linear modes split into two in-phase modes (f+ 0.88, f− 1.58 Hz) and one anti-phase
mode (f_d 1.26 Hz). The anti-phase mode lives only in the children, and its frequency is simply
sqrt(G2/M2): a child swinging from a fixed pivot. Release the parent with the children hanging still and
the anti-phase mode starts at zero, so it can only get energy from the other two. The authors build this
as a model of a turbulent cascade you can actually hold: shell models can't be built as an apparatus, and
a pendulum can. They derive an energy-transfer function W = {H_d, H_m} (Eqs. 12–13) and evaluate it on
video data. In the experiment θ_D grows for θ1(0) ≥ 47.8° and doesn't for ≤ 47.6°. The critical angle
is left for "numerical and theoretical studies ... planned."

## What I checked (all runs; see this directory)
- **System ID holds.** From their lengths, masses and fitted I1, I2 (g = 9.80): f+ 0.873, f− 1.582,
  f_d 1.258 Hz, against their 0.88 / 1.58 / 1.26 (`params.py`). My nonlinear equations of motion
  (Euler–Lagrange from their Eq. 3–5) conserve energy to 6e-11 over 20 s, and at small amplitude their
  spectra land on the same three frequencies (`validate.py`).
- **Eq. 13 is not the leading-order transfer.** Test: shrink the initial disagreement ε from 1e-2 to 1e-5.
  A correct leading-order W must converge to the exact dH_d/dt. Theirs doesn't. The relative rms error
  stays fixed at 71% (θ1(0) = 20°) and 52% (47.8°) for every ε, with best-fit scale 0.67–0.68 every time
  (`wcheck.py`). So the error isn't truncation; a term is missing.
  - **Why:** their expansion of K (Eq. 11) keeps the term linear in θ_D and files the rest as
    "O(p²θ_D²)." But p1 and p_M are O(1), and only p_D is O(θ_D). So p1²θ_D² is the same order as H_d
    itself, and its θ_D-derivative feeds W at leading order. The bracket also has the cross term's
    **other half**: they kept −(∂H_d/∂p_D)(∂K_x/∂θ_D) (that's W_K) and dropped +(∂H_d/∂θ_D)(∂K_x/∂p_D).
  - **Closed form (mine), a Schur-complement expansion of M̃⁻¹:**
    W = W_P + W_K + W_X + W_I, with
    W_X = 2 G2 μ sinΔ_m θ̇1 θ_D² / M2 and
    W_I = −(p_D θ_D / M2)[μ cosΔ_m θ̇1 θ̇_M + (μ²/M2) sin²Δ_m θ̇1²].
    Its error falls as ε² (1.4e-4 → 1.4e-6, then a 1e-7 floor from the numerical derivative)
    (`wclosed.py`). The missing pieces are big: at 47.8°, rms W_X ≈ 0.34 and W_I ≈ 0.21 of the true
    transfer, against W_P 1.05 and W_K 0.42.
  - **What survives:** their qualitative story. W is positive while H_d grows, and W_P is the largest
    single piece at 47.8°, "dominant although not overwhelming," as they say. But the W curve in their
    Fig. 3b/4 is not dH_d/dt, even early on, when θ_D is small.
  - Minor: the text below Eq. 10 says "G_1 sin²θ_D"; the equation (correctly) has G_2.
- **Where the anti-phase mode switches on, in the ideal model** (`rate.py`, ε = 1e-9, 40 s windows):
  the growth rate is 0.000 ± 0.001 /s from 30° to 42°. It turns positive near 42.5–44° (≈0.1–0.25 /s)
  and jumps to 1–5 /s between 46.5° and 48.25°. Fast enough to show up within seconds from a 0.01 rad
  start is **≈46.5°**, near their experimental 47.6/47.8°. The rate above threshold is jagged, not a
  smooth curve.
- **Their 47.8° run, in the model** (ε = 0.01): θ_D first swings to ≈0.4 rad at t ≈ 2.5–3.3 s, and
  H_d reaches ≈7e-3 J. They measured 0.4 rad at 3.8 s and ≈8e-3 J. Same size, about a second earlier.
  (Their children started from video-noise disagreement, not a clean 0.01 rad, and the real pendulum
  has friction.)
- **Mechanism: the onset coincides with chaos in the in-phase motion.** First I validated the ruler: a
  finite-time Lyapunov exponent over 60 s reads +0.05 to +0.13 /s even on regular orbits (linear shear,
  ~ln t / t), so it can't tell them apart. Over 240 s the regular cases fall (5°: 0.050→0.021; 40°: 0.106→0.031)
  and the chaotic ones hold or rise. In-phase-only motion (θ_D ≡ 0) is regular up to 42° (0.037) and chaotic
  from 43° (0.59), with 42.5° on the edge (0.071). Anti-phase growth (ε = 1e-9, 40 s) is 0.005–0.008 /s up to
  42.25° and 0.083 /s at 42.5°. **Both switch on between 42.25° and 42.5°.** Reading: a regular in-phase
  motion has a line spectrum that misses the 2f_d pumping resonance (their f+ + f− = 2.46 < 2.52), while a
  chaotic one is broadband and always contains it. That's a sharper form of the authors' "low-frequency
  components fill the mismatch." It's a coincidence of two thresholds along one family of initial conditions,
  not a proof. (`chaos.py`, `chaos2.py`, `lyap_*.log`)

## Taste
The best thing in the paper is its opening move: turbulence you can hold in your hands. The weak thing is
the one they called the "manageable expression." It became manageable by dropping terms of the same size.

## Files
- `params.py`: parameters and linear frequencies · `sim.py`: nonlinear EOM, H_d and the paper's W
- `validate.py`: energy and small-amplitude checks · `wcheck.py`: the ε-convergence test of Eq. 13
- `wclosed.py`: the corrected leading-order W · `rate.py`, `chaos.py`, `chaos2.py`: thresholds
- `figure.py` → `branched_pendulum_rebuild.png` · `reported.md`: what I sent the authors
