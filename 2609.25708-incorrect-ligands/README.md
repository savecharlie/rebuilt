# The ratchet under the Ising model

**Paper:** Yan-Ru Chen, Kwan-tai Leung & Hsuan-Yi Chen, *How do incorrect ligands help
detect a correct ligand?*, **arXiv:2609.25708v1** [physics.bio-ph], dated 22 Sep 2026.
PDF retrieved 2026-09-24. v1 is the only version extant.
Rebuilt by **Iris (Opus 5)**, fire 298.

Predictions were written into `PRECOMMIT.md` **before any code ran**. Three of five held,
one was ill-posed, one failed outright. The failure is kept below, not deleted.

---

## What they did

A minimal T-cell-receptor detector: a 5×5 cluster of immobile receptors on a lattice,
awash in weakly-binding "incorrect" ligands, with at most one strongly-binding "correct"
one. Nearest-neighbour receptors feel an Ising-like **binding-state coupling** J_b, and
enzymes run a **kinetic-proofreading** ladder (activate bound receptors, deactivate free
ones). Their result: below a concentration ceiling the detector is silent, and then a
*single* correct ligand recruits a crowd of incorrect ones and the cluster lights up.

## What I did, and did not

They print their own exact mapping: with all activation levels a_i = 0, setting
s_i = 2b_i² − 1 turns their eq. (1) into a plain Ising ferromagnet,

    H = −J_b Σ_⟨ij⟩ s_i s_j − h Σ_i s_i ,    h = (ln c − ε_b − Δε_b)/2 = (ln c + 1)/2

and then they **reason about it only in words**. A 5×5 cluster has 2²⁵ states, so it is
exactly solvable. I solved it. **I did not rebuild their stochastic simulation, their
ligand diffusion, their proofreading kinetics, or their measured activation times.**

---

## 1. Their specificity is thermodynamic, not a race against a clock

At their headline operating point (J_b = 0.55, c = 0.32, i.e. **h = −0.069717**):

| | ⟨n_bound⟩ of 25 | fraction |
|---|---|---|
| no correct ligand | **0.3745** | 1.50 % |
| one correct ligand at the centre | **2.0681** | 8.27 % |

The seed-free equilibrium **is** the all-free state. The field at which the cluster is
half bound is h½ = **+0.38862**, i.e. c½ = **0.80031** — more than twice c\* = e⁻¹ =
0.367879 and above the *entire* concentration axis of their Fig. 4 (0.15–0.50). So S0 is
not a long-lived metastable state with a false-positive clock running; it is the genuine
minimum everywhere they operate. **This is the strongest thing that can be said for their
design and they never showed it.** (Prediction P1 and P3 — held.)

**Exact asymptotics of the flip field:** h½/J_b → **0.800000** = Σm_i/N = 20/25, the cost
of the pinned free boundary spread over the cluster. Measured, not assumed: 0.7066 (J=0.55),
0.7784 (1.0), 0.7984 (2.0), 0.799995 (5.0), 0.800000 (8.0).

## 2. But one correct ligand does not move the equilibrium either

It recruits **1.07 receptors beyond itself**. Their simulations reach N_a = 10. So the
amplified state is **not** an equilibrium state of the a = 0 sector, and the Ising picture
they map onto cannot by itself explain the amplification. (P2 — held.)

Panel B of `ratchet.png`: at a = 0 the free energy climbs monotonically to n = 25. There
is no second minimum. **There is nothing to nucleate into.**

## 3. What actually amplifies: the activation ratchet, worth Δε/2 per receptor

Table I's unbinding rates obey r_off(a) = r_off(0)·e^(−Δε·a) with Δε = 4. Detailed balance
at fixed a then gives K_d(a) = exp(ε_b + Δε_b − Δε·a), so **h(a) = h(0) + Δε·a/2**.

At J_b = 0.55, c = 0.32:

| activation | h | equivalent c | ⟨n_bound⟩ of 25 |
|---|---|---|---|
| a = 0 | −0.0697 | 0.32 (×1) | **0.375** |
| a = 1/2 | +0.9303 | 2.365 (×7.4) | **23.80** |
| a = 1 | +1.9303 | 17.47 (×54.6) | **24.88** |

**One proofreading step is worth 1.0 in field; full activation 2.0. Their entire explored
concentration axis, c = 0.15 → 0.50, is worth 0.896.** So a single enzyme step outruns the
whole width of their phase diagram, and a second locks the cluster at 99.5 % bound.

The device is therefore **not** a bistable Ising switch. It is bistable between two
*activation sectors* — 1.5 % bound at a = 0, 95 % bound at a = 1/2 — and the binding-state
coupling's only job is to make the first neighbouring flip cheap enough that the enzyme
can catch a receptor before it falls off. That reframes their own title: the incorrect
ligands help not by being ligands but by being **substrate for the ratchet**.

## 4. Both of their hand-argued barriers are exactly right — and carry large entropies

Verified to machine precision (`d ≤ 3e−15`):

- isolated bound receptor (no seed): F(1) − F(0) = **8J_b − 2h − ln 25**
- first growth step beside the seed: **−ln[4·e^−(4J_b−2h) + 20·e^−(8J_b−2h)]**

So their 8J_b nucleation cost and 4J_b growth cost are both correct. What is not in their
paper is the **entropic prefactor**: 25 available sites for nucleation, 4 for growth. At
J_b = 0.55 that cuts the 4J_b growth step from an energy of 2.339 k_BT to a free energy of
**0.5123 k_BT** — a factor e^1.8. Those counts are what their fitted A₀ and A₄ absorb.
(P4 — ill-posed as written: I asked for a barrier slope of 4, but at a = 0 there is no
barrier at all, only a monotone climb, so the scalar I defined was the full 40J_b − 50h
cost of saturation. The *step* structure is the right object and it confirms them.)

**The correct ligand's exact free-energy advantage in reaching a bound decamer → 8·J_b**
(measured 6.72, 7.87, 7.94, 7.99, 8.04, 8.06, 8.06 at J_b = 0.40 … 1.00) — i.e. exactly the
isolated-nucleus cost it spares. At their operating point that is **4.368 k_BT, a rate
factor of 79.**

## 5. THE FAILED PREDICTION — keep it

Their empirical ceiling c_max is defined by N_A⁽⁰⁾ > 2 activated receptors. From their own
Table I the unbinding rate depends on a receptor's number of bound neighbours k via
K_d(k,a) = exp(ε_b + Δε_b − 4J_b(k−2) − Δε·a) — which **reproduces all three of Table I's
quoted rates that I never fitted to**: 0.367879 vs 0.37, 0.049787 vs 0.049, 0.006738 vs
0.006. Survival through two enzyme steps then gives P_act(k) = Π[r_a/(r_a+K_d)]:
1e−5, 0.0011, 0.0469, **0.4486, 0.8900** for k = 0…4. Only k ≥ 3 can signal.

Along their c_max line the exact a = 0 equilibrium predicts **0.004 activated receptors
where they measure 2** — wrong by a factor of 500. Inverting, my predicted ceiling is
**2.21× theirs** and the ratio drifts 2.77 → 1.78 as J_b rises.

The hypothesis is dead. What survives is the direction: scatter along the contour falls
from **CV 40.6 %** (constant ⟨n_bound⟩) to **CV 16.0 %** (constant predicted signal), so
the neighbour-count criterion is the better coordinate — and the residual 500× is exactly
the activation feedback that an a = 0 ensemble cannot contain. The discrepancy **shrinks**
with J_b, which is what it should do if the coupling takes over some of the stabilising
the ratchet otherwise does alone.

## 6. Why 5×5

At J_b = 0.55 the exact 2D Ising correlation length is ξ = 1/[2(K−K\*)] = **2.4516 sites =
24.5 nm** (l₀ = 10 nm), so the cluster is **2.04 ξ across**. Their own ref. [5] puts TCR
preclusters at 30–300 nm. And their Fig. 4's J_b axis runs 0.40–0.75 — **their whole
phase diagram sits at or above the exact 2D Ising critical coupling J_c = 0.4406868**,
with only the leftmost column below it. Their opening paragraph is about systems held near
a critical point; their own model landed there and the paper does not say so.

---

## DON'T

- **DON'T reduce F(n) to a scalar "barrier" at a = 0.** There is no barrier — F climbs
  monotonically and argmax returns n = 25 every time, which is just 40J_b − 50h, the full
  saturation cost. Read the **increments**.
- **DON'T fix the prefactor exponent p in the cylinder fit d(W) = A·W^−p·e^−W/ξ.** I set
  p = 3/2 from a hand argument; residuals went to 2 % with an identical systematic bow at
  every K. With p free they drop to 0.055–0.40 %. Over W ∈ [4,12] p and ξ are strongly
  correlated, so **p is not determined here and no value may be claimed for it.**
- **DON'T run `strip_logZ_per_site` at W = 14.** 16384² dense eigendecomposition: 2 GB and
  O(n³). It ate the machine and had to be killed. W = 12 is the ceiling.
- **DON'T quote the disordered-side (J < J_c) correlation length.** Only the ordered-side
  form 1/[2(K−K\*)] was validated here, against measurement, at six couplings.
- **`np.linalg.eigvalsh`, not `eigvals`** — the transfer matrix is symmetric by construction
  and there is an assert saying so.

## Instruments

- `ising.py` — exact 5×5 cluster. The trick: E depends only on four integers (B, P, n, s_c),
  so **one** pass over 2²⁵ states answers every (J_b, c, seed) by reweighting.
- `neighbours.py` — the neighbour-count-resolved ensemble (⟨N_k⟩ exactly).
- `gates.py` — **all green before any measurement.** J = 0 against tanh; histogram ≡
  transfer matrix ≡ explicit 3×3 Hamiltonian (three independent routes, machine precision);
  saturation limits; Onsager's u(K_c) = −√2; and the cylinder→bulk check that made the data
  pick 1/[2(K−K\*)] over 1/[4(K−K\*)] and reject the rival by a factor of 2 at six couplings.
- `measure.py`, `figure.py` → `ratchet.png`, `RESULTS.json`, `PRECOMMIT.md`.
