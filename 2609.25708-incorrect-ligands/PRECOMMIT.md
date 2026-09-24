# Pre-commitment — written BEFORE any code was run

Paper: arXiv:**2609.25708v1** (Chen, Leung & Chen, dated 23 Sep 2026), PDF retrieved
2026-09-24. Version pinned; v1 is the only version extant.

## The convention I am fixing (law 11 — a prediction without its convention is not scoreable)

The paper prints its own exact mapping. With all activation levels a_i = 0, set
s_i = 2 b_i^2 - 1 (free -> -1, bound -> +1). Their eq. (1) becomes

    H = -J_b * sum_<ij> s_i s_j  -  h * sum_i s_i ,     h = (ln c - eps_b - d_eps_b)/2

with eps_b = -3.5, d_eps_b = +2.5 (Table I), so **h = (ln c + 1)/2** and h = 0 at
c* = e^-1 = 0.367879.

Cluster: **5x5, free spins.** End Matter: each receptor on the cluster boundary has its
absent outside neighbours treated as b = 0, i.e. **s = -1 pinned**. So a corner site
carries 2 pinned -1 neighbours and an edge site 1. Energies in k_B T, beta = 1.

Correct ligand: b = -1 at the seeded site, so s = +1 there as well; its extra stability
is d_eps_b = 2.5 k_B T relative to an incorrect complex. I model it as the centre site
**pinned to s = +1** (the strong-binding limit, K_d = e^-6 vs e^-1), and separately as a
finite extra field, and I will report both.

Observation window for their specificity criterion: N_A^(0) > 2 at t ~ 220 s.

## Predictions (scoreable, numbered, no hedging)

**P1.** At their headline operating point J_b = 0.55, c = 0.32 (h = -0.0697), the exact
equilibrium of the seed-free 5x5 cluster is strongly FREE: mean bound fraction < 0.10.

**P2.** With the centre pinned bound, the exact equilibrium bound fraction at
(J_b, c) = (0.55, 0.32) is still LOW — below 0.30. i.e. their tenfold amplification is a
**kinetic transient**, not an equilibrium switch, and the activated state is metastable.

**P3.** The equilibrium field at which the seed-free cluster flips (bound fraction = 1/2)
sits at h > 0, i.e. at c > e^-1 = 0.368, because the pinned -1 boundary costs surface
free energy. It should RISE with J_b (a stiffer interface needs a bigger field).

**P4.** The exact projected free-energy barrier dF(n) for the seeded cluster grows
LINEARLY in J_b with slope **4** over their range 0.3 <= J_b <= 0.7 — which is what
would justify their hand-argued t_c ~ (A0 + A4 e^{4 J_b}) t0 / c. I predict the slope is
4 to within +-0.5 and that the barrier is NOT 8 J_b (the isolated-spin value).

**P5 (the one I am least sure of).** The seed-free barrier is larger than the seeded one
by an amount of order 4 J_b, not by an order of magnitude — so specificity comes from
kinetic proofreading plus the observation window, not from the Ising barrier alone.

## What would make each fail
P1 fails if equilibrium is already bound at c = 0.32 (then their S0 is metastable from
the start and the specificity clock is running the whole time). P2 fails if one seed is
enough to move the exact equilibrium past 30% bound (then the amplification is a genuine
thermodynamic first-order switch triggered by a single molecule, which is a stronger and
more beautiful claim than the one they made). P3 fails if the flip field is negative.
P4 fails if the fitted slope is 8, or if ln(barrier) is not linear in J_b at all.

## Scope — said up front (law 10)
I am NOT rebuilding their stochastic simulation, their diffusion, their kinetic
proofreading cascade, or their measured t_a. I am rebuilding the **exact equilibrium
statistical mechanics of the cluster they themselves mapped their model onto**, which
they solved only by words, and asking what it says about the device.
