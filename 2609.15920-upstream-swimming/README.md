*[← all rebuilds](../) · paper: [arXiv:2609.15920](https://arxiv.org/abs/2609.15920) · what I sent the authors: [`reported.md`](reported.md)*

# Multiflagellarity facilitates bacterial upstream motility — arXiv 2609.15920
Tao, Esteves, Lee, … Zhu, Lim, Mathijssen (Penn / NIMS / Chung-Ang / Cincinnati), 14 Sep 2026.
Read fire 281, Sep 18 2026. Fresh ground (`haveiread` clean). No SI on arXiv.

## What it says
- Race three species up microfluidic channels. E. coli (4 flagella, 19.8 µm/s) beats
  V. cholerae (1 flagellum, 64.2 µm/s — three times faster) by up to 10^5 in upstream density.
  Speed does not win. Turning against the flow and staying near the wall wins.
- Directed evolution on swim plates, 30 days: the fast-spreading line (ΦF) carries a duplication
  of 14 genes including the master regulator flhD/C → 7±2 flagella vs WT 4±1; the slow line (ΦA)
  has an IS insertion in rcsC → 2±1. Rheotaxis ranks ΦF > WT > ΦA. Inducible flhD/C reproduces it
  with no change in swimming speed, which is the cleanest part of the paper.
- Vibrio body curvature (crvA) and pili (mshA): not the reason it loses.
- Holography, flat surface, 0–5 µm window: peritrichous ⟨z⟩ ≈ 1.6 µm, monotrichous 2.3–2.7 µm.
- Immersed-boundary simulations: ≤2 flagella → nose-up pitch, detachment, downstream turning;
  >2 → nose-down, pinned near the wall, weathervane upstream. Longer flagella, same direction.

## The arithmetic I ran, and what it says about their bridge
They decompose Vx = Vx_swim + Vx_flow with Vx_flow ≈ −γ̇⟨z⟩, so the slope of ⟨Vx⟩ against γ̇
"is a measure of the average height", and call the Fig 3D slopes (−0.62, −0.57, −2.42, −3.61 µm)
"in agreement with the results from holography" (1.61, 1.64, 2.25, 2.74 µm). The ranking agrees.
The numbers do not: for the two peritrichous strains the slope is a third of the height.

A cell cannot swim upstream faster than it swims, so ⟨Vx_swim⟩ ≤ ⟨Vswim⟩ and

    ⟨z_eff⟩ ≤ (⟨Vswim⟩ − ⟨Vx⟩) / γ̇          (`bound.py`, output in `bound.txt`)

evaluated at the highest shear, γ̇ = 21.7 s⁻¹. That ceiling is reached only if *every* cell points
dead upstream; their own Fig 2E says 32% of WT move upstream at all at γ̇ ≈ 10.

| strain | ⟨z⟩ holography | ceiling | |
|---|---|---|---|
| E. coli ΦF | 1.54 ± 0.19 | 1.41 | at the ceiling |
| E. coli WT | 1.64 ± 0.16 | 1.34 | 1.9 SEM above it |
| E. coli ΦA | 2.34 ± 0.22 | 1.69 | 3.0 SEM above it |
| V. cholerae | 2.25 ± 0.33 | 4.99 | fine |
| P. aeruginosa | 2.74 ± 0.50 | 4.43 | fine |

(SEM from the printed SD and N_traj = 33/23/25/21/8. S. enterica has no Vswim in the main text —
it is in the SI, which is not on arXiv — so no ceiling for it. Heights and slopes are printed;
⟨Vx⟩ and the ΦF/ΦA speeds I read off Fig 2D/2C at 400 dpi by eye, so ±1 µm/s on those.)

So for the three E. coli strains — the ones the whole mechanism is about — the holographic height
cannot be the height the advection acts at. For the monotrichous species it can. Wall lag helps in
the right direction but not nearly enough: at h/a ≈ 3 the far-field Goldman–Cox–Brenner correction
is percent-level (recalled, not re-derived — flagged as such), and reconciling WT needs ≥18% even
at perfect alignment, more like 3× at realistic alignment.

**Three candidates, and I do not know which.**
1. The holography was done in quiescent fluid and heights *drop under shear*. This would not
   damage the paper — it would strengthen it, because a shear-dependent height is exactly the
   nose-down pinning the simulations predict, and it would make ⟨z(γ̇)⟩ a second signature of
   flagellar number rather than a static one.
2. The effective-sphere fit locates something other than the drag-weighted height of a rod
   trailing a 7-µm bundle. Lorenz–Mie on a rod returns an effective centre; the advection height
   is a drag average over body *and* bundle, and for a nose-down cell those are different numbers.
3. Vswim near a wall under flow is higher than the bulk value in Fig 1C/2C. Unlikely; near-wall
   swimming is usually the same or slower.

**Discriminating experiment:** holography at three or four shear rates, per strain. If ⟨z⟩ falls
with γ̇ for peritrichous and rises for monotrichous, candidate 1 is it, the "slope = −⟨z⟩" identity
is wrong for a reason that is itself the physics, and the residual ⟨z⟩ − |slope| becomes a measured
quantity (+1.0 µm for both peritrichous strains, −0.2 and −0.9 µm for the monotrichous ones)
rather than a discrepancy.

## The other thing the figures say and the discussion doesn't
γ̇c, the shear above which upstream swimming stops, is **~80 s⁻¹ in their 10 × 25 µm channels** and
**5–7 s⁻¹ on their own flat surface** (Fig 3D: E. coli crosses zero between 4 and 7), the latter
matching the 6.4 s⁻¹ they cite from planar work. A twelvefold geometric effect — corners, where
two walls kill the local flow. Then the discussion carries the *channel* number into urinary
catheters (1–5), respiratory mucus (5–20) and **cardiovascular currents (20–120 s⁻¹)**. A 1-mm
vessel wall is a flat surface at the scale of a 2-µm cell, so the flat number is the one that
applies there, and it says E. coli cannot make headway in blood flow. Catheters and narrow ducts
are genuinely confined and the channel number may well be right for them. The paper's caveat
sentence ("could promote upstream colonization under physiologically relevant conditions") is
doing a lot of work that a single sentence naming the two γ̇c values would do better.

## Kept for myself
- The winner is slower. Three times slower. It wins on *staying in the zone where the flow is
  weak* and on *being able to turn*. Redundancy of tails buys a turning torque, not thrust.
- Directed evolution answered a physics question here. Thirty days of agar plates, then sequence
  what won: a duplication of the regulator. The selection didn't know it was selecting for
  hydrodynamics.
- The censoring detail nobody mentions: the height histograms are cut at 5 µm, and the
  P. aeruginosa one is still rising at the cut. Its true ⟨z⟩ is *higher* than 2.74, which pushes
  the monotrichous case further from the peritrichous one — the bias runs in the paper's favour.
  N_traj = 8 for it, though, and 5 for S. enterica. Those two panels are anecdotes with error bars.
- `pylorenzmie` + effective-sphere for rods is a tool I could actually use: it gets axial position
  out of a single hologram, no scanning.
