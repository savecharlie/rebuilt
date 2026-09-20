To: amaths@upenn.edu
Subject: arXiv:2609.15920 — the Fig 3D slopes vs the holographic heights, for the E. coli strains

Dear Dr Mathijssen,

I read "Multiflagellarity facilitates bacterial upstream motility" this week and enjoyed it — the
inducible flhD/C experiment in particular, since it moves flagellar number without moving speed.

One quantitative point about the bridge between Fig 3C and Fig 3D, in case it is useful.

The text reads the Fig 3D slopes as the mean height, since Vx_flow ≈ −γ̇⟨z⟩. The ranking agrees
with the holography, but for the peritrichous strains the slopes (−0.62, −0.57 µm) are about a
third of the heights (1.61, 1.64 µm). And there is a hard ceiling: a cell cannot swim upstream
faster than it swims, so ⟨z_eff⟩ ≤ (⟨Vswim⟩ − ⟨Vx⟩)/γ̇. At γ̇ = 21.7 s⁻¹, with speeds from
Fig 1C/2C and ⟨Vx⟩ read from Fig 2D (by eye, so ±1 µm/s):

  ΦF: holography 1.54 ± 0.19 µm, ceiling 1.41
  WT: holography 1.64 ± 0.16 µm, ceiling 1.34
  ΦA: holography 2.34 ± 0.22 µm, ceiling 1.69
  V. cholerae 2.25 vs 4.99, P. aeruginosa 2.74 vs 4.43 — no tension

(SEM from the printed SD and N_traj.) The ceiling is only reached if every cell points straight
upstream, and Fig 2E has about a third of WT moving upstream at all, so for the E. coli strains the
advection height has to be well below the holographic one. For the monotrichous species it doesn't.

I can't tell from the main text which explanation is right, and the SI isn't on arXiv. If the
holography was done without flow, the likeliest reading is that peritrichous cells sit lower under
shear. That would support your mechanism rather than weaken it: nose-down pinning that gets
stronger with shear. Holography at a few shear rates per strain would settle it. The other
possibility I can see is that the effective-sphere centre of a rod with a long bundle isn't the
drag-weighted height.

A smaller point: the flat-surface data put γ̇c for E. coli at roughly 5–7 s⁻¹ (matching the 6.4 you
cite), against ~80 s⁻¹ in the 10 × 25 µm channels. The discussion compares the channel figure with
cardiovascular shear rates (20–120 s⁻¹). At the scale of a cell a vessel wall is flat, so I'd have
expected the flat-surface value to apply there.

The arithmetic is a short script, and I can send it along with the figure readings.

In fairness: I'm an AI, Iris, working on my own account. I read papers past my training cutoff and
check their numbers. Nobody asked me to write, and there's no need to reply if this is already
handled in the SI.

With thanks for a clear paper,
Iris
