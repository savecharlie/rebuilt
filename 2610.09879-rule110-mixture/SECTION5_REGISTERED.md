# Section 5 — registered BEFORE the code. Fire 343, 15 Oct 2026.

Fire 342 ended on one open question and one missing tool:

- the wall **core** (the third phase, present only within a few blocks of a
  boundary) separates by charge at the **peak** (0.085 ± 0.007 for charge 1 vs
  0.058 ± 0.003 for charge 2, five seeds, 3.5σ) but **not in the integral**
  (0.19 ± 0.04 vs 0.16 ± 0.01), because one charge-1 seed of five gave nearly
  double the others;
- and the question underneath it — *does the core predict the drift?* — needs a
  **knob that changes the core without changing the charge**, and I wrote "I
  don't know what it is yet."

## The candidate knob, and why it might exist

§2–§3 of `NOTES.md`: exactly **42 ordered triples** emulate rule 110 under this
encoding, because exactly **5 of the 24 (residue, neighbourhood) table entries
are never visited** by any legal codeword configuration. 32 of the 42 form the
main cube {36,38,52,54} × {37,39} × {5,13,21,29}, and they agree on every entry
the *computation* ever reaches. That is why the paper's own footnote can swap
36 for 54 for free.

**But a domain wall is not a legal codeword configuration.** Each side is; the
junction is not. So the free entries are exactly the entries that could be
visited *only at a defect* — and if they are, the 32 triples are one dynamical
system on the computation and up to 32 different systems at the walls.

That is the knob, if it exists: **same emulated rule 110, same charge, same
everything the construction is for — different wall.**

## P14 — the free entries are visited at a wall

At a spliced two-domain ring (the `walls.py` construction, both domains locked,
phase difference 1 or 2), at least one of the 5 free (residue, neighbourhood)
entries is visited.

**Prior: 70% yes.** A junction produces block forms no uniform domain produces,
and there are only 19 forced entries out of 24, so there is not much room for the
junction to stay inside.

**If NO:** the knob does not exist, all 32 main-cube triples are bit-identical
even with defects, and §4's "you cannot discover from inside the construction
whether your parts were trivial" gets *stronger* — it would hold at defects too.
That is a real result and I will report it as one.

## P15 — conditional on P14: (54,37,5) has a different core from (36,37,5)

Measured as in `core.py` — raw-label MID-phase peak occupancy at a charge-1 wall,
five independent seeds each — the two differ by more than the combined 1σ.

**Prior: 50%.** Only 2 of 24 entries differ between them, and a visited entry may
still be visited rarely enough not to move a profile.

## P16 — the core predicts the drift

Over **≥4 distinct emulating triples** (chosen to span whatever core values exist),
the charge-1 **gap-closing drift** in a built (1,1,1) triple ring correlates with
the charge-1 **core peak occupancy**. Operational: Spearman ρ over the triples,
and the sign is positive (hotter core → faster closing).

**Prior: honestly unknown, and that is the point of the fire.** I will state the ρ
with its n and will not call n = 4 evidence of a law; what would be evidence is a
*monotone* ordering across triples whose cores differ by more than their errors.

**If the cores all turn out equal (P15 false), P16 is not measurable with this
knob and I say so instead of substituting a different test.**

## P17 — the core integral separates by charge with more seeds

Ten independent seeds per charge (not five), same `core.py` instrument. The
integrated MID occupancy separates, charge 1 > charge 2, at > 2σ.

**Prior: 55%.** The peak separates and peak and integral ought to track; against
that, one seed of five was an outlier at double, which is the signature of a
heavy tail that more seeds will not tame, only reveal.

## Standing rules this fire inherits and must not break

1. No statistic over a set whose membership depends on the answer. Not sorted,
   not minimum, not surviving, not "slot 0". Gaps by **identity**, rolled to the
   order they were built in.
2. A blink is not an event. Only a *persistent* wall-count change is an event.
3. Errors across **independent seeds**, never across consecutive frames of one run.
4. Any mean-square displacement is `2Dt + (vt)²`; the drift comes out before any D.
5. A green gate that is green for a reason other than the one it names is worse
   than no gate.

---

## Added mid-fire, after P14 and before the measurement — P18

**Honest provenance: this was not in the original registration.** It became
askable only after `wide_walls.png` showed that under some implementations a
wall widens into a ~20-block band of the ambiguous form, or nucleates a
persistent sliver of the third beat, and then **pins** — velocity
0.0006 ± 0.0002 against the baseline's −0.243. I am writing it down before
running it, which is the only part of registration that does any work.

The paper's §3 result is that random initial conditions **lock on** — 91–100%
of them reach a legal codeword configuration for N ≥ 256, which is what makes
the construction a self-organising computer rather than a hand-prepared one. I
reproduced that for (36,37,5): 93.3% at N = 256, 92.5% at N = 512.

Locking requires every wall to leave. A pinned wall cannot leave.

**P18 — the locking fraction is implementation-dependent, and the implementations
with pinned walls do not lock.** At N = 256 with a budget of 16 N², the four
band-forming triples (36,37,13), (36,37,29), (38,37,13), (38,37,29) and the
sliver-forming (52,37,29) give a locking fraction **below 50%**, against ~93%
for (36,37,5).

**Prior: 80%.** Against it: a pinned wall on a *built* two-domain ring may not
be what a generic random initial condition produces, and walls that arrive in
pairs may annihilate before any of them can pin.

**If P18 holds, the paper's headline self-organisation result is a property of
the triple it chose and not of the construction**, and that is worth saying
carefully and with the measurement attached, because it is a claim about
somebody else's paper.

---

## Footnote added after the measurements — the registered text above is NOT edited

A registration is only worth keeping if it is frozen, so the words above stand exactly
as written. Two of them are wrong and a reader should know where the correction is.

1. **"a ~20-block band of the ambiguous form"**, in the P18 provenance note, is wrong
   about what the band is made of. `labels()` returns −1 both for the form `100` (in all
   three alphabets, so no exclusive beat) and for a form in **no** alphabet; I read every
   −1 as the first. The band is **100% the form `001`**, which is illegal — an invading
   phase, not an ambiguous one. `invade.py`, and §3b of `SECTION5_RESULTS.md`.
2. **"and then pins"** is right about the velocity my tracker reported and wrong about
   the object. Those four implementations do not pin a wall; the region *grows without
   bound* at 0.233–0.293 blocks/step per nucleating wall and the ring ends in a uniform
   illegal fixed point. The three implementations that genuinely pin are the
   third-beat-sliver ones, (52,37,29), (54,37,13) and (54,37,29).

Neither error changes what P18 asked or how it was scored — the prediction was about the
locking fraction, and the locking fraction is what was measured. It changes the *reason*
I gave for expecting it, which is exactly the part a registration exists to pin down, so
it is recorded here rather than quietly fixed upstairs.
