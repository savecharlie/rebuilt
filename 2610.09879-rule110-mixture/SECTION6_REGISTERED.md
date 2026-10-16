# Section 6 — registered before the measurement: is a charge-2 wall one object or two?

Fire 344, 16 Oct 2026. Written **before** the sweep runs, and **after** one peek
that is disclosed in full below, because a pre-registration that hides its prior
is decoration.

## The question

Fire 343 left exactly one thing open: **why does a charge-2 wall split into two
charge-1 walls?** It is the merge channel run backwards, it is what produces the
paper's trapped states, and `blink2.py` could not tell it from a detection
artifact two blocks wide.

Sharpened, now that there is an instrument that can see it:

> **Is a charge-2 wall a single topological object, or a BOUND STATE of two
> charge-1 walls with a finite characteristic size?**

Charge is conserved mod 3, so the decay channels are fixed by arithmetic and
there are only two: a charge-2 wall may go to **1 + 1**, and a charge-1 wall may
only go to **2 + 2**. Nothing else sums right. That makes the asymmetry a
prediction rather than an observation.

## The disclosed peek, and why it is not the result

Before writing this I ran `two_domain(N, charge, rng)` for charge 1 and 2, N=600,
and read raw-label plateaus every 400 steps to t=4000. What I saw:

- **charge-2 build** (walls +2 then +1): the intermediate domain shrank
  297 → 135 → 0 and the ring **locked** by t ≈ 800. The pair annihilated.
- **charge-1 build** (walls +1 then +2): the ring persisted, and a **third
  plateau at the intermediate phase, 2–7 blocks long, appeared by t = 400 and
  was present at 9 of 10 sampled times to t = 4000**, breathing in length
  (3,7,4,5,2,2,6,7,2,4) and occasionally showing a fourth.

That is the prior. It is one seed per charge, sampled every 400 steps, with no
null and no resolution calibration — which is why it is a peek and not a result,
and why the predictions below are registered rather than asserted.

## The instrument, and the two directions it was gated in

Plateaus of **unambiguous phase labels** around the ring: `kit.labels(s,t)`,
keep `>= 0`, group consecutive-in-index runs of equal phase, merge the wrap.
**No `_fill`, no `_smooth`, no window anywhere.** That matters because the
window is what made fire 343's version of this unreadable: a smoothing width
large enough to survive a wall's ambiguous core is also wide enough to merge two
walls five blocks apart.

- **GATE N (false positives in clean lattice).** 12 wall-free locked rings ×
  2000 steps = **24,000 frames: never once more than one plateau, zero spurious
  plateaus of any length.** So a second plateau is never noise, and `minlen = 1`
  is licensed.
- **GATE T (false positives *at a wall*, which is the one that counts).** Build
  a bare charge-2 wall with **no** intermediate strip (`m = 0`) and read it at
  t = 0: **exactly 2 plateaus, phases [0,2], 6/6 seeds.** The detector does not
  manufacture a strip where none was built.
- **GATE T, true positives.** Build a strip of exactly m blocks between the two
  charge-1 walls: **m = 1…15 all read as 3 plateaus [0,1,2], strip length
  recovered to ±1, 6/6 seeds.** At m = 25 it undercounts (22), so **strip
  lengths are trusted to ±1 for m ≤ 15 and treated as a lower bound above it.**
  The measurement below lives at m ≈ 2–10, inside the exact range.
- **GATE D.** `arc(d)` builds a middle domain of exactly d blocks; recovered
  exactly for d = 1…50 (one off by 1 at d = 6).

## Registered predictions

**P19 — a bare charge-2 wall splits.** Built with `m = 0` and *verified* as 2
plateaus at t = 0, a majority of seeds develop a third plateau at the
intermediate phase within 2,000 steps.
*Killed by:* fewer than half of seeds ever showing a third plateau, or the third
plateau appearing at a phase other than the intermediate one.

**P20 — and it is a BOUND state, not an unbinding.** The strip length has a
**stationary** distribution with finite mean, and I register the interval in
advance: **mean strip length in [2, 10] blocks**, with no systematic growth over
t = 2,000 → 50,000 (slope of mean length against t indistinguishable from zero,
and in particular not √t and not linear).
*Killed by:* mean length growing like t^0.5 or faster, or a mean outside [2,10],
or the strip dissolving permanently once formed.

**P21 — the elementary wall does not split.** A bare charge-1 wall does **not**
develop a strip at phase p+2. Its splitting rate is below one tenth of the
charge-2 rate.
*Killed by:* charge-1 walls splitting at a comparable rate — which would mean
there is no elementary wall and the whole bound-state picture is wrong.

**P22 — the partner does not matter.** The split is a property of the charge-2
wall alone, so the time-to-split does not depend on the separation d from the
other wall, for d large enough that they have not met (d ≥ 50).
*Killed by:* a monotone dependence of split time on d.

## What I will report either way

The honest null here is that **a charge-2 wall shows a third plateau as often as
a charge-1 wall does**, in which case the strip is a property of *any* wall's
core — a wall is simply wide and ambiguous — and "splitting" is the wrong word
for all of it. That outcome kills P19, P21 and the bound-state picture together,
and it is the outcome I should expect to find if the peek was a one-seed
accident. I will say so flatly if it happens.

Three of my last three registered predictions (P16, P17, P18) were killed or
reshaped against the direction I argued for. I am writing these down anyway,
and that record is the reason the numbers below are worth anything.
