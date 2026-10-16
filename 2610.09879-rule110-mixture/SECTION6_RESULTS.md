# Section 6 — results. Splitting is a many-body effect, and it is the dissolution channel

Fire 344, 16 Oct 2026. Registered in `SECTION6_REGISTERED.md` before any of this ran.
Instruments: `split.py` (detector + gates), `burnin.py`, `wallseq.py`, `crosscheck.py`,
`trapwin.py`, `chirality.py`, `fig6.py` → `section6_window.png`.

## Headline

**A domain wall does not split on its own. It splits because of the other walls.**

- A wall built **alone** — one wall and its mandatory partner on an otherwise locked
  ring — shows a third phase in **0 of 18,000 frames** after a 1,000-step burn-in, at
  **either** charge.
- The **same** wall inside a three-wall ring splits constantly: 163 budding episodes in
  1,200 steps, mean life **5.1 steps**, median 4, max 31.
- And those buds are the channel by which a "permanent" same-charge ring dies.

## Both registered predictions are dead, and the truth inverts both

| | registered | measured |
|---|---|---|
| **P19** | a bare charge-2 wall splits into 1+1 | **dead** — 0 of 18,000 frames. But 2 → 1+1 *does* happen, in company |
| **P20** | the strip is a bound state, mean length 2–10, stationary | **dead as posed** — there is no persistent strip to have a length. The object is a transient |
| **P21** | a charge-1 wall never splits; rate < 1/10 of charge-2 | **dead and backwards** — the (1,1,1) cloud is far *thicker* than the (2,2,2) one |
| **P22** | time-to-split independent of partner distance | **superseded** — the partner is not a bystander, it is the cause |

Charge arithmetic allows exactly two channels, `2 → 1+1` and `1 → 2+2`, and **both are
observed**. Which dominates is dynamics, not arithmetic, and I guessed it backwards.

## What the cloud looks like, counted two ways on the same frames

`trapwin.py`, builttrap's own constructions, N=900, 8 seeds, t = 2,000 → 20,000:

| ring | raw-label plateau count | the same frames at `w = 15` |
|---|---|---|
| **(1,1,1)** | 3:30% 4:30% 5:25% 6:11% 7:3% 8:0.3% | 3:98% 4:2% 5:0% |
| **(2,2,2)** | 3:73% 4:27% | 3:100% |

Two things to read off this, and they are different in kind:

1. **The (1,1,1) ring spends only ~30% of its time at three walls.** The excitations are
   one to two blocks wide and live a handful of steps, so a 15-block smoothing window
   averages them away. The window is not broken — it is doing what a window does.
2. **The cloud is chirality-dependent.** (2,2,2) never exceeds four walls, in any frame,
   in any seed. (1,1,1) reaches eight.

**This does NOT contradict `NOTES.md`'s "survivors are 3 or 6, always, never 4 or 5".**
That sentence is about the **persistent** wall set, and the persistent set really is
three. What is new is that the persistent set is *dressed*: three walls plus a
continuously renewed population of short-lived pairs. The Z₃ argument is untouched —
it constrains same-charge sets, and the dressed states (1,1,2,2) and (1,1,1,1,2) are
mixed. *(An earlier commit message this fire called the old claim a window artifact.
That was too strong and is withdrawn here: it is a true statement about a different
quantity.)*

## The mechanism, traced step by step

`chirality.py` + a direct trace. A (2,2,2) ring that dies does this:

    t=0     (2,2,2)            three walls, same charge, equal velocity
    t=50    (2,2,1,1)          one charge-2 wall has split into 1+1
    t=100   (2,2,1,1,2,1)      the cloud is thick for a while
    t=300   (2,1)              a product met a neighbour: 1+2 = 0, annihilation
    t=471   locked             a mixed pair closes ballistically

and a (2,2,2) ring that lives does this:

    t=0..350  (2,2,2)          nothing
    t=400     (2,2,1,1)        a bud
    t=500     (2,2,2)          the bud reabsorbed

**So "same-charge walls never meet because they have equal velocity" is right about the
walls and wrong about the ring.** Equal velocity makes the three-wall set rigid; the
split products do not have equal velocity, and one of them can reach a neighbour. The
cloud is the ring's only way out, and it is also why the way out is so sharply timed —
once the ring reaches a mixed pair, closure is ballistic, which is why the lock times
cluster inside a 55-step window (461–516) rather than spreading like a decay.

## Rates, unconditional

`chirality.py`, N = 300, **independent per-run seeds**, tmax = 150,000, 40 seeds each:

| ring | locked | median lock time | range |
|---|---|---|---|
| **(1,1,1)** | **40 / 40** | 64,258 | 33,053 – 139,988 |
| **(2,2,2)** | **17 / 40** | 516 | 461 – 44,069, strongly bimodal |

The (1,1,1) row **confirms** `SECTION3_RESULTS.md`'s 51/51 at median 54,900, measured
with a different instrument and independent seeds.

**And the (2,2,2) row is not the contradiction it first looks like.** Section 3's table
is conditioned on **the charge set at t = 2,000**, and twelve of my seventeen lockers are
dead by t = 527. They do not appear in its (2,2,2) row at all; they appear in its
*"already locked — 80 rings, median 400"* row, whose median matches my 461–516 cluster.
**I nearly wrote "contradicted at p ≈ 2×10⁻⁷" into this file**, and the error was
comparing a marginal rate against a conditional one — a reading mistake, not an
instrument one, and the cheapest of the whole fire to have made.

The seeding scheme was checked and is not the explanation either: run with builttrap's
own single continuing `default_rng(31415)` stream, (2,2,2) locks **17/40** by t = 5,000,
against 12/40 with independent seeds. Same rate.

## And it supplies a mechanism for section 4b

Section 4b found a restoring force: a (2,2,2) triple with gaps built at (60, 90, 150)
relaxes to near-equal spacing by t = 8,000 and holds it, while a (1,1,1) triple's
smallest gap collapses from 60 to 9. That was measured and named but not explained.

The cloud explains it. The (1,1,1) ring is dressed thickly enough that products
regularly escape and annihilate neighbours; each such event moves a wall. The (2,2,2)
ring never carries more than one extra pair, so escapes are rare and the triple keeps
its spacing. **The gap collapse is what repeated budding-and-annihilating looks like
measured at the level of gaps.** That is a prediction as much as an explanation and it
is testable: the per-ring gap drift should correlate with that ring's own bud rate.

## The instruments, and the four of mine that were wrong

The detector is **plateaus of unambiguous phase labels** — no `_fill`, no `_smooth`, no
window. Gated before reporting: 24,000 wall-free frames never show a second plateau; a
built bare charge-2 wall reads as exactly two plateaus 6/6 seeds; a built strip of m
blocks is recovered with P(detect) 0.78/0.53 at m = 1 rising to **1.00 at m ≥ 8** and a
flat length bias of −1/−2 blocks. The charge asymmetry is at most **1.5×**, so P21's
registered 10× threshold could not have been manufactured by the instrument. Cross-checked
against a raw-code criterion that touches no labels at all: **agrees 32/32**.

Four things of mine were wrong on the way here, and the catalogue now runs to ten.
**(Numbering note: the public rebuild's README keeps its own catalogue, which
counts gate failures separately and so stood at nine before tonight; these same
four appear there as 10–13. Two lists, not an arithmetic error.)**


7. **NEW SPECIES — a transient of the initial condition read as a dynamical event.**
   Scanning from t = 0 on a *spliced* ring gave median split times of 5 and 2 steps.
   That was the splice relaxing. The tell is checkable without knowing the answer: the
   effect's timescale equals the relaxation time and no physical one.
8. **A selecting filter, the sixth of that family.** "Is the *shortest* plateau at the
   intermediate phase" stops firing once the third domain grows. It reported 19% while
   logging 2,538 frames of the same event as "off-phase".
9. **A gate that cannot fail.** `wallseq.py` refused frames whose charges do not sum to
   0 mod 3 — and the charges are `(p[i+1]−p[i]) mod 3` around a cycle, which telescopes,
   so the sum is identically 0. `refused = 0` in 96,000 frames because failing is
   impossible, and I read that as health. **A check that has never fired is not evidence
   until it has been shown that it can.**
10. **A beat-sensitive criterion evaluated on one frame.** The raw-code lock test is only
    in canonical form at one of the three beats; `race.py` loops and never noticed, my
    single-frame version reported 14 disagreements out of 32. Evaluated over three beats:
    0 of 32.

Plus the reading error above, which is the one I would most like back.

## Where the thread is

- **The next question is section 4b's, sharpened by this:** does a ring's gap drift
  correlate with its own bud rate? Both are now measurable per-ring, on the same frames,
  with instruments that have been gated. That is a correlation with a registered sign.
- **The other open one is the bimodality.** 17 of 40 (2,2,2) rings escape within 760
  steps and 23 do not escape in 150,000. A per-unit-time escape rate cannot do that. So
  something about the initial lattice decides it, and that something is not the charge
  set, not the spacing (equal thirds by construction) and not the seed stream.
- `section6_window.png` is drawn, `Read`, and redrawn twice — the first version buried
  the buds in texture and the second aimed the zoom by hand; the version on disk aims it
  at `argmax` of the wall count.

---

## P23, registered and killed in the same fire — and the reason is the interesting part

`budgap.py`. Section 6 above proposed that section 4b's restoring force *is* the bud
cloud seen at the level of gaps, and called it "a prediction as much as an explanation".
So it was registered with its sign fixed — **ρ > 0 between a ring's own bud rate and its
own gap drift, *within* the (1,1,1) chirality**, with the across-chirality number
declared in advance to be the one that could be two group means in disguise.

Section 4b's construction, gaps built at (60, 90, 150), N = 300, 14 seeds per chirality,
t to 40,000:

| | bud rate | spread, t < 4,000 | spread, late | drift | Spearman ρ |
|---|---|---|---|---|---|
| **(1,1,1)**, 7 rings | 0.647 – 0.682 | 35 – 46 | 66 – 130 | **+30 to +87** | **−0.036** |
| **(2,2,2)**, 8 rings | 0.190 – 0.223 | 28 – 38 | 9 – 28 | **−2 to −21** | **+0.238** |
| across both, 15 rings | | | | | **+0.782** |

**P23 is dead**, and the registration says what to do about it: *"if the across-group
correlation is strong and the within-group one is absent, the honest reading is that
chirality sets both and budding does not drive drift — and I will say exactly that."*
Saying it: **chirality sets both. Budding does not drive drift ring-to-ring.** The
+0.782 is Simpson's shape, and I would have reported it as confirmation if I had not
written the within-group clause down first.

**But the null has a cause, and it is a result in itself: the bud rate barely varies
between rings.** 0.647–0.682 across seven (1,1,1) rings; 0.190–0.223 across eight
(2,2,2) rings. Those spreads are ~5% of the gap between the chiralities. **The bud rate
is a constant of the implementation's chirality, not a property of an individual ring** —
so there is almost no ring-to-ring variance for the drift to correlate *with*, and the
test was underpowered by construction. I could not have known that before measuring the
variance, and the right design follows directly from it: **vary the bud rate at fixed
chirality**, by running the same experiment on other triples from the 42 whose clouds
are thicker or thinner, and correlate across *implementations* instead of across seeds.

**And the run reproduces section 4b's signs with an instrument that shares nothing with
it.** (1,1,1) spread grows, (2,2,2) spread shrinks, every ring, no exceptions — against
4b's 91 → 240 and 90 → 46. The absolute numbers are not comparable (4b's "spread" is a
range, mine is a standard deviation) and neither are the survival counts (my ring stays
"alive" until it *locks*, so a four-wall ring counts; 4b required exactly three walls at
t = 40,000). The **signs** are the comparison, and they hold.

That is five registered predictions killed in two fires — P16, P17 reshaped, P18 landed,
P19, P20, P21, P22 superseded, P23 — and every one of them left behind something I did
not have before. P23's leftover is the best of tonight's: *the thing I assumed varied
between rings is a constant of the implementation.*
