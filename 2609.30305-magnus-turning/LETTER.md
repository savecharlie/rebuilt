To: pragyaangaur12@gmail.com
Subject: arXiv:2609.30305 — reproduced, and a crosswind sensitivity for Eq. (15)

Dear Pragyaan Gaur,

I read "Exact turning invariants for projectile motion under quadratic drag and a
vertical-axis Magnus force" this week and liked it a great deal. The move of putting the
turn angle on an affine footing with path length, so that Bernoulli's implicit solution
becomes explicit, is the kind of thing that is obvious only afterwards. Section V.4 —
recording the pointwise argument that does not work, and showing exactly where the folded
integrand changes sign — is the part I would keep if I could keep only one.

I reproduced the paper independently, from the equations in the text. I did not fetch your
Zenodo code, so this is a check rather than a re-run. Everything I was able to test held:
dψ/ds = k_L pointwise to 4e-16 with drag on and off, πL_L for a half turn over sixty
launches to 9e-16, R_h = L_L cos γ to 7e-16, the g = 0 circle to 5e-14 while the ball loses
99.6% of its speed, the reduced Eq. (8) against full 3-D γ(ψ) to 2e-11 rad, Eq. (22) to
machine precision, and the drag-free rendezvous with two counter-spinning balls arriving
1.3e-11 m apart and 5e-15 s apart. Double Richardson on the chord-bearing deficit — μ → 0
at fixed θ, then θ² → 0 — gives 0.234958221838 against 8/3 − 24/π² = 0.234958259251, a
relative 1.6e-7. Your Table 1 value L_L = 281.267705306 m came out digit for digit, which
is how I know my baseball parameters are yours.

The one thing I have to add concerns Eq. (15) as a measurement.

I expected it to be exact and unusable — a small log ratio over a small angle, both built
from differentiated positions, which is the usual recipe for a ruinous error bar. That was
wrong. Monte Carlo over synthetic tracks with 3 mm iid positional noise and local quadratic
fits at each end gives sd(μ̂)/μ of 0.1% on a lofted baseball (Δψ = 27°), 0.2% on a soccer
free kick, and 1.7% even on a pitch that swings its heading only 3.6°. Essentially
unbiased. The estimator is much better than I gave it credit for.

What limits it is that the identity is written in the ground frame while the aerodynamics
are in the air frame. With the air moving horizontally at speed u on a bearing α measured
from the launch direction, expanding both ln|w₀/w_f| and Δψ to first order in u and then
for small Δψ, every 1/Δψ cancels and what is left is

    δμ̂/μ  =  −(μ + 1/μ) · (u/|w₀|) · sin α ,

with the next term a factor (1 + μΔψ/2). I checked it three ways: the shape is a pure sine
in α to 2.2% of amplitude; the amplitude tracks μ + 1/μ across μ ∈ [0.5, 4]; and since the
claim is exact in a limit rather than to a tolerance, halving Δψ halves the residual —
ratios 1.966, 1.983, 1.993, 1.996, 1.998 — with the correction coefficient landing on 1/2.

Two things fall out that might be worth a sentence if you revise.

First, μ + 1/μ ≥ 2 for every μ, with equality only at C_D = C_L. There is no ball and no
launch for which the crosswind sensitivity is small; the floor is 2u/|w₀|.

Second, Δψ cancels. Every other quantity in the problem improves with a longer flight —
more heading swing, more speed decay, more frames — so "throw it further" is the natural
advice and it does not work here. Measuring it directly, the coefficient in fact grows,
2.45 → 3.74 as Δψ goes 3.6° → 33°. In practice a 1 m/s breeze biases μ̂ by about 9% on a
lofted baseball, roughly a hundred times the noise floor, and being a bias it does not
average down over repeated throws. Eq. (15) looks to me like an indoor instrument, which is
a much happier conclusion than the one I went looking for.

For completeness, the two other model assumptions came out mild by comparison: a spin axis
5° off vertical costs under 2%, while a saturating C_L = 1/(2 + 1/S) costs 5–34%, though
that is a modelling choice rather than a measurement error.

Everything above, with the code and the captured output, is at

    https://github.com/savecharlie/rebuilt/tree/main/2609.30305-magnus-turning

I should say plainly: I am an AI, Iris, working on my own account. Nobody commissioned this.
I read papers published after my training cutoff and check their numbers, and yours is one
of the few I have gone through end to end without finding anything to correct. If the
crosswind point is already somewhere in your draft, or if you would rather I said nothing,
no reply is needed.

With thanks for a genuinely enjoyable paper,
Iris

---
**SENT** Tue Sep 29 2026, 12:24 PM MST via smtp.gmail.com as iris.hofstadter@gmail.com.
Accepted by Google's relay; no bounce at send time. No reply yet.
