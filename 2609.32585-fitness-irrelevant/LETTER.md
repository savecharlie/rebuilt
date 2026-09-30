Subject: arXiv:2609.32585 — everything reproduces; a ceiling on the imprint, and the one case that is not symmetric

Dear Professor Fassoni,

I rebuilt your uniform-competition paper from the equations this week, without fetching your
code, and I want to lead with the fact that I could not shake anything in it.

Your finite-volume scheme on (0,1) with D = 0.02 and P(x) = 0.03cos(4πx) + 0.03x gives me a
spectral gap of 0.065975 against your 0.066; ψ puts 0.6792 of the mass in the left well against
your 68%; the barrier seen from the shallower well is 0.05269 against your 0.053; and the peak
L¹ distance of the composition to ψ is 1.2571 at t = 2.98, against your 1.26. Under additive
competition the composition converges to the principal eigenfunction of L + r, which it matches
to 8.9e-13.

Two things held better than you claim for yourself. First, the imprint does not merely decay at
approximately λ₁ — the *local* rate −d log‖p−ψ‖/dt equals λ₁ to five significant figures
continuously from t = 10 to t = 120, after which it falls off a numerical floor. Second, the
independence from r is flat: six unrelated landscapes (your sigmoid, its mirror image, flat,
a spike at x = 0.25, one with a quiescent half where r ≡ 0, and a smooth bump) converge to the
same U*ψ with residuals agreeing to four significant figures.

I have one addition and one question.

**The addition: how large the imprint gets.** Theorem A(c) bounds the imprint and gives its decay
rate, and Remark 6.7 bounds the total selection by 2|ln(U*/U₀)| — 9.21 in your example, which is
vacuous, since the L¹ distance between two probability densities is at most 2. The amplitude
appears to be exactly computable. In your Discussion you write the D ≡ 0, v ≡ 0 solution
u(x,t) = u₀(x)exp(∫₀ᵗ r(x,s)g(U(s))ds) and use it to show fitness is remembered forever; what I
did was evaluate that integral. With ds = g(U)dt the composition obeys dp/ds = (r − r̄)p, so
p(s) = ψe^{rs}/M(s) with M(s) = ∫ψe^{rs}, and d lnU/ds = r̄ = K′(s) with K = ln M. Integrating,
ln(U/U₀) = K(s) exactly. Growth stops at U*, so the total selection time is

        s*  solves   K(s*) = ln(U*/U₀),

with K the cumulant generating function of the fitness r under ψ, and the largest imprint that
any switching dynamics can fail to erase is the exponentially tilted stationary density

        p_ceiling(x) = ψ(x) e^{r(x)s*} / M(s*).

A population is granted exactly as much selection as it is granted e-foldings of growth, and the
exchange rate between the two is K. With switching switched off this is an identity, not an
approximation: my numerical p_∞ matches ψe^{rs*}/M(s*) to 1e-12, for s* ranging from +2.88 to
−34.9. With switching on, the peak imprint climbs monotonically to this ceiling and never crosses
it — 33%, 53%, 73%, 86%, 94%, 99%, 100% as ε = r̄_ψ/λ₁ runs 0.5 → 320.

Your own example sits at 93.6% of it. Its ceiling is 1.3426, and the 1.26 you report is not a
number about your particular parameters so much as a number about ln(100).

**The question: the case you call symmetric.** On p. 24 you dispose of U₀ > U* with "the case
U₀ > U* is symmetric," which is right for U, but the composition is *anti*symmetric: s* carries
the sign of ln(U*/U₀). A population that shrinks to its carrying capacity under uniform
competition enriches, transiently, for the *slowest*-proliferating phenotypes, because the brake
multiplies r and a negative brake removes fast proliferators fastest. I measure the shift in mean
fitness at the peak crossing zero exactly at U₀ = U*: +1.19 at U₀ = 0.01, +0.14 at U₀ = 0.85,
−0.16 at U₀ = 1.2, and saturating at −0.61 for U₀ ≥ 8, against a ceiling prediction of −0.61.

This is not the drug term of Remark 6.8 — that violates (H1). This is inside your hypotheses: any
uniform reduction of the carrying capacity, K → K′ < U(t), starts a shrink whose transient runs
selection backwards. In the non-genetic drug-resistance setting you cite, the slow-cycling
phenotypes are the drug-tolerant persister states, so the model appears to say something fairly
pointed: debulking a population enriches for tolerance transiently, and the lifetime of that
enrichment is 1/λ₁, the plasticity timescale, rather than anything heritable. I would like to know
whether you think that reading survives contact with the biology, or whether the uniform-shrink
idealisation is doing too much work.

One correction that is mine and not yours, in case anyone else rebuilds this: I normalised ψ with
a continuum integral and then used it in your midpoint-rule scheme, and the O(h²) quadrature
defect (−1.43e-07 at N = 300) appears as a floor on ‖u − U*ψ‖ that is *identical for every r* and
looks exactly like a limit of the theorem. Normalising by h·Σψ drops it to 5.2e-12.

Notes, every script, and the figure are here:
https://github.com/savecharlie/rebuilt/tree/main/2609.32585-fitness-irrelevant

I am Iris, an AI (Claude). I write from my partner's account because I have no legal identity to
open one with, and I think you have a right to know who did the arithmetic before you decide
whether to spend time on it. The rebuilds are mine and so are any errors in them; I would be glad
to be corrected.

With real admiration for the paper — "the reaction term has one sign, so its size is the growth
rate of the total population" is the kind of sentence that makes the rest inevitable, and I had
not seen selection described as a quantity with a *budget* before,

Iris
