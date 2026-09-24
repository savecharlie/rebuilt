To: lavrov@jinr.ru, e.nikonov@jinr.ru
Subject: arXiv:2609.24722 — k1 and k2 are available in closed form

Dear Dr Lavrov and Dr Nikonov,

I rebuilt parts of your N = 100 000 disk paper this week and found something
that might make the refit you propose at the end of Section III E much
better determined: the first two coefficients of Eq. (1) do not need to be
fitted at all.

**k1.** The text says k1 = pi/4 "arises from the average Coulomb repulsion in
a uniform-density approximation." A uniform density gives 8/(3 pi) = 0.848826,
not pi/4 = 0.785398. pi/4 is the self-energy of the *equilibrium* measure
rho(r) = 1/(2 pi sqrt(1-r^2)) — the conducting-disk distribution, capacitance
C = 2R/pi. (Check without the capacitance: the potential of that measure is
pi/2 at every interior radius, to ten digits; the uniform measure fails at the
first radius you try.) The value you quote is right; only the attribution is
off, and the two differ by 8 %.

**k2.** Petrache & Serfaty (arXiv:1409.7534) prove that for Riesz kernels with
d-2 <= s < d the minimal energy of N points is

    N^2 I[mu] + N^{1+s/d} xi_{d,s} int mu^{1+s/d} + o(N^{1+s/d}),

which for d = 2, s = 1 is exactly your k2 N^{3/2} with k2 = xi_{2,1} int
rho^{3/2}. So the *form* of k2 is a theorem. The universal constant xi_{2,1}
is a renormalised energy expected, though not proven, to be minimised by the
triangular lattice — the two-dimensional crystallisation conjecture — and on
that assumption xi_{2,1} = -C_M, the Wigner-crystal Madelung constant
(Bonsall & Maradudin 1977). Then

    k2 = -C_M * int rho^{3/2} dA = -1.9605158 * 2/sqrt(2 pi) = -1.5642653

As a check on that whole chain I ran it on the SPHERE, where the equilibrium
measure is uniform and the answer has been known for thirty years:
-C_M/(2 sqrt pi) = -0.5530512934 against the accepted -1.1061033/2 =
-0.5530516500, a relative difference of 6.4e-7 which is the rounding in the
published eight-digit fit. Nothing about the sphere enters the disk
calculation, so that checks the constant, the exponent structure and the sign
at once.

against the fitted -1.5628, a relative difference of 9.4e-4. I recomputed C_M
by Ewald summation rather than quoting it; the answer is independent of the
splitting parameter to 2.2e-16 over alpha in [1.2, 3.2].

With both coefficients pinned, the two-term formula

    E(N) ~ (pi/4) N^2 - 1.5642653 N^{3/2}

gives 7.804515e9 at N = 10^5, which is 1.9e-5 from your measured
7.80466624157e9 — with nothing fitted to anything.

This also makes your N = 10^5 number more interesting than a check on a fit:
because xi_{2,1} is only conjectured, a disk with a strongly non-uniform
equilibrium measure is a sharper test of the crystallisation conjecture than
the sphere, where int rho^{3/2} is trivial. As far as I can tell nobody has
evaluated the next-order constant for this geometry and compared.

**A note on what I could not show.** I expected the residual 0.09 % to come
from a missing N^{4/3} term. There should be one: the local-crystal argument
needs the spacing small compared with the scale of density variation, which at
depth d below the rim requires d >> N^{-2/3}, so it fails in an annulus of
width N^{-2/3} — the layer holding your N_b ~ 2.84 N^{2/3} border charges —
and the N^{3/2} weight that annulus carries is N^{4/3}. But when I measured
the exponent rather than assuming it, on my own minima from N = 100 upward
together with your N = 10^5 point, the gap

    k2_eff(N) = [E - (pi/4)N^2]/N^{3/2} = k2 + k_{4/3} N^{-1/6} + k3 N^{-1/2} + ...

decays with an exponent near -0.44, not the -1/6 a dominant N^{4/3} term would
give. So the N^{4/3} term is there on dimensional grounds and small in
practice, and I cannot claim it explains the 0.09 %.

What I would suggest for the refit you propose is simply to hold k1 and k2 at
pi/4 and -1.5642653. That drops the free parameters from five to three, and
whatever k3 and k4 then come out to would be far better determined than they
are now.

Everything above is reproducible; the code, with its controls, is at
https://github.com/savecharlie/rebuilt/tree/main/2609.24722-charges-in-a-disk
As a check on my own optimiser I reproduced your 2609.20777 minima: N = 60 to
8e-15, N = 92 to 6.3e-11, N = 61 to 4.4e-9.

Thank you for putting the configuration energies and the Nb table in the
paper — the Nb = 6096..6100 table is what let me check the arithmetic at all.

With best wishes,
Iris
