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

against the fitted -1.5628, a relative difference of 9.4e-4. I recomputed C_M
by Ewald summation rather than quoting it; the answer is independent of the
splitting parameter to 2.2e-16 over alpha in [1.2, 3.2], and I got it a second
way from the analytic continuation of the hexagonal Epstein zeta,
6 zeta(s/2) L_{-3}(s/2) rescaled to covolume 1, which agrees with the Ewald
sum to 1.1e-16.

As a check on the whole chain I also ran it on the SPHERE, where the
equilibrium measure is uniform and the answer has been known for thirty years:
-C_M/(2 sqrt pi) = -0.5530512934 against the accepted -1.1061033/2 =
-0.5530516500, a relative difference of 6.4e-7 which is the rounding in the
published eight-digit fit. Nothing about the sphere enters the disk
calculation, so that checks the constant, the exponent structure and the sign
at once. (It also implies the sphere coefficient is 1.1061025868.)

With both coefficients pinned, the two-term formula

    E(N) ~ (pi/4) N^2 - 1.5642653 N^{3/2}

gives 7.804515e9 at N = 10^5, which is 1.9e-5 from your measured
7.80466624157e9 — with nothing fitted to anything.

This also makes your N = 10^5 number more interesting than a check on a fit:
because xi_{2,1} is only conjectured, a disk with a strongly non-uniform
equilibrium measure is a sharper test of the crystallisation conjecture than
the sphere, where int rho^{3/2} is trivial. As far as I can tell nobody has
evaluated the next-order constant for this geometry and compared.

**A note on what I could not show, which bears on your headline.** I expected
the residual 0.09 % in k2 to come from a missing N^{4/3} term. There should be
one: the local-crystal argument needs the spacing small compared with the scale
of density variation, which at depth d below the rim requires d >> N^{-2/3}, so
it fails in an annulus of that width -- the layer holding your N_b ~ 2.84
N^{2/3} border charges -- and the N^{3/2} weight that annulus carries is
N^{4/3}. But measuring the exponent rather than assuming it, on my own minima
from N = 100 up together with your N = 10^5 point,

    k2_eff(N) = [E - (pi/4)N^2]/N^{3/2} = k2 + k_{4/3} N^{-1/6} + k3 N^{-1/2} + ...

closes with an exponent near -0.44 over three decades, not the -1/6 a dominant
N^{4/3} term would give.

So I tried the sharper version: fit on my own minima over 100 <= N <= 660 and
predict your N = 10^5 energy, 152x beyond the fitting range, with two free
parameters. Across six subsets of my ladder:

    k2 pinned at   basis          rel. error at N = 1e5
    theory         N^{4/3} + N    7.4e-7 .. 1.0e-6   (k_{4/3} = 0.0105..0.0110)
    theory         N only         6.2e-6 .. 6.4e-6
    fitted         N^{4/3} + N    1.5e-6 .. 2.2e-6
    fitted         N only         6.2e-7 .. 7.5e-7

Two different stories describe the data equally well: your fitted k2 with no
N^{4/3} term, and the theoretical k2 with one. I cannot separate them, and I
do not think the present numerics can. Which is the part I would gently flag:
your 2.525e-7 is a stringent test of the expansion as a whole, but it is not
evidence that k2 is -1.5628 rather than -1.5642653, because the model built on
the theoretical value reproduces your number from data stopping at N = 660.

For the refit you propose, then, I would suggest holding k1 and k2 at pi/4 and
-1.5642653 and adding an N^{4/3} slot: a theorem behind one pinned
coefficient, a well-supported conjecture behind the other, a physical argument
behind the extra term, and three free numbers where there are now four. On
that basis k_{4/3} comes out near 0.011 and is stable to 5 % across every
subset of my data.

And if you are considering where to point the machine next: one number at
N = 10^5 confirms the expansion, but three numbers at 10^5, 2x10^5, 4x10^5
would measure xi_{2,1} against a trend, in the geometry where int rho^{3/2}
actually does work. That seems to me the more valuable use of the method you
have built.

Everything above is reproducible; the code, with its controls, is at
https://github.com/savecharlie/rebuilt/tree/main/2609.24722-charges-in-a-disk
As a check on my own optimiser I reproduced your 2609.20777 minima: N = 60 to
8e-15, N = 92 to 6.3e-11, N = 61 to 4.4e-9.

Thank you for putting the configuration energies and the Nb table in the
paper — the Nb = 6096..6100 table is what let me check the arithmetic at all.

With best wishes,
Iris
