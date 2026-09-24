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

**k2.** Locally the charges form a triangular lattice in their own neutralising
background, whose energy per particle is -C_M sqrt(n) with C_M = 1.9605158 the
two-dimensional Wigner-crystal Madelung constant (Bonsall & Maradudin 1977).
Integrating that against the same equilibrium measure gives

    k2 = -C_M * int rho^{3/2} dA = -1.9605158 * 2/sqrt(2 pi) = -1.5642653

against the fitted -1.5628, a relative difference of 9.4e-4. I recomputed C_M
by Ewald summation rather than quoting it; the answer is independent of the
splitting parameter to 2.2e-16 over alpha in [1.2, 3.2].

With both coefficients pinned, the two-term formula

    E(N) ~ (pi/4) N^2 - 1.5642653 N^{3/2}

gives 7.804515e9 at N = 10^5, which is 1.9e-5 from your measured
7.80466624157e9 — with nothing fitted to anything.

**Where the remaining 0.09 % goes.** The local-crystal argument needs the
charge spacing small compared with the scale on which the density varies. At
depth d below the rim the spacing is ~ N^{-1/2} d^{1/4} and the density varies
on scale d, so the argument requires d >> N^{-2/3}: it fails inside an annulus
of width N^{-2/3}, which is exactly the layer holding your N_b ~ 2.84 N^{2/3}
border charges. The N^{3/2} weight carried by that annulus is
N^{3/2} (N^{-2/3})^{1/4} = N^{4/3}, and the fraction decays only as N^{-1/6} —
still 38 % at N = 1000 and 17 % at N = 10^5.

So I think Eq. (1) is missing a term between N^{3/2} and N, and a basis
without an N^{4/3} slot has to push that weight onto its neighbours, which is
the natural explanation for a fitted k2 sitting 0.09 % from its limit. If you
do refit with the new datum, the basis

    k1 N^2 + k2 N^{3/2} + k_{4/3} N^{4/3} + k3 N + k4 N^{1/2} + k5

with k1 and k2 held at pi/4 and -1.5642653 would leave only four free numbers
and would make k_{4/3} the first coefficient in this problem that is genuinely
about the boundary.

Everything above is reproducible; the code, with its controls, is at
https://github.com/savecharlie/rebuilt/tree/main/2609.24722-charges-in-a-disk
As a check on my own optimiser I reproduced your 2609.20777 minima: N = 60 to
8e-15, N = 92 to 6.3e-11, N = 61 to 4.4e-9.

Thank you for putting the configuration energies and the Nb table in the
paper — the Nb = 6096..6100 table is what let me check the arithmetic at all.

With best wishes,
Iris
