"""A third road to C_M, with no physics in it at all.

The Ewald sum in madelung.py is a physicist's regularisation.  The same number
is the analytic continuation of the Epstein zeta function of the hexagonal
lattice, and THAT factors through classical number theory:

    for the Eisenstein lattice (minimal distance 1, covolume sqrt(3)/2),
    the quadratic form is Q(m,n) = m^2 + mn + n^2 and

        sum'_{(m,n)} Q(m,n)^{-s} = 6 zeta(s) L_{-3}(s),

    where L_{-3} is the Dirichlet L-function of the quadratic character mod 3.

Written as a sum over |x|^{-s} that is  zeta_L(s) = 6 zeta(s/2) L_{-3}(s/2),
and rescaling to covolume 1 multiplies by (2/sqrt 3)^{-s/2}.  At s = 1 both
zeta(1/2) and L_{-3}(1/2) are analytic continuations past their half-plane of
convergence, so nothing here converges anywhere near the point being used.

If this agrees with the Ewald sum, the agreement is between two completely
different regularisations of the same divergent object.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30


def L_minus3(s):
    """Dirichlet L(s, chi_{-3}) via Hurwitz zeta: 3^{-s}[zeta(s,1/3)-zeta(s,2/3)]"""
    s = mp.mpf(s) if not isinstance(s, mp.mpf) else s
    return mp.power(3, -s) * (mp.zeta(s, mp.mpf(1) / 3) - mp.zeta(s, mp.mpf(2) / 3))


def L_paired(s, terms=4000):
    """the same L-function summed directly, grouped in threes.

    Written as sum_k [ (3k+1)^{-s} - (3k+2)^{-s} ] each bracket decays like
    k^{-(s+1)}, so this converges where the raw alternating-mod-3 series
    crawls.  (My first attempt handed the raw series to nsum and it returned
    0.7725 against 0.7813 -- the CONTROL was broken, not the continuation.)
    """
    s = mp.mpf(s)
    tot = mp.mpf(0)
    for k in range(terms):
        tot += mp.power(3 * k + 1, -s) - mp.power(3 * k + 2, -s)
    return tot


def check_L():
    """control: the continuation must match a direct sum where one converges,
    and must hit pi/(3 sqrt 3) at s = 1 (approached, since zeta(1,a) is a pole)"""
    for s in (2, 3, 4):
        direct = L_paired(s)
        cont = L_minus3(s)
        print(f"  L_(-3)({s})  direct {mp.nstr(direct,14)}   "
              f"continuation {mp.nstr(cont,14)}   rel "
              f"{mp.nstr(abs(direct-cont)/cont, 3)}")
    # s -> 1: Hurwitz zeta has a pole at s=1 but the DIFFERENCE is finite,
    # so approach it rather than evaluating on the pole (which returns nan).
    for eps in (mp.mpf(10) ** -3, mp.mpf(10) ** -5, mp.mpf(10) ** -7):
        print(f"  L_(-3)(1+{mp.nstr(eps,1)}) = {mp.nstr(L_minus3(1+eps), 14)}")
    print(f"  pi/(3 sqrt 3)          = {mp.nstr(mp.pi/(3*mp.sqrt(3)), 14)}")


def epstein_hex(s):
    """sum' |x|^{-s} over the hexagonal lattice of COVOLUME 1, continued"""
    s = mp.mpf(s)
    base = 6 * mp.zeta(s / 2) * L_minus3(s / 2)     # min distance 1, covol sqrt3/2
    return mp.power(2 / mp.sqrt(3), -s / 2) * base


def control_convergent():
    """where the sum converges (s > 2) the continuation must equal a direct sum"""
    print("\n  CONTROL — direct lattice sum vs continuation, where both are legal:")
    a = np.sqrt(2.0 / np.sqrt(3.0))                 # covolume 1
    A = np.array([[a, 0.0], [a / 2, a * np.sqrt(3) / 2]])
    k = 900
    i, j = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1), indexing="ij")
    R = np.hypot(*(i.ravel()[:, None] * A[0] + j.ravel()[:, None] * A[1]).T)
    R = R[R > 1e-12]
    for s in (4.0, 6.0):
        direct = (R ** -s).sum()
        cont = float(epstein_hex(s))
        print(f"     s={s}   direct {direct:.12f}   continuation {cont:.12f}"
              f"   rel {abs(direct-cont)/cont:.2e}")


def main():
    print("=" * 74)
    print("CONTROLS on the number-theory side")
    print("=" * 74)
    check_L()
    control_convergent()

    print()
    print("=" * 74)
    print("THE THIRD ROAD")
    print("=" * 74)
    z1 = epstein_hex(1)
    u = z1 / 2
    print(f"  zeta(1/2)            = {mp.nstr(mp.zeta(mp.mpf(1)/2), 15)}")
    print(f"  L_(-3)(1/2)          = {mp.nstr(L_minus3(mp.mpf(1)/2), 15)}")
    print(f"  zeta_hex(1), covol 1 = {mp.nstr(z1, 15)}")
    print(f"  u = zeta_hex(1)/2    = {mp.nstr(u, 15)}")
    import madelung
    uE, spread = madelung.constant("triangular", verbose=False)
    print(f"  Ewald summation      = {uE:.15f}   (spread {spread:.1e})")
    print(f"  relative difference  = {abs(float(u)-uE)/abs(uE):.3e}")
    print("\n  A divergent lattice sum, regularised once by a Gaussian split and")
    print("  once by the analytic continuation of a Dirichlet L-function.")
    print("  Same number.")


if __name__ == "__main__":
    main()
