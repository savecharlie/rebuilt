"""Which powers belong in the expansion, and what that does to k2.

Amore & Zarate fit  E = k1 N^2 + k2 N^{3/2} + k3 N + k4 N^{1/2} + k5.
There is a reason to expect a term between N^{3/2} and N, and if it is there
and has no slot in the fit, its weight has to go somewhere — mostly into k2.

Where the local-density argument breaks.  At distance d from the rim the
equilibrium density is rho ~ 1/(2 pi sqrt(2 d)), so the spacing between
charges is a ~ (N rho)^{-1/2} ~ N^{-1/2} d^{1/4}.  A local crystal argument
needs the spacing to be small compared with the distance over which the
density changes, which here is d itself:

        N^{-1/2} d^{1/4} << d     <=>     d >> N^{-2/3}.

So the argument fails in a rim annulus of width N^{-2/3} — which is exactly
the layer that holds the Nb ~ N^{2/3} rim charges Amore & Zarate found
empirically.  The whole N^{3/2} weight carried by that annulus is

        N^{3/2} * int_{rim layer} rho^{3/2} dA  ~  N^{3/2} * (N^{-2/3})^{1/4}
                                               =  N^{3/2 - 1/6}  =  N^{4/3}.

That is the size of the error in k2's term, and it says two things:
  (1) k2 = -C_M int rho^{3/2} is exact as the N^{3/2} coefficient;
  (2) the next term is N^{4/3}, not N.
"""
import numpy as np
from scipy import integrate

RHO = lambda s: 1.0 / (2.0 * np.pi * np.sqrt(max(1.0 - s * s, 1e-300)))


def rim_weight(N):
    """fraction of int rho^{3/2} that lives inside d < N^{-2/3} of the rim"""
    d = N ** (-2 / 3)
    f = lambda t: RHO(np.sin(t)) ** 0.5 * np.sin(t)      # rho^{3/2} dA = rho^{1/2} sin t dt
    tot, _ = integrate.quad(f, 0, np.pi / 2, limit=200)
    tail, _ = integrate.quad(f, np.arcsin(max(0.0, 1 - d)), np.pi / 2, limit=200)
    return tail, tot


def design(N, powers):
    return np.stack([np.asarray(N, float) ** p for p in powers], axis=1)


def fit(N, E, powers, pin=None):
    """least squares in the given powers; `pin` maps power -> fixed coefficient"""
    N = np.asarray(N, float); E = np.asarray(E, float)
    free = [p for p in powers if pin is None or p not in pin]
    y = E.copy()
    if pin:
        for p, c in pin.items():
            y = y - c * N ** p
    if free:
        A = design(N, free)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        out = dict(zip(free, c))
    else:
        out = {}                      # every power pinned: nothing to fit
    if pin:
        out.update(pin)
    res = E - sum(out[p] * N ** p for p in powers)
    return out, res


def main():
    print("=" * 74)
    print("How much of the N^{3/2} weight sits in the layer where LDA fails")
    print("=" * 74)
    print(f"  {'N':>8}  {'rim width':>11}  {'tail/total':>11}   implied N^{{4/3}}-size")
    for N in (1e2, 1e3, 1e4, 1e5, 1e6, 1e8):
        tail, tot = rim_weight(N)
        print(f"  {N:8.0e}  {N**(-2/3):11.3e}  {tail/tot:11.3e}"
              f"   {tail/tot * 1.5643 * N**1.5:12.4e}")
    rs = [rim_weight(N)[0] / rim_weight(N)[1] for N in (1e2, 1e3, 1e4, 1e5, 1e6)]
    ratios = [b / a for a, b in zip(rs[:-1], rs[1:])]
    print(f"\n  measured ratio per decade: "
          f"{', '.join(f'{r:.4f}' for r in ratios)}")
    print(f"  10^(-1/6) = {10**(-1/6):.4f}   <- so the fraction really does decay")
    print("  as N^{-1/6}, which is the N^{4/3} term hiding inside the N^{3/2} one.")
    print(f"\n  READ THE COLUMN, NOT THE EXPONENT: at N = 1000, {rs[1]*100:.0f} % of the")
    print(f"  N^{{3/2}} weight sits in the annulus where the argument is not valid;")
    print(f"  at N = 1e6 it is still {rs[4]*100:.0f} %.  N^{{-1/6}} is a very slow decay.")


if __name__ == "__main__":
    main()
