"""Road 1: the two continuum constants, computed from scratch.

The paper writes E_GS(N) ~ k1 N^2 + k2 N^{3/2} + ... and says
"the leading term k1 = pi/4 arises from the average Coulomb repulsion in a
uniform-density approximation."

Two different measures are in play and they give different numbers:
  * the UNIFORM measure on the disk  ->  8/(3 pi) = 0.848826...
  * the EQUILIBRIUM (conducting / arcsine) measure -> pi/4 = 0.785398...
Only the second is pi/4.  This file computes both without assuming either.

Every constant here is checked two ways: a closed form and a quadrature that
knows nothing about it.
"""
import numpy as np
from scipy import integrate, special

# --- potential of an axisymmetric surface density sigma(s) on the unit disk,
#     evaluated at in-plane radius r.  Angular integral done in closed form:
#       int_0^{2pi} dtheta / sqrt(r^2+s^2-2 r s cos theta) = (4/(r+s)) K(m),
#       m = 4 r s / (r+s)^2        (K = complete elliptic integral, m-convention)
def potential(sigma, r, n=4000):
    """phi(r) = int sigma(s) / |x - y| dA(y)."""
    def integrand(s):
        if s <= 0.0:
            return 0.0
        m = 4.0 * r * s / (r + s) ** 2
        # K has a LOG singularity at s = r (m -> 1).  It is integrable, but in
        # floating point 4rs can exceed (r+s)^2 by an ulp and ellipk returns
        # nan; clamping keeps the integrand finite and the log integrable.
        if m >= 1.0:
            m = 1.0 - 1e-15
        return sigma(s) * (4.0 / (r + s)) * special.ellipk(m) * s
    # split at s = r so the quadrature sees the singularity as an endpoint
    pts = [0.0, r, 1.0] if 0.0 < r < 1.0 else [0.0, 1.0]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b > a:
            v, _ = integrate.quad(integrand, a, b, limit=n, epsabs=1e-12,
                                  epsrel=1e-12, weight='alg-loga' if False else None)
            tot += v
    return tot


def potential_arcsine(r):
    """phi(r) for the equilibrium measure, with the rim singularity removed.

    rho_eq(s) s ds = sin(t)/(2 pi) dt under s = sin t, so the ONLY singularity
    left in the integrand is the log of K at s = r, which quad handles once
    the interval is split there.  Written separately because the generic
    `potential` above integrates in s and blows up when the evaluation point
    sits on the rim.
    """
    ts = arcsin_clip(r)

    def integrand(t):
        s = np.sin(t)
        if r + s == 0.0:
            return 0.0
        m = 4.0 * r * s / (r + s) ** 2
        if m >= 1.0:
            m = 1.0 - 1e-15
        return (4.0 / (r + s)) * special.ellipk(m) * np.sin(t) / (2.0 * np.pi)

    pts = sorted({0.0, ts, np.pi / 2})
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b > a:
            v, _ = integrate.quad(integrand, a, b, limit=400,
                                  epsabs=1e-13, epsrel=1e-13)
            tot += v
    return tot


def arcsin_clip(r):
    return np.arcsin(min(max(r, 0.0), 1.0))


def self_energy(sigma, arcsine=False):
    """E = (1/2) int sigma(x) phi(x) dA(x).

    For the equilibrium measure the radial weight rho*2*pi*s diverges like
    (1-s^2)^{-1/2} at the rim.  The substitution s = sin t kills it EXACTLY:
        rho_eq(s) 2 pi s ds = [1/(2 pi cos t)] 2 pi sin t cos t dt = sin t dt,
    so the integrand becomes phi(sin t) sin t with no singularity at all.
    """
    if arcsine:
        f = lambda t: potential_arcsine(np.sin(t)) * np.sin(t)
        v, err = integrate.quad(f, 0.0, np.pi / 2, limit=200,
                                epsabs=1e-11, epsrel=1e-11)
    else:
        f = lambda s: sigma(s) * potential(sigma, s) * 2.0 * np.pi * s
        v, err = integrate.quad(f, 0.0, 1.0, limit=400,
                                epsabs=1e-11, epsrel=1e-11)
    return 0.5 * v, 0.5 * err


def main():
    print("=" * 74)
    print("CONTROL 1 — the equilibrium measure is the one that flattens the potential")
    print("=" * 74)
    # arcsine / conducting-disk measure, total charge 1
    rho_eq = lambda s: 1.0 / (2.0 * np.pi * np.sqrt(max(1.0 - s * s, 0.0)))
    norm, _ = integrate.quad(lambda s: rho_eq(s) * 2 * np.pi * s, 0, 1,
                             points=[1.0], limit=200)
    print(f"  normalisation  int rho_eq dA = {norm:.12f}   (want 1)")
    print("  phi(r) inside the disk, which must be CONSTANT for a conductor:")
    for r in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0):
        print(f"     r={r:4.2f}   phi = {potential_arcsine(r):.10f}" + ("   <- ON the rim" if r >= 1.0 else ""))
    print(f"  pi/2 = {np.pi/2:.10f}     <- Q/C with C = 2R/pi (Thomson 1867)")

    E_eq, err = self_energy(rho_eq, arcsine=True)
    print(f"\n  E[rho_eq]  = {E_eq:.12f}  +-{err:.1e}   (quadrature)")
    print(f"  pi/4       = {np.pi/4:.12f}              (closed form)")
    print(f"  rel diff   = {abs(E_eq-np.pi/4)/(np.pi/4):.3e}")

    print()
    print("=" * 74)
    print("CONTROL 2 — the UNIFORM measure gives a different number")
    print("=" * 74)
    rho_u = lambda s: 1.0 / np.pi
    print("  phi(r) for the uniform disk is NOT constant:")
    for r in (0.0, 0.3, 0.6, 0.9):
        print(f"     r={r:4.2f}   phi = {potential(rho_u, r):.10f}")
    print(f"  phi(0) closed form = 2                      -> {2.0:.10f}")
    E_u, err_u = self_energy(rho_u)
    print(f"\n  E[uniform] = {E_u:.12f}  +-{err_u:.1e}   (quadrature)")
    print(f"  8/(3 pi)   = {8/(3*np.pi):.12f}              (closed form)")
    print(f"  rel diff   = {abs(E_u-8/(3*np.pi))/(8/(3*np.pi)):.3e}")
    print(f"\n  uniform exceeds equilibrium by "
          f"{100*(E_u-np.pi/4)/(np.pi/4):.2f} %  — and pi/4 is the LOWER one,")
    print("  as it must be: the equilibrium measure is the minimiser.")

    print()
    print("=" * 74)
    print("CONTROL 3 — the shape functional that sets k2")
    print("=" * 74)
    f = lambda s: rho_eq(s) ** 1.5 * 2 * np.pi * s
    v, e = integrate.quad(f, 0, 1, points=[1.0], limit=400,
                          epsabs=1e-13, epsrel=1e-13)
    closed = 2.0 / np.sqrt(2.0 * np.pi)
    print(f"  int rho_eq^{{3/2}} dA = {v:.12f} +-{e:.1e}   (quadrature)")
    print(f"  2/sqrt(2 pi)        = {closed:.12f}              (closed form)")
    print(f"  rel diff            = {abs(v-closed)/closed:.3e}")
    return dict(E_eq=E_eq, E_uniform=E_u, shape=v)


if __name__ == "__main__":
    main()
