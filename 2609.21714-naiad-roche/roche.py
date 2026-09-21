"""Independent check of Agrusa, Cuk, Nesvorny & Minton, arXiv:2609.21714, Section 2
("Naiad's proximity to the Roche limit").

What they did: the net acceleration at Naiad's sub-Neptune point using MacCullagh's
formula for the self-gravity (Eq. 1), and a minimum Drucker-Prager cohesion (Fig. 2,
friction angle 35 deg) after Holsapple (2006, 2008).

What this does differently: EXACT gravity of a homogeneous ellipsoid (Carlson R_D
integrals), not a far-field expansion evaluated at the surface, and Holsapple's
volume-averaged stresses derived here from the tensor virial identity
    integral sigma_ij dV = integral x_j f_i dV   (traction-free surface),
so for a body force f_i = rho k_i x_i the average stress is rho k_i a_i^2 / 5.

Validated BEFORE it says anything about Naiad:
  1. sphere: A_i = 2/3 exactly, surface g = (4/3) pi G rho r.
  2. the classical fluid Roche ellipsoid: the maximum of mu = Omega^2/(pi G rho) over
     shapes whose averaged stress is isotropic must be 0.09009 at axis ratios ~(0.50, 0.47),
     i.e. the fluid Roche limit d/R = 2.455 (rho_p/rho)^(1/3)  (Chandrasekhar 1969).

Iris, fire 290.
"""
import numpy as np
from scipy.special import elliprd
from scipy.optimize import brentq, minimize

G = 6.67430e-11


def Acoef(a, b, c):
    """Chandrasekhar's index-symbol A_i for semi-axes a>=b>=c (sum = 2).
    Interior potential U = -pi G rho (I - sum A_i x_i^2); g_i = -2 pi G rho A_i x_i."""
    A1 = (2 / 3) * a * b * c * elliprd(b * b, c * c, a * a)
    A2 = (2 / 3) * a * b * c * elliprd(c * c, a * a, b * b)
    A3 = (2 / 3) * a * b * c * elliprd(a * a, b * b, c * c)
    return np.array([A1, A2, A3])


def surface_g_exact(rho, a, b, c, GMp, d):
    """Net outward acceleration at the sub-planet point (x=a), synchronous rotation.
    tidal (+centrifugal) along the planet line: 3 Omega^2 x in the Hill approximation;
    written out exactly like their Eq. 1 for the planet terms."""
    A = Acoef(a, b, c)
    g_self = -2 * np.pi * G * rho * A[0] * a
    m = rho * 4 / 3 * np.pi * a * b * c
    tide = GMp / (d - a) ** 2 - GMp / d ** 2
    cent = G * (GMp / G + m) / d ** 3 * a
    return g_self + tide + cent


def surface_g_maccullagh(rho, a, b, c, GMp, d):
    """Their Eq. 1 verbatim: MacCullagh far-field self-gravity evaluated at r=a."""
    m = rho * 4 / 3 * np.pi * a * b * c
    Am = m * (b * b + c * c) / 5
    Bm = m * (a * a + c * c) / 5
    Cm = m * (a * a + b * b) / 5
    g_self = -G * m / a ** 2 - 1.5 * G / a ** 4 * (Bm + Cm - 2 * Am)
    tide = GMp / (d - a) ** 2 - GMp / d ** 2
    cent = G * (GMp / G + m) / d ** 3 * a
    return g_self + tide + cent


def avg_stress(rho, a, b, c, Om2):
    """Volume-averaged stresses (tension +), synchronous satellite, Hill tides:
    body-force coefficients k = (3 Om2, 0, -Om2)
    (tide 2,-1,-1 plus centrifugal 1,1,0; y cancels -- my first draft missed that and
    the Roche validation caught it at mu=0.0946) - 2 pi G rho A_i."""
    A = Acoef(a, b, c)
    k = np.array([3 * Om2, 0.0, -Om2]) - 2 * np.pi * G * rho * A
    return rho * k * np.array([a * a, b * b, c * c]) / 5


def dp_min_cohesion(sig, phi_deg=35.0):
    """Holsapple (2007/2008) Drucker-Prager: sqrt(J2) <= k - s I1 (tension +).
    Returns (k_min, c_MohrCoulomb) in Pa, 0 if cohesionless material holds."""
    phi = np.radians(phi_deg)
    s = 2 * np.sin(phi) / (np.sqrt(3) * (3 - np.sin(phi)))
    I1 = sig.sum()
    dev = sig - I1 / 3
    J2 = 0.5 * (dev ** 2).sum()
    k = np.sqrt(J2) + s * I1
    k = max(k, 0.0)
    c = k * np.sqrt(3) * (3 - np.sin(phi)) / (6 * np.cos(phi))
    return k, c


# ---------------------------------------------------------------- validation
def validate():
    A = Acoef(1.0, 1.0, 1.0)
    assert np.allclose(A, 2 / 3, atol=1e-14), A
    assert abs(Acoef(3.0, 2.0, 1.0).sum() - 2) < 1e-13
    # sphere surface gravity
    rho, r = 1000.0, 1e4
    g = -2 * np.pi * G * rho * Acoef(r, r, r)[0] * r
    assert abs(g / (-(4 / 3) * np.pi * G * rho * r) - 1) < 1e-13

    # Roche fluid ellipsoid: isotropic averaged stress <=> k_i a_i^2 equal.
    # For given (b/a, c/a) the two equalities fix mu; take the family, find max mu.
    def resid(q):
        b, c = q
        Ai = Acoef(1.0, b, c)
        # (3mu - 2A1) = (mu - 2A2) b^2 = (-mu - 2A3) c^2 with mu = Om2/(pi G rho)
        mu1 = (2 * Ai[0] - 2 * Ai[1] * b * b) / 3                 # x = y
        mu2 = (2 * Ai[0] - 2 * Ai[2] * c * c) / (3 + c * c)       # x = z
        return mu1, mu2

    # along the equilibrium curve mu1 = mu2; parametrize by b, solve c, maximize mu
    best = None
    for b in np.linspace(0.30, 0.95, 1301):
        f = lambda c: resid((b, c))[0] - resid((b, c))[1]
        cs = np.linspace(0.05, b - 1e-6, 400)
        v = [f(c) for c in cs]
        for i in range(len(cs) - 1):
            if np.sign(v[i]) != np.sign(v[i + 1]):
                c = brentq(f, cs[i], cs[i + 1], xtol=1e-14)
                mu = resid((b, c))[0]
                if best is None or mu > best[0]:
                    best = (mu, b, c)
    mu, b, c = best
    dR = (4 / (3 * mu)) ** (1 / 3)
    print(f"[validate] Roche fluid ellipsoid: mu_max = {mu:.5f} at b/a={b:.3f}, c/a={c:.3f};"
          f" d/R = {dR:.4f} (rho_p/rho)^(1/3)")
    assert abs(mu - 0.09009) < 2e-5 and abs(dR - 2.455) < 1e-3
    print("[validate] sphere A_i=2/3, sum A_i = 2, sphere g exact: OK")


# ---------------------------------------------------------------- Naiad
GMp = 6835100e9            # their Table 3, m^3 s^-2
Rp = 24764e3               # their Table 3
d = 48227e3                # Naiad semi-major axis (their Table 1: 0.482e5 km)
a0, b0, c0 = 48e3, 30e3, 26e3   # Karkoschka 2003, as quoted
GMn = 0.008e9


def naiad():
    V = 4 / 3 * np.pi * a0 * b0 * c0
    rho = GMn / G / V
    print(f"density from GM and shape: {rho:.0f} kg/m^3 (paper: 800 +- 480); "
          f"R_eff = {(a0*b0*c0)**(1/3)/1e3:.1f} km (paper ~33)")
    rhoP = GMp / G / (4 / 3 * np.pi * Rp ** 3)
    for r in (800.0, rho):
        print(f"  delta = d/Rp (rho/rhoP)^(1/3) at rho={r:.0f}: {d/Rp*(r/rhoP)**(1/3):.3f} (paper 1.54)")

    Om2 = (GMp + GMn) / d ** 3
    print("\nnet outward acceleration at sub-Neptune point, nominal shape [cm/s^2]:")
    for r in (800.0, 1280.0):
        ge = surface_g_exact(r, a0, b0, c0, GMp, d) * 100
        gm = surface_g_maccullagh(r, a0, b0, c0, GMp, d) * 100
        print(f"  rho={r/1000:.2f}: exact {ge:+.4f}   MacCullagh (their Eq.1) {gm:+.4f}")
    rc_e = brentq(lambda r: surface_g_exact(r, a0, b0, c0, GMp, d), 100, 5000)
    rc_m = brentq(lambda r: surface_g_maccullagh(r, a0, b0, c0, GMp, d), 100, 5000)
    print(f"  density to hold loose material at sub-Neptune point: exact {rc_e/1000:.3f}, "
          f"MacCullagh {rc_m/1000:.3f} g/cc")

    print("\nDrucker-Prager minimum cohesion, phi=35, nominal shape (Fig. 2):")
    rows = []
    for r in (300, 570, 800, 1000, 1100, 1150, 1200, 1300):
        sig = avg_stress(r, a0, b0, c0, Om2)
        k, c = dp_min_cohesion(sig)
        rows.append((r, k, c))
        print(f"  rho={r/1000:.2f}: sig=({sig[0]/1e3:+.2f},{sig[1]/1e3:+.2f},{sig[2]/1e3:+.2f}) kPa"
              f"  k_min={k/1e3:7.3f} kPa   c_MC={c/1e3:7.3f} kPa")
    rz = brentq(lambda r: dp_min_cohesion(avg_stress(r, a0, b0, c0, Om2))[0] - 1e-3, 500, 3000)
    print(f"  cohesion reaches zero at rho = {rz/1000:.3f} g/cc (Fig. 2 nominal curve: ~1.15)")
    return rows


if __name__ == "__main__":
    validate()
    naiad()
