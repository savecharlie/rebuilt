"""Figures of rotating bodies, and the local surface metric on them.

Rebuild support for arXiv:2609.22474 (Lock 2026), "Tectonics on early Earth
driven by Earth's changing shape during tidal recession of the Moon".

CONVENTION.  Everything here uses GEOCENTRIC LATITUDE phi, measured from the
equator, in radians, 0 <= phi <= pi/2.  The paper uses theta measured from the
pole (its Eq. 4 is r sin(theta) dphi).  So paper-theta = pi/2 - our-phi.  This
file never uses colatitude.  Mixing the two is exactly the kind of dull fault
tools/dull.py exists for, so it is stated once and obeyed.

The surface model is the paper's own (its Sec. 2.1, last paragraph): the crust
has negligible strength, so the surface IS a gravitational equipotential, and
strain is accommodated LOCALLY -- a material point keeps its latitude and the
strain at that latitude is the change in the local line elements

    dl_lat = sqrt(r^2 + (dr/dphi)^2) dphi        (paper Eq. 3)
    dl_lon = r cos(phi) dlambda                  (paper Eq. 4)
    dA     = dl_lat dl_lon                       (paper Eq. 5)

Written by Iris (Opus 5, 1M), fire 297, Sep 23 2026.
"""
import numpy as np
from scipy import optimize, special

# ---------------------------------------------------------------- oblate spheroid

def spheroid_radius(phi, a, c):
    """Geocentric radius of the oblate spheroid x^2/a^2 + z^2/c^2 = 1."""
    return 1.0 / np.sqrt((np.cos(phi) / a) ** 2 + (np.sin(phi) / c) ** 2)


def spheroid_drdphi(phi, a, c):
    """d r / d phi, analytic."""
    s, co = np.sin(phi), np.cos(phi)
    q = (co / a) ** 2 + (s / c) ** 2
    dq = 2 * s * co * (1.0 / c ** 2 - 1.0 / a ** 2)
    return -0.5 * q ** -1.5 * dq


def spheroid_axes(f, R=1.0):
    """Equatorial a and polar c of a spheroid of flattening f=(a-c)/a and
    volume equal to a sphere of radius R."""
    u = 1.0 - f                      # c/a
    a = R * u ** (-1.0 / 3.0)        # a^2 c = R^3
    return a, a * u

# ---------------------------------------------------------------- local metric

def line_elements(phi, r, drdphi):
    """Paper Eqs. 3 and 4, per unit dphi and per unit dlambda."""
    return np.sqrt(r ** 2 + drdphi ** 2), r * np.cos(phi)

# ---------------------------------------------- exact answers, for validation only

def spheroid_volume(a, c):
    return 4.0 / 3.0 * np.pi * a ** 2 * c


def spheroid_area(a, c):
    """Exact surface area of an oblate spheroid (a > c)."""
    e = np.sqrt(1.0 - (c / a) ** 2)
    return 2 * np.pi * a ** 2 * (1.0 + (1.0 - e ** 2) / e * np.arctanh(e))


def meridian_quarter_arc(a, c):
    """Exact pole-to-equator meridian arc of the generating ellipse = a E(e^2)."""
    e2 = 1.0 - (c / a) ** 2
    return a * special.ellipe(e2)          # scipy takes m = e^2

# ------------------------------------------------------- Maclaurin spheroid (exact)

def maclaurin_omega2(e):
    """omega^2 / (pi G rho) for a homogeneous Maclaurin spheroid of eccentricity e.
    Chandrasekhar (1969), Eq. 32: 2 (1-e^2)^{1/2} (3-2e^2) arcsin(e)/e^3 - 6(1-e^2)/e^2.
    """
    e = np.asarray(e, dtype=float)
    return (2 * np.sqrt(1 - e ** 2) * (3 - 2 * e ** 2) * np.arcsin(e) / e ** 3
            - 6 * (1 - e ** 2) / e ** 2)

# --------------------------------------------------- Roche model (point mass + spin)

def roche_surface(phi, q):
    """Equipotential of a point mass plus rotation: GM/r + omega^2 r^2 cos^2(phi)/2 = const.

    q = omega^2 R_p^3 / (GM) with R_p the POLAR radius (where the centrifugal term
    vanishes, so the polar radius fixes the constant).  Returns r/R_p.
    Exact: 1/x + q x^2 cos^2(phi) / 2 = 1  ->  cubic in x.
    """
    cc = np.cos(phi) ** 2
    out = np.empty_like(np.atleast_1d(phi), dtype=float)
    for i, c2 in enumerate(np.atleast_1d(cc)):
        if c2 * q < 1e-15:
            out[i] = 1.0
            continue
        # solve  q c2 x^3 / 2 - x + 1 = 0 for the small positive root
        roots = np.roots([q * c2 / 2.0, 0.0, -1.0, 1.0])
        real = np.sort([z.real for z in roots if abs(z.imag) < 1e-9 and z.real > 0])
        out[i] = real[0]
    return out if np.ndim(phi) else out[0]

# --------------------------------------------------------- Radau-Darwin (validation)

def radau_darwin_f(q, moi_factor):
    """Flattening from the Darwin-Radau relation.

        f = lam * q,   lam = (5/2) / [1 + (5/2)^2 (1 - (3/2) C/MR^2)^2]

    with q = omega^2 a^3 / GM.  Checked two ways in validate.py: the homogeneous
    limit C/MR^2 = 2/5 must give exactly f = 5q/4 (the small-e Maclaurin result),
    and Earth's C/MR^2 = 0.3307 must land within a fraction of a per cent of the
    observed 1/298.257.  Both are asserted there; neither is taken on trust here.
    """
    lam = 2.5 / (1.0 + (2.5 * (1.0 - 1.5 * moi_factor)) ** 2)
    return q * lam
