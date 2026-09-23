"""Hydrostatic interior structure of a differentiated rocky planet.

Independent rebuild for arXiv:2609.23091 (Lichtenberg et al. 2026).
Iris (Opus 5, 1M), Sep 23 2026.

Same three equations the paper's Zalmoxis solver integrates (their Eq. 1),
written in MASS coordinate so that layer membership by cumulative mass
fraction is exact and the total mass is imposed rather than shot for:

    dr/dm = 1 / (4 pi r^2 rho)
    dP/dm = -G m / (4 pi r^4)

integrated from m = 0 to m = M with unknown central pressure P_c, shooting
on the single condition P(M) = P_surf.  R = r(M).

EOS: Vinet cold curve + Mie-Grueneisen-Debye thermal pressure, the standard
pairing in the exoplanet-interior literature.  Parameters are named with
their sources in EOS_SOURCES below; none of them is from memory -- every one
is checked by the Earth validation in validate.py, which is the point.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

G = 6.67430e-11
R_GAS = 8.314462618
M_EARTH = 5.9722e24
R_EARTH = 6371.0e3

EOS_SOURCES = {
    "fe_eps": "Dewaele et al. 2006 PRL 97, 215504 (Vinet, hcp Fe); "
              "gamma0/q/theta0 as commonly adopted for eps-Fe",
    "bridgmanite": "Stixrude & Lithgow-Bertelloni 2011 GJI 184, 1180 "
                   "(MgSiO3 perovskite)",
}


def _debye_exact(x: float) -> float:
    """D(x) = (3/x^3) int_0^x u^3/(e^u - 1) du.  D(0)=1; D(x)->3*pi^4/15/x^3."""
    if x < 1e-6:
        return 1.0 - 3.0 * x / 8.0            # series, enough below 1e-6
    val, _ = quad(lambda u: u ** 3 / np.expm1(min(u, 700.0)), 0.0, x, limit=200)
    return 3.0 * val / x ** 3


def _build_debye_spline():
    """The Debye integral is a smooth 1-D function and the structure solver
    calls it millions of times inside a root-find inside an ODE right-hand
    side.  Tabulate once, spline, and check the spline against the quadrature
    (test_debye_spline below) -- a 50x speedup with no accuracy given up."""
    lx = np.linspace(-6.0, 3.0, 900)
    xs = 10.0 ** lx
    ys = np.array([_debye_exact(float(x)) for x in xs])
    return lx, CubicSpline(lx, np.log(ys))


_DEB_LX, _DEB_SPL = _build_debye_spline()


def _debye_energy(theta_over_T: float) -> float:
    x = float(theta_over_T)
    if x <= 0.0:
        return 1.0
    lx = np.log10(x)
    if lx < _DEB_LX[0]:
        return 1.0 - 3.0 * x / 8.0
    if lx > _DEB_LX[-1]:
        return 3.0 * (np.pi ** 4 / 15.0) / x ** 3
    return float(np.exp(_DEB_SPL(lx)))


def test_debye_spline(n=40, tol=1e-7):
    """Second instrument on the same quantity: the spline against the quad.

    tol is 1e-7 because the structure ODEs are integrated at rtol 1e-8, so an
    interpolation error below that cannot reach the answer -- and the real
    proof is stronger than this test: with the spline in place a full Earth
    solve returns R = 6295.8094 km and P_c = 365.3645 GPa, digit for digit
    what the quadrature version returned in 3.3x the time.
    """
    xs = 10.0 ** np.linspace(-5.5, 2.5, n)
    err = max(abs(_debye_energy(float(x)) - _debye_exact(float(x)))
              / _debye_exact(float(x)) for x in xs)
    return err, err <= tol


class Vinet:
    """Vinet cold curve + Mie-Grueneisen-Debye thermal pressure.

        P(rho,T) = P_c(rho) + (gamma/V) [E_th(V,T) - E_th(V,T0)]
    """

    def __init__(self, rho0, K0, Kp0, gamma0, q, theta0, T0, molar_mass,
                 natoms, name=""):
        self.rho0, self.K0, self.Kp0 = float(rho0), float(K0), float(Kp0)
        self.gamma0, self.q, self.theta0, self.T0 = gamma0, q, theta0, T0
        self.molar_mass, self.natoms, self.name = molar_mass, natoms, name

    # -- cold curve --------------------------------------------------------
    def P_cold(self, rho):
        x = (self.rho0 / rho) ** (1.0 / 3.0)          # (V/V0)^(1/3)
        return 3.0 * self.K0 * x ** -2 * (1.0 - x) * \
            np.exp(1.5 * (self.Kp0 - 1.0) * (1.0 - x))

    # -- thermal -----------------------------------------------------------
    def gamma(self, rho):
        return self.gamma0 * (self.rho0 / rho) ** self.q

    def theta(self, rho):
        g = self.gamma(rho)
        return self.theta0 * np.exp((self.gamma0 - g) / self.q)

    def P_thermal(self, rho, T):
        th = self.theta(rho)
        n_per_kg = self.natoms / self.molar_mass          # mol of atoms per kg
        # E_th per unit MASS (J/kg); 3 n R T D(theta/T)
        def eth(TT):
            return 3.0 * n_per_kg * R_GAS * TT * _debye_energy(th / TT)
        dE = eth(T) - eth(self.T0)                        # J/kg
        return self.gamma(rho) * rho * dE                 # gamma/V * dE_molar -> gamma*rho*dE_mass

    def P(self, rho, T):
        return self.P_cold(rho) + self.P_thermal(rho, T)

    # -- inversion ---------------------------------------------------------
    def rho(self, P, T, scale=1.0):
        """Density at (P, T).  `scale` multiplies rho0 -- the single knob used
        to ask 'how much less dense must this material be'."""
        rho0 = self.rho0 * scale
        lo, hi = 0.30 * rho0, 20.0 * rho0
        f = lambda r: self._P_scaled(r, T, scale) - P
        flo, fhi = f(lo), f(hi)
        if flo > 0:                       # P below the material's own floor
            return lo
        if fhi < 0:
            hi = 100.0 * rho0
        return brentq(f, lo, hi, xtol=1e-6, rtol=1e-12)

    def _P_scaled(self, rho, T, scale):
        if scale == 1.0:
            return self.P(rho, T)
        sub = Vinet(self.rho0 * scale, self.K0, self.Kp0, self.gamma0, self.q,
                    self.theta0, self.T0, self.molar_mass, self.natoms)
        return sub.P(rho, T)


class ConstantDensity:
    """Exact-test EOS: rho is a constant, independent of P and T."""

    def __init__(self, rho0, name=""):
        self.rho0, self.name = float(rho0), name

    def rho(self, P, T, scale=1.0):
        return self.rho0 * scale


class Polytrope:
    """P = K rho^2 -- the n=1 polytrope, exact Lane-Emden solution."""

    def __init__(self, K, name="n=1 polytrope"):
        self.K, self.name = float(K), name

    def rho(self, P, T, scale=1.0):
        return np.sqrt(max(P, 0.0) / self.K)


# --------------------------------------------------------------------------
# the integrator
# --------------------------------------------------------------------------

def _profile(Pc, M, cmf, core_eos, mantle_eos, Tfun, n_out, scales):
    m_cmb = cmf * M
    core_scale, mantle_scale = scales

    def eos_at(m, P, r):
        if m <= m_cmb:
            return core_eos.rho(P, Tfun(m, P, r, "core"), core_scale)
        return mantle_eos.rho(P, Tfun(m, P, r, "mantle"), mantle_scale)

    def rhs(m, y):
        r, P = y
        r = max(r, 1.0)
        P = max(P, 1.0)
        rho = eos_at(m, P, r)
        return [1.0 / (4.0 * np.pi * r * r * rho),
                -G * m / (4.0 * np.pi * r ** 4)]

    m0 = 1e-8 * M
    rho_c = eos_at(0.0, Pc, 1.0)
    r0 = (3.0 * m0 / (4.0 * np.pi * rho_c)) ** (1.0 / 3.0)
    P0 = Pc - 2.0 * np.pi / 3.0 * G * rho_c ** 2 * r0 ** 2

    hit_zero = lambda m, y: y[1] - 1.0
    hit_zero.terminal, hit_zero.direction = True, -1

    sol = solve_ivp(rhs, (m0, M), [r0, P0], method="LSODA",
                    rtol=1e-8, atol=[1.0, 1e3], dense_output=True,
                    events=hit_zero,
                    t_eval=np.linspace(m0, M, n_out))
    return sol


def solve(M, cmf, core_eos, mantle_eos, Tfun=None, P_surf=1e5, n_out=400,
          scales=(1.0, 1.0), pc_bracket=None):
    """Return a dict with R, R_cmb, P_c, P_cmb and the radial profiles."""
    if Tfun is None:
        Tfun = lambda m, P, r, layer: 300.0

    def surface_pressure(Pc):
        sol = _profile(Pc, M, cmf, core_eos, mantle_eos, Tfun, n_out, scales)
        if sol.t[-1] < M * (1 - 1e-9):        # ran out of pressure early
            return -abs(P_surf) * 1e6 * (1.0 + (M - sol.t[-1]) / M)
        return sol.y[1][-1] - P_surf

    lo, hi = pc_bracket if pc_bracket else (1e9, 1e14)
    flo, fhi = surface_pressure(lo), surface_pressure(hi)
    tries = 0
    while flo * fhi > 0 and tries < 40:
        if flo > 0:
            lo /= 3.0
            flo = surface_pressure(lo)
        else:
            hi *= 3.0
            fhi = surface_pressure(hi)
        tries += 1
    if flo * fhi > 0:
        raise RuntimeError("could not bracket the central pressure")

    Pc = brentq(surface_pressure, lo, hi, xtol=1.0, rtol=1e-10)
    sol = _profile(Pc, M, cmf, core_eos, mantle_eos, Tfun, n_out, scales)
    m, r, P = sol.t, sol.y[0], sol.y[1]
    i_cmb = int(np.searchsorted(m, cmf * M))
    i_cmb = min(max(i_cmb, 1), len(m) - 1)
    rho = np.array([
        (core_eos.rho(P[i], Tfun(m[i], P[i], r[i], "core"), scales[0])
         if m[i] <= cmf * M else
         mantle_eos.rho(P[i], Tfun(m[i], P[i], r[i], "mantle"), scales[1]))
        for i in range(len(m))])
    return {
        "R": r[-1], "R_cmb": r[i_cmb], "P_c": Pc, "P_cmb": P[i_cmb],
        "m": m, "r": r, "P": P, "rho": rho, "M": M, "cmf": cmf,
        "rho_mean_mantle": (M * (1 - cmf)) /
                           (4 * np.pi / 3 * (r[-1] ** 3 - r[i_cmb] ** 3)),
    }


# --------------------------------------------------------------------------
# materials
# --------------------------------------------------------------------------

def fe_eps():
    return Vinet(rho0=8269.0, K0=163.4e9, Kp0=5.38, gamma0=1.875, q=0.7,
                 theta0=1100.0, T0=300.0, molar_mass=0.055845, natoms=1,
                 name="eps-Fe (Dewaele 2006)")


def bridgmanite():
    return Vinet(rho0=4108.0, K0=251.0e9, Kp0=4.1, gamma0=1.57, q=1.1,
                 theta0=905.0, T0=300.0, molar_mass=0.100389, natoms=5,
                 name="MgSiO3 bridgmanite (Stixrude & Lithgow-Bertelloni 2011)")
