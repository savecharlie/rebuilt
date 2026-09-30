"""
Finite-volume switching dynamics + uniform-competition growth, as specified in
Fassoni, arXiv:2609.32585 Section 6.4 (subsec:numerics).

    d u_i/dt = (J_{i+1/2} - J_{i-1/2})/h + g(U) r_i u_i
    J_{i+1/2} = a_{i+1/2} (u_{i+1}/psi_{i+1} - u_i/psi_i)/h,   a = D*psi

Properties the scheme is supposed to have, and which are ASSERTED here rather
than assumed: columns of the generator sum to zero (mass conservation), psi is
an exact stationary state, positivity is preserved.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig


class Switching:
    def __init__(self, ell, N, D, v=None, P=None):
        """D, v, P are callables of x.  Give either v or P (with v = -P')."""
        self.ell, self.N = ell, N
        self.h = h = ell / N
        self.x = (np.arange(N) + 0.5) * h
        self.xf = np.arange(1, N) * h          # interior faces
        self.D = D

        # psi propto exp(int_0^x v/D).  Build on a fine grid, sample at centres
        # and faces so that a = D*psi is consistent with psi.
        M = 200_000
        xs = np.linspace(0.0, ell, M + 1)
        if v is None:
            eps = ell * 1e-7
            vv = -(P(xs + eps) - P(xs - eps)) / (2 * eps)
        else:
            vv = v(xs)
        integ = np.concatenate([[0.0], np.cumsum(0.5 * (vv[1:] / D(xs[1:]) + vv[:-1] / D(xs[:-1])) * np.diff(xs))])
        integ -= integ.max()
        w = np.exp(integ)
        Z = np.trapz(w, xs)
        self.psi_of = lambda q: np.interp(q, xs, w) / Z
        self.psi = self.psi_of(self.x)
        # Normalise psi in the DISCRETE measure (h*sum), not the continuum one.
        # The continuum normalisation leaves an O(h^2) midpoint-vs-trapezoid
        # defect -- 1.43e-7 at N=300 -- and that defect shows up as a spurious
        # floor on ||u(t)-U* psi||, identical for every r, which looks exactly
        # like "the theorem is only true to 1e-7".  It is the ruler, not the world.
        self.quad_defect = self.h*self.psi.sum() - 1.0
        self.psi = self.psi/(self.h*self.psi.sum())
        self.psif = self.psi_of(self.xf)/(1.0+self.quad_defect)
        self.a = D(self.xf) * self.psif

        # generator
        A = np.zeros((N, N))
        c = self.a / h**2
        idx = np.arange(N - 1)
        A[idx, idx + 1] += c / self.psi[idx + 1]
        A[idx, idx] -= c / self.psi[idx]
        A[idx + 1, idx] += c / self.psi[idx]
        A[idx + 1, idx + 1] -= c / self.psi[idx + 1]
        self.A = A

    # ---- invariants of the discretisation ------------------------------
    def col_sums(self):
        return np.abs(self.A.sum(axis=0)).max()

    def psi_stationarity(self):
        return np.abs(self.A @ self.psi).max()

    def mass(self, u):
        return self.h * np.sum(u, axis=-1)

    def spectral_gap(self):
        mu = eig(self.A, right=False)
        mu = np.sort(mu.real)          # all real here (self-adjoint in weighted inner product)
        return -mu[-2], mu

    # ---- evolution ------------------------------------------------------
    def run(self, u0, r, g, T, n_out=1200, mode="uniform", gp=None, rtol=1e-10):
        """mode 'uniform': du = Au + g(U) r u.   mode 'additive': du = Au + (r - U) u.
        gp is g'(U); supplied so BDF gets an ANALYTIC Jacobian (a numerical one
        costs N extra RHS evaluations per step and made this run take minutes)."""
        h = self.h
        rv = r(self.x) if callable(r) else r
        ones = np.ones(self.N)

        if mode == "uniform":
            def rhs(t, u):
                return self.A @ u + g(h*u.sum()) * rv * u
            def jac(t, u):
                U = h*u.sum()
                J = self.A + np.diag(g(U)*rv)
                if gp is not None:
                    J = J + h*gp(U)*np.outer(rv*u, ones)
                return J
        else:
            def rhs(t, u):
                return self.A @ u + (rv - h*u.sum()) * u
            def jac(t, u):
                U = h*u.sum()
                return self.A + np.diag(rv - U) - h*np.outer(u, ones)

        ts = np.linspace(0, T, n_out)
        sol = solve_ivp(rhs, (0, T), u0, t_eval=ts, method="BDF",
                        jac=jac, rtol=rtol, atol=1e-14)
        assert sol.success, sol.message
        return sol.t, sol.y.T
