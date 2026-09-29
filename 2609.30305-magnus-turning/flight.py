"""Independent integrator for a ball under gravity, quadratic drag and a
vertical-axis Magnus force -- Eq. (1) of Gaur, arXiv:2609.30305.

    dv/dt = -g zhat - kD |v| v + kL |v| (what x v)

Nothing in this module knows about any of the paper's results. The heading psi
and the path length s are carried as state variables and their derivatives are
computed from the ACCELERATION, so that testing dpsi/ds = kL is a genuine test
and not a restatement of the integrator.

Iris, fire 308 (29 Sep 2026).
"""
import numpy as np
from scipy.integrate import solve_ivp

G = 9.80665          # m/s^2, CODATA standard gravity (the paper's value: checked
                     # against its v0 = 86.89 m/s at theta = 45 deg for a baseball)
RHO = 1.225          # kg/m^3, ISA sea level


def coeffs(m, r, CD, CL, rho=RHO):
    """Return (kD, kL) in 1/m from ball mass, radius and coefficients."""
    A = np.pi * r * r
    return rho * CD * A / (2 * m), rho * CL * A / (2 * m)


def rhs(t, y, kD, kL, spin, g=G):
    """y = [x, y, z, vx, vy, vz, psi, s]. spin = +1 or -1 (sign of omega_z)."""
    vx, vy, vz = y[3], y[4], y[5]
    V = np.sqrt(vx * vx + vy * vy + vz * vz)
    # zhat x v = (-vy, vx, 0); the vertical component of v contributes nothing.
    mx, my = -spin * vy, spin * vx
    ax = -kD * V * vx + kL * V * mx
    ay = -kD * V * vy + kL * V * my
    az = -g - kD * V * vz
    W2 = vx * vx + vy * vy
    # psidot from the acceleration, independently of anything the paper claims.
    psidot = (vx * ay - vy * ax) / W2 if W2 > 0 else 0.0
    return [vx, vy, vz, ax, ay, az, psidot, V]


def fly(v0, theta, kD, kL, spin=+1, bearing=0.0, z0=0.0,
        t_max=400.0, rtol=1e-12, atol=1e-12, events=None, dense=True, g=G):
    """Integrate a launch. theta in radians (elevation), bearing in radians."""
    w0 = v0 * np.cos(theta)
    y0 = [0.0, 0.0, z0,
          w0 * np.cos(bearing), w0 * np.sin(bearing), v0 * np.sin(theta),
          bearing, 0.0]
    return solve_ivp(rhs, (0.0, t_max), y0, args=(kD, kL, spin, g),
                     method="DOP853", rtol=rtol, atol=atol,
                     dense_output=dense, events=events, max_step=np.inf)


def ground_event(z_target=0.0, rising_ok=False):
    """Terminal event: crossing back down through z = z_target."""
    def ev(t, y, *a):
        return y[2] - z_target
    ev.terminal = True
    ev.direction = 0.0 if rising_ok else -1.0
    return ev


def turn_event(psi_target):
    """Terminal event: heading reaches psi_target (unwrapped)."""
    def ev(t, y, *a):
        return y[6] - psi_target
    ev.terminal = True
    ev.direction = 1.0
    return ev


# Representative balls. Coefficients are the paper's "representative constants",
# not measurements -- kept here only so the reproduction uses the same inputs.
BALLS = {
    #            m (kg)   r (m)    CD     CL
    "baseball":  (0.145,  0.0366, 0.35,  0.20),
    "soccer":    (0.430,  0.110,  0.25,  0.25),
    "tabletennis": (0.0027, 0.020, 0.45, 0.30),
    "golf":      (0.0459, 0.0213, 0.25,  0.25),
    "tennis":    (0.0577, 0.0335, 0.55,  0.25),
    "frisbee":   (0.175,  0.135,  0.08,  0.24),
}
