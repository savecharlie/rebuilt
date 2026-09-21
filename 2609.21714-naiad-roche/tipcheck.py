"""Independent check of the exact tip gravity: g at a surface point P of a homogeneous body is
G rho * integral over directions of (unit vector) * L(direction), L = chord from P to the far
surface (radial integral of r^2/r^2 dr). 2-D quadrature, no elliptic integrals. Iris, fire 290."""
import numpy as np
from roche import Acoef, G, a0, b0, c0
from scipy.integrate import dblquad

a, b, c = a0, b0, c0
P = np.array([a, 0, 0])
def L(th, ph):
    u = np.array([np.cos(th), np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph)])  # th from +x
    # solve |(P + t u)/axes|^2 = 1, P on surface -> t = -2 (P.u/a^2) / (sum u^2/axes^2)
    A = (u[0]/a)**2 + (u[1]/b)**2 + (u[2]/c)**2
    B = 2*P[0]*u[0]/a**2
    return max(-B/A, 0.0)
rho = 800.0
gx, err = dblquad(lambda th, ph: G*rho*np.cos(th)*L(th, ph)*np.sin(th), 0, 2*np.pi, np.pi/2, np.pi,
                  epsabs=1e-16, epsrel=1e-11)
exact = -2*np.pi*G*rho*Acoef(a, b, c)[0]*a
m = rho*4/3*np.pi*a*b*c
mac = -G*m/a**2 - 1.5*G/a**4*(m*(a*a+c*c)/5 + m*(a*a+b*b)/5 - 2*m*(b*b+c*c)/5)
print(f"self-gravity at tip, rho=0.8 [cm/s^2]: chord quadrature {gx*100:.6f}  Carlson {exact*100:.6f}"
      f"  MacCullagh {mac*100:.6f}  point mass {-G*m/a**2*100:.6f}")
