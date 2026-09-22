"""
A closed form for the critical rank fraction.

(Q_0 * Q_0)(t) can be done in closed form, because
    Qt_0(s)^2 = (1 - 2 e^{-d sqrt(s/D)} + e^{-2 d sqrt(s/D)}) / s^2
and  L^{-1}{ e^{-a sqrt(s)} / s^2 } = (t + a^2/2) erfc(a/(2 sqrt t))
                                       - a sqrt(t/pi) e^{-a^2/(4t)} =: h_a(t).

With tau = D t / d^2 (so d = D = 1 without loss), g(tau) = 0 becomes

    (tau+1) erfc(1/(2 sqrt tau)) - (tau+2) erfc(1/sqrt tau)
        = 2 sqrt(tau/pi) ( e^{-1/(4 tau)} - e^{-1/tau} )

and  phi_c = 1 - erf( 1 / (2 sqrt(tau_c)) ).
"""
import numpy as np
from scipy.special import erf, erfc
from scipy.optimize import brentq
import phic   # the quadrature version, for cross-checking

def h(a, t):
    return (t + a*a/2)*erfc(a/(2*np.sqrt(t))) - a*np.sqrt(t/np.pi)*np.exp(-a*a/(4*t))

def conv_closed(t):
    return t - 2*h(1.0, t) + h(2.0, t)

def g_closed(t):
    return (t+1)*erfc(1/(2*np.sqrt(t))) - (t+2)*erfc(1/np.sqrt(t)) \
           - 2*np.sqrt(t/np.pi)*(np.exp(-1/(4*t)) - np.exp(-1/t))

if __name__ == "__main__":
    print("closed form vs quadrature -- two roads, same door")
    print(f"{'t':>8} {'conv quad':>15} {'conv closed':>15} {'g quad':>15} {'g closed':>15}")
    for t in [0.05, 0.2, 0.7439, 1.6, 10.0, 100.0]:
        cq, _ = phic.conv(t)
        print(f"{t:8.4f} {cq:15.10f} {conv_closed(t):15.10f} {phic.g(t):+15.10f} {g_closed(t):+15.10f}")

    tc = brentq(g_closed, 0.3, 1.5, xtol=1e-15, rtol=8.9e-16)
    phic_val = 1 - erf(1/(2*np.sqrt(tc)))
    print()
    print(f"  tau_c = {tc:.12f}   d^2/D")
    print(f"  phi_c = {phic_val:.12f}")
    print(f"  paper : 0.412   |  quadrature route: 0.4123101755")
    print(f"  agreement between routes: {abs(phic_val - 0.4123101755):.2e}")
