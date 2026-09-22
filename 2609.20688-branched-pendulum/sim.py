"""Full nonlinear branched double pendulum (Toda & Ooshida Eqs. 3-5), energy check, exact H_d and W.
Iris (Opus 5)."""
import numpy as np
from scipy.integrate import solve_ivp
from params import M1, M2, mu, G1, G2

def rhs(t, y):
    t1, t2, t3, w1, w2, w3 = y
    d2, d3 = t2 - t1, t3 - t1
    M = np.array([[M1, mu*np.cos(d2), mu*np.cos(d3)],
                  [mu*np.cos(d2), M2, 0.0],
                  [mu*np.cos(d3), 0.0, M2]])
    f = np.array([mu*np.sin(d2)*w2**2 + mu*np.sin(d3)*w3**2 - G1*np.sin(t1),
                  -mu*np.sin(d2)*w1**2 - G2*np.sin(t2),
                  -mu*np.sin(d3)*w1**2 - G2*np.sin(t3)])
    return np.concatenate([[w1, w2, w3], np.linalg.solve(M, f)])

def energy(Y):
    t1, t2, t3, w1, w2, w3 = Y
    K = 0.5*(M1*w1**2 + M2*(w2**2 + w3**2)) + mu*np.cos(t2-t1)*w1*w2 + mu*np.cos(t3-t1)*w1*w3
    U = G1*(1-np.cos(t1)) + G2*(2-np.cos(t2)-np.cos(t3))
    return K + U

def run(th1_deg, eps=0.01, T=8.0, n=8001, eps_mode="angle"):
    th1 = np.radians(th1_deg)
    y0 = [th1, eps, -eps, 0, 0, 0]     # children displaced antisymmetrically by +-eps (theta_D = eps)
    t = np.linspace(0, T, n)
    s = solve_ivp(rhs, (0, T), y0, t_eval=t, rtol=1e-11, atol=1e-12, method="DOP853")
    return s.t, s.y

def mode_vars(Y):
    t1, t2, t3, w1, w2, w3 = Y
    tM, tD = (t2+t3)/2, (t2-t3)/2
    wM, wD = (w2+w3)/2, (w2-w3)/2
    dM = tM - t1
    p1 = M1*w1 + 2*mu*np.cos(dM)*np.cos(tD)*wM - 2*mu*np.sin(dM)*np.sin(tD)*wD
    pM = 2*mu*np.cos(dM)*np.cos(tD)*w1 + 2*M2*wM
    pD = -2*mu*np.sin(dM)*np.sin(tD)*w1 + 2*M2*wD
    Hd = pD**2/(4*M2) + G2*np.sin(tD)**2
    D = 2*M1*M2 - 4*mu**2*np.cos(dM)**2
    WK = mu*np.sin(dM)*(-M2*p1 + mu*np.cos(dM)*pM)*pD**2/(M2**2*D)
    WP = G2/M2*(1-np.cos(tM))*pD*tD
    return dict(tM=tM, tD=tD, dM=dM, p1=p1, pM=pM, pD=pD, Hd=Hd, WK=WK, WP=WP, W=WK+WP)

if __name__ == "__main__":
    t, Y = run(47.8, T=20)
    E = energy(Y)
    print("energy drift over 20 s: %.2e relative" % ((E.max()-E.min())/E[0]))
