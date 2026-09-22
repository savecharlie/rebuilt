"""Branched double pendulum, Toda & Ooshida arXiv:2609.20688 — parameters and linear modes.
Rebuilt by Iris (Opus 5)."""
import numpy as np
g = 9.80
l1, l2, l12 = 0.084, 0.043, 0.270
m1, m2 = 0.12880, 0.10531
I1, I2 = 1.88e-3, 0.516e-3
M1 = I1 + m1*l1**2 + 2*m2*l12**2
M2 = I2 + m2*l2**2
mu = m2*l12*l2
G1 = m1*g*l1 + 2*m2*g*l12
G2 = m2*g*l2

def linear_freqs(M1=M1, M2=M2, mu=mu, G1=G1, G2=G2):
    Mi = np.array([[M1, 2*mu], [2*mu, 2*M2]])
    Ki = np.diag([G1, 2*G2])
    w2 = np.sort(np.linalg.eigvals(np.linalg.solve(Mi, Ki)).real)
    fp, fm = np.sqrt(w2)/(2*np.pi)
    fd = np.sqrt(G2/M2)/(2*np.pi)
    return fp, fm, fd

if __name__ == "__main__":
    print(f"M1={M1:.5g} M2={M2:.5g} mu={mu:.5g} G1={G1:.5g} G2={G2:.5g}")
    print("f+ f- fd = %.3f %.3f %.3f  (paper 0.88 1.58 1.26)" % linear_freqs())
