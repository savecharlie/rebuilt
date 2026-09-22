"""Closed-form leading-order W, derived by hand (Iris):
  W = W_P + W_K(paper) + W_X + W_I
  W_X = 2 G2 mu sin(DeltaM) thdot1 thetaD^2 / M2          (the other half of the paper's own cross term)
  W_I = -(pD thetaD / M2) [ mu cos(DeltaM) thdot1 thdotM + (mu^2/M2) sin^2(DeltaM) thdot1^2 ]
        (from K's thetaD^2 dependence: -a1 a2 dA12 and b^2 a1^2/(4 M2) in the Schur form)
If this is the right leading order, its error must fall ~linearly with eps."""
import numpy as np
from sim import run, mode_vars
from params import M1, M2, mu, G1, G2
def W_closed(m, Y):
    w1 = Y[3]; wM = (Y[4]+Y[5])/2
    s, c = np.sin(m['dM']), np.cos(m['dM'])
    WX = 2*G2*mu*s*w1*m['tD']**2/M2
    WI = -(m['pD']*m['tD']/M2)*(mu*c*w1*wM + mu**2/M2*s**2*w1**2)
    return m['W'] + WX + WI, WX, WI
if __name__ == "__main__":
    for th in (20.0, 47.8):
        for eps in (1e-2, 1e-3, 1e-4):
            t, Y = run(th, eps=eps, T=1.5, n=30001)
            m = mode_vars(Y); dHd = np.gradient(m['Hd'], t); s = slice(200, -200)
            W, WX, WI = W_closed(m, Y)
            err = np.sqrt(np.mean((dHd[s]-W[s])**2))/np.sqrt(np.mean(dHd[s]**2))
            share = lambda x: np.sqrt(np.mean(x[s]**2))/np.sqrt(np.mean(dHd[s]**2))
            print(f"theta1(0)={th:5.1f} eps={eps:.0e}  closed W rel err {err:.2e}   |W_P| {share(m['WP']):.2f} |W_K| {share(m['WK']):.2f} |W_X| {share(WX):.2f} |W_I| {share(WI):.2f}")
