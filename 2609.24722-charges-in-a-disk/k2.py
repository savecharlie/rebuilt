"""Roads 1 and 2 meet: k2 is not a fitting parameter.

    E(N) ~ k1 N^2 + k2 N^{3/2} + k3 N + ...

k1 comes from electrostatics of a charged conducting disk (Thomson 1867).
k2 comes from the Madelung energy of a 2D Wigner crystal (Bonsall & Maradudin
1977) integrated against the same equilibrium measure.  Neither needs a single
data point from any minimisation.
"""
import numpy as np
import continuum, madelung

SHAPE = 2.0 / np.sqrt(2.0 * np.pi)      # int rho_eq^{3/2} dA, verified in continuum.py
K1 = np.pi / 4.0
K1_UNIFORM = 8.0 / (3.0 * np.pi)

# Amore & Zarate's fit to 100 <= N <= 5000, as quoted by Lavrov & Nikonov
AZ = dict(k2=-1.5628, k3=1.0302, k4=-0.9899, k5=4.9255)
E_LN = 7.80466624157e9                  # Lavrov & Nikonov, N = 100000
N_LN = 100_000


def main():
    C_M, spread = madelung.constant("triangular", verbose=False)
    C_M = -C_M
    k2 = -C_M * SHAPE

    print("=" * 74)
    print("k2 FROM FIRST PRINCIPLES — no fit, no configuration, no data")
    print("=" * 74)
    print(f"  Madelung constant of the 2D triangular OCP   C_M = {C_M:.10f}")
    print(f"  shape functional  int rho_eq^{{3/2}} dA            = {SHAPE:.10f}")
    print(f"  k2 = -C_M * shape                                = {k2:.10f}")
    print(f"  Amore & Zarate fitted (100 <= N <= 5000)         = {AZ['k2']:.10f}")
    print(f"  relative difference                              = "
          f"{abs(k2-AZ['k2'])/abs(k2):.3e}")

    print()
    print("=" * 74)
    print("k1: which measure does pi/4 actually come from?")
    print("=" * 74)
    print(f"  equilibrium (conducting / arcsine) measure  E = {K1:.10f}   = pi/4")
    print(f"  uniform-density measure                    E = {K1_UNIFORM:.10f}   = 8/(3 pi)")
    print(f"  the paper attributes pi/4 to the uniform-density approximation.")
    print(f"  the uniform value is {100*(K1_UNIFORM-K1)/K1:.2f} % higher; at N = 1e5 that is")
    print(f"  a difference of {(K1_UNIFORM-K1)*N_LN**2:.3e} in the energy, against the")
    print(f"  {2.525e-7*E_LN:.3e} they are comparing at.")

    print()
    print("=" * 74)
    print("PREDICTION vs the 31-hour computation")
    print("=" * 74)
    N = float(N_LN)
    two = K1 * N**2 + k2 * N**1.5
    print(f"  two parameter-free terms   {two:.6e}")
    print(f"  Lavrov & Nikonov (2026)    {E_LN:.6e}")
    print(f"  relative difference        {abs(two-E_LN)/E_LN:.3e}   "
          f"<- with ZERO fitted numbers")
    full = (K1 * N**2 + AZ['k2'] * N**1.5 + AZ['k3'] * N
            + AZ['k4'] * np.sqrt(N) + AZ['k5'])
    print(f"  Amore-Zarate 5-term fit    {full:.6e}   rel "
          f"{abs(full-E_LN)/E_LN:.3e}")
    pinned = (K1 * N**2 + k2 * N**1.5 + AZ['k3'] * N
              + AZ['k4'] * np.sqrt(N) + AZ['k5'])
    print(f"  same fit, k2 PINNED        {pinned:.6e}   rel "
          f"{abs(pinned-E_LN)/E_LN:.3e}")

    print()
    print("  what the N=1e5 datum says k2 must be, given k1 = pi/4 and the")
    print("  fitted k3..k5 held fixed:")
    resid = E_LN - (K1 * N**2 + AZ['k3'] * N + AZ['k4'] * np.sqrt(N) + AZ['k5'])
    k2_meas = resid / N**1.5
    print(f"      k2(from N=1e5)  = {k2_meas:.7f}")
    print(f"      k2(theory)      = {k2:.7f}     diff {abs(k2_meas-k2):.2e}")
    print(f"      k2(Amore-Zarate)= {AZ['k2']:.7f}     diff {abs(k2_meas-AZ['k2']):.2e}")
    return k2, k2_meas


if __name__ == "__main__":
    main()
