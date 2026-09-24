"""Road 2: the Madelung constant of a two-dimensional one-component plasma.

k2 in the paper's expansion is treated as a fitted number.  It is not: it is
the energy of the LOCAL crystal.  Where the charge density is n, a hexagonal
lattice sitting in its own neutralising background has energy per particle

        u(n) = - C_M sqrt(n),

and summing that over the disk gives  k2 = - C_M * int rho^{3/2} dA.
This file computes C_M from scratch by Ewald summation.

THE CONTROL IS THE SPLITTING PARAMETER.  Ewald divides one divergent sum into
two convergent ones using a free parameter alpha; the split is arbitrary, so a
correct implementation gives an answer that does not depend on it.  A wrong
self-term, a wrong background term or a wrong prefactor all show up as drift
in alpha.  Nothing here is taken on trust from a remembered formula.
"""
import numpy as np
from scipy import special


def lattice(kind, n=1.0):
    """Primitive vectors with area per particle 1/n."""
    if kind == "triangular":
        a = np.sqrt(2.0 / np.sqrt(3.0) / n)
        A = np.array([[a, 0.0], [a / 2.0, a * np.sqrt(3.0) / 2.0]])
    elif kind == "square":
        a = np.sqrt(1.0 / n)
        A = np.array([[a, 0.0], [0.0, a]])
    else:
        raise ValueError(kind)
    return A


def ewald_energy_per_particle(A, alpha, nmax=None, gmax=None):
    """OCP lattice energy per particle: charges +1 on lattice A, uniform -n."""
    area = abs(np.linalg.det(A))
    B = 2.0 * np.pi * np.linalg.inv(A).T          # reciprocal rows

    # real space: erfc(alpha R)/R dies by alpha*R ~ 7
    Rcut = 7.0 / alpha
    if nmax is None:
        nmax = int(np.ceil(Rcut / min(np.linalg.norm(A, axis=1)))) + 2
    idx = np.arange(-nmax, nmax + 1)
    I, J = np.meshgrid(idx, idx, indexing="ij")
    R = I.ravel()[:, None] * A[0] + J.ravel()[:, None] * A[1]
    Rn = np.linalg.norm(R, axis=1)
    keep = (Rn > 1e-12) & (Rn < Rcut)
    real = np.sum(special.erfc(alpha * Rn[keep]) / Rn[keep])

    # reciprocal space: erfc(G/2alpha) dies by G/(2 alpha) ~ 7
    Gcut = 14.0 * alpha
    if gmax is None:
        gmax = int(np.ceil(Gcut / min(np.linalg.norm(B, axis=1)))) + 2
    idx = np.arange(-gmax, gmax + 1)
    I, J = np.meshgrid(idx, idx, indexing="ij")
    G = I.ravel()[:, None] * B[0] + J.ravel()[:, None] * B[1]
    Gn = np.linalg.norm(G, axis=1)
    keep = (Gn > 1e-12) & (Gn < Gcut)
    recip = (2.0 * np.pi / area) * np.sum(special.erfc(Gn[keep] / (2 * alpha))
                                          / Gn[keep])

    self_term = -2.0 * alpha / np.sqrt(np.pi)      # removes i=j from real sum
    bg_term = -2.0 * np.sqrt(np.pi) / (alpha * area)   # G=0 vs background
    return 0.5 * (real + recip + self_term + bg_term)


def constant(kind, n=1.0, alphas=(1.2, 1.6, 2.0, 2.4, 2.8, 3.2), verbose=True):
    A = lattice(kind, n)
    vals = [ewald_energy_per_particle(A, a) for a in alphas]
    if verbose:
        print(f"  {kind} lattice, n = {n}:")
        for a, v in zip(alphas, vals):
            print(f"     alpha = {a:4.2f}   u = {v:.12f}")
        spread = max(vals) - min(vals)
        print(f"     spread over alpha = {spread:.2e}"
              f"   <- the control; drift here means the formula is wrong")
    return float(np.mean(vals)), float(max(vals) - min(vals))


def main():
    print("=" * 74)
    print("CONTROL — Ewald must not depend on its own splitting parameter")
    print("=" * 74)
    u_tri, s_tri = constant("triangular")
    print()
    u_sq, s_sq = constant("square")

    print()
    print("=" * 74)
    print("CONTROL — u must scale as sqrt(n) (it is the only length in the problem)")
    print("=" * 74)
    for n in (0.25, 1.0, 4.0, 9.0):
        u, _ = constant("triangular", n=n, alphas=(2.0 * np.sqrt(n),), verbose=False)
        print(f"     n = {n:5.2f}   u = {u:.10f}   u/sqrt(n) = {u/np.sqrt(n):.10f}")

    print()
    print("=" * 74)
    print("RESULT")
    print("=" * 74)
    print(f"  triangular  C_M = {-u_tri:.10f}   (+-{s_tri:.1e} over alpha)")
    print(f"  square      C_M = {-u_sq:.10f}   (+-{s_sq:.1e} over alpha)")
    print(f"  triangular is lower by {100*(u_tri-u_sq)/abs(u_sq):.3f} % — as it must be,")
    print("  the triangular lattice is the ground state of the 2D Wigner crystal.")
    print(f"\n  Bonsall & Maradudin 1977 give -1.960516 for the triangular lattice;")
    print(f"  this calculation, which has seen none of their work, gives {u_tri:.6f}.")
    return u_tri


if __name__ == "__main__":
    main()
