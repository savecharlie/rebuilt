"""Effective surface gravity of a rotating figure, pole vs equator -- and what
that allows for crustal thickness.

Lock (2026) argues the first crust was thicker at the equator than at the poles
because the depth of melting scales as 1/g, and puts the extreme ratio at "a
factor of ten" (his Sec. 4.1 and plain-language summary).  A tenfold thickness
ratio under h ~ 1/g requires g_pole / g_equator = 10.

That ratio is computable exactly for the two limiting structures, so it can be
bracketed without any equation of state.

MACLAURIN (uniform).  Chandrasekhar's index symbols; the interior field is
linear, g_i = 2 pi G rho A_i x_i, with

    A_i = a1 a2 a3 Int_0^inf du / [(a_i^2 + u) Delta],   sum A_i = 2

computed here by quadrature and checked against the sum rule and against the
sphere (A_i = 2/3).

ROCHE (point mass).  g = |grad(GM/r + w^2 r^2 cos^2 phi / 2)| on the figure.

Run:  python3 gravity.py
"""
import sys, pathlib
import numpy as np
from scipy import integrate, optimize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "tools"))
from figures import maclaurin_omega2
import dull

def index_symbols(a, c):
    """A_1 = A_2, A_3 for an oblate spheroid, by quadrature."""
    def A(ai2):
        f = lambda u: 1.0 / ((ai2 + u) * np.sqrt((a * a + u) ** 2 * (c * c + u)))
        return a * a * c * integrate.quad(f, 0, np.inf, limit=400)[0]
    return A(a * a), A(c * c)


def maclaurin_gravity(e):
    """(g_eq_eff, g_pole) / (2 pi G rho R) for a Maclaurin spheroid of eccentricity e,
    volume-normalised to a sphere of radius R=1."""
    c_over_a = np.sqrt(1 - e * e)
    a = (1.0 / c_over_a) ** (1.0 / 3.0)       # a^2 c = 1
    c = a * c_over_a
    A1, A3 = index_symbols(a, c)
    # omega^2 / (pi G rho) from the Maclaurin relation -> omega^2 / (2 pi G rho) = w2/2
    w2 = maclaurin_omega2(e) / 2.0
    return A1 * a - w2 * a, A3 * c


def roche_gravity(q):
    """(g_eq_eff, g_pole) / (GM/R_pol^2) for the Roche figure."""
    x = optimize.brentq(lambda z: q * z ** 3 / 2 - z + 1, 1.0, 1.5, xtol=1e-15)
    return 1.0 / x ** 2 - q * x, 1.0


if __name__ == "__main__":
    print(__doc__.split("Run:")[0])
    with dull.watch() as w:
        print("=" * 72)
        print("0. CONTROL: index symbols")
        A1, A3 = index_symbols(1.0, 1.0)
        print(f"   sphere: A1={A1:.10f} A3={A3:.10f}  (exact 2/3 = {2/3:.10f})")
        assert abs(A1 - 2 / 3) < 1e-9 and abs(A3 - 2 / 3) < 1e-9
        for e in (0.3, 0.6, 0.866):
            c_over_a = np.sqrt(1 - e * e); a = (1 / c_over_a) ** (1 / 3)
            A1, A3 = index_symbols(a, a * c_over_a)
            print(f"   e={e:.3f}: 2A1+A3 = {2*A1+A3:.12f}  (sum rule = 2)")
            assert abs(2 * A1 + A3 - 2) < 1e-9

        print()
        print("=" * 72)
        print("1. MACLAURIN: effective gravity ratio pole/equator")
        print(f"{'a/c':>7} {'e':>7} | {'g_pole/g_eq':>12}   {'(h_eq/h_pol under h~1/g)':>26}")
        for ratio in (1.0001, 1.2, 1.5, 1.7161, 2.0, 2.5, 2.7198, 3.0):
            e = np.sqrt(1 - 1 / ratio ** 2)
            ge, gp = maclaurin_gravity(e)
            print(f"{ratio:7.4f} {e:7.4f} | {gp/ge:12.4f}   {gp/ge:26.2f}")

        print()
        print("   The Maclaurin sequence's own landmarks:")
        res = optimize.minimize_scalar(lambda e: -maclaurin_omega2(e), bounds=(0.5, 0.999),
                                       method="bounded", options={"xatol": 1e-13})
        print(f"     max omega^2 at e={res.x:.6f}, a/c={1/np.sqrt(1-res.x**2):.4f}")
        print(f"     Jacobi bifurcation at e=0.812670, a/c={1/np.sqrt(1-0.812670**2):.4f}")
        # CAREFUL.  A first pass took argmin|g_eq| on a grid and printed "g_eff -> 0
        # at a/c = 22.4".  That is a small number read as a zero.  g_eq = 2 pi G rho
        # A3 c^2/a is strictly positive at every finite a/c, so a UNIFORM body has no
        # corotation limit in the g_eq -> 0 sense at all; it has a maximum-omega
        # point instead.  The decay is what is real, so measure its exponent.
        # And a second catch, in the same run: I first wrote "g_eq ~ (a/c)^-1" here
        # from a head calculation that forgot volume conservation moves a as well as
        # c.  dull.scaling() measured -1.56 against my claimed -1.0 and said believe
        # the data.  The truth is g_eq = 2 pi G rho A3 c^2/a with a^2 c fixed, so
        # c^2/a = (a/c)^-5/3, and A3 -> 2: the exponent is -5/3, reached as -1.6607
        # by a/c = 2236.  Written down because the prose was wrong and the tool was
        # the only thing that knew.
        ee = np.array([0.99, 0.995, 0.998, 0.999, 0.9995])
        ac = 1 / np.sqrt(1 - ee ** 2)
        ge = np.array([maclaurin_gravity(x)[0] for x in ee])
        print(f"     uniform body: g_eq at a/c = {np.array2string(ac, precision=1)}")
        print(f"                        = {np.array2string(ge, precision=5)}")
        print(f"     pre-asymptotic exponent here: {dull.scaling(ac, ge, 'g_eq (near)'):.4f}"
              f"  (A3 is still only {index_symbols((1/np.sqrt(1-0.999**2))**(1/3)*np.sqrt(1-0.999**2)**0, np.sqrt(1-0.999**2))[1]:.3f} of its limit 2)")
        ef = np.array([0.9999, 0.99999, 0.999999, 0.9999999])
        acf = 1 / np.sqrt(1 - ef ** 2)
        gef = np.array([maclaurin_gravity(x)[0] for x in ef])
        print(f"     asymptotic:  {dull.scaling(acf, gef, 'g_eq (far)', expect=-5/3)}"
              f"  <- g_eq = 2 pi G rho A3 c^2/a ~ (a/c)^-5/3; never zero, so NO g_eq=0 limit")

        print()
        print("=" * 72)
        print("2. ROCHE: effective gravity ratio pole/equator")
        print(f"{'q':>10} {'a/c':>7} | {'g_pole/g_eq':>12}")
        for q in (1e-6, 0.05, 0.1, 0.2, 0.28, 0.2960, 8 / 27 - 1e-9):
            ge, gp = roche_gravity(q)
            x = optimize.brentq(lambda z: q * z ** 3 / 2 - z + 1, 1.0, 1.5, xtol=1e-15)
            print(f"{q:10.6f} {x:7.4f} | {gp/ge:12.4f}")
        print("   -> the point-mass limit reaches ANY gravity ratio, because g_eq -> 0")
        print("      at break-up -- but only at a/c = 3/2, where the figure ceases to exist.")
        print("      So the two limiting structures disagree about the ceiling by 15x in a/c")
        print("      and by infinity in gravity ratio.  Nothing here is geometry.")

        print()
        print("=" * 72)
        print("3. The paper's stated numbers")
        print("   'equatorial surface gravity as low as ~4 m/s^2' at ~2.7 L_EM")
        print("   'equatorial gravity of only ~8 m/s^2' in the canonical case")
        print("   polar gravity 'somewhat increased' over the non-rotating 9.8 m/s^2")
        for gp in (10.0, 11.0, 12.0):
            print(f"     with g_pole = {gp:4.1f}:  ratio at g_eq=4 is {gp/4:.2f},"
                  f"  at g_eq=8 is {gp/8:.2f}")
    print(f"\ndull.watch: {len(w.findings)} FP events")
    for ev in w.findings[:6]:
        print("   ", ev)
