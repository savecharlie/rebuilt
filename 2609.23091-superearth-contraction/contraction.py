"""What melt-solid density contrast does an 11 % radius contraction need?

geometry.py showed, with no equation of state at all, that the paper's quoted
11.3 % contraction at 1 Mearth is equivalent to the mantle's MEAN density
rising by about 50 %.  The published isochemical melt-solid contrast is
nothing like that large: it falls from roughly 20 % in the shallow mantle to
about 4 % at core-mantle-boundary pressure (Karki et al. 2018 GRL; Petitgirard
et al. 2015 PNAS put MgSiO3 glass 1.6 % below bridgmanite at 133 GPa).

The two are not the same quantity.  A less dense mantle inflates the planet,
which lowers the interior pressure, which decompresses the mantle further.
The mean-density change is the LOCAL contrast times that feedback.  This file
measures the feedback instead of arguing about it.

The experiment: take the validated solid planet, make the mantle a uniform
fraction delta less dense AT THE SAME PRESSURE AND TEMPERATURE -- which is
exactly the quantity the melt literature reports -- and re-solve the structure
self-consistently.

    contraction(delta) = 1 - R(0) / R(delta)

Iris (Opus 5, 1M), Sep 23 2026.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path.home() / "iris-the-maker" / "tools"))
from dull import watch, sane, scaling                           # noqa: E402
import structure as st                                          # noqa: E402

CORE_SCALE = 0.90          # 10 wt% light elements; validated against PREM
CMF = 0.325


class DensityDeficit:
    """rho(P,T) = (1 - delta) * rho_base(P,T).

    Deliberately crude and deliberately explicit: it is the published
    melt-solid contrast imposed at every pressure, not a thermodynamically
    self-consistent melt EoS.  It answers 'if the melt is delta less dense
    than the solid everywhere, how much bigger is the planet', which is the
    question a radius contraction actually poses.
    """

    def __init__(self, base, delta):
        self.base, self.delta = base, delta
        self.name = f"{base.name} minus {delta*100:.1f}%"

    def rho(self, P, T, scale=1.0):
        return (1.0 - self.delta) * self.base.rho(P, T, scale)


def geotherm_solid(m, P, r, layer):
    if layer == "mantle":
        return 300.0 + 2200.0 * (max(P, 0.0) / 136e9) ** 0.35
    return 4000.0 + 1500.0 * float(np.clip((P - 136e9) / (364e9 - 136e9), 0, 1))


def solve_planet(M, delta, n_out=260):
    mantle = st.bridgmanite() if delta == 0.0 else \
        DensityDeficit(st.bridgmanite(), delta)
    out = st.solve(M, CMF, st.fe_eps(), mantle, Tfun=geotherm_solid,
                   n_out=n_out, scales=(CORE_SCALE, 1.0))
    return out


if __name__ == "__main__":
    masses = [1.0, 3.0, 5.0, 10.0]
    deltas = np.array([0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])

    with watch() as w:
        print("Radius contraction produced by a uniform mantle density deficit")
        print(f"(core = eps-Fe x {CORE_SCALE} [PREM-validated], CMF = {CMF})\n")
        base = {}
        for M in masses:
            base[M] = solve_planet(M * st.M_EARTH, 0.0)
        print("  solid reference planets:")
        for M in masses:
            b = base[M]
            print(f"    {M:4.0f} Me   R = {b['R']/st.R_EARTH:6.4f} Re "
                  f"({b['R']/1e3:7.1f} km)   R_cmb = {b['R_cmb']/1e3:6.1f} km"
                  f"   P_c = {b['P_c']/1e9:7.1f} GPa")

        print(f"\n  {'delta':>7} | " +
              " | ".join(f"{M:>4.0f} Me contraction" for M in masses))
        table = {M: [] for M in masses}
        meanrho = {M: [] for M in masses}
        for d in deltas:
            row = []
            for M in masses:
                o = solve_planet(M * st.M_EARTH, float(d))
                c = 1.0 - base[M]["R"] / o["R"]
                table[M].append(c)
                meanrho[M].append(base[M]["rho_mean_mantle"] /
                                  o["rho_mean_mantle"])
                row.append(f"{c*100:16.2f} %")
            print(f"  {d*100:6.1f} % | " + " | ".join(row))

        for M in masses:
            sane(np.array(table[M]), f"contraction({M} Me)")

        # -- the inversion, by interpolation then ONE confirming solve ------
        # Inverting with brentq would cost ~15 full structure solves per
        # target.  The sweep above already contains the curve, so interpolate
        # it and then spend a single solve checking the answer -- which also
        # makes the confirming solve a second instrument on the interpolation.
        print("\n  Inverting for the paper's quoted contraction:")
        targets = {1.0: 0.113, 10.0: 0.095}
        found = {}
        for M, tgt in targets.items():
            cs_ = np.array(table[M])
            spl = CubicSpline(cs_, deltas)          # contraction -> delta
            d_star = float(spl(tgt))
            o = solve_planet(M * st.M_EARTH, d_star)
            c_check = 1.0 - base[M]["R"] / o["R"]
            # rho_mean_mantle of the SOLID over the MOLTEN, i.e. how much
            # denser the mantle gets.  Written the other way up on the first
            # run and it printed a negative densification -- a formula error,
            # invisible to dull.py by construction, caught only by the sign
            # looking wrong and by geometry.py disagreeing.
            amp_ratio = base[M]["rho_mean_mantle"] / o["rho_mean_mantle"]
            dens = amp_ratio - 1.0
            found[M] = d_star
            print(f"    {M:4.0f} Me  target {tgt*100:5.1f} %  ->  "
                  f"delta = {d_star*100:5.2f} %")
            print(f"              confirming solve gives {c_check*100:5.2f} % "
                  f"(interpolation error {abs(c_check-tgt)*100:.3f} points)")
            print(f"              mean mantle densification {dens*100:5.1f} %"
                  f"   amplification {dens/d_star:4.2f}x"
                  f"   R_molten {o['R']/1e3:.0f} km")

        # -- the mass trend, their claim, tested independently ---------------
        print("\n  Their claim: contraction falls from ~11 % at 1 Me to "
              "~9.5 % at 10 Me.")
        d_fix = found[1.0]
        cs2 = []
        for M in masses:
            spl = CubicSpline(deltas, np.array(table[M]))
            cs2.append(float(spl(d_fix)))
        print(f"    at a FIXED delta = {d_fix*100:.2f} % (the value that "
              f"reproduces their 1 Me number):")
        for M, c in zip(masses, cs2):
            print(f"      {M:4.0f} Me  contraction {c*100:6.2f} %")
        cs2 = np.array(cs2)
        sane(cs2, "contraction vs mass")
        pfit = scaling(np.array(masses), cs2, "contraction(mass)")
        print(f"    measured mass exponent {pfit:+.4f}")
        print(f"    paper 11.3 -> 9.5 % over 1 -> 10 Me is exponent "
              f"{np.log(9.5/11.3)/np.log(10):+.4f}")

        np.savez("sweep.npz", deltas=deltas, masses=np.array(masses),
                 **{f"c{M:g}": np.array(table[M]) for M in masses},
                 **{f"R{M:g}": base[M]["R"] for M in masses})

    print()
    w.report()
