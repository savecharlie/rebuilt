"""What mean mantle density change does an 11 % radius contraction REQUIRE?

No equation of state anywhere in this file.  Mantle mass is conserved, so the
mantle's mean density is just its mass over its shell volume, and the ratio
between the molten and solidified states follows from three numbers the paper
states in its own Results section:

    R_int  : 7150 km (molten, Phi = 0.99)  ->  6340 km (solid, Phi = 0)
    core   : radius decreases by about 3 % over the same interval
    mantle : shell thins by roughly 20 %

    rho_solid / rho_molten = V_molten / V_solid
                           = (R_m^3 - Rc_m^3) / (R_s^3 - Rc_s^3)

The answer is the quantity the paper's headline is equivalent to, and it is
not stated anywhere in the paper.  Iris (Opus 5, 1M), Sep 23 2026.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path.home() / "iris-the-maker" / "tools"))
from dull import watch, sane, scaling                           # noqa: E402

R_M, R_S = 7150.0, 6340.0          # km, paper Sec. 3.1
CORE_SHRINK = 0.03                 # paper: "core radius decreases by only ~3 %"


def shell_volume(R, Rc):
    return 4 * np.pi / 3 * (R ** 3 - Rc ** 3)


def density_ratio(Rc_solid, R_m=R_M, R_s=R_S, core_shrink=CORE_SHRINK):
    Rc_molten = Rc_solid / (1.0 - core_shrink)
    return shell_volume(R_m, Rc_molten) / shell_volume(R_s, Rc_solid)


def shell_thinning(Rc_solid, R_m=R_M, R_s=R_S, core_shrink=CORE_SHRINK):
    Rc_molten = Rc_solid / (1.0 - core_shrink)
    h_m, h_s = R_m - Rc_molten, R_s - Rc_solid
    return (h_m - h_s) / h_m


if __name__ == "__main__":
    with watch() as w:
        print("Implied mean mantle densification, from the paper's own numbers")
        print(f"  R_int  {R_M:.0f} km (molten)  ->  {R_S:.0f} km (solid)"
              f"   =  {(R_M - R_S) / R_M * 100:.1f} % contraction")
        print(f"  core radius decreases by {CORE_SHRINK*100:.0f} %\n")
        print(f"  {'Rc_solid [km]':>14} {'Rc_molten':>10} "
              f"{'shell thins':>12} {'rho_s/rho_m':>12}")
        Rcs = np.arange(2800.0, 4001.0, 200.0)
        ratios = []
        for Rc in Rcs:
            r = density_ratio(Rc)
            t = shell_thinning(Rc)
            ratios.append(r)
            flag = "  <- matches their '~20 %' thinning" if 0.185 <= t <= 0.215 else ""
            print(f"  {Rc:14.0f} {Rc/(1-CORE_SHRINK):10.0f} "
                  f"{t*100:11.1f} % {r:12.3f}{flag}")
        ratios = np.array(ratios)
        sane(ratios, "rho_solid/rho_molten over core radius")

        print(f"\n  Over a core radius from {Rcs[0]:.0f} to {Rcs[-1]:.0f} km the "
              f"implied densification\n  ranges only from {ratios.min():.3f} to "
              f"{ratios.max():.3f} -- the conclusion does not depend on the\n"
              f"  core radius, which is why it can be drawn without an EoS.")

        # Earth-like core: the value their '~20 % shell thinning' selects
        Rc_earth = brentq(lambda Rc: shell_thinning(Rc) - 0.20, 2500.0, 4500.0)
        r_earth = density_ratio(Rc_earth)
        print(f"\n  Their stated ~20 % shell thinning is consistent with a solid "
              f"core radius of\n  {Rc_earth:.0f} km (Earth's is 3480 km), and at "
              f"that core radius the mantle's mean\n  density must rise by a "
              f"factor of {r_earth:.3f}  -- i.e. by {100*(r_earth-1):.0f} %.")

        # what the Bower 2019 comparison amounts to
        print("\n  The same arithmetic applied to the comparison the paper draws:")
        for label, contraction in (("this paper, 1 Me", 0.113),
                                   ("this paper, 10 Me", 0.095),
                                   ("Bower et al. 2019", 0.05)):
            Rs = 6340.0
            Rm = Rs / (1.0 - contraction)
            r = density_ratio(3480.0, R_m=Rm, R_s=Rs)
            print(f"    {label:<20} {contraction*100:5.1f} % contraction"
                  f"  ->  mean mantle densification {100*(r-1):5.1f} %")

        # is the contraction linear in the densification?  measure, do not assert.
        # NOTE: the thin-shell formula below assumes a FIXED core, so this scan
        # must fix the core too.  Comparing a shrinking-core scan against a
        # fixed-core formula is a definition mismatch, and it read as a factor
        # of two before I noticed -- dull.py cannot see that kind of fault.
        print("\n  Sensitivity, measured rather than asserted "
              "(core held FIXED at Rc = 3480 km):")
        dens = np.linspace(0.02, 0.60, 30)
        contr = np.array([
            1.0 - R_S / brentq(
                lambda Rm: density_ratio(3480.0, R_m=Rm, R_s=R_S,
                                         core_shrink=0.0) - (1.0 + d),
                R_S, 3.0 * R_S)
            for d in dens])
        sane(contr, "contraction vs densification")
        p_fit = scaling(dens, contr, "contraction(densification)")
        thin = (1 - (3480.0 / R_S) ** 3) / 3
        print(f"    small-strain slope d(contraction)/d(densification) = "
              f"{contr[0]/dens[0]:.4f}")
        print(f"    thin-shell prediction (1 - (Rc/R)^3)/3            = {thin:.4f}")
        print(f"    agreement {100*abs(contr[0]/dens[0] - thin)/thin:.2f} % "
              f"-- two roads, one door")
        print(f"    over 2-60 % the relation bends: "
              f"contraction ~ densification^{p_fit:.3f}")
        print(f"\n    so, at fixed core, a 50 % densification gives "
              f"{100*(1 - R_S/brentq(lambda Rm: density_ratio(3480.0, R_m=Rm, R_s=R_S, core_shrink=0.0) - 1.5, R_S, 3*R_S)):.1f} % "
              f"contraction;\n    the paper's extra ~1 point comes from its "
              f"core shrinking 3 % as well.")

    print()
    w.report()
