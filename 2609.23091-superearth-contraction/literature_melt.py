"""What contraction do the PUBLISHED melt-solid contrasts actually produce?

contraction.py inverts the paper's 11.3 % to a required uniform density
deficit.  This file does the forward problem with a pressure-dependent deficit
built from published anchors instead of a constant:

    delta(P) = delta_0 * exp(-P / P*)

anchored on the two numbers the melt literature actually reports:
    delta(0)       ~ 0.20   shallow mantle melt-solid contrast
    delta(136 GPa) ~ 0.04   Karki et al. 2018, GRL 45 (CMB pressure)
                            Petitgirard et al. 2015, PNAS 112 give 0.016 for
                            MgSiO3 glass at 133 GPa, so 0.04 is generous.

This is an interpolation between two published points, not a derived equation
of state, and it is deliberately generous at both ends.  The structure is then
solved with it directly, so no weighting argument is needed -- the pressure
re-adjustment is handled by the solver exactly as in the uniform-delta case.

Thermal contraction of the superheated melt is NOT included here; the paper
accounts for it separately at 1-2 percentage points of radius, and that number
is used as theirs rather than re-derived.

Iris (Opus 5, 1M), Sep 23 2026.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path.home() / "iris-the-maker" / "tools"))
from dull import watch, sane                                    # noqa: E402
import structure as st                                          # noqa: E402
from contraction import CORE_SCALE, CMF, geotherm_solid, solve_planet  # noqa


class PressureDependentDeficit:
    def __init__(self, base, delta0, p_star):
        self.base, self.delta0, self.p_star = base, delta0, p_star
        self.name = f"{base.name} melt, delta0={delta0}, P*={p_star/1e9:.1f} GPa"

    def delta(self, P):
        return self.delta0 * np.exp(-max(P, 0.0) / self.p_star)

    def rho(self, P, T, scale=1.0):
        return (1.0 - self.delta(P)) * self.base.rho(P, T, scale)


def p_star_from(delta0, delta_cmb, P_cmb=136e9):
    return P_cmb / np.log(delta0 / delta_cmb)


if __name__ == "__main__":
    M = st.M_EARTH
    with watch() as w:
        base = solve_planet(M, 0.0)
        print(f"solid 1 Me reference: R = {base['R']/1e3:.1f} km, "
              f"P_cmb = {base['P_cmb']/1e9:.1f} GPa\n")
        print("Forward: contraction from a published, pressure-dependent contrast")
        print(f"  {'delta(0)':>9} {'delta(136GPa)':>14} {'P* [GPa]':>9} "
              f"{'R_molten [km]':>14} {'contraction':>12}")
        rows = []
        for d0, dcmb in ((0.20, 0.040), (0.20, 0.016), (0.25, 0.040),
                         (0.30, 0.060), (0.15, 0.040), (0.35, 0.080)):
            ps = p_star_from(d0, dcmb)
            eos = PressureDependentDeficit(st.bridgmanite(), d0, ps)
            o = st.solve(M, CMF, st.fe_eps(), eos, Tfun=geotherm_solid,
                         n_out=220, scales=(CORE_SCALE, 1.0))
            c = 1.0 - base["R"] / o["R"]
            rows.append(c)
            print(f"  {d0:9.2f} {dcmb:14.3f} {ps/1e9:9.1f} "
                  f"{o['R']/1e3:14.1f} {c*100:11.2f} %")
        rows = np.array(rows)
        sane(rows, "contraction from literature contrasts")
        print(f"\n  range {rows.min()*100:.2f} to {rows.max()*100:.2f} %"
              f"   against the paper's 11.3 %")
        print(f"  the most generous case here reaches "
              f"{rows.max()/0.113*100:.0f} % of their number,")
        print(f"  and the paper attributes a further 1-2 points of its 11.3 %"
              f" to thermal\n  contraction of the superheated melt, which is "
              f"not included above.")
    print()
    w.report()
