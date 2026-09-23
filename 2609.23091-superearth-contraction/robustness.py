"""Try to break Finding 3 before reporting it.

contraction.py finds that at a FIXED, pressure-independent mantle density
deficit, the fractional radius contraction RISES with planet mass -- roughly
+30 % from 1 to 10 Mearth.  The paper reports it FALLING (11.3 -> 9.5 %) and
explains that as "a more strongly compressed massive interior returns a smaller
fractional radius change for the same melt-to-solid density contrast".

Same controlled experiment, opposite sign.  So before saying so, try to make my
sign flip with every knob I am unsure about:

  A  core geotherm -- mine is anchored on Earth (4000 -> 5500 K) and is surely
     too cold for a 10 Mearth core.  Make it much hotter.
  B  mantle geotherm exponent -- the adiabat shape.
  C  core-mass fraction -- 0.2 and 0.5 as well as Earth's 0.325.
  D  mantle stiffness K0' -- 4.1 is bridgmanite; try 3.5 and 5.0 for a mantle
     that stiffens more or less quickly under compression.

If the sign survives all four, it is a property of the structure equations and
not of my parameter choices.  Iris (Opus 5, 1M), Sep 23 2026.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path.home() / "iris-the-maker" / "tools"))
from dull import watch, sane                                    # noqa: E402
import structure as st                                          # noqa: E402
from contraction import DensityDeficit, CORE_SCALE              # noqa: E402

DELTA = 0.10


def make_geotherm(mantle_exp=0.35, core_hot=1500.0, core_base=4000.0):
    def g(m, P, r, layer):
        if layer == "mantle":
            return 300.0 + 2200.0 * (max(P, 0.0) / 136e9) ** mantle_exp
        frac = float(np.clip((P - 136e9) / (364e9 - 136e9), 0.0, 1.0))
        return core_base + core_hot * frac
    return g


def trend(label, cmf=0.325, Kp=4.1, **gt):
    geo = make_geotherm(**gt)
    out = []
    for M in (1.0, 10.0):
        mantle = st.bridgmanite()
        mantle.Kp0 = Kp
        b = st.solve(M * st.M_EARTH, cmf, st.fe_eps(), mantle, Tfun=geo,
                     n_out=220, scales=(CORE_SCALE, 1.0))
        d = st.solve(M * st.M_EARTH, cmf, st.fe_eps(),
                     DensityDeficit(mantle, DELTA), Tfun=geo,
                     n_out=220, scales=(CORE_SCALE, 1.0))
        out.append((1.0 - b["R"] / d["R"], b["R"], b["R_cmb"]))
    c1, c10 = out[0][0], out[1][0]
    ratio = c10 / c1
    sign = "RISES" if ratio > 1 else "falls"
    print(f"  {label:<34} 1 Me {c1*100:5.2f} %   10 Me {c10*100:5.2f} %   "
          f"ratio {ratio:5.3f}  {sign}")
    return ratio


if __name__ == "__main__":
    print(f"Contraction at a fixed mantle density deficit delta = {DELTA*100:.0f} %")
    print("  the paper's trend, for reference:  11.3 % -> 9.5 %   ratio 0.841  falls\n")
    with watch() as w:
        ratios = [
            trend("baseline (as in contraction.py)"),
            trend("core 3x hotter gradient", core_hot=4500.0),
            trend("core base 6000 K, +4500 K", core_base=6000.0, core_hot=4500.0),
            trend("mantle adiabat exponent 0.25", mantle_exp=0.25),
            trend("mantle adiabat exponent 0.45", mantle_exp=0.45),
            trend("CMF 0.20", cmf=0.20),
            trend("CMF 0.50", cmf=0.50),
            trend("mantle K0' = 3.5", Kp=3.5),
            trend("mantle K0' = 5.0", Kp=5.0),
        ]
        r = np.array(ratios)
        sane(r, "mass-trend ratios")
    print()
    w.report()
    print(f"\n  {int((r > 1).sum())}/{len(r)} variants give a RISING trend; "
          f"range {r.min():.3f} to {r.max():.3f}.")
    print(f"  the paper's 0.841 is outside that range by "
          f"{(r.min()-0.841)/0.841*100:.0f} % at the closest variant.")
