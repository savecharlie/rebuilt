"""The three traditional AMOC fingerprints of Emirzade et al. 2026 (2609.18551),
computed from NOAA ERSSTv5 -- an observational dataset independent of the
HadISST they used.  Definitions verbatim from their Methods:

  SST_SG-G  = mean SST over 46-61N, 55-20W, Nov-Mar, minus the GLOBAL mean SST
              over the same time steps.
  SST_SG-NH = same box minus the NORTHERN HEMISPHERE mean.
  SST_DP    = annual mean SST over 45-80N, 70W-30E, minus annual mean SST over
              45S-0 in the same longitude range.
"""
import numpy as np, ersst

WINTER = [11, 12, 1, 2, 3]
SG_BOX = dict(s=46, n=61, w=-55, e=-20)
DP_N   = dict(s=45, n=80, w=-70, e=30)
DP_S   = dict(s=-45, n=0,  w=-70, e=30)

def build(e=None):
    e = e or ersst.ERSST()
    out = {}
    yr, sg   = e.season(e.box(**SG_BOX),        WINTER)
    _,  glob = e.season(e.box(-90, 90, 0, 360), WINTER)
    _,  nh   = e.season(e.box(0, 90, 0, 360),   WINTER)
    ya, dn   = e.annual(e.box(**DP_N))
    _,  ds   = e.annual(e.box(**DP_S))
    out['years']  = yr
    out['SG']     = sg
    out['GLOBAL'] = glob
    out['NH']     = nh
    out['SG-G']   = sg - glob
    out['SG-NH']  = sg - nh
    out['years_a'] = ya
    out['DP_N'], out['DP_S'] = dn, ds
    out['DP']      = dn - ds
    return e, out

if __name__ == '__main__':
    e, f = build()
    from ersst import trend
    print(f"years {f['years'][0]}-{f['years'][-1]}  (annual {f['years_a'][0]}-{f['years_a'][-1]})")
    print(f"\n{'series':10s} {'period':12s} {'trend degC/century':>20s} {'+/- 1sigma':>11s}")
    for period in [(1870,2025),(1900,2025),(1950,2025),(1980,2025)]:
        for k in ['SG','GLOBAL','NH','SG-G','SG-NH']:
            y = f['years']; m = (y>=period[0])&(y<=period[1])
            b,se = trend(y[m], f[k][m])
            print(f"{k:10s} {str(period):12s} {b:+20.4f} {se:11.4f}")
        y = f['years_a']; m=(y>=period[0])&(y<=period[1])
        for k in ['DP_N','DP_S','DP']:
            b,se = trend(y[m], f[k][m])
            print(f"{k:10s} {str(period):12s} {b:+20.4f} {se:11.4f}")
        print()
