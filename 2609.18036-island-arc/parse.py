"""Parse the GCMT ndk catalogue. 5 lines per event.
Validated below against two events whose GCMT solutions are published:
2004-12-26 Sumatra-Andaman and 2011-03-11 Tohoku."""
import numpy as np, datetime

def read_ndk(path):
    L = open(path).read().splitlines()
    ev = []
    for i in range(0, len(L) - 4, 5):
        l1, l2, l3, l4, l5 = L[i:i+5]
        try:
            date = l1[5:15]; lat = float(l1[27:33]); lon = float(l1[34:41])
            dep_h = float(l1[42:47])
            # line 3: centroid depth is the one to use
            p3 = l3.split()
            cdep = float(p3[7]); dtype = p3[9] if len(p3) > 9 else "?"
            p4 = l4.split()
            expo = int(p4[0])
            Mrr = float(p4[1]) * 10.0**expo          # dyne-cm
            p5 = l5.split()
            # scalar moment: exponent from line 4, value is p5[10] (after 3x(eig,plunge,az))
            M0 = float(p5[10]) * 10.0**expo          # dyne-cm
            ev.append((date, lat, lon, cdep, Mrr, M0))
        except Exception:
            continue
    dt = np.dtype([("date","U10"),("lat","f8"),("lon","f8"),("dep","f8"),
                   ("Mrr","f8"),("M0","f8")])
    return np.array(ev, dtype=dt)

if __name__ == "__main__":
    e = read_ndk("jan76_dec20.ndk")
    print(f"events parsed: {len(e)}   {e['date'][0]} .. {e['date'][-1]}")
    DYNECM_TO_NM = 1e-7
    for want, label, pub in [("2004/12/26","Sumatra-Andaman", 3.95e22),
                             ("2011/03/11","Tohoku",          5.31e22)]:
        sel = e[(e["date"] == want)]
        big = sel[np.argmax(sel["M0"])]
        m0 = big["M0"] * DYNECM_TO_NM
        Mw = (np.log10(big["M0"]) - 16.1) / 1.5
        print(f"  {label:16s} M0 = {m0:.3e} N m   Mw = {Mw:.2f}   "
              f"(published GCMT {pub:.2e})  ratio {m0/pub:.3f}")
