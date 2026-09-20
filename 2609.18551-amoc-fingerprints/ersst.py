"""Minimal, validated reader for NOAA ERSSTv5 monthly SST (2 deg grid).

Every number in this directory comes through here.  The box-mean is
area-weighted by cos(latitude) and masks land/ice via the file's own
missing-value mask, recomputed per month (the mask moves with sea ice).
"""
import numpy as np
from netCDF4 import Dataset

PATH = 'data/sst.mnmean.nc'

class ERSST:
    def __init__(self, path=PATH):
        d = Dataset(path)
        self.lat = d.variables['lat'][:].astype(float)      # 88 .. -88, 2 deg
        self.lon = d.variables['lon'][:].astype(float)      # 0 .. 358, 2 deg
        sst = d.variables['sst'][:]                          # masked array
        self.sst = np.ma.masked_invalid(sst.astype(float))
        t = d.variables['time'][:]
        # days since 1800-01-01; monthly, so reconstruct year/month exactly
        import datetime as dt
        base = dt.date(1800, 1, 1)
        dates = [base + dt.timedelta(days=int(x)) for x in t]
        self.year = np.array([x.year for x in dates])
        self.month = np.array([x.month for x in dates])
        self.w = np.cos(np.deg2rad(self.lat))                # area weight

    def lon_mask(self, w, e):
        """Longitude mask on the 0-360 grid; w,e in degrees east (may be negative).

        BUG FOUND AND FIXED Sep 20 2026: the first version wrapped BOTH ends with
        mod 360, so a full circle (0,360) collapsed to w==e==0 and selected a
        single meridian.  On ERSST that quietly returned a plausible-looking
        global mean (15-17 degC) computed along the Greenwich line alone; on
        HadISST, whose grid has no cell at exactly 0, it returned an all-masked
        array, which is what exposed it.  A plausible number is the dangerous
        failure.  Full-circle spans are now detected before any wrapping.
        """
        if (e - w) >= 360.0 or np.isclose((e - w) % 360.0, 0.0) and e != w:
            return np.ones(self.lon.shape, bool)
        lo = np.mod(self.lon, 360.0)
        w, e = np.mod(w, 360.0), np.mod(e, 360.0)
        return (lo >= w) & (lo <= e) if w <= e else (lo >= w) | (lo <= e)

    def box(self, s, n, w, e):
        """Area-weighted monthly mean SST over a lat/lon box -> (2072,) array."""
        jm = (self.lat >= s) & (self.lat <= n)
        im = self.lon_mask(w, e)
        sub = self.sst[:, jm][:, :, im]                       # (t, ny, nx)
        wt = np.broadcast_to(self.w[jm][None, :, None], sub.shape)
        wt = np.ma.array(wt, mask=np.ma.getmaskarray(sub))
        return (sub * wt).sum((1, 2)) / wt.sum((1, 2))

    def season(self, series, months):
        """Collapse a monthly series to one value per season-year.

        For a wrapped season (Nov-Mar) the year is that of the January,
        which is the Caesar/Rahmstorf convention.  Returns (years, values);
        incomplete seasons at the ends are dropped.
        """
        wrap = months[0] > months[-1]
        yr = self.year.copy()
        if wrap:
            yr = np.where(np.isin(self.month, [m for m in months if m >= months[0]]),
                          yr + 1, yr)
        sel = np.isin(self.month, months)
        out_y, out_v = [], []
        for y in range(yr[sel].min(), yr[sel].max() + 1):
            m = sel & (yr == y)
            if m.sum() == len(months):
                out_y.append(y); out_v.append(series[m].mean())
        return np.array(out_y), np.array(out_v)

    def annual(self, series):
        return self.season(series, list(range(1, 13)))

def anomaly(years, vals, base=(1981, 2010)):
    m = (years >= base[0]) & (years <= base[1])
    return vals - vals[m].mean()

def trend(years, vals):
    """OLS slope per century and its OLS standard error."""
    g = np.isfinite(vals)
    x, y = years[g].astype(float), np.asarray(vals)[g]
    n = len(x)
    b, a = np.polyfit(x, y, 1)
    resid = y - (b * x + a)
    se = np.sqrt((resid**2).sum() / (n - 2) / ((x - x.mean())**2).sum())
    return b * 100, se * 100
