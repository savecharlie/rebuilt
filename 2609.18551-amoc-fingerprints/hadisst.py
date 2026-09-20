"""HadISST1 reader with the same interface as ersst.ERSST.

TRAP, checked and not assumed: HadISST marks SEA ICE with the sentinel -1000.0
(198 cells in the last month, 1079 in 1870), separately from the land fill of
-1e30.  Leaving those in destroys every hemispheric and global mean.  They are
masked here.  Values in [-1.8, -1.79) are real freezing-point water and stay.
"""
import datetime as dt
import numpy as np
from netCDF4 import Dataset
from ersst import ERSST

class HadISST(ERSST):
    def __init__(self, path='data/HadISST_sst.nc'):
        d = Dataset(path)
        self.lat = d.variables['latitude'][:].astype(float)     # 89.5 .. -89.5
        self.lon = np.mod(d.variables['longitude'][:].astype(float), 360.0)
        sst = np.ma.masked_invalid(d.variables['sst'][:].astype(float))
        sst = np.ma.masked_less(sst, -100.0)                    # the -1000 ice flag
        # reorder longitude to ascending 0..359 so lon_mask works unchanged
        order = np.argsort(self.lon)
        self.lon = self.lon[order]
        self.sst = sst[:, :, order]
        base = dt.date(1870, 1, 1)
        dates = [base + dt.timedelta(days=int(x)) for x in d.variables['time'][:]]
        self.year = np.array([x.year for x in dates])
        self.month = np.array([x.month for x in dates])
        self.w = np.cos(np.deg2rad(self.lat))
