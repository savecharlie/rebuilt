# arXiv:2609.18551 — AMOC fingerprints, rebuilt from two SST datasets

Emirzade, Ajagun-Brauns, Ben-Yami, Bathiany, Shin, Kug & Boers, *Disentangled Fingerprints
suggest no historical weakening of Atlantic Overturning and Subpolar Gyre* (16 Sep 2026).

The paper says the observational AMOC "fingerprints" in common use — an SST box south of
Greenland minus a reference field — correlate only weakly with the true overturning inside CMIP6
models, and that a ridge-regression fingerprint trained on piControl does better and shows no
recent weakening.

I can't touch the model half without the CMIP6 fields. I rebuilt the observational half, on
**HadISST (theirs) and NOAA ERSSTv5 (independent)**.

## What came out

- **Their observational result is dataset-independent.** Every sign and nearly every magnitude
  survives the swap from HadISST to ERSSTv5. Not shown in the paper; worth a sentence.
- **The subpolar box's own trend is zero after 1900** (+0.011 ± 0.094 °C/century, HadISST), so
  the index's −0.611 is carried entirely by the global mean it subtracts.
- **…and the control says that is still remarkable.** Sliding the same-sized box over 868
  fully-covered ocean regions gives a median index trend of −0.048; only 7 of 868 reach −0.546.
  The warming hole is a bottom-1% outlier and the traditional index is not measuring nothing.
  **This killed my first hypothesis.**
- **The three traditional fingerprints agree with each other**, r = 0.89–0.99 detrended at
  30-year smoothing, in both datasets. **This killed my second hypothesis**, and sharpens the
  paper: they are not three noisy shots at the AMOC, they are one consistent measurement of
  something else.
- **The box is a near-cancellation:** 63% of its area cools at −0.255 °C/century and 37% warms
  at +0.266, leaving a mean of −0.064. The index's *level* is a box-placement decision.

Full write-up, including the two bugs I found in my own code: [`notes.md`](notes.md).
The letter sent to the corresponding author: [`reported.md`](reported.md).
The readable version: [`../essays/2609.18551-amoc-fingerprints.md`](../essays/2609.18551-amoc-fingerprints.md).

## Instrument validation, run before anything new was computed

| check | against | result |
|---|---|---|
| Niño3.4, 918 months, ERSSTv5 | CPC `ersst5.nino.mth.91-20.ascii` | rms **0.0090 °C**, r = 0.99996 |
| global ocean annual anomaly, 146 yr | NOAA NCEI published series | rms **0.024 °C**, r = 0.997; trend 0.600 vs 0.612 °C/century |

## Scripts

| file | what it makes |
|---|---|
| `ersst.py` | the area-weighted box-mean reader, seasonal collapse, OLS trend + σ |
| `hadisst.py` | same interface, HadISST loader (masks the −1000 sea-ice sentinel) |
| `fingerprints.py` | the three indices and their term-by-term decomposition, both datasets |
| `warming_hole.py` | per-cell North Atlantic Nov–Mar trend, coldest cells, box composition |
| `null_boxes.py` | **the control** — the same index for 868 ocean boxes (~3 min) |
| `agreement.py` | do the three fingerprints agree with each other |
| `figure.py` | `warming_hole_decomposition.png` |

## Data

Not committed — 600 MB. Fetch into `data/`:

```sh
mkdir -p data && cd data
curl -LO https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
curl -L https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz | gunzip > HadISST_sst.nc
curl -o cpc_nino.txt https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii
curl -o noaa_ocean.csv "https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/global/time-series/globe/ocean/12/12/1880-2026/data.csv"
```

Then `python3 fingerprints.py`, etc. The scripts write `.npy` intermediates beside themselves.
