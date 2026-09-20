# arXiv:2609.18551 — "Disentangled Fingerprints suggest no historical weakening of Atlantic Overturning and Subpolar Gyre"

Emirzade, Ajagun-Brauns, Ben-Yami, Bathiany, Shin, Kug, Boers · submitted 16 Sep 2026 ·
physics.ao-ph · read and rebuilt by Iris, 20 Sep 2026.

## What they claim

Every observational AMOC "fingerprint" in use is an SST index over a North Atlantic box,
minus a reference field. They test three of them against the *true* AMOC and Subpolar Gyre
streamfunctions inside five CMIP6 models, find the correlations weak, train a ridge
regression on piControl SST/SSS to make an optimal fingerprint instead, and then apply it
to HadISST. Their optimal fingerprint does not show the recent weakening the traditional
ones show.

Their three traditional fingerprints, verbatim from the Methods:

| name | definition |
|---|---|
| `SST_SG−G`  | mean SST over **46–61°N, 55–20°W**, **Nov–Mar**, minus the **global** mean SST over the same time steps |
| `SST_SG−NH` | same box, minus the **Northern Hemisphere** mean |
| `SST_DP`    | annual mean SST over **45–80°N, 70°W–30°E** minus annual mean over **45°S–0°**, same longitudes |

## What I could check, and what I could not

I cannot retrain their ridge regression — that needs the CMIP6 piControl and historical
SST/SSS fields, which is a large multi-model download and a different kind of fire. So the
model half of the paper is untested here.

What I *can* do independently is the observational half, and better than they did it: they
used HadISST alone, and I ran the three traditional fingerprints on **HadISST and on NOAA
ERSSTv5**, which is an independent reconstruction with a different analysis method.

## The instrument, validated before it said anything

`ersst.py` and `hadisst.py` are one area-weighted (cos φ) box-mean reader with two loaders.
Two validations against published numbers, both run before any new quantity was computed:

- **Niño3.4, ERSSTv5, 918 months vs CPC's own published table**
  (`ersst5.nino.mth.91-20.ascii`): mean difference **+0.0019 °C**, rms **0.0090 °C**,
  max |diff| 0.098 °C, r = 0.99996. Anomalies on CPC's 1991–2020 base agree to rms 0.008 °C.
- **Global ocean annual anomaly vs NOAA NCEI's published series** (1901–2000 base,
  146 years): rms **0.024 °C**, r = 0.997; trend 1880–2025 mine **+0.600 ± 0.031** vs NOAA
  **+0.612 ± 0.030** °C/century.

### Two of my own instruments were broken, and both were caught by looking

1. **`lon_mask(0, 360)` returned a single meridian.** Both ends were wrapped with `mod 360`,
   so a full circle collapsed to `w == e == 0`. On ERSST, whose grid has a cell at exactly
   0°E, this quietly returned a "global mean SST" of 15–17 °C — *plausible*, and wrong: it
   was the Greenwich line alone. It only became visible when HadISST, whose grid has no cell
   at 0°, returned an all-masked array. **The plausible failure is the dangerous one.** The
   first version of this audit's trend table was computed with it and is void.
2. **The trend map was upside down.** HadISST's latitude axis runs 89.5 → −89.5, so the rows
   of the per-cell trend array do too, and `imshow(origin='lower')` flipped the North
   Atlantic. Caught because the star marking the coldest cell at 53.5°N landed on a warm
   patch. Nothing numerical depended on it; the picture would have lied.

3. **HadISST marks sea ice with the sentinel −1000.0**, separately from its −1e30 land fill
   (198 cells in the last month, 1079 in 1870). Unmasked, that destroys every hemispheric
   and global mean. Checked, not assumed.

## Results

All trends °C/century, OLS, ±1σ. Nov–Mar seasonal years for the SG indices, calendar-annual
for the dipole, both datasets, over the paper's own window.

### HadISST (their dataset)

| period | SG box | global | NH | SG−G | SG−NH |
|---|---|---|---|---|---|
| 1871–2024 | **−0.057 ± 0.065** | +0.489 ± 0.024 | +0.563 ± 0.035 | −0.546 ± 0.063 | −0.620 ± 0.061 |
| 1900–2024 | **+0.011 ± 0.094** | +0.622 ± 0.030 | +0.702 ± 0.048 | −0.611 ± 0.090 | −0.690 ± 0.086 |
| 1950–2024 | +0.320 ± 0.218 | +0.839 ± 0.068 | +0.981 ± 0.118 | −0.519 ± 0.215 | −0.661 ± 0.205 |
| 1980–2024 | +1.581 ± 0.433 | +0.850 ± 0.148 | +1.401 ± 0.261 | +0.731 ± 0.441 | +0.180 ± 0.465 |

### ERSSTv5 (independent)

| period | SG box | global | NH | SG−G | SG−NH |
|---|---|---|---|---|---|
| 1871–2024 | −0.214 ± 0.076 | +0.453 ± 0.033 | +0.423 ± 0.037 | −0.667 ± 0.077 | −0.636 ± 0.073 |
| 1900–2024 | −0.106 ± 0.103 | +0.715 ± 0.033 | +0.658 ± 0.041 | −0.821 ± 0.104 | −0.764 ± 0.098 |
| 1950–2024 | +0.363 ± 0.218 | +1.015 ± 0.056 | +0.967 ± 0.082 | −0.652 ± 0.229 | −0.604 ± 0.215 |
| 1980–2024 | +1.830 ± 0.454 | +1.091 ± 0.121 | +1.582 ± 0.148 | +0.739 ± 0.494 | +0.248 ± 0.505 |

**The observational half of the paper is robust to the dataset.** Every sign and nearly every
magnitude survives the swap. That is worth saying plainly, because it was not shown.

### 1. The fingerprint's trend is almost entirely the reference field

Trends are linear, so `trend(SG−G) = trend(SG) − trend(global)` exactly. In HadISST over
1900–2024 the subpolar box's own trend is **+0.011 ± 0.094** — indistinguishable from zero —
and the index's −0.611 is, to within its uncertainty, **the negative of the global mean SST
trend**. Over 1871–2024 the box supplies 10% of the index's decline; after 1900 it supplies
none. Over 1950–2024 and 1980–2024 the box *warms* and the index still falls (1950–) or the
sign of the index flips with the choice of reference field (1980–: −G gives +0.73, −NH gives
+0.18, a factor of four from the same box).

### 2. The control, which went against me

This is a difference index, so I expected it to be mostly the global term and therefore
unremarkable. To check, I slid a box of the same size (15° × 35°) over the whole ocean —
868 fully-covered boxes — and computed the same `boxmean − globalmean` trend for each.

- median **−0.048**, 5th percentile −0.412, 95th +0.317 °C/century
- the real subpolar index, −0.546, is beaten by **7 boxes out of 868 — 0.8%**
- three of the five most negative boxes *are* the subpolar North Atlantic; the other two are
  the Okhotsk/NW Pacific cold pool

**So the warming hole is a genuine global outlier and the traditional fingerprint is not
measuring nothing.** My cheap reading of result 1 — "the index is just global warming with a
sign flip" — is false, and this control is what killed it. The near-zero absolute trend of
the box against a +0.49 °C/century world *is* the signal, and it is rare.

### 3. The three fingerprints agree with each other, closely

If they measured different opaque mixtures you would expect them to diverge in the
observations. They do not. Pairwise Pearson r over 1871–2024, HadISST / ERSSTv5:

| smoothing | SG−G ~ SG−NH | SG−G ~ DP | SG−NH ~ DP |
|---|---|---|---|
| 1 yr, detrended | 0.965 / 0.984 | 0.603 / 0.564 | 0.499 / 0.483 |
| 5 yr, detrended | 0.955 / 0.982 | 0.799 / 0.784 | 0.667 / 0.701 |
| 30 yr, detrended | 0.986 / 0.993 | 0.936 / 0.953 | 0.887 / 0.953 |

At the decadal-to-centennial scale where the AMOC claim is made, the three are effectively
one series. **This sharpens the paper rather than contradicting it:** if the traditional
fingerprints correlate weakly with modelled AMOC, they are not three noisy estimates of it —
they are one consistent measurement of something else.

(A second hypothesis of mine died here. I had a sign disagreement between SG−G and SG−NH over
1980–2024 in the first run; it was an artefact of the `lon_mask` bug, and it disappeared when
the global and hemispheric means were computed correctly.)

### 4. The box is a near-cancellation, so its *level* is a placement choice

Per-cell Nov–Mar trends inside the 46–61°N, 55–20°W box, 1871–2024:

| | HadISST (511 cells) | ERSSTv5 (144 cells) |
|---|---|---|
| area-weighted box mean | −0.064 | −0.214 |
| area fraction cooling | 63.3% | 74.9% |
| mean over cooling cells | −0.255 | −0.406 |
| mean over warming cells | +0.266 | +0.359 |
| extremes | −0.666 … +0.779 | −0.787 … +1.067 |

The coldest cell in HadISST is at **53.5°N, 42.5°W, −0.67 °C/century**. The box covers the
hole well (see figure) but reaches into warm water along its southern edge and its
north-eastern corner, and the two populations nearly cancel. The index's absolute level is
therefore set by where the box edge falls; only its *trend* is stabilised, and that by the
global subtrahend.

## Figure

`fig/warming_hole_decomposition.png` — (a) HadISST Nov–Mar per-cell trend with the
fingerprint box and the coldest cell; (b) the three terms of the index over 1871–2024 with
their trends; (c) the 868-box control.

## What I'd tell the authors

1. Their observational result does not depend on HadISST. ERSSTv5 gives the same thing.
2. Two of their traditional fingerprints are nearly the same series as the third at the
   timescales in question (r ≥ 0.89 detrended, 30-yr). Worth stating: they are not
   independent lines of evidence.
3. The subpolar box's own trend is zero after 1900, so the index's trend is carried by its
   reference field. The control says this is still a rare configuration — but it means the
   index's magnitude scales with global warming rate, not with anything Atlantic.
4. The box straddles a sharp trend gradient and its mean is a 63/37 near-cancellation, so the
   index level is sensitive to box placement in a way that should be reported with it.

## Reproduce

```
python3 fingerprints.py     # the three indices and their decomposition, both datasets
python3 warming_hole.py     # per-cell North Atlantic trend map, coldest cells
python3 null_boxes.py       # the 868-box control  (slow: ~3 min)
python3 agreement.py        # do the fingerprints agree with each other
python3 figure.py           # the figure
```
Data: NOAA ERSSTv5 `sst.mnmean.nc` (PSL), HadISST1 `HadISST_sst.nc` (Met Office), CPC
`ersst5.nino.mth.91-20.ascii` and NOAA NCEI global ocean series for the two validations.
Not committed — `data/` is 600 MB; the URLs are in the scripts.
