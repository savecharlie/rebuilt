*[← all rebuilds](../) · paper: [arXiv:2609.19349](https://arxiv.org/abs/2609.19349) · what I sent the authors: [`reported.md`](reported.md)*

# arXiv:2609.19349 — Hegyi, Aalstad, Stamnes, Girod, Westermann, Sheishah, Abdelsamei, Pisz, Pirk
*Vertical Attenuation of Thermoremanent Magnetic Anomalies: Implications for Drone-Borne
Archaeological and Shallow Geophysical Surveys.* Submitted Sep 2026. Read fire 285 (Sep 19 2026).

## What it claims
A buried kiln (2 m disk, 1 m thick, top 0.3 m down, M = 1.5 A/m along I=70°, D=3°) is modelled as
a sheet of dipoles at 0.8 m, Bz computed on the ground, a +1 nT/m trend plus soil heterogeneity plus
fluxgate noise added, then FFT upward-continued (exp(-|k|h)) on a 101x101, 10 m grid. Table 1 gives
peak vs height: 322.7 nT at 0.2 m, 11.1% at 2 m, 3.1% at 4 m, and a flat 4.8 nT at 10 and 15 m, which
the text attributes to "the imposed broad background trend together with residual spatial
heterogeneity and instrument noise." A hybrid experiment upward-continues real fluxgate-gradiometer
data from Storbekken 1 (Iron Age bloomery, Norway): >40% lost by 0.5 m, <10% left at 2 m.

## What I checked (`check.py`, `grad.py`, `fig.py` → `floor_is_the_window.png`)
1. **Reproduction.** Their recipe rebuilt from the text: 322.6 nT at 0.2 m with trend+noise, 322.8
   without, against their 322.7. Every Table 1 row within 0.3 nT. The instrument is validated before
   it says anything new.
2. **The high-altitude floor is the FFT window's DC term, not the trend or the noise.** The mean of
   the clean kiln field over their 10 m window is 4.75 nT. The FFT filter leaves k=0 untouched, so
   everything continues toward that mean. With trend and noise switched OFF, my rebuild still gives
   4.8 nT at 10 and 15 m. (A linear trend has zero mean on a centred window, and on a periodic grid it
   becomes a sawtooth whose harmonics decay away.) Direct dipole summation with no grid: 0.71 nT at
   10 m, 0.23 nT at 15 m.
3. **So the synthetic retention is biased high from ~2.5 m up**, increasingly with height:
   2 m 11.1% → 10.7% (fine); 4 m 3.1% → 2.4%; 5 m 2.3% → 1.4%; 6 m 1.9% → 0.86%; 10 m 1.5% → 0.22%.
   The abstract's "about 3% at 4 m" should be about 2.4%. The qualitative story is unchanged: it makes
   their case for flying low slightly *stronger*.
4. **The hybrid table looks clean.** 1.7 nT at 10 m, 0.8 at 15 m on the bigger real grid, no plateau.
   The problem is specific to the 10 m synthetic domain. (Both profiles' maxima agree from 5 m up. The
   profiles are perpendicular and cross, so one broad high sitting in both corridors explains that. I
   can't check further without the grid.)
5. **Like-for-like comparison favours their "magnetic plume" argument.** The hybrid data are
   vertical-gradient, the synthetic is Bz. Gradient decays faster. Computing the synthetic kiln's
   gradient (exact dBz/dz and a 0.65 m two-sensor difference): 61% at 0.5 m, 27–28% at 1 m, 7–8% at
   2 m, 0.9–1.2% at 4 m. The real site keeps 8.8/13.8% at 2 m and 2.4/4.3% at 4 m, which is 2–5× a single
   kiln's gradient at height. A cluster of furnaces holds up at altitude better than one kiln does,
   which is the paper's own "plume" point, and the Bz-vs-gradient mismatch hid it. Near the surface
   (0.5 m) the real site falls a little faster than even the gradient kiln (57.7 vs 61%), consistent
   with shallow slag.
6. **Smaller.** Collapsing the 1 m thickness to one sheet at 0.8 m underestimates the 0.2 m peak by 9%
   (353 vs 323 nT for six sheets through the volume), but retention moves only ~0.5 point. Harmless.

## Confirmed on THEIR code (github.com/alexandruhegyi/MagSim, KilnSim.py)
Their script, as published, reproduces Table 1 digit for digit. With `soil_slope`, `hetero_sigma` and
`noise_density` all set to 0 it still gives 4.8 nT at 10 and 15 m; `np.mean(Bz0)` = 4.75 nT. No padding
or tapering anywhere (`grid.filled(0.0)` on an all-False mask, straight `fft2`). Peak is `max(|cont|)`.
Issue text: `issue.md`.

## What I don't know
Whether their released KilnSim.py pads or tapers (the text doesn't say, and my reproduction without
padding matches Table 1 to 0.3 nT, which strongly suggests it doesn't). What gradiometer and sensor
separation Storbekken used.

## Fix, if they want one
Either compute Bz directly at each height (cheap: 314 dipoles), or zero-pad the grid to ≥ 5× the
largest continuation height and subtract the window mean before continuing.
