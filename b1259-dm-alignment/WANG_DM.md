# Whole-profile DM from the two outer edges

Wang, Johnston & Manchester (2004, MNRAS 351, 599) timed B1259−63 against a separate high-S/N standard at each frequency. Their common phase reference lay halfway between the midpoints of the steep outer edges of the two components. This script adapts that principle for a folded wideband archive. The smoothing, half-height edge finder, trial-DM iteration and error calculation below are implementation choices. [Paper, Sections 2–3 and Fig. 1](https://doi.org/10.1111/j.1365-2966.2004.07806.x)

## Inputs

Keep `wang_dm.py` and `DM_comparison_align_v2.py` from this branch in the same directory. Use Python with PSRCHIVE, NumPy, SciPy and Matplotlib. The comparison additionally uses the OzSTAR `pdmp` executable and Ghostscript for embedding PostScript. Use a separate set of templates per receiver/band.

Supply one high-S/N, frequency-resolved, full-period **two-component** archive from a quiet epoch, or a phase-coherent stack of independent quiet epochs. Inspect the dispersion correction in the reference: substantial smearing within its frequency bins or sub-bands will degrade the templates. Retain the original separation between the two components. If a component straddles phase zero, rotate the entire reference once before choosing windows. Exclusive component files or separately centred pulse1/pulse2 profiles lose the edge separation required here. Use the original whole-profile archives for the targets too. Do not include epochs that contributed to a stacked reference in an independent validation sample.

## Run

Start in the script directory and set `ref` to the actual quiet archive.

```bash
ref=/path/to/quiet_full_profile.ar
python wang_dm.py --reference "$ref" -o preview
```

Open `preview_templates.pdf`. Choose one phase interval enclosing the whole leading component and another enclosing the whole trailing component. Include baseline outside both outer edges. These example windows are illustrative and must be replaced for your archive.

```bash
python wang_dm.py --reference "$ref" --windows 0.22 0.42 0.55 0.76 -o reference_check
```

Inspect **each panel** of `reference_check_templates.pdf`. Red lines show the outer half-height crossings, and the black dashed line shows their midpoint. Confirm the same outer features are followed at every frequency, that neither crossing is a noise spike or an internal shoulder, and that template smoothing preserves the edges. The per-band crossings are saved in `reference_check_edges.csv`. Weak or unmeasurable reference bands are omitted. The programme does not assess component identities automatically.

```bash
python wang_dm.py --reference "$ref" --windows 0.22 0.42 0.55 0.76 \
  --input /path/to/whole_profile_epochs --suffix .clfd --nsub 16 -o p116_wang
```

Use a file suffix or a dedicated directory that selects **only the desired whole-profile archives**. In a mixed directory `.clfd` can also match gated components. The programme searches one directory level and excludes the reference path. Add `--no-pdmp` for a quick initial run. To check a strongly dispersed epoch for phase-wrap ambiguity, repeat it from another independently plausible absolute `--start-dm` and compare both the DM and the residuals. Separate input scripts and Slurm submission are optional; run on an allocated compute node.

## What the script measures

For each fixed sub-band, the smoothed reference profile is aligned by shifting the midpoint of the two outer half-height edge positions to phase 0.5. This preserves its frequency-specific widths, asymmetry and component amplitudes. The target is dedispersed at the trial DM; each band is cross-correlated against its own standard. The script fits a common DM correction plus a free phase offset to the sub-band lags. Fractional lags use parabolic interpolation of the circular correlation peak. Lags are unwrapped in frequency using the nearest phase turn. Both wrap errors and competing peaks are possible at low S/N.

Conditional lag errors use the residual noise, fitted template amplitude and template phase derivative. The DM standard error is derived from the inverse-variance-weighted line fit and enlarged by the square root of reduced chi-square when greater than one. Iteration stops at a remaining correction below 2% of that error (absolute floor 10⁻⁷ DM units); after 20 unsuccessful updates the fit is flagged. The formal error assumes high-S/N fixed templates, independent sub-band noise, a correct phase connection and an adequate profile model. It omits uncertainty from reference selection, outer-edge placement and time-variable scattering. Near periastron inspect the residuals and try an independent quiet template and narrower matched frequency coverage before physically interpreting a DM excursion. To compare across bands or with the older component estimates, use the same epochs and accepted frequencies wherever possible.

## Outputs

- `p116_wang_templates.pdf`, `p116_wang_edges.csv` and `p116_wang_settings.json` document the templates, edge locations and run arguments.
- `p116_wang.csv` contains one row per target with fit status, whole-profile DM, conditional formal error, pdmp comparison, retained sub-band count, reduced chi-square and convergence diagnostics.
- `p116_wang.pdf` shows template-fit and pdmp waterfalls, integrated profiles, **white crosses at fitted fiducial phases**, and relative sub-band delays with formal errors. A constant offset is removed in the delay panel. Original pdmp diagnostic pages follow, and the whole-sample comparison is on the last page.
- `p116_wang_pdmp/` holds raw logs and available PostScript plots. pdmp uses the full target and a per-archive offset to centre its grid at 157.5 pc cm⁻³, with half-range 17.5 and step 0.001. The template fit requires overlap with the reference band and sufficient S/N.

The estimator and plotting path were checked with two-component Gaussian profiles of frequency-dependent widths and amplitude ratios, independent target noise and nonzero achromatic phase offset. Injected DMs of 144 and 149 were recovered to within 0.001 pc cm⁻³ in these two example runs. This checks basic recovery and sign; calibrating interval coverage or scattering bias for a periastron epoch requires further data-specific tests.