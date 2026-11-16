# B1259 component rotation measures on OzSTAR

[RMSYNTH_v2_PSRSALSA.ipynb](RMSYNTH_v2_PSRSALSA.ipynb) measures rotation measure (RM) from folded, full-Stokes pulsar archives. It runs two analyses on each archive and writes separate results for each detected pulse component.

## What happens with the current settings?

1. PSRCHIVE converts the archive to Stokes parameters, dedisperses at the header DM and averages the subintegrations. The notebook sorts channels by frequency and identifies usable channels from their weights and finite data.
2. The notebook smooths the band-averaged Stokes I profile and searches for up to two significant peaks. A gate extends to the two-sigma baseline crossing or the valley between peaks. Two detected components are labelled leading and trailing in archive phase order. The original analysis plots these gates for inspection.
3. The original analysis subtracts the off-pulse baseline and sums Q + iU across each gate before measuring its Faraday spectrum. It refines the sampled maximum with a three-point parabola and estimates uncertainty from 500 draws of measured off-pulse spectra. Its trial grid represents total RM and accounts for an RM already removed from the archive.
4. The PSRSALSA-style analysis forms a Faraday spectrum for each selected longitude bin and adds their amplitudes within each gate. It fits a shifted instrumental response, with free scale and offset, to the collapsed spectrum. This keeps emission with opposite polarisation vectors in different bins from cancelling in the gate sum. This branch requires an archive before Faraday correction and assumes the baseline has already been removed.
5. The PSRSALSA-style branch calculates an analytic error from the fitted-RM derotated profile, debiased linear polarisation and nominal channel frequencies. Its Gaussian bootstrap is disabled by default. When enabled, each draw adds independent Q and U noise using the channel’s off-pulse Stokes I standard deviation and refits the spectrum. The reported uncertainty is the population standard deviation of those fitted RMs.

## Run the notebook

Use a Jupyter Python kernel with PSRCHIVE, NumPy, Matplotlib and a working C compiler available as `cc`. The notebook compiles its bundled upstream simplex optimiser in a temporary directory. Set `INPUT` to one archive or an archive directory and run the cells in order. In directory mode the default extension is `.center`. Both analysis branches are enabled. For the synthetic check alone, set `RUN_ORIGINAL = RUN_COMPARISON = False`.

| Setting | Current value or behaviour |
| --- | --- |
| `FREQ_MHZ` | Full band |
| `RM_LIMITS`, `GRID_POINTS` | −10,000 to +10,000 rad/m², 20,000 samples including endpoints |
| `REGIONS` | Automatic detection. Override the comparison gates with e.g. `{"leading": [(80, 91)], "trailing": [(104, 115)]}`. Bin indices are zero-based and inclusive. |
| `ONPULSE2` | Union of comparison gates. Its complement supplies the noise samples. Override independently to match the C run. The original branch uses its own padded off-pulse mask. |
| `PS_GRID`, `FTOL` | Same initial grid as the original branch, fit tolerance 0.001 |
| `PS_DRAWS`, `SEED` | 0 draws, seed 42. Use at least 2 draws to report a Gaussian-bootstrap error. |
| `PS_LONG` | False. Per-bin spectra are saved. Set True for individual bin RM fits and uncertainties. |
| `SUBTRACT_BASELINE` | False for the comparison. Set True to subtract the per-channel mean outside `ONPULSE2`. |

## Read the outputs

The default output directory is `files_close_to_peri_2e4/`, relative to the notebook’s working directory.

- `component_rm.csv`, `*.stokes.png` and `*.rm.png` contain the original results and diagnostics. Check the gates and search-boundary flags.
- `psrsalsa_rm.csv` contains one collapsed fit per region and optional per-bin fits. Compare `rm` with `grid_peak` and inspect the fitted curve for sidelobe solutions. A blank or NaN `sigma_bootstrap` is expected with zero draws.
- `*.ps.npz` contains the trial grid, archive bin indices, RM-by-bin amplitude array, collapsed spectrum, fitted curve, response frequencies, noise estimates and bootstrap RM draws. The comparison plots appear inline.

## Reproducibility

The reference is [PSRSALSA commit 2efbaf4, 22 February 2022](https://github.com/weltevrede/psrsalsa/tree/2efbaf4592e08fdb76002bbc9dc25e8e4441813b). Match the archive, preprocessing, frequency grid, bin windows, RM grid, fit tolerance and draw count to the C run. The synthetic control checks recovery of an injected RM of +7.5 rad/m², the collapsed spectrum against an independent channel sum, and repeatability at fixed seed. Agreement with the installed OzSTAR binary on a B1259 archive remains unverified. Differences can arise from channel handling, floating-point reductions and the NumPy and GSL random generators.
