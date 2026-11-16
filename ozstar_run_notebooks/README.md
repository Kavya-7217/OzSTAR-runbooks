# Python RM synthesis for automatically detected pulse components

[RMSYNTH_v2_PSRSALSA.ipynb](RMSYNTH_v2_PSRSALSA.ipynb) implements the default instrumental-response technique from [PSRSALSA commit 2efbaf4](https://github.com/weltevrede/psrsalsa/tree/2efbaf4592e08fdb76002bbc9dc25e8e4441813b). Python and NumPy perform the synthesis, fit and uncertainty calculations. It runs without the rmsynth executable, a PSRSALSA checkout or a C compiler.

## Run

Use a Jupyter kernel with NumPy and Matplotlib. Archive processing also requires the PSRCHIVE Python bindings. Set `INPUT` to one baseline-subtracted, single-subintegration, full-Stokes archive before Faraday correction, or to a directory of archives. Directory mode selects `.center` files by default. Run the cells in order. To run just the injected-signal test, set `RUN_ARCHIVES = False`.

## What the notebook measures

- PSRCHIVE dedisperses the archive. A smoothed, weighted Stokes I profile identifies up to two pulse components. The notebook saves their inclusive bin ranges and uses their union as the on-pulse window for uncertainty estimates.
- At each trial RM, the code derotates Q+iU in every channel and sums the channel vectors separately for each selected longitude. The component spectrum is the sum of those **bin amplitudes**. This preserves signal from bins with different polarisation angles.
- The fitted curve is a shifted instrumental response with free scale and offset. The fit uses the entire spectrum, the grid midpoint as its starting RM, and the pinned simplex update and convergence rules. Compare the fitted RM with the sampled peak to identify possible sidelobe solutions.
- Each component gets an analytic error. Set `BOOTSTRAP` to at least 2 (e.g. 500) for the population spread of Gaussian-noise refits. The noise scale comes from each channel’s off-pulse Stokes I. Set `LONGITUDE_FITS = True` to fit the individual phase bins as well.

## Settings and outputs

The initial search covers −10,000 to +10,000 rad/m² at 20,000 grid points, with `FTOL = 0.001` and bootstrapping disabled. For comparisons, match the RM grid, phase bins and signal window. Each additional bootstrap draw repeats the search and fit. For timing on your machine, use one archive first.

Each run creates a fresh folder under `python_psrsalsa_rm/`. The archive folder contains `gates.png` and `gates.json`, plus a CSV table, spectrum/map PNG and compressed NumPy archive per component. The CSV records the RM, grid peak, analytic and bootstrap uncertainties, and fit diagnostics. The NPZ contains the RM samples, selected bin indices, amplitude map, collapsed spectrum, fitted curve, channel frequencies and bootstrap estimates.

## Reproducibility

The synthetic control checks a known positive RM, independent channel sums, response shape, opposite polarisation in adjacent bins, and repeatability of Gaussian draws. A comparison on a real B1259 archive remains outstanding. Differences from the C program can arise from archive preprocessing, channel flags, floating-point accumulation and NumPy versus GSL random draws. The notebook zeros unusable channels and retains their positions in the nominal instrumental response. It derives nominal channel centres from the signed archive bandwidth and centre frequency. Check those conventions, especially for irregular channel grids, before interpreting a discrepancy.
