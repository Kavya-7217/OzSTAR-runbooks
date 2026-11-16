# Component RM on OzSTAR with native PSRSALSA

[RMSYNTH_v2_PSRSALSA.ipynb](RMSYNTH_v2_PSRSALSA.ipynb) identifies pulse components from a Parkes archive and runs the installed PSRSALSA `rmsynth` program separately for each component. Python selects the phase bins. PSRSALSA reads the original archive and computes every RM result.

## Run on OzSTAR

Start Jupyter in the environment used for your PSRSALSA build. The kernel needs PSRCHIVE’s Python bindings, NumPy and Matplotlib. The executable needs its usual shared libraries and PGPLOT colour PostScript driver.

```bash
module load gcc/11.3.0 fftw/3.3.10 gsl/2.7 cfitsio/4.2.0
conda activate psrsalsa_env
```

Set `INPUT` to a single archive or a directory. Directory mode selects files ending in `EXTENSION` (default `.center`). Use baseline-subtracted, full-Stokes data with one integrated subintegration and with Faraday rotation still present. Run the notebook cells in order.

## What the notebook measures

1. **Find the gates.** PSRCHIVE supplies a dedispersed, channel-weighted Stokes I profile. The notebook smooths it and selects up to two significant peaks. A valley separates resolved components. It labels two peaks leading and trailing in archive phase order. If it finds no peak, it saves the profile and skips the RM fit.
2. **Pass the gates to PSRSALSA.** Each component becomes its own `-onpulse` selection. The union of the detected components becomes `-onpulse2` in every call. PSRSALSA uses this second window for integrated signal and its complement for noise. Two disjoint components therefore share the same signal window for analytic errors. Wrapped gates use repeated inclusive bin ranges.
3. **Measure RM.** Native `rmsynth` computes the longitude-resolved spectra, sums their amplitudes over the component and fits its instrumental response. Its ASCII table gives fitted RM and analytic error. Set `BOOTSTRAP` to a positive draw count (e.g. 500) to add its Gaussian-bootstrap uncertainty. The default is 0. Set `LONGITUDE_FITS = True` to also fit each selected longitude bin with `-long`.

The default RM grid spans −10,000 to +10,000 rad/m² with 20,000 samples. Set `RM_LOW`, `RM_HIGH`, `RM_STEPS` and `FTOL` to match the comparison run. The default fit tolerance is 0.001. With `FIX_SEED = True`, the notebook passes the native `-fixseed` option.

## Inspect and reproduce a result

Each execution creates a new directory under `OUTPUT`. Open the archive’s `gates.png` to check the selected emission and `gates.json` to see the exact bins. Each component directory contains the archive symlink, exact invocation in `command.txt`, native `rmsynth.log`, `input.RMtable`, and native spectrum and RM-map plots (`spectrum.ps` and `map.ps`). Replay the recorded invocation from that directory. A collapsed table has one result, and a longitude table has one row per selected bin.

The notebook checks the source checkout against [PSRSALSA commit 2efbaf4](https://github.com/weltevrede/psrsalsa/tree/2efbaf4592e08fdb76002bbc9dc25e8e4441813b) and records the executable path and SHA-256 in `binary.json`. A checkout revision alone cannot establish how an untracked executable was built. To reproduce an OzSTAR result, match the executable hash, input archive, gates, native options and runtime environment.
