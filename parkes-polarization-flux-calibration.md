# Polarization and Flux Calibration for Parkes UWL Data

Steps for calibrating folded pulsar archives from Parkes UWL data, for getting peak flux and phase-averaged flux values. Based on the procedure used for MTP0096 single-pulse calibration, adapted for folded profiles.

## Do you need polarization calibration first?

Not strictly, if all you want is peak flux and phase-averaged flux for a pulsar — polarization calibration can be skipped by passing the `-X` flag to `pac`.

That said, it's worth running anyway: it only adds a small amount of processing time to the calibration step, unless there are no noise diode observations available for that session.

## Steps

### 1. Fold the noise diode observation

```bash
dspsr -c 0.089903802930864 -L 20 -A --scloffs <noise-diode-files> -e calib
```

(period `-c` is pulsar-specific — this value is just an example)

### 2. Change the archive type to PolnCal

```bash
psredit -c type=PolnCal -m <filename>.calib
```

### 3. Clean the noise diode square wave

Use `pazi` interactively on the polcal archive to remove RFI/artifacts from the square wave.

### 4. Build the calibration database

```bash
pac -w -u pcm -u cal -u fluxcal
```

This creates `database.txt`. The `fluxcal` file itself comes from the [Parkes Calibration and Data Processing Files page](https://www.parkes.atnf.csiro.au/observing/Calibration_and_Data_Processing_Files.html) — download the one matching the relevant calibrator/epoch.

### 5. Apply the calibration to the target archive

```bash
pac -cTS -d database.txt target_archive
```

Note the `-c` flag here — it's needed alongside `-T` and `-S` when applying from this database.

## Things to double-check if flux values look off

A few points that came up debugging another folded dataset where flux values came out roughly 5x lower than expected:

- **Signal disappearing when `--scloffs` is used on the target fold.** If `--scloffs` was used when folding the noise diode but *not* when folding the target, and adding `--scloffs` to the target fold makes the pulsar signal vanish, that's worth isolating on its own — try folding the target with and without `--scloffs` and compare the two profiles directly before doing anything else, since a downstream flux mismatch is hard to diagnose while this is still unresolved.
- **Unexpected file extensions from `pac -w`.** Watch what `-u` flags actually produce — e.g. `pac -w -u fluxcal -u pcm -u dzT` produces a `.dzT` file, which isn't part of the flow above. Stick to `-u pcm -u cal -u fluxcal` when building the database for this procedure.
- **Which flags to apply the calibration with.** Use `pac -cTS -d database.txt target_archive` (as above) rather than `pac -d database.txt -TFSb <archive>` — the extra `-c` flag matters for how the flux calibration solution gets applied.
- **Order of operations.** Polarization-calibrate first (`pac -cT` or similar against the polcal database), then flux-calibrate the polarization-calibrated archive — not the other way around.

