# Polarization and Flux Calibration for Parkes UWL Data

Steps for calibrating folded pulsar archives from Parkes UWL data, for getting peak flux and phase-averaged flux values. 

## Steps

### 1. Fold the noise diode observation

```bash
dspsr -c 0.089903802930864 -L 20 -A --scloffs <noise-diode-files> -e cal
```


### 2. Change the archive type to PolnCal

```bash
psredit -c type=PolnCal -m <filename>.cal
```

### 3. Clean the noise diode square wave

Use `pazi` interactively on the polcal archive to remove RFI/artifacts from the square wave.

### 4. Build the calibration database

```bash
pac -w -u pcm -u cal -u fluxcal
```

This creates `database.txt`. The file with extension `fluxcal` should be downloaded from the [Parkes Calibration and Data Processing Files page](https://www.parkes.atnf.csiro.au/observing/Calibration_and_Data_Processing_Files.html) — download the one closer to the relevant epoch. The file with name `pcm.fits` can be shared upon request. File with extension `.cal` is the polarisation calibrator built from step 2. Make sure all these files are present inside the current working directory.

### 5. Apply the calibration to the target archive

```bash
pac -cTS -d database.txt target_archive
```

Note the `-c` flag here — it's needed alongside `-T` and `-S` when applying from this database.
