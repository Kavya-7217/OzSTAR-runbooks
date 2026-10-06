# Flagging the First N Integrations of a Scan

This workflow flags the first **N integrations** of a selected scan in a Measurement Set (MS).

## 1. Activate the environment

```bash
conda activate /fred/oz002/skavya/miniconda3/envs/ms_flags
```

## 2. Run the flagging script

Example:

```bash
/fred/oz341/skavya/scripts/flag_timerange.py 1773330393_2026-09-09T23-59-12_fkb.ms -s 3 -n 1
```

Here:

- `1773330393_2026-09-09T23-59-12_fkb.ms` is the input Measurement Set.
- `-s 3` selects scan 3.
- `-n 1` flags the first integration of that scan.
- `-l` can optionally be used to specify the output log file.

For example:

```bash
/fred/oz341/skavya/scripts/flag_timerange.py 1773330393_2026-09-09T23-59-12_fkb.ms -s 3 -n 5 -l flag_scan3.log
```

## 3. Script

```python
#!/usr/bin/env python3

import argparse, logging, datetime
import numpy as np
from casatools import table as tb
from casatasks import flagdata

def mjds_to_casa(t):
    unix = (t / 86400.0 - 40587.0) * 86400.0
    return datetime.datetime.fromtimestamp(
        unix, tz=datetime.timezone.utc
    ).strftime("%Y/%m/%d/%H:%M:%S.%f")[:-3]

parser = argparse.ArgumentParser(
    description="Flag first N integrations of a selected scan."
)
parser.add_argument("ms", help="Input Measurement Set")
parser.add_argument("-s", "--scan", type=int, required=True, help="Scan number")
parser.add_argument("-n", "--nint", type=int, required=True, help="Number of integrations to flag")
parser.add_argument("-l", "--log", default="flag_integrations.log", help="Output log file")
args = parser.parse_args()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(args.log), logging.StreamHandler()]
)
log = logging.getLogger(__name__)

t = tb()
t.open(args.ms)
time_col = t.getcol("TIME")
scan_col = t.getcol("SCAN_NUMBER")
interval_col = t.getcol("INTERVAL")
t.close()

mask = scan_col == args.scan

if not np.any(mask):
    raise ValueError(
        f"Scan {args.scan} not found. Available scans: {np.unique(scan_col)}"
    )

if args.nint <= 0:
    raise ValueError("Number of integrations must be greater than zero.")

times = time_col[mask]
intervals = interval_col[mask]
unique_times = np.sort(np.unique(times))

if args.nint > len(unique_times):
    raise ValueError(
        f"Requested {args.nint} integrations, but scan "
        f"{args.scan} has only {len(unique_times)}."
    )

first_time, last_time = unique_times[0], unique_times[args.nint - 1]

first_int = np.median(intervals[np.isclose(times, first_time)])
last_int = np.median(intervals[np.isclose(times, last_time)])

t_start = first_time - first_int / 2.0
t_end = last_time + last_int / 2.0

timerange = f"{mjds_to_casa(t_start)}~{mjds_to_casa(t_end)}"

log.info(
    f"MS: {args.ms}, scan: {args.scan}, "
    f"integrations to flag: {args.nint}/{len(unique_times)}"
)
log.info(f"Flagging timerange: {timerange}")

flagdata(
    vis=args.ms,
    scan=str(args.scan),
    mode="manual",
    timerange=timerange,
    action="apply",
    flagbackup=True
)

log.info("Flagging complete.")
```

## Notes

- The script uses the MS `TIME` column to identify the unique integrations in the selected scan.
- `TIME` corresponds to the midpoint of an integration, so the start and end of the selected integrations are calculated using the `INTERVAL` column.
- CASA `flagdata` is then used to flag the complete timerange containing the requested integrations.
- `flagbackup=True` creates a CASA flag backup before applying the new flags.
- The default log file is `flag_integrations.log`.
