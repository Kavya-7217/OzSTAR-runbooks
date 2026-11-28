"""Inspect the sub-band 3 lag switch in the 12 May trailing-component archive."""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def normalise(p):
    p = p - np.median(p)
    return p / max(np.max(abs(p)), 1e-30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--script", required=True, type=Path, help="DM script used for the original fit")
    parser.add_argument("-o", default="may12_pulse2_band3_peaks.pdf")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("dm_fit", args.script)
    dm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dm)
    archive = dm.load(str(args.archive))
    order = np.argsort(archive["freq"])
    labels = np.empty(len(order), int)
    for i, group in enumerate(np.array_split(order, 16)):
        labels[group] = i
    omit = order[::8]  # jackknife group 1
    keep = np.ones(len(order), bool)
    keep[omit] = False
    d = dm.subset(dict(archive, subband=labels), keep)
    group = d["subband"] == 2  # sub-band 3
    freqs = d["freq"][group]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), constrained_layout=True)
    overlays = []
    peaks = []
    trials = (153.643497, 153.645765)
    for row, trial in enumerate(trials):
        ddm = trial - d["dm0"]
        wf = dm.shift(d, ddm)
        n = wf.shape[1]
        sub = wf[group].sum(axis=0)
        ref = np.fft.rfft(wf.sum(axis=0) - sub)
        ref[0] = 0
        width = min(7, len(ref))
        power = np.convolve(abs(ref)**2, np.ones(width)/width, mode="same")
        ref *= np.maximum(0., 1. - 1.5*n*(~group).sum()/np.maximum(power, 1e-30))
        up = 16
        ccf = np.fft.irfft(np.fft.rfft(sub)*ref.conj(), n=up*n)*up
        if not ccf.max() > 0:
            raise ValueError("No measurable correlation peak")
        ccf /= ccf.max()
        lag = dm.subband_delays(d, ddm, 16)["lag"][2]
        grid_peak = np.argmax(ccf)/up
        mismatch = abs((grid_peak-lag+n/2) % n - n/2)
        if mismatch > 1/up:
            raise RuntimeError("Reconstructed peak differs from the production script; check its template filter.")
        peaks.append(lag)
        lag_grid = (np.arange(up*n)-up*n//2)/up
        ax, profile_ax = axes[row]
        ax.plot(lag_grid, np.fft.fftshift(ccf), color="0.2")
        ax.set(xlabel="Template lag (phase bins)", ylabel="Correlation / maximum",
               title=f"Trial DM {trial:.6f}")
        # Show the same pulse region for both candidate template alignments.
        centre = np.argmax(np.fft.irfft(ref, n=n))
        x = (np.arange(n)-centre+n//2) % n - n//2
        idx = np.argsort(x)
        profile_ax.plot(x[idx], normalise(sub)[idx], color="0.65", label="Sub-band profile")
        profile_ax.set(xlabel="Bins from reference peak", ylabel="Relative amplitude")
        # Store the current template so both lag candidates can be drawn on each row.
        overlays.append((ref, n, x, idx))
    for row, (ax, profile_ax) in enumerate(axes):
        ref, n, x, idx = overlays[row]
        for i, lag in enumerate(peaks):
            color = ("tab:blue", "tab:orange")[i]
            label = f"{lag:+.2f} bins" + (" (selected)" if i == row else "")
            ax.axvline(lag, color=color, ls="--", label=label)
            shifted = np.fft.irfft(ref*np.exp(-2j*np.pi*np.arange(len(ref))*lag/n), n=n)
            profile_ax.plot(x[idx], normalise(shifted)[idx], color=color, ls="--", label=label)
        margin = max(8., abs(peaks[0]-peaks[1]))
        ax.set_xlim(min(peaks)-margin, max(peaks)+margin)
        profile_ax.set_xlim(min(0., *peaks)-2*margin, max(0., *peaks)+2*margin)
        ax.legend(fontsize=8)
        profile_ax.legend(fontsize=8)
    fig.suptitle(f"12 May trailing component | jackknife group 1 omitted | sub-band 3: {freqs.min():.0f}–{freqs.max():.0f} MHz")
    fig.savefig(args.o)
    plt.close(fig)
    print(f"Selected lags: {peaks[0]:+.3f}, {peaks[1]:+.3f} bins. Wrote {args.o}")


if __name__ == "__main__":
    main()
