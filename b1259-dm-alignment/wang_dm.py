"""Frequency-specific, outer-edge-referenced DMs for folded two-component profiles."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import gaussian_filter1d
from scipy.signal import resample

import DM_comparison_align_v2 as existing


def roll_fraction(profile, bins):
    k = np.fft.rfftfreq(len(profile))
    return np.fft.irfft(np.fft.rfft(profile) * np.exp(-2j*np.pi*k*bins), n=len(profile))


def groups(freq, edges):
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (freq >= lo) & (freq < hi if i < len(edges)-2 else freq <= hi)
        yield np.flatnonzero(mask)


def template(profile, windows, smooth, min_snr):
    n = len(profile)
    phase = np.arange(n)/n
    inside = [(phase >= lo) & (phase <= hi) for lo, hi in windows]
    off = ~(inside[0] | inside[1])
    if off.sum() < 8 or any(mask.sum() < 4 for mask in inside):
        raise ValueError("component windows need more bins or a larger off-pulse region")
    noise = np.std(profile[off])
    model = gaussian_filter1d(profile - np.median(profile[off]), smooth, mode="wrap")
    points = []
    for side, mask in enumerate(inside):
        bins = np.flatnonzero(mask)
        peak = bins[np.argmax(model[bins])]
        if noise <= 0 or model[peak] < min_snr * noise:
            raise ValueError("weak template component")
        half = model[peak]/2
        # Follow the edge outwards from the component peak until half height.
        steps = np.arange(peak, bins[0]-1, -1) if side == 0 else np.arange(peak, bins[-1]+1)
        below = steps[model[steps] <= half]
        if not len(below):
            raise ValueError("outer half-height crossing outside window")
        a = below[0]
        b = a + (1 if side == 0 else -1)
        points.append(a + (half-model[a])/(model[b]-model[a])*(b-a))
    centre = np.mean(points)
    return model, points, centre, roll_fraction(model, n/2-centre)


def measure(d, standards, edges, dm, sn_min):
    wf = existing.shift(d, dm - d["dm0"])
    n = wf.shape[1]
    rows = []
    for i, ch in enumerate(groups(d["freq"], edges)):
        if not len(ch) or standards[i] is None:
            continue
        ref = resample(standards[i], n)
        ref -= ref.mean()
        y = wf[ch].mean(axis=0)
        y -= y.mean()
        # Circular correlation and local fractional-bin peak interpolation.
        corr = np.fft.irfft(np.fft.rfft(y)*np.fft.rfft(ref).conj(), n=n)
        j = np.argmax(corr)
        left, peak, right = corr[j-1], corr[j], corr[(j+1)%n]
        curvature = 2*peak-left-right
        lag = j + ((right-left)/(2*curvature) if curvature > 0 else 0.)
        lag = (lag + n/2) % n - n/2
        shifted = roll_fraction(ref, lag)
        amp = np.dot(y, shifted)/np.dot(ref, ref)
        residual = y - amp*shifted
        sigma = 1.4826*np.median(abs(residual-np.median(residual)))
        deriv = np.fft.irfft(2j*np.pi*np.fft.rfftfreq(n)*np.fft.rfft(ref), n=n)
        if sigma <= 0 or amp <= 0 or np.linalg.norm(deriv) == 0:
            continue
        sn = amp*np.linalg.norm(ref)/sigma
        if sn >= sn_min:
            x = existing.KDM * (np.mean(d["freq"][ch]**-2)-d["fref"]**-2)*n/d["period"]
            rows.append((i, np.mean(d["freq"][ch]), x, lag, sigma/(amp*np.linalg.norm(deriv)), sn))
    if not rows:
        return np.empty((0, 6))
    return np.array(rows)


def fit(d, standards, edges, start, sn_min):
    dm = start
    n = d["data"].shape[1]
    for iteration in range(21):
        rows = measure(d, standards, edges, dm, sn_min)
        if len(rows) < 3:
            raise ValueError("fewer than three usable frequency templates")
        x, lag, sigma = rows[:, 2:5].T
        lag = np.unwrap(lag*2*np.pi/n)*n/(2*np.pi)
        design = np.column_stack([np.ones(len(x)), x])
        weighted = design/sigma[:, None]
        pars, _, rank, _ = np.linalg.lstsq(weighted, lag/sigma, rcond=None)
        if rank < 2:
            raise ValueError("insufficient dispersion leverage")
        residual = lag-design@pars
        redchi = np.sum((residual/sigma)**2)/(len(rows)-2)
        err = np.sqrt(np.linalg.inv(weighted.T@weighted)[1, 1]*max(redchi, 1.))
        if abs(pars[1]) <= max(1e-7, .02*err):
            return dm, err, rows, redchi, iteration, pars[1]
        if iteration < 20:
            dm += pars[1]
    raise ValueError(f"iteration limit at DM {dm:.5f}, remaining correction {pars[1]:+.5f}")


def plot_epoch(pdf, d, edges, standards, row, result, sn):
    dm, err, rows, chi, _, _ = result
    n = d["data"].shape[1]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), gridspec_kw={"height_ratios":[1, 3]},
                             constrained_layout=True)
    axes[0,2].axis("off")
    axes[0,2].text(0, .4, f"{len(rows)} sub-bands\nReduced χ² = {chi:.2f}")
    for col, (label, value) in enumerate([("Outer-edge DM", dm), ("pdmp DM", row["dm_pdmp"])]):
        top, ax = axes[:,col]
        if not np.isfinite(value):
            top.set_title("pdmp unavailable"); ax.axis("off"); continue
        image, prof, roll, px, fy = existing.small_waterfall(d, value-d["dm0"])
        top.plot((px[:-1]+px[1:])/2, prof, lw=1)
        top.set(xlim=(0,1), title=f"{label} = {value:.4f}")
        ax.pcolormesh(px, fy, image, rasterized=True, cmap="viridis",
                      vmin=np.nanpercentile(image,1), vmax=np.nanpercentile(image,99.5))
        bands = rows if col == 0 else measure(d, standards, edges, value, sn)
        if len(bands):
            loc = ((n/2 + bands[:,3] + roll + .5)/n) % 1
            ax.plot(loc, bands[:,1], "wx", ms=5)
            lag = np.unwrap(bands[:,3]*2*np.pi/n)*n/(2*np.pi)
            lag -= np.average(lag, weights=bands[:,4]**-2)
            ms = 1000*d["period"]/n
            axes[1,2].errorbar(lag*ms, bands[:,1], xerr=bands[:,4]*ms,
                                fmt="o" if col == 0 else "s", ms=3, label=label)
        ax.set(xlabel="Pulse phase", ylabel="Frequency (MHz)")
    axes[1,2].axvline(0, color=".6", lw=.7)
    axes[1,2].set(xlabel="Relative arrival delay (ms)", ylabel="Frequency (MHz)")
    axes[1,2].legend(fontsize=8)
    fig.suptitle(f"{row['file']}\nWhole-profile DM {dm:.4f} ± {err:.4f} pc cm⁻³ | × outer-edge midpoint")
    pdf.savefig(fig); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--input", type=Path)
    p.add_argument("--suffix", default=".ar")
    p.add_argument("--windows", type=float, nargs=4, metavar=("L0","L1","T0","T1"))
    p.add_argument("--nsub", type=int, default=16)
    p.add_argument("--smooth", type=float, default=2., help="reference Gaussian smoothing, bins")
    p.add_argument("--min-snr", type=float, default=5.)
    p.add_argument("--start-dm", type=float)
    p.add_argument("--no-pdmp", action="store_true")
    p.add_argument("-o", default="wang_dm")
    a = p.parse_args()
    if a.nsub < 3 or a.smooth <= 0 or a.min_snr <= 0:
        p.error("require nsub >= 3 and positive smoothing and S/N")
    r = existing.load(str(a.reference))
    n = r["data"].shape[1]
    phase = np.arange(n)/n
    edges = np.linspace(min(r["freq"]), max(r["freq"]), a.nsub+1)
    prefix = Path(a.o)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(str(prefix)+"_settings.json", "w") as fp:
        json.dump(vars(a), fp, default=str, indent=2)
    pdfpath = str(prefix)+"_templates.pdf"
    if a.windows:
        l0,l1,t0,t1 = a.windows
        if not 0 <= l0 < l1 < t0 < t1 <= 1:
            p.error("supply ordered non-overlapping windows in phase 0..1")
    standards = []
    edge_rows = []
    with PdfPages(pdfpath) as pdf:
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.plot(phase, r["data"].mean(axis=0))
        ax.set(xlabel="Reference pulse phase", ylabel="Mean normalised intensity",
               title="Quiet reference | choose a window around each full component", xlim=(0,1))
        ax.grid(alpha=.25)
        pdf.savefig(fig); plt.close(fig)
        if a.windows:
            fig, axs = plt.subplots((a.nsub+3)//4, 4, squeeze=False,
                                    figsize=(12, 2.2*((a.nsub+3)//4)), constrained_layout=True)
            for i, ch in enumerate(groups(r["freq"], edges)):
                ax = axs.flat[i]
                ax.set_title(f"{edges[i]:.0f}–{edges[i+1]:.0f} MHz", fontsize=9)
                try:
                    if not len(ch):
                        raise ValueError("empty band")
                    profile = r["data"][ch].mean(axis=0)
                    ax.plot(phase, profile, lw=.5, color=".6")
                    model, points, centre, standard = template(profile, [(l0,l1),(t0,t1)], a.smooth, a.min_snr)
                    edge_rows.append((i+1, edges[i], edges[i+1], points[0]/n, points[1]/n, centre/n))
                    ax.plot(phase, model, lw=.8)
                    for edge in points:
                        ax.axvline(edge/n, color="C3", lw=.9)
                    ax.axvline(centre/n, color="k", ls="--", lw=.8)
                except ValueError as exc:
                    standard = None
                    ax.text(.04,.85,str(exc),transform=ax.transAxes, fontsize=6)
                standards.append(standard)
            for ax in list(axs.flat)[a.nsub:]:
                ax.axis("off")
            pdf.savefig(fig); plt.close(fig)
    if a.windows:
        with open(str(prefix)+"_edges.csv", "w", newline="") as fp:
            writer = csv.writer(fp)
            writer.writerow(["subband","low_mhz","high_mhz","leading_edge_phase","trailing_edge_phase","fiducial_phase"])
            writer.writerows(edge_rows)
    print("Reference diagnostic:", pdfpath)
    if not a.input:
        return
    if not a.windows or sum(t is not None for t in standards) < 3:
        p.error("inspect reference diagnostic and supply windows yielding at least three templates")
    files = sorted(f for f in a.input.glob("*"+a.suffix) if f.is_file() and f.resolve()!=a.reference.resolve())
    if not files:
        p.error("no target archives matched")
    pd_dir = Path(str(prefix)+"_pdmp"); pd_dir.mkdir(exist_ok=True)
    summary = []
    if not a.no_pdmp:
        print("pdmp cross-check uses the whole input archive on the 140–175 grid.")
    with PdfPages(str(prefix)+".pdf") as pdf, open(str(prefix)+".csv","w",newline="") as out:
        fields = "file status dm_wang err_formal dm_pdmp err_pdmp nsub reduced_chi2 updates remaining_dm".split()
        writer = csv.DictWriter(out, fieldnames=fields); writer.writeheader()
        for path in files:
            row = {"file":path.name}
            try:
                d = existing.load(str(path))
                if d["source"] != r["source"]:
                    raise ValueError("reference and target source names differ")
                result = fit(d, standards, edges, d["dm0"] if a.start_dm is None else a.start_dm, a.min_snr)
                dm, err, rows, chi, updates, correction = result
                pd = dict(dm=np.nan, err=np.nan, img=None)
                if not a.no_pdmp:
                    pd = existing.run_pdmp(str(path), str(pd_dir), path.name,
                        ["-f","-do",str(157.5-d["dm0"]),"-dr","17.5","-ds","0.001"])
                row.update(status="converged",dm_wang=dm,err_formal=err,dm_pdmp=pd["dm"],
                           err_pdmp=pd["err"],nsub=len(rows),reduced_chi2=chi,updates=updates,
                           remaining_dm=correction)
                plot_epoch(pdf,d,edges,standards,row,result,a.min_snr)
                if not a.no_pdmp:
                    existing.pdmp_page(pdf,dict(file=path.name,pdmp_img=pd["img"]))
                summary.append(row.copy())
            except (RuntimeError, ValueError, OSError, np.linalg.LinAlgError) as exc:
                row["status"] = str(exc)
            writer.writerow(row); out.flush()
            print(path.name, row["status"], row.get("dm_wang",""))
        fig, ax = plt.subplots()
        valid = [q for q in summary if np.isfinite(q["dm_pdmp"])]
        if valid:
            x = np.array([q["dm_pdmp"] for q in valid])
            y = np.array([q["dm_wang"] for q in valid])
            ax.errorbar(x,y,yerr=[q["err_formal"] for q in valid],fmt="o",ms=3)
            ends = [min(x.min(),y.min()),max(x.max(),y.max())]
            ax.plot(ends,ends,"k--",lw=.8)
        else:
            ax.text(.5,.5,"No paired measurements",ha="center",transform=ax.transAxes)
        ax.set(xlabel="pdmp DM (pc cm⁻³)",ylabel="Outer-edge template DM (pc cm⁻³)",
               title=f"{len(summary)} successful fits | {len(files)} input archives")
        pdf.savefig(fig); plt.close(fig)
    if not summary:
        raise SystemExit("No successful fits; examine the status column in the CSV.")

if __name__ == "__main__":
    main()
