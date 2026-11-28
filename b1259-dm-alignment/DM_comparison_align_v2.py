#!/usr/bin/env python
"""Measure component DMs by sub-band alignment and compare with pdmp.

Run on one folded component per input archive. See README.md for statistical
assumptions, frequency-cut checks, quality flags and output definitions.
"""
import json
import shlex
import argparse, csv, glob, os, re, shutil, subprocess, sys, tempfile, traceback
from multiprocessing import Pool
import numpy as np
KDM = 4.148808e3  # s MHz^2 / (pc cm^-3)
# ---------------------------------------------------------------- read archive
def load(path):
    import psrchive
    ar = psrchive.Archive_load(path)
    ar.pscrunch()
    ar.dedisperse()
    ar.tscrunch()
    ar.remove_baseline()
    data = ar.get_data()[0, 0].astype(float)          # (nchan, nbin)
    w = ar.get_weights()[0]
    try:
        f = np.asarray(ar.get_frequencies(), float)
    except Exception:
        integ = ar.get_Integration(0)
        f = np.array([integ.get_centre_frequency(c) for c in range(ar.get_nchan())])
    f = f.reshape(-1)
    good = (w > 0) & np.isfinite(f) & (f > 0) & np.isfinite(data).all(axis=1) & (data.std(axis=1) > 0)
    data, f = data[good], f[good]
    # normalise every channel to zero median, unit noise
    med = np.median(data, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(data - med), axis=1, keepdims=True)
    data = (data - med) / np.where(mad > 0, mad, data.std(axis=1, keepdims=True))
    return dict(data=data, freq=f, period=ar.get_Integration(0).get_folding_period(),
                dm0=ar.get_dispersion_measure(), fref=ar.get_centre_frequency(),
                source=ar.get_source())
def shift(d, ddm):
    """Dedisperse the (already header-DM-dedispersed) data by an extra ddm."""
    nbin = d["data"].shape[1]
    k = np.arange(nbin // 2 + 1)
    delay_bins = KDM * ddm * (d["freq"] ** -2 - d["fref"] ** -2) / d["period"] * nbin
    X = np.fft.rfft(d["data"], axis=1) * np.exp(2j * np.pi * np.outer(delay_bins, k) / nbin)
    return np.fft.irfft(X, n=nbin, axis=1)
# ---------------------------------------------------------- alignment method
def subband_delays(d, ddm, nsub, up=16):
    """Delay (bins) and S/N of each sub-band relative to the rest of the band."""
    wf = shift(d, ddm)
    nbin = wf.shape[1]
    order = np.argsort(d["freq"])
    total = wf.sum(axis=0)
    res = []
    subbands = (np.split(order, np.flatnonzero(np.diff(d["subband"][order])) + 1)
                if "subband" in d else np.array_split(order, min(nsub, len(order))))
    for g in subbands:
        sub = wf[g].sum(axis=0)
        tmpl = total - sub                                   # rest of the band
        # low-pass the template: keep harmonics clearly above its noise level
        T = np.fft.rfft(tmpl); T[0] = 0
        noise_power = nbin * (len(order) - len(g))
        width = min(7, len(T))
        power = np.convolve(np.abs(T)**2, np.ones(width)/width, mode="same")
        T *= np.maximum(0., 1. - 1.5 * noise_power / np.maximum(power, 1e-30))
        tmpl = np.fft.irfft(T, n=nbin)
        ccf = np.fft.irfft(np.fft.rfft(sub) * np.conj(T), n=up * nbin) * up
        j = int(np.argmax(ccf))
        left, peak, right = ccf[j-1], ccf[j], ccf[(j+1) % ccf.size]
        curvature = 2 * peak - left - right
        offset = (right - left) / (2 * curvature) if curvature > 0 else 0.
        lag = ((j + offset) / up + nbin / 2) % nbin - nbin / 2
        snr = ccf[j] / (np.sqrt(len(g)) * np.linalg.norm(tmpl) + 1e-30)   # matched-filter S/N
        a = KDM * (d["freq"][g] ** -2 - d["fref"] ** -2) / d["period"] * nbin  # bins per unit DM
        peak_bin = np.argmax(tmpl) + lag if np.any(T) else np.nan
        res.append((d["freq"][g].mean(), a.mean(), lag, snr, peak_bin))
    f, a, lag, snr, peak_bin = map(np.array, zip(*res))
    return dict(f=f, a=a, lag=lag, snr=snr, peak_bin=peak_bin)
def line_fit(x, y, sig):
    """Weighted straight-line fit; errors scaled by the actual scatter (chi2/dof)."""
    w = 1 / sig ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ (A.T @ (w * y))
    chi2 = np.sum(w * (y - A @ p) ** 2); dof = max(len(x) - 2, 1)
    return p, np.sqrt(np.diag(cov) * chi2 / dof), chi2 / dof
def align_dm(d, nsub, min_snr, niter=20, start_offset=0.):
    """Evaluate every diagnostic at the returned DM, including the last slope."""
    if len(d["freq"]) < 3 or np.ptp(d["freq"]) == 0:
        raise ValueError("Need at least three usable channels spanning frequency")
    ddm = start_offset
    for iteration in range(niter + 1):
        sb = subband_delays(d, ddm, nsub)
        use = np.isfinite(sb["snr"]) & (sb["snr"] >= min_snr)
        if use.sum() < 3:
            raise ValueError(f"Only {use.sum()} usable sub-bands above the S/N cut")
        a, lag, snr = (sb[k][use] for k in ("a", "lag", "snr"))
        p, e, _ = line_fit(a, lag, 1 / snr)
        converged = abs(p[1]) <= max(1e-7, 0.05 * e[1])
        if converged or iteration == niter:
            break
        ddm += p[1]
    resid = lag - (p[0] + p[1] * a)
    rms_ms = np.sqrt(np.average(resid**2, weights=snr**2)) * d["period"] * 1e3 / d["data"].shape[1]
    return dict(dm=d["dm0"] + ddm, err=e[1], rms_ms=rms_ms, nused=int(use.sum()),
                use=use, sb=sb, updates=iteration, converged=converged, correction=p[1])


def subset(d, keep):
    out = {**d, "freq": d["freq"][keep], "data": d["data"][keep]}
    if "subband" in d:
        out["subband"] = d["subband"][keep]
    return out


def measure(d, opt):
    labels = np.empty(len(d["freq"]), int)
    for i, group in enumerate(np.array_split(np.argsort(d["freq"]), opt["nsub"])):
        labels[group] = i
    d = dict(d, subband=labels)
    start = opt.get("start_dm")
    offset = 0. if start is None else start - d["dm0"]
    fit = align_dm(d, opt["nsub"], opt["min_snr"], start_offset=offset)
    status = "ok" if fit["converged"] else "iteration_limit"
    count = opt["jackknife"]
    # Distribute each omitted group across frequency, preserving band coverage.
    groups = np.empty(len(d["freq"]), int)
    groups[np.argsort(d["freq"])] = np.arange(len(groups)) % count
    estimates = []
    if fit["converged"] and len(groups) >= max(count, opt["nsub"] * 2):
        for group in range(count):
            try:
                part = align_dm(subset(d, groups != group), opt["nsub"], opt["min_snr"],
                                start_offset=fit["dm"] - d["dm0"])
            except (ValueError, np.linalg.LinAlgError):
                continue
            if part["converged"]:
                estimates.append(part["dm"])
    error = np.nan
    if len(estimates) == count:
        fractions = np.bincount(groups, minlength=count) / len(groups)
        pseudo = fit["dm"] + (1 / fractions - 1) * (fit["dm"] - np.asarray(estimates))
        centre = np.average(pseudo, weights=fractions)
        error = np.sqrt(np.average((pseudo - centre)**2, weights=fractions) / (count - 1))
    if status == "ok" and not (np.isfinite(error) and error > 0):
        status = "uncertainty_unresolved"
    return {**fit, "status": status, "jack_err": error, "jack_ok": len(estimates)}

# ------------------------------------------------------------------- pdmp
def run_pdmp(path, outdir, base, extra):
    tmp = tempfile.mkdtemp()
    try:
        ps = os.path.join(tmp, "pdmp.ps")
        cmd = ["pdmp", *extra, "-g", ps + "/vcps", os.path.abspath(path)]
        try:
            result = subprocess.run(cmd, cwd=tmp, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
            out, code = result.stdout, result.returncode
        except FileNotFoundError as exc:
            out, code = str(exc), 127
        with open(os.path.join(outdir, base + "_pdmp.txt"), "w") as log:
            log.write(shlex.join(cmd) + f"\nexit={code}\n" + out)
        if code:
            return dict(dm=np.nan, err=np.nan, img=None)
        img = ps_to_png(ps)
        if os.path.exists(ps):
            shutil.copy(ps, os.path.join(outdir, base + "_pdmp.ps"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    num = r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)"
    m = re.search(r"Best DM\s*=\s*" + num, out)
    # e.g. "Best DM =  145.6256  Correction = -0.1789  Error = 0.06807"
    e = re.search(r"Best DM\s*=[^\n]*?Error\s*=\s*" + num, out)
    return dict(dm=float(m.group(1)) if m else np.nan,
                err=float(e.group(1)) if e else np.nan, img=img)
def ps_to_png(ps):
    """PostScript -> image array with ghostscript (None if gs is missing)."""
    if not os.path.exists(ps):
        return None
    png = ps + ".png"
    try:
        subprocess.run(["gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
                        "-r100", "-dFirstPage=1", "-dLastPage=1", "-sOutputFile=" + png, ps],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        import matplotlib.image as mpimg
        return (mpimg.imread(png)[..., :3] * 255).astype(np.uint8)
    except Exception:
        return None
# ------------------------------------------------------------- per archive
def pulse_label(name):
    m = re.search(r"pulse(\d+)", name)
    return f"pulse{m.group(1)}" if m else "other"
def small_waterfall(d, ddm, nch=64, nb=256):
    wf = shift(d, ddm)
    roll = wf.shape[1] // 2 - int(np.argmax(wf.sum(axis=0)))
    wf = np.roll(wf, roll, axis=1)
    phase_groups = np.array_split(np.arange(wf.shape[1]), min(nb, wf.shape[1]))
    profile = wf.sum(axis=0)
    prof = np.array([profile[g].mean() for g in phase_groups])
    edges = np.linspace(d["freq"].min(), d["freq"].max(), min(nch, len(wf)) + 1)
    idx = np.clip(np.searchsorted(edges, d["freq"], side="right") - 1, 0, len(edges)-2)
    image = np.full((len(edges)-1, len(prof)), np.nan)
    for i in np.unique(idx):
        row = wf[idx == i].mean(axis=0)
        image[i] = [row[g].mean() for g in phase_groups]
    phase_edges = np.r_[[g[0] for g in phase_groups], wf.shape[1]] / wf.shape[1]
    return image, prof, roll, phase_edges, edges
def process(job):
    path, opt = job
    base = os.path.basename(path)[: -len(opt["ext"])]
    try:
        d = load(path)
        al = measure(d, opt)
        cuts = []
        for low in opt["cuts"]:
            band = subset(d, d["freq"] >= low)
            row = dict(cut_mhz=low, fmin=np.nan, fmax=np.nan, dm=np.nan, error=np.nan,
                       shift=np.nan, status="no_channels")
            if len(band["freq"]):
                row.update(fmin=band["freq"].min(), fmax=band["freq"].max())
                try:
                    cut = measure(band, opt)
                    row.update(dm=cut["dm"], error=cut["jack_err"], shift=cut["dm"]-al["dm"],
                               status=cut["status"])
                except (ValueError, np.linalg.LinAlgError) as exc:
                    row["status"] = str(exc)
            cuts.append(row)
        extra = opt["pdmp_args"]
        if opt["pdmp_centre"] is not None:
            extra = ["-do", format(opt["pdmp_centre"] - d["dm0"], ".15g"), *extra]
        pd = run_pdmp(path, opt["outdir"], base, extra)
        nbin = d["data"].shape[1]
        sb_al = al["sb"]
        sb_pd = subband_delays(d, pd["dm"] - d["dm0"], opt["nsub"]) if np.isfinite(pd["dm"]) else None
        wf_al = small_waterfall(d, al["dm"] - d["dm0"])
        wf_pd = small_waterfall(d, pd["dm"] - d["dm0"]) if np.isfinite(pd["dm"]) else None
        return dict(ok=True, file=base, pulse=pulse_label(base), source=d["source"],
                    period=d["period"], nbin=nbin, fmin=d["freq"].min(), fmax=d["freq"].max(),
                    dm_al=al["dm"], err_al=al["jack_err"], fit_err=al["err"], quality=al["status"], updates=al["updates"],
                    correction=al["correction"], jack_ok=al["jack_ok"], cuts=cuts, header_dm=d["dm0"], rms_ms=al["rms_ms"], nused=al["nused"],
                    use=al["use"], dm_pd=pd["dm"], err_pd=pd["err"], pdmp_img=pd["img"],
                    sb_al=sb_al, sb_pd=sb_pd, wf_al=wf_al, wf_pd=wf_pd)
    except Exception as e:
        if opt["debug"]:
            raise
        return dict(ok=False, file=base, error=f"{type(e).__name__}: {e!r}", tb=traceback.format_exc())
# ------------------------------------------------------------------- plots
COL = {"pulse1": "C0", "pulse2": "C3", "other": "C2"}
def first_page(pdf, res):
    import matplotlib.pyplot as plt
    groups = sorted({r["pulse"] for r in res})
    fig, axs = plt.subplots(1, len(groups), figsize=(11.69, 6), squeeze=False)
    fig.suptitle("Alignment DM vs pdmp DM", fontsize=14, weight="bold")
    for ax, g in zip(axs[0], groups):
        rr = [r for r in res if r["pulse"] == g and np.isfinite(r["dm_pd"]) and r["quality"] == "ok"]
        if not rr:
            ax.set_title(g); ax.text(0.1, 0.5, "No valid paired measurements"); continue
        x = np.array([r["dm_pd"] for r in rr]); xe = np.array([r["err_pd"] for r in rr])
        y = np.array([r["dm_al"] for r in rr]); ye = np.array([r["err_al"] for r in rr])
        ax.errorbar(x, y, xerr=np.nan_to_num(xe), yerr=ye, fmt="o", ms=4, color=COL.get(g, "k"),
                    ecolor="0.5", capsize=2)
        lo, hi = np.nanmin([x, y]), np.nanmax([x, y]); pad = 0.05 * (hi - lo + 1e-3)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8, label="1:1")
        txt = f"N = {len(rr)}"
        if len(rr) >= 3:
            dd = y - x
            txt += f"\nmean(align − pdmp) = {dd.mean():+.3f}\nstd = {dd.std():.3f}"
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_title(g); ax.set_xlabel("pdmp DM (pc cm\\(^{-3}\\))")
        ax.set_ylabel("Alignment DM (pc cm\\(^{-3}\\))"); ax.legend(loc="lower right", fontsize=8)
        ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    pdf.savefig(fig); plt.close(fig)
def archive_page(pdf, r):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11.69, 8.27))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 3], hspace=0.05)
    axs = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[:, 2])]
    tops = [fig.add_subplot(gs[0, 0], sharex=axs[0]), fig.add_subplot(gs[0, 1], sharex=axs[1])]
    fig.suptitle(f"{r['file']}  ({r['pulse']})\n"
                 f"Alignment DM ({r['quality']}) = {r['dm_al']:.3f} ± {r['err_al']:.3f}      "
                 f"pdmp DM = {r['dm_pd']:.3f} ± {r['err_pd']:.3f}      "
                 f"difference = {r['dm_al'] - r['dm_pd']:+.3f}", fontsize=11)
    u = r["use"]
    for ax, top, W, sb, lab, c in [(axs[0], tops[0], r["wf_al"], r["sb_al"], "Alignment DM", "C0"),
                                   (axs[1], tops[1], r["wf_pd"], r["sb_pd"], "pdmp DM", "C3")]:
        if W is None:
            top.set_title("pdmp DM not found"); continue
        wf, prof, roll, phase_edges, freq_edges = W
        top.plot((phase_edges[:-1] + phase_edges[1:]) / 2, prof, color=c, lw=0.8)
        top.set_title(lab); top.set_yticks([]); top.tick_params(labelbottom=False)
        ax.pcolormesh(phase_edges, freq_edges, wf, cmap="viridis", shading="flat", rasterized=True,
                      vmin=np.nanpercentile(wf, 1), vmax=np.nanpercentile(wf, 99.5))
        phase = ((sb["peak_bin"] + roll + 0.5) / r["nbin"]) % 1
        show = u & np.isfinite(phase)
        ax.plot(phase[show], sb["f"][show], "wx", ms=5, mew=1.2,
                label="Sub-band arrival phase")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Pulse phase"); ax.set_ylabel("Frequency (MHz)")
    ms = r["period"] / r["nbin"] * 1e3
    ax = axs[2]
    ax.plot(r["sb_al"]["lag"][u] * ms, r["sb_al"]["f"][u], "o", color="C0", label="at alignment DM")
    if r["sb_pd"] is not None:
        sb = r["sb_pd"]
        ax.plot(sb["lag"][u] * ms, sb["f"][u], "s", color="C3", label="at pdmp DM")
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("Sub-band delay (ms)"); ax.set_ylabel("Frequency (MHz)")
    ax.set_title(f"Sub-band delays\n({r['nused']} used, rms resid = {r['rms_ms']:.3f} ms)")
    ax.legend(fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.86, wspace=0.3)
    pdf.savefig(fig); plt.close(fig)
def pdmp_page(pdf, r):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.suptitle(f"{r['file']}  --  pdmp", fontsize=10)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.93]); ax.axis("off")
    if r["pdmp_img"] is None:
        ax.text(0.5, 0.5, "pdmp plot not available (is ghostscript 'gs' loaded?)", ha="center")
    else:
        ax.imshow(r["pdmp_img"])
    pdf.savefig(fig, dpi=150); plt.close(fig)
# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory"); ap.add_argument("extension")
    ap.add_argument("-j", "--nproc", type=int, default=os.cpu_count())
    ap.add_argument("--nsub", type=int, default=16, help="number of sub-bands (default 16)")
    ap.add_argument("--min-snr", type=float, default=5, help="min sub-band S/N used in fit (default 5)")
    ap.add_argument("-o", "--out", default="dm_align")
    ap.add_argument("--debug", action="store_true", help="run first file only and show the full error")
    ap.add_argument("--start-dm", type=float, help="optional absolute alignment seed (default: archive header)")
    ap.add_argument("--check-fmin", type=float, action="append", default=[],
                    help="repeat alignment above this frequency in MHz (repeatable)")
    ap.add_argument("--jackknife", type=int, default=8, help="interleaved channel groups")
    ap.add_argument("--pdmp-centre", type=float, help="absolute DM grid centre converted to a per-archive pdmp -do offset")
    ap.add_argument("--pdmp-args", default="", help="quoted extra options passed to pdmp")
    a = ap.parse_args()
    if a.pdmp_centre is not None and (not np.isfinite(a.pdmp_centre) or
                                      {"-dm", "-do"}.intersection(shlex.split(a.pdmp_args))):
        ap.error("--pdmp-centre requires a finite DM and supplies the -do option itself")
    if a.jackknife < 3 or a.nsub < 3 or a.min_snr <= 0 or a.nproc < 1:
        ap.error("require jackknife and nsub >= 3, positive S/N and worker count")
    if a.start_dm is not None and not np.isfinite(a.start_dm):
        ap.error("start DM must be finite")
    if any(not np.isfinite(f) or f <= 0 for f in a.check_fmin):
        ap.error("frequency cuts must be positive finite MHz")
    files = sorted(glob.glob(os.path.join(a.directory, "*" + a.extension)))
    if not files:
        sys.exit("No files found")
    outdir = a.out + "_pdmp"; os.makedirs(outdir, exist_ok=True)
    opt = dict(ext=a.extension, nsub=a.nsub, min_snr=a.min_snr, outdir=outdir, debug=a.debug,
               cuts=a.check_fmin, start_dm=a.start_dm, jackknife=a.jackknife, pdmp_centre=a.pdmp_centre, pdmp_args=shlex.split(a.pdmp_args))
    with open(a.out + "_settings.json", "w") as fh:
        json.dump(vars(a), fh, indent=2)
    if a.debug:
        r = process((files[0], opt))
        print({k: v for k, v in r.items() if k in ("file", "dm_al", "err_al", "rms_ms", "nused", "dm_pd", "err_pd", "quality", "updates", "correction", "cuts")})
        return
    res = []
    failures = []
    with Pool(min(a.nproc, len(files))) as pool:
        for i, r in enumerate(pool.imap(process, [(f, opt) for f in files]), 1):
            if r["ok"]:
                res.append(r)
                print(f"[{i}/{len(files)}] {r['file']}: align {r['dm_al']:.3f}±{r['err_al']:.3f}  pdmp {r['dm_pd']:.3f}")
            else:
                failures.append(r)
                print(f"[{i}/{len(files)}] {r['file']}: FAILED {r['error']}\n{r['tb']}")
    with open(a.out + "_failures.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["file", "error"])
        writer.writerows((r["file"], r["error"]) for r in failures)
    if not res:
        sys.exit("Nothing processed.")
    with open(a.out + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "pulse", "dm_align", "err_align", "rms_resid_ms", "nsub_used", "dm_pdmp", "err_pdmp", "diff",
                    "quality", "updates", "remaining_dm", "fit_err_internal", "jack_ok", "jack_requested",
                    "header_dm", "fmin_mhz", "fmax_mhz"])
        for r in res:
            w.writerow([r["file"], r["pulse"], r["dm_al"], r["err_al"], f"{r['rms_ms']:.4f}",
                        r["nused"], r["dm_pd"], r["err_pd"], f"{r['dm_al'] - r['dm_pd']:+.4f}", r["quality"], r["updates"],
                        r["correction"], r["fit_err"], r["jack_ok"], a.jackknife, r["header_dm"],
                        r["fmin"], r["fmax"]])
    if a.check_fmin:
        with open(a.out + "_cuts.csv", "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["file", "pulse", "cut_mhz", "fmin", "fmax",
                                                    "dm", "error", "shift", "status"])
            writer.writeheader()
            for r in res:
                for row in r["cuts"]:
                    writer.writerow(dict(file=r["file"], pulse=r["pulse"], **row))
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(a.out + ".pdf") as pdf:
        first_page(pdf, res)
        for r in res:
            archive_page(pdf, r)
            pdmp_page(pdf, r)
    print(f"Wrote {a.out}.pdf, {a.out}.csv (raw pdmp output + PS in {outdir}/)")
if __name__ == "__main__":
    main()
