#!/usr/bin/env python
"""
dm_align_simple.py -- measure DM by aligning the pulse in frequency vs time, compare with pdmp.
Method (per archive):
  1. Load with PSRCHIVE: pscrunch, dedisperse at header DM, tscrunch.
  2. Split the band into sub-bands; cross-correlate each sub-band profile with the
     rest of the band to get its arrival-time delay.
  3. Fit  delay = c + dDM * 4.15e3 * (f^-2 - fref^-2)  ->  DM = DM_header + dDM.
     Repeat until converged (the pulse is then vertical in freq vs time).
  4. Run pdmp on the same archive.
PDF output:
  page 1   : alignment DM vs pdmp DM (with error bars), pulse1 and pulse2 fitted separately
  then per archive: (a) waterfalls at both DMs + sub-band delays, (b) pdmp's own plot
Usage:
  python dm_align_simple.py <directory> <extension> [-j 16] [--nsub 16] [--min-snr 5]
  python dm_align_simple.py DM_latest/P116/ .T --debug      # first file only, full error
"""
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
    good = (w > 0) & (data.std(axis=1) > 0)
    data, f = data[good], f[good]
    # normalise every channel to zero median, unit noise
    med = np.median(data, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(data - med), axis=1, keepdims=True)
    data = (data - med) / np.where(mad > 0, mad, 1)
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
    for g in np.array_split(order, nsub):
        sub = wf[g].sum(axis=0)
        tmpl = total - sub                                   # rest of the band
        # low-pass the template: keep harmonics clearly above its noise level
        T = np.fft.rfft(tmpl); T[0] = 0
        T[np.abs(T) ** 2 < 3 * nbin * (len(order) - len(g))] = 0
        tmpl = np.fft.irfft(T, n=nbin)
        ccf = np.fft.irfft(np.fft.rfft(sub) * np.conj(T), n=up * nbin) * up
        j = int(np.argmax(ccf))
        lag = j / up if j < up * nbin / 2 else j / up - nbin
        snr = ccf[j] / (np.sqrt(len(g)) * np.linalg.norm(tmpl) + 1e-30)   # matched-filter S/N
        a = KDM * (d["freq"][g] ** -2 - d["fref"] ** -2) / d["period"] * nbin  # bins per unit DM
        res.append((d["freq"][g].mean(), a.mean(), lag, snr))
    f, a, lag, snr = map(np.array, zip(*res))
    return dict(f=f, a=a, lag=lag, snr=snr)
def line_fit(x, y, sig):
    """Weighted straight-line fit; errors scaled by the actual scatter (chi2/dof)."""
    w = 1 / sig ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ (A.T @ (w * y))
    chi2 = np.sum(w * (y - A @ p) ** 2); dof = max(len(x) - 2, 1)
    return p, np.sqrt(np.diag(cov) * chi2 / dof), chi2 / dof
def align_dm(d, nsub, min_snr, niter=10):
    ddm = 0.0
    for _ in range(niter):
        sb = subband_delays(d, ddm, nsub)
        use = sb["snr"] >= min_snr
        if use.sum() < 3:
            raise RuntimeError(f"only {use.sum()} sub-bands with S/N >= {min_snr}")
        p, e, red = line_fit(sb["a"][use], sb["lag"][use], 1 / sb["snr"][use])
        ddm += p[1]                       # slope = remaining DM error
        if abs(p[1]) < 0.05 * e[1]:
            break
    resid = sb["lag"][use] - (p[0] + p[1] * sb["a"][use])           # bins
    w = sb["snr"][use] ** 2
    rms_ms = np.sqrt(np.sum(w * resid ** 2) / w.sum()) * d["period"] / d["data"].shape[1] * 1e3
    return dict(dm=d["dm0"] + ddm, err=e[1], rms_ms=rms_ms, nused=int(use.sum()), use=use)
# ------------------------------------------------------------------- pdmp
def run_pdmp(path, outdir, base):
    tmp = tempfile.mkdtemp()
    try:
        ps = os.path.join(tmp, "pdmp.ps")
        out = subprocess.run(["pdmp", "-g", ps + "/vcps", os.path.abspath(path)], cwd=tmp,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
        open(os.path.join(outdir, base + "_pdmp.txt"), "w").write(out)     # raw output kept
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
    wf = shift(d, ddm)[np.argsort(d["freq"])[::-1]]          # high freq at top
    roll = wf.shape[1] // 2 - int(np.argmax(wf.sum(axis=0)))  # peak at phase 0.5
    wf = np.roll(wf, roll, axis=1)
    bn0 = np.array_split(np.arange(wf.shape[1]), min(nb, wf.shape[1]))
    prof = np.array([wf.sum(axis=0)[b].mean() for b in bn0])    # frequency-averaged profile
    ch = np.array_split(np.arange(wf.shape[0]), min(nch, wf.shape[0]))
    bn = np.array_split(np.arange(wf.shape[1]), min(nb, wf.shape[1]))
    wf = np.array([[wf[np.ix_(c, b)].mean() for b in bn] for c in ch])
    return wf, prof
def process(job):
    path, opt = job
    base = os.path.basename(path)[: -len(opt["ext"])]
    try:
        d = load(path)
        al = align_dm(d, opt["nsub"], opt["min_snr"])
        pd = run_pdmp(path, opt["outdir"], base)
        nbin = d["data"].shape[1]
        sb_al = subband_delays(d, al["dm"] - d["dm0"], opt["nsub"])
        sb_pd = subband_delays(d, pd["dm"] - d["dm0"], opt["nsub"]) if np.isfinite(pd["dm"]) else None
        wf_al = small_waterfall(d, al["dm"] - d["dm0"])
        wf_pd = small_waterfall(d, pd["dm"] - d["dm0"]) if np.isfinite(pd["dm"]) else None
        return dict(ok=True, file=base, pulse=pulse_label(base), source=d["source"],
                    period=d["period"], nbin=nbin, fmin=d["freq"].min(), fmax=d["freq"].max(),
                    dm_al=al["dm"], err_al=al["err"], rms_ms=al["rms_ms"], nused=al["nused"],
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
        rr = [r for r in res if r["pulse"] == g and np.isfinite(r["dm_pd"])]
        x = np.array([r["dm_pd"] for r in rr]); xe = np.array([r["err_pd"] for r in rr])
        y = np.array([r["dm_al"] for r in rr]); ye = np.array([r["err_al"] for r in rr])
        ax.errorbar(x, y, xerr=np.nan_to_num(xe), yerr=ye, fmt="o", ms=4, color=COL.get(g, "k"),
                    ecolor="0.5", capsize=2)
        lo, hi = np.nanmin([x, y]), np.nanmax([x, y]); pad = 0.05 * (hi - lo + 1e-3)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.8, label="1:1")
        txt = f"N = {len(rr)}"
        if len(rr) >= 3:
            p, e, _ = line_fit(x, y, ye)
            xx = np.array([lo - pad, hi + pad])
            ax.plot(xx, p[0] + p[1] * xx, color=COL.get(g, "k"), lw=1,
                    label=f"fit: slope {p[1]:.3f} ± {e[1]:.3f}")
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
                 f"Alignment DM = {r['dm_al']:.3f} ± {r['err_al']:.3f}      "
                 f"pdmp DM = {r['dm_pd']:.3f} ± {r['err_pd']:.3f}      "
                 f"difference = {r['dm_al'] - r['dm_pd']:+.3f}", fontsize=11)
    u = r["use"]
    for ax, top, W, sb, lab, c in [(axs[0], tops[0], r["wf_al"], r["sb_al"], "Alignment DM", "C0"),
                                   (axs[1], tops[1], r["wf_pd"], r["sb_pd"], "pdmp DM", "C3")]:
        if W is None:
            top.set_title("pdmp DM not found"); continue
        wf, prof = W
        top.plot(np.linspace(0, 1, prof.size, endpoint=False), prof, color=c, lw=0.8)
        top.set_title(lab); top.set_yticks([]); top.tick_params(labelbottom=False)
        ax.imshow(wf, aspect="auto", cmap="viridis", extent=[0, 1, r["fmin"], r["fmax"]],
                  vmin=np.percentile(wf, 1), vmax=np.percentile(wf, 99.5))
        ax.plot(0.5 + sb["lag"][u] / r["nbin"], sb["f"][u], "wx", ms=5)   # sub-band arrival times
        ax.set_xlabel("Pulse phase"); ax.set_ylabel("Frequency (MHz)")
    ms = r["period"] / r["nbin"] * 1e3
    ax = axs[2]
    ax.plot(r["sb_al"]["lag"][u] * ms, r["sb_al"]["f"][u], "o", color="C0", label="at alignment DM")
    if r["sb_pd"] is not None:
        sb = r["sb_pd"]
        ax.plot(sb["lag"][u] * ms, sb["f"][u], "s", color="C3", label="at pdmp DM")
        p, _, _ = line_fit(sb["a"][u], sb["lag"][u], 1 / sb["snr"][u])
        o = np.argsort(sb["f"]); o = o[u[o]]
        ax.plot((p[0] + p[1] * sb["a"][o]) * ms, sb["f"][o], "C3--", lw=1, label=r"\\(\nu^{-2}\\) fit")
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
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.directory, "*" + a.extension)))
    if not files:
        sys.exit("No files found")
    outdir = a.out + "_pdmp"; os.makedirs(outdir, exist_ok=True)
    opt = dict(ext=a.extension, nsub=a.nsub, min_snr=a.min_snr, outdir=outdir, debug=a.debug)
    if a.debug:
        r = process((files[0], opt))
        print({k: v for k, v in r.items() if k in ("file", "dm_al", "err_al", "rms_ms", "nused", "dm_pd", "err_pd")})
        return
    res = []
    with Pool(min(a.nproc, len(files))) as pool:
        for i, r in enumerate(pool.imap(process, [(f, opt) for f in files]), 1):
            if r["ok"]:
                res.append(r)
                print(f"[{i}/{len(files)}] {r['file']}: align {r['dm_al']:.3f}±{r['err_al']:.3f}  pdmp {r['dm_pd']:.3f}")
            else:
                print(f"[{i}/{len(files)}] {r['file']}: FAILED {r['error']}\n{r['tb']}")
    if not res:
        sys.exit("Nothing processed.")
    with open(a.out + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "pulse", "dm_align", "err_align", "rms_resid_ms", "nsub_used", "dm_pdmp", "err_pdmp", "diff"])
        for r in res:
            w.writerow([r["file"], r["pulse"], f"{r['dm_al']:.4f}", f"{r['err_al']:.4f}", f"{r['rms_ms']:.4f}",
                        r["nused"], f"{r['dm_pd']:.4f}", f"{r['err_pd']:.4f}", f"{r['dm_al'] - r['dm_pd']:+.4f}"])
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
