"""Exercise the v2 alignment core on synthetic folded profiles."""
import numpy as np
from DM_comparison_align_v2 import align_dm, measure, subset

rng = np.random.default_rng(1259)
freq = np.linspace(704., 4032., 64)
bins = np.arange(256)
period = .04776
dm0 = 146.7
options = {"nsub": 16, "min_snr": 5, "jackknife": 8}
for dm, tau in [(149., 0.), (144., 0.), (149., .00037)]:
    errors, reported, cut_bias = [], [], []
    rejected = 0
    for trial in range(60):
        delay = 4148.808 * (dm-dm0) * (freq**-2 - 2000.**-2) * len(bins) / period
        pulse = 4*np.exp(-0.5*((bins[None,:]-100-delay[:,None])/4)**2)
        if tau:
            tau_bins = tau * (freq/1400.)**-4 * len(bins) / period
            kernel = np.exp(-bins[None,:]/tau_bins[:,None])
            kernel /= kernel.sum(axis=1, keepdims=True)
            pulse = np.fft.irfft(np.fft.rfft(pulse)*np.fft.rfft(kernel))
        data = pulse + rng.normal(size=pulse.shape)
        data -= np.median(data, axis=1, keepdims=True)
        data /= 1.4826*np.median(abs(data), axis=1, keepdims=True)
        archive = dict(data=data, freq=freq, fref=2000., period=period, dm0=dm0)
        fit = measure(archive, options)
        if trial == 0:
            assert not align_dm(archive, 16, 5, niter=0)["converged"]
        if fit["status"] != "ok":
            rejected += 1
            continue
        if tau:
            high = align_dm(subset(archive, freq >= 1400.), 16, 5)
            if high["converged"]:
                cut_bias.append(high["dm"]-dm)
        errors.append(fit["dm"]-dm)
        reported.append(fit["jack_err"])
    covered = np.count_nonzero(abs(np.array(errors)) <= reported)
    print(f"DM={dm}, tau={tau}: valid={len(errors)}/60, bias={np.mean(errors):.5f}, "
          f"median_error={np.median(reported):.5f}, covered={covered}/{len(errors)}")
    assert len(errors) >= 50
    if tau:
        print(f"Mean bias above 1400 MHz: {np.mean(cut_bias):.5f}")
        assert abs(np.mean(cut_bias)) < abs(np.mean(errors))
    else:
        assert abs(np.mean(errors)) < 0.005
        assert 0.4 < covered / len(errors) < 0.9
