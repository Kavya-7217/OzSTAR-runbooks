"""Exercise the v2 alignment core on synthetic folded profiles."""
import numpy as np
from DM_comparison_align_v2 import align_dm

rng = np.random.default_rng(1259)
freq = np.linspace(704., 4032., 64)
bins = np.arange(256)
period = .04776
dm0 = 146.7
for dm, tau in [(149., 0.), (144., 0.), (149., .00037)]:
    errors, reported = [], []
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
        fit = align_dm(dict(data=data, freq=freq, fref=2000., period=period, dm0=dm0),16,5)
        errors.append(fit["dm"]-dm)
        reported.append(fit["err"])
    print(dm,tau,np.mean(errors),np.median(reported),sum(abs(np.array(errors)) <= reported))
