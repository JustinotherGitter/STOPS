import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev
from scipy import signal


def continuum(w, spec, deg=11, std=1.6, steps=5, pos=False, plot=False) -> np.array:
    if plot:
        fig, axs = plt.subplots(2, 1, True)
        axs[0].plot(w, spec, label='data')

    nw = w.copy()
    nspec = spec.copy()

    for i in range(steps):
    
        p = chebyshev.chebfit(nw, nspec, deg=deg)
        ch = chebyshev.chebval(nw, p)
        diff = nspec - ch
        sigma = np.std(diff)
        
        ok = np.where( np.abs(diff) > std*sigma) if pos else np.where( diff > -1*std*sigma)
        
        nw = nw[ok]
        nspec = nspec[ok]

        if plot:
            axs[0].plot(w, chebyshev.chebval(w, p), label=f"fit {i}")
            
    if plot:
        axs[1].plot(w, spec/chebyshev.chebval(w,p), label="normalised data")
        for ax in axs: ax.legend()
        plt.show()
        
    return p


def filtered_continuum(spec, cutoff: float = 0.005, std: float = 1.6, steps: int = 5, plot: bool = False) -> np.ndarray:
    """
    Define the continuum as the low frequency signal of the spectrum
    Filtering as shown in https://swharden.com/blog/2020-09-23-signal-filtering-in-python/
    """
    if plot:
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=[20, 8])
        axs[0].plot(spec, label="data")
        
    data = np.ma.masked_where(spec == 0, spec, copy=True)
        
    for i in range(steps):
        b, a = signal.butter(1, cutoff, "lowpass")
        filt_cont = signal.filtfilt(b, a, data, method="gust")
        
        if plot:
            axs[0].plot(filt_cont, label=f"{i}")
        
        diff = data - filt_cont
        data = np.ma.masked_where(np.abs(diff) > std * data.std(), data)
        
    if plot:
        axs[1].plot(spec / filt_cont, label="normalised data")
        for ax in axs: ax.legend()
    
    return filt_cont


def grow(maskedarray: np.ma.masked_array, growth: int = 1) -> np.ma.masked_array:
    """
    Accepts a masked array and grows the mask by a specified amount
    """
    mArr = maskedarray.copy()
    for i, val in enumerate(maskedarray.mask):
        if not val: continue

        mArr.mask[max(0, i - growth): i] = True
        mArr.mask[i: min(i + growth + 1, len(mArr.mask))] = True
            
    return mArr