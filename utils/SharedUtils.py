# Shared helper functions for convenience
import os
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev
from astropy.io import fits as pyfits
from scipy import signal

# TODO@JustinotherGitter: Add pathlike (filelike?) typing


def get_files(data_dir: str, filenames: list[str] | None, prefix: str="m", extention: str="fits") -> list[os.PathLike]:
    nfiles = []

    if filenames == None:
        # Handle finding valid files
        nfiles = Path(data_dir).glob(f"{prefix}*{'.' if extention else ''}{extention}")
        nfiles = list(sorted(nfiles))

    else:
        # Handle recieving list of files
        for file in sorted(filenames):
            # Validate
            if Path(data_dir, file).is_file():
                nfiles.append(file)
            else:
                errMsg = f"{file} not found in the data directory {data_dir}, file unlisted."
                logging.warn(errMsg)

    # raise error if nfiles is empty
    if nfiles == []:
        if filenames == None:
            errMsg = f"Wildcard search, '{data_dir}/{prefix}*{'.' if extention else ''}{extention}', "
        else:
            errMsg = f"Filenames provided, {filenames}, "
        errMsg += "returned no matches."
        logging.error(errMsg)
        raise FileNotFoundError(errMsg)
    
    logging.debug(f"Files returned: {nfiles}")
    return nfiles


def get_arc(filenames: str, exclude_arc: bool=False) -> str:
    # No files provided
    if filenames == []:
        errMsg = f"No files to search for the arc in"
        logging.error(errMsg)
        raise FileNotFoundError(errMsg)
    
    # Handle exclusion of arc
    if exclude_arc:
        return ''

    # Handle inclusion of arc
    for file in filenames:
        with pyfits.open(file) as hdu:
            if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                return file

    # Handle arc not found
    errMsg = f"No arc file found within provided files, {filenames}."
    logging.error(errMsg)
    raise FileNotFoundError(errMsg)


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