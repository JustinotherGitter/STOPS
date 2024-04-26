"""Utility functions for modules"""

# MARK: Imports
import os
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev
from astropy.io import fits as pyfits
from scipy import signal


# MARK: Find files
def find_files(
    data_dir: Path,
    filenames: list[str] | None = None,
    prefix: str = "mxgbp",
    ext: str = "fits",
) -> list[Path]:
    """
    Checks if `filenames` in `data_dir` are valid,
    otherwise, finds `prefix`*.`ext` files.

    Parameters
    ----------
    data_dir : str
        Directory path where the FITS files are located.
    filenames : List[str], optional
        List of filenames to search for. If provided, only these files will be searched for.
        (Default is None)
    prefix : str, optional
        Filename prefix to search for.
        (Default is "mxgbp")
    ext : str, optional
        File extension to search for.
        (Default is ".fits")

    Returns
    -------
    List[Path]
        List of ``pathlib.Path`` objects representing the valid or found files.

    Raises
    ------
    FileNotFoundError
        If `filenames` is provided and any of the files are not found
        in the `data_dir` directory.
    FileNotFoundError
        If `filenames` is not provided and no files matching the
        `prefix`*.`ext` regex search are found in the directory.
    """
    valid = []
    errMsg = ""
    # Search for files in `data_dir`
    if filenames is None:
        valid = [
            fl
            for fl in os.listdir(data_dir)
            if fl.startswith(prefix) and fl.endswith(ext)
        ]
        errMsg = f"No files matching '{prefix}*.{ext}' found in '{data_dir}'."

    else:
        # Check files are valid
        for fl in sorted(filenames):
            if not os.path.exists(data_dir / fl):
                logging.warning(
                    f"File '{fl}' not found in '{data_dir}'. Dropped from filenames."
                )
            else:
                valid.append(fl)

        errMsg = f"`{data_dir}` contains no files listed in `filenames`."

    # No files found or valid, raise error
    if not valid:
        logging.error(errMsg)
        raise FileNotFoundError(errMsg)

    logging.debug(f"find_files will parse and return: {valid}")
    return [Path(data_dir) / fl for fl in valid]


# MARK: Get Arc File
def find_arc(filenames: list[Path]) -> Path:
    for fl in filenames:
        with pyfits.open(fl) as hdu:
            if hdu["PRIMARY"].header["OBJECT"] == "ARC":
                logging.debug(f"find_arc returns arc: {fl}")
                return fl

    # Handle arc not found
    errMsg = f"No arc file found within provided files: {filenames}."
    logging.error(errMsg)
    raise FileNotFoundError(errMsg)


# MARK: Continuum
def continuum(w, spec, deg=11, std=1.6, steps=5, pos=False, plot=False) -> np.array:
    if plot:
        fig, axs = plt.subplots(2, 1, sharex=True)
        axs[0].plot(w, spec, label="data")

    nw = w.copy()
    nspec = spec.copy()

    for i in range(steps):
        p = chebyshev.chebfit(nw, nspec, deg=deg)
        ch = chebyshev.chebval(nw, p)
        diff = nspec - ch
        sigma = np.std(diff)

        ok = (
            np.where(np.abs(diff) > std * sigma)
            if pos
            else np.where(diff > -1 * std * sigma)
        )

        nw = nw[ok]
        nspec = nspec[ok]

        if plot:
            axs[0].plot(w, chebyshev.chebval(w, p), label=f"fit {i}")

    if plot:
        axs[1].plot(w, spec / chebyshev.chebval(w, p), label="normalised data")
        for ax in axs:
            ax.legend()
        plt.show()

    return p


# MARK: Filtered Continuum
def filtered_continuum(
    spec, cutoff: float = 0.005, std: float = 1.6, steps: int = 5, plot: bool = False
) -> np.ndarray:
    """
    Define the continuum as the low frequency signal of the spectrum.

    Parameters
    ----------
    spec
        An ``arrayLike`` list of the spectrum
    cutoff: float, optional
        The low frequency cut-off
        (The default is 0.005)
    std: float, optional
        The standard deviation
        (The default is 1.6)
    steps: int, optional
        The iterations for lowpass filtering
        (The default is 5)
    plot: bool, optional
        Parameter determining whether a plot should be returned
        (The default is False)

    Returns
    -------
    numpy.ndarray
        The spectrum `spec` after lowpass filtering

    See Also
    --------
    signal filtering:
        https://swharden.com/blog/2020-09-23-signal-filtering-in-python/
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
        for ax in axs:
            ax.legend()

    return filt_cont


# MARK: Grow
def grow(maskedarray: np.ma.masked_array, growth: int = 1) -> np.ma.masked_array:
    """
    Accepts a masked array and grows the mask by a specified amount
    """
    mArr = maskedarray.copy()
    for i, val in enumerate(maskedarray.mask):
        if not val:
            continue

        mArr.mask[max(0, i - growth) : i] = True
        mArr.mask[i : min(i + growth + 1, len(mArr.mask))] = True

    return mArr
