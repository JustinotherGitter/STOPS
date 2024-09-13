"""Utility functions for STOPS modules and classes."""

# MARK: Imports
import os
import logging
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev
from astropy.io import fits as pyfits
from scipy import signal


# MARK: Is Arc
def is_arc(file: Path) -> bool:
    """
    Check if the FITS file is an `arc` file.

    Parameters
    ----------
    file : Path
        A Path object representing the FITS file to check.

    Returns
    -------
    bool
        Whether the file is an `arc` file or not.
    """
    is_arc: bool = False

    # Open the FITS file and check the OBJECT keyword
    with pyfits.open(file) as hdul:
        is_arc = hdul["PRIMARY"].header["OBJECT"] == "ARC"
    
    return is_arc


# MARK: Get Arc Lamp
def get_arc_lamp(file: Path) -> str:
    """
    Get the arc lamp from the FITS file.

    Parameters
    ----------
    file : Path
        The FITS file to attempt to extract the arc lamp from.

    Returns
    -------
    str
        The arc lamp name, formatted as a `.txt` file.

    Raises
    ------
    ValueError
        If the file provided is not an arc file.
    """
    lamp: str = ""

    # Check if the file is an arc file
    if not is_arc(file):
        errMsg = f"File '{file}' is not an arc file."
        logging.error(errMsg)
        raise ValueError(errMsg)
    
    # Open the FITS file and get the LAMPID keyword
    with pyfits.open(file) as hdul:
        lamp = hdul["PRIMARY"].header["LAMPID"]

    # Append a file extension, assumed a `.txt` file (sourced from polsalt).
    lamp += '.txt'

    return lamp


# MARK: Find files
def find_files(
    data_dir: Path,
    filenames: list[str] | None = None,
    prefix: str = "mxgbp",
    ext: str = "fits",
    sep_arc: bool = False,
) -> list[Path] | tuple[list[Path], list[Path]]:
    """
    Checks if `filenames` in `data_dir` are valid, otherwise,
    finds `prefix`*.`ext` files.

    Parameters
    ----------
    data_dir : Path
        Directory path where the FITS files are located.
    filenames : list[str] | None, optional
        List of filenames to search for. If provided,
        only these files will be searched for,
        by default None.
    prefix : str, optional
        Filename prefix to search for,
        by default "mxgbp".
    ext : str, optional
        File extension to search for,
        by default "fits".
    sep_arc : bool, optional
        Additional return for arc files,
        by default False.

    Returns
    -------
    list[Path] | tuple[list[Path], list[Path]]
        List of Path objects representing the valid or found files, or a
        tuple of lists representing the valid files and arc files,
        respectively.

    Raises
    ------
    FileNotFoundError
        If `filenames` is provided and any of the files are not found in the
        `data_dir` directory, or
        if `filenames` is **not** provided and **no** files matching the
        `prefix`*.`ext` regex search are found in the directory.
    """
    valid = []
    valid_arcs = []
    errMsg = ""

    # Search for files in `data_dir`
    if filenames is None:
        for fl in os.listdir(data_dir):
            if fl.startswith(prefix) and fl.endswith(ext):
                if is_arc(data_dir / fl) and sep_arc:
                    valid_arcs.append(fl)
                else:
                    valid.append(fl)

        errMsg = f"No files matching '{prefix}*.{ext}' found in '{data_dir}'."

    else:
        # Check files are valid
        for fl in sorted(filenames):
            if not os.path.exists(data_dir / fl):
                logging.warning(
                    f"File '{fl}' not found in '{data_dir}'. Dropped from filenames."
                )
            else:
                if is_arc(data_dir / fl) and sep_arc:
                    valid_arcs.append(fl)
                else:
                    valid.append(fl)

        errMsg = f"`{data_dir}` contains no files listed in `filenames`."

    # No files found or valid, raise error
    if not valid:
        logging.error(errMsg)
        raise FileNotFoundError(errMsg)

    # Return valid files, with the arc files
    if sep_arc:
        logging.debug(f"find_files - {valid} and {valid_arcs}")
        return [
            Path(data_dir) / fl
            for fl in valid
        ], [
            Path(data_dir) / fl
            for fl in valid_arcs
        ]

    logging.debug(f"find_files will parse and return: {valid}")

    return [Path(data_dir) / fl for fl in valid]


# MARK: Get Arc File
def find_arc(filenames: list[Path]) -> Path:
    """
    Find the arc file from a list of files.

    Parameters
    ----------
    filenames : list[Path]
        The list of files to search for the arc file.

    Returns
    -------
    Path
        The Path object representing the arc file.

    Raises
    ------
    FileNotFoundError
        If no arc file is found within the provided files.
    """

    # Check if any of the files are arc files
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
def continuum(
    wav: np.ndarray,
    spec: np.ndarray,
    deg: int = 11,
    std: float = 1.6,
    steps: int = 5,
    pos: bool = False,
    plot: bool = False
) -> np.array:
    """
    Define the continuum using a Chebyshev polynomial fit.

    Parameters
    ----------
    wav : np.ndarray
        The wavelength array related to the spectrum.
    spec : np.ndarray
        The one-dimensional spectrum array.
    deg : int, optional
        The polynomial degree,
        by default 11.
    std : float, optional
        The standard deviation to sigma clip the spectrum by when iterating,
        by default 1.6.
    steps : int, optional
        The amount of iterations to perform for sigma clipping,
        by default 5.
    pos : bool, optional
        The boolean deciding whether the absolute difference should
        be considered or not,
        by default False.
    plot : bool, optional
        The boolean deciding whether to display additional information
        as plots,
        by default False.

    Returns
    -------
    np.array
        The Chebyshev polynomial fit, found using `chebfit`, to the spectrum.

    References
    ----------
    Continuum fitting:
        van Soelen, B., 2019,
        "PHYS6854 - Computational Physics",
        University of the Free State.
    """
    # Ignore RankWarning, obvious to the user when Rank is poorly conditioned
    warnings.simplefilter('ignore', np.exceptions.RankWarning)

    if plot:
        fig, axs = plt.subplots(2, 1, sharex=True)
        axs[0].plot(wav, spec, label="data")

    # Copy the arrays, such that the originals are not modified
    nw = wav.copy()
    nspec = spec.copy()

    # Iteratively fit the Chebyshev polynomial
    for i in range(steps):
        p = chebyshev.chebfit(nw, nspec, deg=deg)
        ch = chebyshev.chebval(nw, p)
        diff = nspec - ch
        sigma = np.std(diff)

        # Sigma-clip the spectrum
        ok = (
            np.where(np.abs(diff) > std * sigma)
            if pos
            else np.where(diff > -1 * std * sigma)
        )

        # Mask the arrays
        nw = nw[ok]
        nspec = nspec[ok]

        if plot:
            axs[0].plot(wav, chebyshev.chebval(wav, p), label=f"fit {i}")

    if plot:
        axs[1].plot(
            wav,
            spec / chebyshev.chebval(wav, p),
            label="normalised data"
        )
        for ax in axs:
            ax.legend()
        plt.show()

    return p


# MARK: Filtered Continuum
def filtered_continuum(
    spec: np.ndarray,
    cutoff: float = 0.005,
    std: float = 1.6,
    steps: int = 5,
    plot: bool = False
) -> np.ndarray:
    """
    Define the continuum as the low frequency signal of the spectrum.

    Parameters
    ----------
    spec : np.ndarray
        The spectrum to filter.
    cutoff : float, optional
        The low frequency cut-off,
        by default 0.005.
    std : float, optional
        The standard deviation,
        by default 1.6.
    steps : int, optional
        The amount of iterations for lowpass filtering,
        by default 5.
    plot : bool, optional
        The boolean deciding whether to display additional information
        as plots,
        by default False.

    Returns
    -------
    np.ndarray
        The filtered spectrum, after lowpass filtering.

    See Also
    --------
    signal filtering:
        https://swharden.com/blog/2020-09-23-signal-filtering-in-python/
    """
    if plot:
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=[20, 8])
        axs[0].plot(spec, label="data")

    # Mask the data where the spectrum is zero
    data = np.ma.masked_where(spec == 0, spec, copy=True)

    # Iteratively filter the spectrum
    for i in range(steps):
        b, a = signal.butter(1, cutoff, "lowpass")
        filt_cont = signal.filtfilt(b, a, data, method="gust")

        if plot:
            axs[0].plot(filt_cont, label=f"{i}")

        # Mask the data where the difference is greater than a desired std.
        diff = data - filt_cont
        data = np.ma.masked_where(np.abs(diff) > std * data.std(), data)

    if plot:
        axs[1].plot(spec / filt_cont, label="normalised data")
        for ax in axs:
            ax.legend()

    return filt_cont


# MARK: Grow
def grow(
    maskedarray: np.ma.masked_array,
    growth: int = 1
) -> np.ma.masked_array:
    """
    Grows the mask of a masked array by a specified amount.

    Parameters
    ----------
    maskedarray : np.ma.masked_array
        The masked array to grow the mask of.
    growth : int, optional
        The amount by which to grow the mask,
        by default 1.

    Returns
    -------
    np.ma.masked_array
        The masked array with the grown mask.
    """
    # Copy the masked array
    mArr = maskedarray.copy()

    # Grow the mask
    for i, val in enumerate(maskedarray.mask):
        if not val:
            continue

        mArr.mask[max(0, i - growth): i] = True
        mArr.mask[i: min(i + growth + 1, len(mArr.mask))] = True

    return mArr
