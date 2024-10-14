#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module for cross correlating polarization beams."""

# MARK: Imports
import sys
import logging
import itertools as iters
from pathlib import Path
from importlib.resources import files
from typing import Callable

import numpy as np
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
import matplotlib.axes
from astropy.io import fits as pyfits
from scipy import signal

from STOPS.utils.SharedUtils import find_files, continuum
from STOPS.utils.Constants import SAVE_CORR, OFFSET
import STOPS.utils

# MARK: Logging init.
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.INFO)


# MARK: Correlate class
class CrossCorrelate:
    """
    Cross correlate allows for comparing the extensions of multiple
    FITS files, or comparing the $O$- and $E$-beams of a single FITS file.

    Parameters
    ----------
    data_dir : str | Path
        The path to the data to be cross correlated
    filenames : list[str]
        The ecwmxgbp*.fits files to be cross correlated.
        If only one filename is defined, correlation is done against the two
        polarization beams.
    split_ccd : bool, optional
        Decides whether the CCD regions should each be individually
        cross correlated.
        (The default is True, which splits the spectrum up into its separate
        CCD regions)
    cont_ord : int, optional
        The degree of a chebyshev to fit to the continuum.
        (The default is 11)
    plot : bool, optional
        Decides whether the continuum fitting should be plotted
        (The default is False, so no continua plots are displayed)
    save_prefix : str, optional
        The name or directory to save the figure produced to.
        "." saves a default name to the current working. A default name is
        also used when save_prefix is a directory.
        (The default is None, I.E. The figure is not saved, only displayed)

    Attributes
    ----------
    data_dir
    fits_list
    beams : str
        The mode of correlation.
        'OE' for same file, and 'O' or 'E' for different files but same
        extension.
    ccds : int
        The number of CCD's in the data.
        Used to split the CCD's if split_ccd is True.
    cont_ord : int
        The degree of the chebyshev to fit to the continuum.
    can_plot : bool
        Decides whether the continuum fitting should be plotted
    offset : int, DEPRECATED
        The amount the spectrum is shifted, mainly to test the effect of the
        cross correlation.
        (The default is 0, I.E. no offset introduced)
    save_prefix
    wav_unit : str
        The units of the wavelength axis.
        (The default is Angstroms)
    wav_cdelt : int
        The wavelength increment.
        (The default is 1)
    alt : Callable
        An alternate method of cross correlating the data.
        (The default is None)

    Methods
    -------
    load_file(filename: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]
        Loads the data from a FITS file.
    get_bounds(bpm: np.ndarray) -> np.ndarray
        Finds the bounds for the CCD regions.
    remove_cont(spec: list, wav: list, bpm: list, plot_cont: bool) -> None
        Removes the continuum from the data.
    correlate(filename1: Path, filename2: Path | None = None) -> None
        Cross correlates the data.
    ftcs(filename1: Path, filename2: Path | None = None) -> None
        Cross correlates the data using the Fourier Transform.
    plot(spec, wav, bpm, corrdb, lagsdb) -> None
        Plots the data.
    process() -> None
        Processes the data.

    Other Parameters
    ----------------
    offset : int, optional
        The amount the spectrum is shifted, mainly to test the effect of the
        cross correlation.
        (The default is 0, I.E. no offset introduced)
    **kwargs : dict
        keyword arguments.
        Allows for passing unpacked dictionary to the class constructor.
        ftcs : bool, optional
            Boolean whether to use Fourier Transform for cross correlation.

    See Also
    --------
    scipy:
        https://docs.scipy.org/doc/scipy/reference/generated/
        correlation:
            scipy.signal.correlate.html
    matplotlib custom style:
        https://matplotlib.org/stable/users/explain/customizing.html

    Notes
    -----
    Constants Imported (See utils.Constants):
        SAVE_CORR:
            The default save name for the correlation plot.
        OFFSET:
            The vertical offset of spectra in the output plot.
    """

    # MARK: Correlate init
    def __init__(
            self,
            data_dir: Path,
            filenames: list[str],
            beams: str = "OE",
            split_ccd: bool = True,
            cont_ord: int = 11,
            plot: bool = False,
            offset: int = 0,
            save_prefix: Path | None = None,
            **kwargs,
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="ecwmxgbp",
            ext="fits",
        )
        self._beams = None
        self.beams = beams

        self.ccds = 1
        if split_ccd:
            # with pyfits.open(self.fits_list[0]) as hdu:
            #     BPM == 2 near center of CCD (extract version != *_sc)
            #     self.ccds = sum(hdu["BPM"].data.sum(axis=1)[0] == 2)
            self.ccds = 3

        self.cont_ord = cont_ord
        self.can_plot = plot

        self.offset = offset
        if offset != 0:
            # logging.warning("'offset' is only for testing.")
            # # Add an offset to the spectra to test cross correlation
            # self.spec1 = np.insert(
            #     self.spec1, [0] * offset, self.spec1[:, :offset], axis=-1
            # )[:, : self.spec1.shape[-1]]

            err_msg = "Offset deprecated after testing finalized."
            logging.error(err_msg)
            raise DeprecationWarning(err_msg)

        self.save_prefix = save_prefix
        # Handle directory save name
        if self.save_prefix and self.save_prefix.is_dir():
            self.save_prefix /= SAVE_CORR
            logging.warning((
                f"Correlation save name resolves to a directory. "
                f"Saving under {self.save_prefix}"
            ))

        self.wav_unit = "\\AA"
        self.wav_cdelt = 1

        self.alt = self.ftcs if kwargs.get("ftcs") else None

        logging.debug(f"__init__ - \n{repr(self)}")

        return

    # MARK: Correlate repr
    def __repr__(self) -> str:
        template = (
            "CrossCorrelate(\n"
            f"\tdata_dir={self.data_dir},\n"
            f"\tfits_list=[\n\t\t{"\n\t\t".join(
                map(str, self.fits_list)
            )}\n\t],\n"
            f"\tbeams={self._beams},\n"
            f"\tsplit_ccd={self.ccds},\n"
            f"\tcont_ord={self.cont_ord},\n"
            f"\tplot={self.can_plot},\n"
            f"\toffset={self.offset},\n"
            f"\tsave_prefix={self.save_prefix},\n"
            f"\twav_unit={self.wav_unit},\n"
            f"\twav_cdelt={self.wav_cdelt},\n"
            f"\talt={self.alt},\n"
            ")\n"
        )

        return template

    # MARK: Beams property
    @property
    def beams(self) -> str:
        return self._beams

    @beams.setter
    def beams(self, mode: str) -> None:
        if mode not in ['O', 'E', 'OE']:
            err_msg = f"Correlation mode '{mode}' not recognized."
            logging.error(err_msg)
            raise ValueError(err_msg)

        self._beams = mode

        return

    # MARK: Load file
    def load_file(
            self,
            filename: Path
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Load the data from a FITS file.

        Parameters
        ----------
        filename : Path
            The name of the FITS file to load.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, np.ndarray]
            The spectrum, wavelength, and bad pixel mask.

        """

        # Open HDU
        with pyfits.open(filename) as hdul:
            spec = hdul["SCI"].data.sum(axis=1)
            wav = (
                    np.arange(spec.shape[-1])
                    * hdul["SCI"].header["CDELT1"]
                    + hdul["SCI"].header["CRVAL1"]
            )
            wav = np.array((wav, wav))
            bpm = hdul["BPM"].data.sum(axis=1)

            self.wav_cdelt = float(hdul["SCI"].header["CDELT1"])

            if hdul["SCI"].header["CTYPE1"] != 'Angstroms':
                self.wav_unit = hdul["SCI"].header["CTYPE1"]

        return spec, wav, bpm

    # MARK: Get bounds
    def get_bounds(self, bpm: np.ndarray, gap_length: int = 30) -> np.ndarray:
        """
        Find the bounds for a file based on the CCD count.

        Parameters
        ----------
        bpm : np.ndarray
            The bad pixel mask.
        gap_length : int, optional
            The minimum length of a gap to be considered a CCD region.
            (Defaults to 30)

        Returns
        -------
        np.ndarray
            The bounds for the CCD regions.

        """
        if not 0 <= gap_length <= 60:
            msg = f"An invalid gap length of {gap_length} was encountered."
            logging.error(msg)
            raise ValueError(msg)

        # Check if get_bounds is needed
        if self.ccds == 1:
            return np.array(
                [[(0, bpm[0].shape[-1])], [(0, bpm[1].shape[-1])]]
            ).astype(int)

        # bounds.shape -> (O|E, CCD's, low.|up. bound)
        bounds = np.zeros((2, self.ccds, 2))

        # Check if `BPM` contains any `2`'s
        if np.any(bpm == 2):
            for ext, ccd in iters.product(range(2), range(self.ccds)):
                mid = np.where(bpm[ext] == 2)[0][ccd]
                ccds = self.ccds * 2
                bounds[ext, ccd] = (
                    max(mid - bpm.shape[-1] // ccds, 0),
                    min(mid + bpm.shape[-1] // ccds, bpm.shape[-1])
                )

            return bounds.astype(int)

        for ext in range(len(self._beams)):

            # Find min and max vals
            min_val, max_val = 0, bpm[ext].shape[-1]

            while True:
                if bpm[ext][min_val] == 0:
                    break
                min_val += 1

            while True:
                if bpm[ext][max_val - 1] == 0:
                    break
                max_val -= 1

            # Find ranges of non zero values
            regions: list[np.ndarray] = np.split(
                np.where(bpm[ext, min_val: max_val] == 1)[0],
                np.where(
                    np.diff(np.where(bpm[ext, min_val: max_val] == 1)[0]) != 1
                )[0] + 1
            )

            # If less than 2 regions, raise error
            if len(regions) < 2:
                msg = "Less than 2 regions found in BPM. Returning bounds."
                logging.error(msg)
                raise ValueError(msg)

            # Find `regions` longer than `gap_length`
            regions = [
                region for region in regions
                if len(region) >= gap_length
            ]

            # Ensure 2 regions are found
            if len(regions) < 2:
                logging.debug(
                    "get_bounds - Less than 2 regions found in BPM." +
                    f"Calling get_bounds with gap_length = {gap_length - 10}"
                )
                return self.get_bounds(bpm, gap_length - 10)

            elif len(regions) > 2:
                logging.debug(
                    "get_bounds - More than 2 regions found in BPM. " +
                    f"Calling get_bounds with gap_length = {gap_length + 5}"
                )
                return self.get_bounds(bpm, gap_length - 10)

            # Ensure region order correct
            if regions[0][0] > regions[1][0]:
                regions = regions[::-1]

            # Assign bounds from regions
            bounds[ext] = np.array([
                (min_val, regions[0][0]),
                (regions[0][-1], regions[1][0]),
                (regions[1][-1], max_val),
            ])

            # Get lower and upper bound for each ccd, save to bounds
            # Lower -> min is zero, Upper -> max is bpm length

            # for ext in range(2):
            #     cedge = ccdcenter(bpm[ext])
            #     bounds[ext] = np.array(cedge)

        logging.debug(f"get_bounds - found bounds at \n{bounds.astype(int)}")

        return bounds.astype(int)

    # MARK: Remove Continua
    def remove_cont(
            self,
            spec: np.ndarray,
            wav: np.ndarray,
            bpm: np.ndarray,
            plot_cont: bool
    ) -> np.ndarray:
        """
        Remove the continuum from the data.

        Parameters
        ----------
        spec : np.ndarray
            The spectrum to remove the continuum from.
        wav : np.ndarray
            The wavelength of the spectrum.
        bpm : np.ndarray
            The bad pixel mask.
        plot_cont : bool
            Decides whether the continuum fitting should be plotted

        Returns
        -------
        spec : np.ndarray

        """
        # Mask out the bad pixels for fitting continua
        okwav = np.where(bpm != 1)

        # Define continua
        ctm = continuum(
            wav[okwav],
            spec[okwav],
            deg=self.cont_ord,
            plot=plot_cont,
        )

        # Normalise spectra
        spec /= chebyshev.chebval(wav, ctm)
        spec -= 1

        return spec

    # MARK: Correlate
    def correlate(
            self,
            filename1: Path,
            filename2: Path | None = None,
            alt: Callable = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[list], list[list]]:
        """
        Cross correlates the data.

        Parameters
        ----------
        filename1 : Path
            The name of the first FITS file to cross correlate.
        filename2 : Path, optional
            The name of the second FITS file to cross correlate.
            (Defaults to None)
        alt : Callable, optional
            An alternate method of cross correlating the data.
            (Defaults to None)

        Returns
        -------
        spec, wav, bpm, lagsdb, corrdb:
            tuple[np.ndarray, np.ndarray, np.ndarray, list[list], list[list]]

        """
        #  mode: OE -> 'O1' & 'E1', O -> 'O1' & 'O2', E -> 'E1' & 'E2'
        # Load data
        spec, wav, bpm = self.load_file(filename1)
        if filename2 and self._beams != 'OE':
            def unpack(exts, *args):
                return [arr[exts] for arr in args]

            if self._beams == 'O':
                spec[-1], wav[-1], bpm[-1] = unpack(
                    0, *self.load_file(filename2)
                )

            else:
                spec[0], wav[0], bpm[0] = spec[-1], wav[-1], bpm[-1]
                spec[-1], wav[-1], bpm[-1] = unpack(
                    -1, *self.load_file(filename2)
                )

        bounds = self.get_bounds(bpm)

        logging.debug(
            f"correlate - data shape:\n\tspec/wav/bpm: {spec.shape}"
        )

        corrdb = [[] for _ in range(self.ccds)]
        lagsdb = [[] for _ in range(self.ccds)]
        for ccd in range(self.ccds):
            sig = []
            for ext in range(2):
                lb, ub = bounds[ext, ccd]

                if self.cont_ord > 0:
                    spec[ext, lb:ub] = self.remove_cont(
                        spec[ext, lb:ub],
                        wav[ext, lb:ub],
                        bpm[ext, lb:ub],
                        self.can_plot
                    )

                # Invert BPM (and account for 2); zero bad pixels
                sig.append((
                    spec[ext, lb:ub]
                    * abs(bpm[ext, lb:ub].astype(np.int8) * -1 + 1)
                ))

            # Finally(!!!) cross correlate signals and scale max -> 1
            corrdb[ccd] = signal.correlate(*sig) if not alt else alt(*sig)
            corrdb[ccd] /= np.max(corrdb[ccd])
            # noinspection PyTypeChecker
            lagsdb[ccd] = signal.correlation_lags(
                sig[0].shape[-1],
                sig[1].shape[-1]
            ) * self.wav_cdelt

        return spec, wav, bpm, corrdb, lagsdb

    # MARK: ftcs alternate
    def ftcs(
            self,
            signal1: np.ndarray,
            signal2: np.ndarray
    ) -> np.ndarray:
        """
        Cross correlates the data using the Fourier Transform.

        Parameters
        ----------
        signal1 : np.ndarray
            The first signal to cross correlate.
        signal2 : np.ndarray
            The second signal to cross correlate.

        Returns
        -------
        np.ndarray
            The correlation data using the Fourier Transform.

        """
        logging.debug(
            f"ftcs - data shape:\n\tspec/wav/bpm: {signal1.shape}"
        )

        # Invert BPM (and account for 2); zero bad pixels
        ft_spec1 = np.fft.fft(signal1)
        ft_spec2 = np.fft.fft(signal2)

        if self.can_plot:
            plt.plot(ft_spec1)
            plt.plot(ft_spec2)
            plt.show()

        # Cross correlate signals
        # ft_spectrum1 * np.conj(ft_spectrum2)
        corr_entry = signal.correlate(ft_spec1, ft_spec2)

        return np.fft.ifft(corr_entry)

    # MARK: Plot
    def plot(self, spec, wav, bpm, corrdb, lagsdb) -> None:
        """
        Plot the data.

        Parameters
        ----------
        spec : np.ndarray
            The spectrum.
        wav : np.ndarray
            The wavelength.
        bpm : np.ndarray
            The bad pixel mask.
        corrdb : np.ndarray
            The cross correlation data.
        lagsdb : np.ndarray
            The `lags` data.

        Returns
        -------
        None

        """
        plt.style.use([
            files(STOPS.utils).joinpath('STOPS.mplstyle'),
            files(STOPS.utils).joinpath('STOPS_correlate.mplstyle'),
        ])
        bounds = self.get_bounds(bpm)

        fig, axs = plt.subplots(2, self.ccds, sharey="row")

        if self.ccds == 1:
            # Convert axs to a 2D array
            # noinspection PyTypeChecker
            axs: np.ndarray[matplotlib.axes.Axes] = np.swapaxes(
                np.atleast_2d(axs), 0, 1
            )

        # for ext, ccd in iters.product(range(2), range(self.ccds)):

        for ccd in range(self.ccds):
            axs[0, ccd].plot(
                lagsdb[ccd],
                corrdb[ccd] * 100,
                color='C4',
                label=f"max lag @ {lagsdb[ccd][corrdb[ccd].argmax()] - (bounds[1, ccd, 0] - bounds[0, ccd, 0])}",
            )

            for ext in range(2):
                lb, ub = bounds[ext, ccd]
                logging.debug(f"fl-{ext}: {wav[ext, lb]}:{wav[ext, ub - 1]}")

                axs[1, ccd].plot(
                    wav[ext, lb:ub],
                    spec[ext, lb:ub]
                    * abs(bpm[ext, lb:ub].astype(np.int8) * -1 + 1)
                    + OFFSET * ext,
                    label=(
                        f"${self._beams if self._beams != 'OE' else self._beams[ext]}"
                        f"_{ext + 1 if self._beams != 'OE' else 1}$"
                        f"{(' (+' + str(OFFSET * ext) + ')') if ext > 0 else ''}"
                    ),
                )

        axs[0, 0].set_ylabel("Normalised Correlation\n($\\%$)")
        for ax in axs[1:, 0]:
            ax.set_ylabel("Normalised Intensity\n(Counts)")
        xcol = int(self.ccds != 1)
        axs[0, xcol].set_xlabel(f"Signal Lag ({self.wav_unit})")
        axs[-1, xcol].set_xlabel(f"Wavelength ({self.wav_unit})")
        for ax in axs.flatten():
            leg = ax.legend()
            leg.set_draggable(True)

        plt.show()

        # Handle do not save
        if not self.save_prefix:
            return

        # Handle save
        fig.savefig(fname=self.save_prefix)

        return

    # MARK: Process all listed images
    def process(self) -> None:
        """
        Process the data.

        Returns
        -------
        None

        """
        if self._beams != 'OE' and len(self.fits_list) == 1:
            # change mode to OE with warning
            logging.warning((
                f"`{self._beams}` correlation not possible for "
                "a single file. correlation `mode` changed to 'OE'."
            ))
            self._beams = 'OE'

        # OE `mode` (same file, diff. ext.)
        if self._beams == 'OE':
            for fl in self.fits_list:
                logging.info(f"'OE' correlation of {fl}.")
                spec, wav, bpm, corr, lags = self.correlate(fl, alt=self.alt)
                self.plot(spec, wav, bpm, corr, lags)

            return

        # O|E `mode` (diff. files, same ext.)
        for fl1, fl2 in iters.combinations(self.fits_list, 2):
            logging.info(f"{self._beams} correlation of {fl1} vs {fl2}.")
            spec, wav, bpm, corr, lags = self.correlate(fl1, fl2, alt=self.alt)
            self.plot(spec, wav, bpm, corr, lags)

        return


# MARK: Main function
def main(argv) -> None:
    return


if __name__ == "__main__":
    main(sys.argv[1:])
