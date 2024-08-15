#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module for cross correlating polarization beams."""

# MARK: Imports
import sys
import logging
import itertools as iters
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
import matplotlib.axes
from astropy.io import fits as pyfits
from scipy import signal

from STOPS.utils.SharedUtils import find_files, continuum
from STOPS.utils.Constants import SAVE_CORR, OFFSET

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.INFO)


# MARK: Correlate class
class CrossCorrelate:

    #----------corr0----------

    """
    Cross correlate allows for comparing the extensions of multiple
    FITS files, or comparing the O and E beams of a single FITS file.

    Parameters
    ----------
    data_dir : str | Path
        The path to the data to be cross correlated
    filenames : list[str]
        The ecwmxgbp*.fits files to be cross correlated.
        If only one filename is defined, correlation is done against the two polarization beams.
    split_ccd : bool, optional
        Decides whether the CCD regions should each be individually cross correlated.
        (The default is True, which splits the spectrum up into its separate CCD regions)
    cont_ord : int, optional
        The degree of a chebyshev to fit to the continuum.
        (The default is 11)
    plot : bool, optional
        Decides whether the continuum fitting should be plotted
        (The default is False, so no continua plots are displayed)
    save_prefix : str, optional
        The name or directory to save the figure produced to.
        "." saves a default name to the current working. A default name is also used when save_prefix is a directory.
        (The default is None, I.E. The figure is not saved, only displayed)

    Attributes
    ----------
    data_dir
    fits_list
    beams : str
        The mode of correlation.
        'OE' for same file, and 'O' or 'E' for different files but same extension.
    ccds : int
        The number of CCD's in the data. Used to split the CCD's if split_ccd is True.
    cont_ord : int
        The degree of the chebyshev to fit to the continuum.
    can_plot : bool
        Decides whether the continuum fitting should be plotted
    offset : int
        The amount the spectrum is shifted, mainly to test the effect of the cross correlation
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
        The amount the spectrum is shifted, mainly to test the effect of the cross correlation
        (The default is 0, I.E. no offset introduced)
    **kwargs : dict
        keyword arguments. Allows for passing unpacked dictionary to the class constructor.
        ftcs : bool, optional
            Decides whether the Fourier Transform should be used for cross correlation.

    See Also
    --------
    scipy
        https://docs.scipy.org/doc/scipy/reference/generated/
        correlation: scipy.signal.correlate.html

        
    Notes
    -----
    Constants Imported (See utils.Constants):
        SAVE_CORR
    
    """

    #----------corr1----------

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
            # BPM == 2 near center of CCD if CCD count varies
            with pyfits.open(self.fits_list[0]) as hdu:
                self.ccds = sum(hdu["BPM"].data.sum(axis=1)[0] == 2)

        self.cont_ord = cont_ord
        self.can_plot = plot
        self.offset = offset
        if offset != 0:
            logging.warning("'offset' is only for testing.")

            err_msg = "Offset removed after finalizing testing."
            logging.error(err_msg)
            raise ValueError(err_msg)
            # # Add an offset to the spectra to test cross correlation
            # self.spec1 = np.insert(
            #     self.spec1, [0] * offset, self.spec1[:, :offset], axis=-1
            # )[:, : self.spec1.shape[-1]]

        self.save_prefix = save_prefix
        # Handle directory save name
        if self.save_prefix and self.save_prefix.is_dir():
            self.save_prefix /= SAVE_CORR
            logging.warning((
                f"Correlation save name resolves to a directory. "
                f"Saving under {self.save_prefix}"
            ))

        self.wav_unit = "$\AA$"
        self.wav_cdelt = 1

        self.alt = self.ftcs if kwargs.get("ftcs") else None

        logging.debug("__init__ - \n", self.__dict__)
        return

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
    def get_bounds(self, bpm: np.ndarray) -> np.ndarray:
        """
        Find the bounds for a file based on the CCD count.
        
        Parameters
        ----------
        bpm : np.ndarray
            The bad pixel mask.
        
        Returns
        -------
        np.ndarray
            The bounds for the CCD regions.
        
        """
        # bounds.shape -> (O|E, CCD's, low.|up. bound)
        if self.ccds == 1:
            return np.array(
                [[(0, bpm[0].shape[-1])], [(0, bpm[1].shape[-1])]]
            ).astype(int)

        bounds = np.zeros((2, self.ccds, 2))

        # Get lower and upper bound for each ccd, save to bounds
        # Lower -> min is zero, Upper -> max is bpm length
        for ext, ccd in iters.product(range(2), range(self.ccds)):
            mid = np.where(bpm[ext] == 2)[0][ccd]
            ccds = self.ccds * 2
            bounds[ext, ccd] = (
                max(mid - bpm.shape[-1] // ccds, 0),
                min(mid + bpm.shape[-1] // ccds, bpm.shape[-1])
            )

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
        spec, wav, bpm, lagsdb, corrdb: tuple[np.ndarray, np.ndarray, np.ndarray, list[list], list[list]]

        """
        #  mode: OE -> 'O1' & 'E1', O -> 'O1' & 'O2', E -> 'E1' & 'E2'
        # Load data
        spec, wav, bpm = self.load_file(filename1)
        if filename2 and self.beams != 'OE':
            def unpack(exts, *args):
                return [arr[exts] for arr in args]

            if self.beams == 'O':
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
                        * abs(bpm[ext, lb:ub] * -1 + 1)
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
            Path(__file__).parent.resolve() / 'utils/STOPS.mplstyle',
            Path(__file__).parent.resolve() / 'utils/STOPS_correlate.mplstyle'
        ])
        bounds = self.get_bounds(bpm)

        fig, axs = plt.subplots(2, self.ccds, sharey="row")

        if self.ccds == 1:
            # Convert axs to a 2D array
            # noinspection PyTypeChecker
            axs: np.ndarray[matplotlib.axes.Axes] = np.swapaxes(np.atleast_2d(axs), 0, 1)

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
                    spec[ext, lb:ub] * abs(bpm[ext, lb:ub] * -1 + 1) + OFFSET * ext,
                    label=(
                        f"${self.beams if self.beams != 'OE' else self.beams[ext]}"
                        f"_{ext + 1 if self.beams != 'OE' else 1}$"
                        f"{(' (+' + str(OFFSET * ext) + ')') if ext > 0 else ''}"
                    ),
                )

        axs[0, 0].set_ylabel("Normalised Correlation\n(\%)")
        for ax in axs[1:, 0]:
            ax.set_ylabel("Normalised Intensity\n(Counts)")
        xcol = int(self.ccds != 1)
        axs[0, xcol].set_xlabel(f"Signal Lag ({self.wav_unit})")
        axs[-1, xcol].set_xlabel(f"Wavelength ({self.wav_unit})")
        for ax in axs.flatten():
            leg = ax.legend()
            leg.set_draggable(True)

        # plt.tight_layout()
        # fig1 = plt.gcf()
        # DPI = fig1.get_dpi()
        # fig1.set_size_inches(700.0/float(DPI), 250.0/float(DPI))
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
        if self.beams != 'OE' and len(self.fits_list) == 1:
            # change mode to OE with warning
            logging.warning((
                f"`{self.beams}` correlation not possible for "
                "a single file. correlation `mode` changed to 'OE'."
            ))
            self.beams = 'OE'

        # OE `mode` (same file, diff. ext.)
        if self.beams == 'OE':
            for fl in self.fits_list:
                logging.info(f"'OE' correlation of {fl}.")
                spec, wav, bpm, corr, lags = self.correlate(fl, alt=self.alt)
                self.plot(spec, wav, bpm, corr, lags)

            return

        # O|E `mode` (diff. files, same ext.)
        for fl1, fl2 in iters.combinations(self.fits_list, 2):
            logging.info(f"{self.beams} correlation of {fl1} vs {fl2}.")
            spec, wav, bpm, corr, lags = self.correlate(fl1, fl2, alt=self.alt)
            self.plot(spec, wav, bpm, corr, lags)

        return


# MARK: Main function
def main(argv) -> None:
    return


if __name__ == "__main__":
    main(sys.argv[1:])
