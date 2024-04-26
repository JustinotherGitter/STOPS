"""Module for cross correlating polarization beams."""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __init__ import __author__, __email__, __version__

# MARK: Imports
import os
import sys
import logging
import itertools as iters
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import signal

from utils.SharedUtils import find_files, continuum
from utils.Constants import SAVE_CORR

OFFSET = 0.3

# MARK: Correlate class
class CrossCorrelate:
    """
    Cross correlate allows for comparing the extensions of multiple
    FITS files, or comparing the O and E beams of a single FITS file.

    Parameters
    ----------
    data_dir : str
        The path to the data to be cross correlated
    filenames : list[str]
        The ecwmxgbp*.fits files to be cross correlated.
        If only one filename is defined, correlation is done against the two polarization beams.
    split_ccd : bool, optional
        Decides whether the CCD regions should each be individually cross correlated.
        (The default is True, which splits the spectrum up into its seperate CCD regions)
    cont_ord : int, optional
        The degree of a chebyshev to fit to the continuum.
        (The default is 11)
    cont_plot : bool, optional
        Decides whether or not the continuum fitting should be plotted
        (The default is False, so no continua plots are displayed)
    save_name : str, optional
        The name or directory to save the figure produced to.
        "." saves a default name to the current working. A default name is also used when save_name is a directory.
        (The default is None, I.E. The figure is not saved, only displayed)

    Attributes
    ----------
    
    
    Methods
    -------

    
    Other Parameters
    ----------------
    offset : int, optional
        The amount the spectrum is shifted, mainly to test the effect of the cross correlation
        (The default is 0, I.E. no offset introduced)
    **kwargs : dict
        keyword arguments. Allows for passing unpacked dictionary to the class constructor.

    See Also
    --------
    scipy.signal.correlate
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html#scipy.signal.correlate

    """

    def __init__(
        self,
        data_dir: Path,
        filenames: list[str],
        beams: str = "OE",
        split_ccd: bool = True,
        cont_ord: int = 11,
        cont_plot: bool = False,
        offset: int = 0,
        save_prefix: Path | None = None,
        **kwargs
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="ecwmxgbp",
            ext="fits",
        )
        self.beams = beams
        self.ccds = 1
        if split_ccd:
            # BPM == 2 near center of CCD if CCD count varies
            with pyfits.open(self.fits_list[0]) as hdu:
                self.ccds = sum(hdu["BPM"].data.sum(axis=1)[0] == 2)

        self.cont_ord = cont_ord
        self.cont_plot = cont_plot
        self.offset = offset
        if offset != 0:
            logging.warning("'offset' is only for testing.")
        # # Add an offset to the spectra to test cross correlation
        # self.spec1 = np.insert(
        #     self.spec1, [0] * offset, self.spec1[:, :offset], axis=-1
        # )[:, : self.spec1.shape[-1]]

        self.save_name = save_prefix
        # Handle directory save name
        if self.save_name and self.save_name.is_dir():
            self.save_name /= SAVE_CORR
            logging.warning((
                f"Correlation save name resolves to a directory. "
                f"Saving under {self.save_name}"
                ))

        self.unit_wav = "$\AA$"

        logging.debug("__init__ - \n", self.__dict__)
        return
    
    def load_file(
        self,
        filename: Path
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        spec, wav, bpm = None, None, None

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

            if hdul["SCI"].header["CTYPE1"] != 'Angstroms':
                self.wav_units = hdul["SCI"].header["CTYPE1"]

        return spec, wav, bpm

    def get_bounds(self, bpm: np.ndarray) -> np.ndarray:
        """Find the bounds for a file based on the CCD count."""
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
            CCDs = self.ccds * 2
            bounds[ext, ccd] = (
                max(mid - bpm.shape[-1] // CCDs, 0),
                min(mid + bpm.shape[-1] // CCDs, bpm.shape[-1])
            )

        return bounds.astype(int)
    
    def split_ccds(
        self,
        spec: np.ndarray,
        wav: np.ndarray,
        bpm: np.ndarray
    ) -> list[list]:
        out = []
        bounds = self.find_bounds(bpm)
        for arr in [spec, wav, bpm]:
            arr = [[arr[ext, bounds[ext, ccd, 0]:bounds[ext, ccd, 1]] for ccd in range(self.ccds)] for ext in range(2)]
            out.append(arr)

        return out

    def remove_cont(
        self,
        spec: list,
        wav: list,
        bpm: list,
        plotCont: bool
    ) -> None:
        # Mask out the bad pixels for fitting continua
        okwav = np.where(bpm != 1)

        # Define continua
        ctm = continuum(
            wav[okwav],
            spec[okwav],
            deg=self.cont_ord,
            plot=plotCont,
        )

        # Normalise spectra
        spec /= chebyshev.chebval(wav, ctm)
        spec -= 1

        return spec

    def correlate(
        self,
        filename1: Path,
        filename2: Path | None = None
    ) -> None:
        #  mode: OE -> 'O1' & 'E1', O -> 'O1' & 'O2', E -> 'E1' & 'E2'
        # Load data
        spec, wav, bpm = self.load_file(filename1)
        if filename2 and self.beams != 'OE':
            unpack = lambda ext, *args: [arr[ext] for arr in args]

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
                        self.cont_plot
                    )

                # Invert BPM (and account for 2); zero bad pixels
                sig.append((
                    spec[ext, lb:ub]
                    * abs(bpm[ext, lb:ub] * -1 + 1)
                ))

            # Finally(!!!) cross correlate signals and scale max -> 1
            corrdb[ccd] = signal.correlate(*sig)
            corrdb[ccd] /= np.max(corrdb[ccd])
            lagsdb[ccd] = signal.correlation_lags(
                sig[0].shape[-1],
                sig[1].shape[-1]
            )

        return (spec, wav, bpm), (corrdb, lagsdb)
    
    def FTCS(
        self,
        filename1: Path,
        filename2: Path | None = None
    ) -> None:
        #  mode: OE -> 'O1' & 'E1', O -> 'O1' & 'O2', E -> 'E1' & 'E2'
        # Load data
        spec, wav, bpm = self.load_file(filename1)
        if filename2 and self.beams != 'OE':
            unpack = lambda ext, *args: [arr[ext] for arr in args]

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
            f"FTSC - data shape:\n\tspec/wav/bpm: {spec.shape}"
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
                        self.cont_plot
                    )

                # Invert BPM (and account for 2); zero bad pixels
                ft_spec = np.fft.fft(spec[ext, lb:ub] * abs(bpm[ext, lb:ub] * -1 + 1))
                sig.append(ft_spec)

            # Finally(!!!) cross correlate signals and scale max -> 1
            corrdb[ccd] = np.fft.ifft(signal.correlate(*sig)) # ft_spectrum1 * np.conj(ft_spectrum2)
            corrdb[ccd] /= np.max(corrdb[ccd])
            lagsdb[ccd] = signal.correlation_lags(
                sig[0].shape[-1],
                sig[1].shape[-1]
            )

        return (spec, wav, bpm), (corrdb, lagsdb)

    def plot(self, spec, wav, bpm, corrdb, lagsdb) -> None:
        plt.style.use(Path(__file__).parent.resolve() / 'utils/STOPS.mplstyle')
        bounds = self.get_bounds(bpm)

        fig, axs = plt.subplots(2, self.ccds, sharey="row")

        if self.ccds == 1:
            # Convert axs to a 2D array
            axs = np.swapaxes(np.atleast_2d(axs), 0, 1)

        # for ext, ccd in iters.product(range(2), range(self.ccds)):

        for ccd in range(self.ccds):
            axs[0, ccd].plot(
                lagsdb[ccd],
                corrdb[ccd] * 100,
                color='C4',
                label=f"max lag @ {lagsdb[ccd][corrdb[ccd].argmax()]}",
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
        for ax in axs[0, :]:
            ax.set_xlabel("Signal Lag")
        for ax in axs[1:, 0]:
            ax.set_ylabel(f"Norm. Intensity\n(Counts)")
        for ax in axs[-1, :]:
            ax.set_xlabel(f"Wavelength ({self.unit_wav})")
        for ax in axs.flatten():
            ax.legend()

        # plt.tight_layout()
        plt.show()

        # Handle do not save
        if not self.save_name:
            return

        # Handle save
        fig.savefig(fname=self.save_name)

        return

    def process(self) -> None:
        if self.beams not in ['O', 'E', 'OE']:
            errMsg = f"Correlation mode '{self.beams}' not recognized."
            logging.error(errMsg)
            raise ValueError(errMsg)
        
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
                (spec, wav, bpm), (corr, lags) = self.correlate(fl)#self.FTCS(fl)
                self.plot(spec, wav, bpm, corr, lags)

            return
        
        # O|E `mode` (diff. files, same ext.)
        for fl1, fl2 in iters.combinations(self.fits_list, 2):
            logging.info(f"{self.beams} correlation of {fl1} vs {fl2}.")
            (spec, wav, bpm), (corr, lags) = self.correlate(fl1, fl2)#self.FTCS(fl1, fl2)
            self.plot(spec, wav, bpm, corr, lags)

        return


def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])
