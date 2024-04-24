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
        if not self.save_name and self.save_name.is_dir():
            self.save_name /= SAVE_CORR
            logging.warning((
                f"Correlation save name resolves to a directory. "
                f"Saving under {self.save_name}"
                ))

        self.unit_wav = "$\AA$"

        # self.bounds = None # bounds1, bounds2
        # # In process:
        # self.bounds = self.setBounds()

        self.corrdb = []
        self.lagsdb = []

        logging.debug(self.__dict__)
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
            bpm = hdul["BPM"].data.sum(axis=1)

            if hdul["SCI"].header["CTYPE1"] != 'Angstroms':
                self.wav_units = hdul["SCI"].header["CTYPE1"]

        return spec, wav, bpm

    def find_bounds(self, bpm: np.ndarray) -> np.ndarray:
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

    def remove_cont(
        self,
        spec: np.ndarray,
        wav: np.ndarray,
        bpm: np.ndarray,
        bounds: np.ndarray,
        plotCont: bool
    ) -> None:
        for ext, ccd in iters.product(range(2), range(self.ccds)):
            # Get the range for current extension, ccd combination
            ccdBound = range(*bounds[ext][ccd])

            # Mask out the bad pixels for fitting continua
            okwav = np.where(bpm[ext][ccdBound] != 1)

            # Define continua
            ctm = continuum(
                wav[ccdBound][okwav],
                spec[ext][ccdBound][okwav],
                deg=self.cont_ord,
                plot=plotCont,
            )

            # Normalise spectra
            spec[ext][ccdBound] /= chebyshev.chebval(wav[ccdBound], ctm)
            spec[ext][ccdBound] -= 1

        return spec

    def correlate(
        self,
        filename1: Path,
        filename2: Path | None = None
    ) -> None:
        # (spec / wav / bpm).shape -> (2, 'CCDs', 'val.')
        #  OE -> 'O1' and 'E1', O -> 'O1' and 'O2', E -> 'E1' and 'E2'
        spec, wav, bpm = self.load_file(filename1)
        if not filename2:
            unpack = lambda a, b, c, ext: (a[ext], b[ext], c[ext])

            if self.beams == 'O':
                spec[-1], wav[-1], bpm[-1] = unpack(
                    *self.load_file(filename2), 0
                )
            else:
                spec[0], wav[0], bpm[0] = spec[-1], wav[-1], bpm[-1]
                spec[-1], wav[-1], bpm[-1] = unpack(
                    *self.load_file(filename2), 1
                )
        
        # bounds.shape -> (2, 'CCDs', 2)
        bounds = self.find_bounds(bpm)

        if self.cont_ord > 0:
            spec = self.remove_cont(
                spec,
                wav,
                bpm,
                bounds,
                self.cont_plot
            )

        corrdb = np.zeros_like(spec[0])
        lagsdb = np.zeros_like(spec[0])
        # for ext, ccd in iters.product(2, range(self.ccds)):
        # Get the range for current ext./CCD combination
        curr_bounds = np.array([
            range(*bounds[0][ccd]),
            range(*bounds[1][ccd])
        ])

        # Invert BPM (and account for 2 in BPM) to zero bad pixels
        sig = spec[:, :, curr_bounds] * abs(bpm[:, :, curr_bounds] * -1 + 1)
        logging.debug(wav[0][curr_bounds[0]][0], wav[1][curr_bounds[1]][0])

        # Finally(!!!) cross correlate signals
        corr = signal.correlate(sig[0], sig[1])
        corr /= np.max(corr)  # Scales array so that the maximum correlation is at 1
        lags = signal.correlation_lags(sig[0].shape[-1], sig[1].shape[-1])

        corrdb = corr
        lagsdb = lags
        # end for loop
        return

    def checkPlot(self) -> None:
        # Plot
        fig, axs = plt.subplots(3, 3, sharey="row")

        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Add cross correlation to plots
            axs[0, ccd].plot(
                self.lagsdb[ext][ccd],
                self.corrdb[ext][ccd] * 100,
                label=f"Ext: {ext + 1}, max lag @ {self.lagsdb[ext][ccd][self.corrdb[ext][ccd].argmax()]}",
            )

            ccdBound1 = range(*self.bounds1[ext][ccd])
            ccdBound2 = range(*self.bounds2[ext][ccd])

            axs[ext + 1, ccd].plot(
                self.wav2[ccdBound2],
                self.spec2[ext][ccdBound2] * abs(self.bpm2[ext][ccdBound2] * -1 + 1),
                label="sig2",
            )
            axs[ext + 1, ccd].plot(
                self.wav1[ccdBound1],
                self.spec1[ext][ccdBound1] * abs(self.bpm1[ext][ccdBound1] * -1 + 1),
                label="sig1",
            )

        axs[0, 0].set_ylabel("Normalised Correlation\n(%)")
        for ax in axs[0, :]:
            ax.set_xlabel("Signal Lag")
        for i, ax in enumerate(axs[1:, 0]):
            ax.set_ylabel(f"Ext. {i + 1} - Norm. Intensity\n(Counts)")
        for ax in axs[-1, :]:
            ax.set_xlabel(f"Wavelength ({self.wavUnits})")
        for ax in axs.flatten():
            ax.legend()

        plt.tight_layout()
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
                logging.debug(f"'OE' correlation of {fl}.")
                self.correlate(fl)
                self.plot()
        
            return
        
        # O|E `mode` (diff. files, same ext.)
        for fl1, fl2 in iters.combinations(self.fits_list):
            logging.debug(f"{self.beams} correlation of {fl1} vs {fl2}.")
            self.correlate(fl1, fl2)
            self.plot()

        return


def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])
