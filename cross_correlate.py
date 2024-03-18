#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__email__ = "justin.jb78+Masters@gmail.com"

import os
import sys
import logging
import itertools as iters

import numpy as np
from numpy.polynomial import chebyshev
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import signal

from utils import SharedUtils as su

# TODO@JustinotherGitter: Update correlate to use [filenames] instead of in1/in2
# TODO@JustinotherGitter: Update correlate to use relevant args:
    # filename <- in1/in2,
    # continuum_order <- cont,
    # continuum_plot <- cont_plot
# TODO@JustinotherGitter: Implement own logging in main()


class CrossCorrelate:
    """
    Cross correlate allows for comparing the extensions of multiple
    FITS files, or comparing the O and E beams of a single FITS file.

    Parameters
    ----------
    in1 : str
        The first ecwmxgbp*.fits file to be cross correlated
    in2 : str, optional
        The second ecwmxgbp*.fits file to be cross correlated.
        Cross correlation against the two extensions occurs if left empty
        (The default is None)
    split_ccd : bool, optional
        Decides whether the CCD regions should each be individually cross correlated.
        (The default is True, which splits the spectrum up into its seperate CCD regions)
    cont : int, optional
        The degree of a chebyshev to fit to the continuum.
        (The default is 11)
    cont_plot : bool, optional
        Decides whether or not the continuum fitting should be plotted
        (The default is False, so no continua plots are displayed)
    offset : int, optional
        The amount the spectrum is shifted, mainly to test the effect of the cross correlation
        (The default is 0, I.E. no offset introduced)
    save_name : str, optional
        The name or directory to save the figure produced to.
        "." saves a default name to the current working. A default name is also used when save_name is a directory.
        (The default is None, I.E. The figure is not saved, only displayed)

    Returns
    -------
    None

    Raises
    ------
    File not found Error
        Raised when input spectra or save_name directories are invalid


    Based on https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html#scipy.signal.correlate

    """

    def __init__(
        self,
        data_dir: os.PathLike,
        fits_list: list[os.PathLike],
        split_ccd: bool = True,
        cont: int = 11,
        cont_plot: bool = False,
        offset: int = 0,
        save_name: str = None,
        **kwargs
    ) -> None:
        # Defined when passed in
        self.data_dir = data_dir
        self.fits_list = fits_list
        self.split_ccd = split_ccd
        self.cont = cont
        self.cont_plot = cont_plot
        self.offset = offset
        self.save_name = save_name

        # Defined when processed
        self.spec = None # spec1, spec2
        self.wav = None # wav1, wav2
        self.bpm = None # bpm1, bpm2

        self.exts = 0
        self.ccds = 1
        self.bounds = None # bounds1, bounds2


        # self.invert = False
        # self.wavUnits = "Å"
        # self.wav1, self.spec1, self.bpm1 = self.checkLoad(in1)
        # self.wav2, self.spec2, self.bpm2 = self.checkLoad(in2, in1)

        # Move to process
        for target in self.fits_list:
            self.spec, self.wav, self.bpm = self.loadFile(target)

            self.exts = self.spec[0].shape[0]

            # Bounds shape [extensions, ccds, lower / upper bound]
            self.bounds = self.setBounds()

            if split_ccd:
                self.splitCCD()

        # self.bounds.append(np.array(
        #     [[[0, self.spec1[0].shape[-1]]], [[0, self.spec1[1].shape[-1]]]], dtype=int
        # ))
        # self.bounds.append(np.array(
        #     [[[0, self.spec2[0].shape[-1]]], [[0, self.spec2[1].shape[-1]]]], dtype=int
        # ))


        # self.exts = self.spec1.shape[0]
        # self.ccds = 1
        # Bounds shape [extensions, ccds, lower / upper bound]
        # self.bounds1 = np.array(
        #     [[[0, self.spec1[0].shape[-1]]], [[0, self.spec1[1].shape[-1]]]], dtype=int
        # )
        # self.bounds2 = np.array(
        #     [[[0, self.spec2[0].shape[-1]]], [[0, self.spec2[1].shape[-1]]]], dtype=int
        # )
        if split_ccd:
            self.splitCCD()

        self.cont = cont
        if cont > 0:
            self.rmvCont(cont_plot)

        # Add an offset to the spectra to test cross correlation
        self.spec1 = np.insert(
            self.spec1, [0] * offset, self.spec1[:, :offset], axis=-1
        )[:, : self.spec1.shape[-1]]

        self.corrdb = []
        self.lagsdb = []
        # self.correlate()

        self.save_name = save_name
        # self.checkPlot()

        return
    
    def loadFile(self, filename: os.PathLike) -> tuple[list, list, list]:
        spec, wav, bpm = None, None, None

        # Open HDU
        with pyfits.open(self.data_dir / self.filename) as hdu:
            #Load spec, wav, and bpm data - indexing [wav, intensity, beam]
            spec = hdu["SCI"].data.sum(axis=1)
            wav  = (
                np.arange(spec.shape[-1]) * hdu["SCI"].header["CDELT1"] + hdu["SCI"].header["CRVAL1"]
            )
            bpm = hdu["BPM"].data.sum(axis=1)

            # Check wavelength units unchanged
            if "Angstroms" not in hdu["SCI"].header["CTYPE1"]:
                self.wavUnits = hdu["SCI"].header["CTYPE1"]

        # TODO@JustinotherGitter: Recheck return of o and e beams.
        return ([spec, spec[::-1]], [wav, wav], [bpm, bpm[::-1]])

    def setBounds(self) -> list[np.ndarray, np.ndarray]:
        bounds = []
        bounds.append(np.array(
            [[[0, self.spec[0].shape[-1]]], [[0, self.spec[1].shape[-1]]]], dtype=int
        ))
        bounds.append(np.array(
            [[[0, self.spec[0].shape[-1]]], [[0, self.spec[1].shape[-1]]]], dtype=int
        ))

        return bounds
    
    # def checkLoad(self, path1: str, path2: str = None) -> np.ndarray:

    #     # If the first path is invalid
    #     if (path1 == None) or (not os.path.isfile(os.path.expanduser(path1))):
    #         # And the second path is not defined, raise an error
    #         if path2 == None:
    #             raise FileNotFoundError(f"{path1} is invalid")

    #         # Use the second path but swap the O and E beams
    #         path1 = path2
    #         self.invert = True

    #     # Load data
    #     with pyfits.open(os.path.expanduser(path1)) as hdu:
    #         spec = hdu["SCI"].data.sum(axis=1)
    #         wav = (
    #             np.arange(spec.shape[-1]) * hdu["SCI"].header["CDELT1"]
    #             + hdu["SCI"].header["CRVAL1"]
    #         )
    #         bpm = hdu["BPM"].data.sum(axis=1)

    #         if "Angstroms" not in hdu["SCI"].header["CTYPE1"]:
    #             self.wavUnits = hdu["SCI"].header["CTYPE1"]

    #     # Return data and implement swap if necessary
    #     return (wav, spec[::-1], bpm[::-1]) if self.invert else (wav, spec, bpm)

    def splitCCD(self) -> None:
        # Assumed BPM has a value of 2 near the center of each CCD (i.e. sum(bpm == 2) = count(ccd))
        self.ccds = sum(self.bpm1[0] == 2)
        
        # update bounds to reflect ccds
        self.bounds1 = np.zeros([self.exts, self.ccds, 2], dtype=int)
        self.bounds2 = np.zeros([self.exts, self.ccds, 2], dtype=int)

        # Get lower and upper bound for each ccd, save to bounds
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            mid1 = np.where(self.bpm1[ext] == 2)[0][ccd]
            mid2 = np.where(self.bpm2[ext] == 2)[0][ccd]

            # Lower bound, min non-zero
            lowb1 = max(mid1 - self.bpm1.shape[-1] // (self.ccds * 2), 0)
            uppb1 = min(
                mid1 + self.bpm1.shape[-1] // (self.ccds * 2), self.bpm1.shape[-1]
            )

            # Upper bound, max bpm length
            lowb2 = max(mid2 - self.bpm2.shape[-1] // (self.ccds * 2), 0)
            uppb2 = min(
                mid2 + self.bpm2.shape[-1] // (self.ccds * 2), self.bpm2.shape[-1]
            )

            self.bounds1[ext, ccd] = (lowb1, uppb1)
            self.bounds2[ext, ccd] = (lowb2, uppb2)

    def rmvCont(self, plotCont) -> None:
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Get the range for current extension, ccd combination
            ccdBound1 = range(*self.bounds1[ext][ccd])
            ccdBound2 = range(*self.bounds2[ext][ccd])

            # Mask out the bad pixels for fitting continua
            okwav1 = np.where(self.bpm1[ext][ccdBound1] != 1)
            okwav2 = np.where(self.bpm2[ext][ccdBound2] != 1)

            # Define continua
            ctm1 = su.continuum(
                self.wav1[ccdBound1][okwav1],
                self.spec1[ext][ccdBound1][okwav1],
                deg=self.cont,
                plot=plotCont,
            )
            ctm2 = su.continuum(
                self.wav2[ccdBound2][okwav2],
                self.spec2[ext][ccdBound2][okwav2],
                deg=self.cont,
                plot=plotCont,
            )

            # Normalise spectra
            self.spec1[ext][ccdBound1] /= chebyshev.chebval(self.wav1[ccdBound1], ctm1)
            self.spec1[ext][ccdBound1] -= 1

            self.spec2[ext][ccdBound2] /= chebyshev.chebval(self.wav2[ccdBound2], ctm2)
            self.spec2[ext][ccdBound2] -= 1

        return

    def correlate(self) -> None:
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Get the range for current extension, ccd combination
            ccdBound1 = range(*self.bounds1[ext][ccd])
            ccdBound2 = range(*self.bounds2[ext][ccd])

            # Add rows/cols for correlation and lags data
            if len(self.corrdb) <= ext:
                self.corrdb.append([])
                self.lagsdb.append([])
            if len(self.corrdb[ext]) <= ccd:
                self.corrdb[ext].append([])
                self.lagsdb[ext].append([])

            # Invert BPM (and account for 2 in BPM) to zero bad pixels
            sig1 = self.spec1[ext][ccdBound1] * abs(self.bpm1[ext][ccdBound1] * -1 + 1)
            sig2 = self.spec2[ext][ccdBound2] * abs(self.bpm2[ext][ccdBound2] * -1 + 1)
            print(self.wav1[ccdBound1][0], self.wav2[ccdBound2][0])

            # Finally(!!!) cross correlate signals
            corr = signal.correlate(sig1, sig2)
            corr /= np.max(corr)  # Scales array so that the maximum correlation is at 1
            lags = signal.correlation_lags(sig1.shape[-1], sig2.shape[-1])

            self.corrdb[ext][ccd] = corr
            self.lagsdb[ext][ccd] = lags

        return

    def checkPlot(self, default_name: str = "OEcorr.pdf") -> None:
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
        if self.save_name == None:
            return

        # Handle lazy save_name
        if self.save_name == ".":
            self.save_name = os.getcwd()

        # Handle save name directory, use a default name (overwrite with warning)
        if self.save_name[-1] == "/" or os.path.isdir(self.save_name):
            self.save_name += default_name
            print(
                f"Save name is a directory. Saving cross correlation results as {default_name}"
            )

        # Check save location valid
        save_dir = os.path.expanduser("/".join(self.save_name.split("/")[:-1]))
        if not os.path.isdir(save_dir):
            raise FileNotFoundError(f"The path ({save_dir}) does not exist")

        # Save
        if self.save_name != None:
            fig.savefig(fname=self.save_name)

        return

    def process(self) -> None:
        for target in self.fits_list:
            logging.debug(f"Processing {target}")
            # self.correlate(target)
            # self.checkPlot()

        return


def main(argv) -> None: # TODO@JustinotherGitter: Handle cross_correlate.py called directly
    return

if __name__ == "__main__":
    main(sys.argv[1:])
