"""Module for analyzing the sky lines of a wavelength calibrated image."""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __init__ import __author__, __email__, __version__

# MARK: Imports
import os
import sys
import logging
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import signal

from utils.SharedUtils import find_files
from utils.Constants import SAVE_SKY

# plt.rcParams['figure.figsize'] = (20, 4)
plt.rcParams['image.origin'] = 'lower'

# MARK: Skylines Class
class Skylines:
    """
    Class representing the Skylines object.

    Parameters
    ----------
    data_dir : Path
        The directory containing the data files.
    filenames : list[str]
        The list of filenames to be processed.
    beam : str, optional
        The beam mode, by default "OE".
    cont_plot : bool, optional
        Flag indicating whether to plot the continuum, by default False.
    save_prefix : Path | None, optional
        The prefix for saving the data, by default None.
    **kwargs
        Additional keyword arguments.

    Attributes
    ----------
    data_dir : Path
        The directory containing the data files.
    fits_list : list[str]
        The list of fits file paths.
    beam : str
        The beam mode.
    cont_plot : bool
        Flag indicating whether to plot the continuum.
    save_prefix : Path | None
        The prefix for saving the data.
    wav_unit : str
        The unit of wavelength.
    rawWav : np.ndarray
        The raw wavelength data.
    rawSpec : np.ndarray
        The raw spectral data.
    rawBpm : np.ndarray
        The raw bad pixel mask data.
    corrWav : np.ndarray
        The corrected wavelength data.
    corrSpec : np.ndarray
        The corrected spectral data.
    spec : np.ndarray
        The median spectrum.
    normSpec : np.ndarray
        The normalized spectrum.

    Methods
    -------
    checkLoad(self, path1: str) -> np.ndarray:
        Checks and loads the data from the given path.
    transform(self, wav_sol: np.ndarray, spec: np.ndarray) -> np.ndarray:
        Transforms the input wavelength and spectral data based on the given wavelength solution.
    rmvCont(self) -> np.ndarray:
        Removes the continuum from the spectrum.
    skylines(self) -> None:
        Placeholder method for processing skylines.
    process(self) -> None:
        Placeholder method for processing the data.
    """

    def __init__(
        self,
        data_dir: Path,
        filenames : list[str],
        beam: str = "OE",
        cont_plot: bool = False,
        save_prefix: Path | None = None,
        **kwargs,
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="t", # t[o|e]beam
            ext="fits",
        )
        self._beam = None
        self.beam = beam

        self.cont_plot = cont_plot
        self.save_prefix = save_prefix
        # Handle directory save name
        if self.save_prefix and self.save_prefix.is_dir():
            self.save_prefix /= SAVE_SKY
            logging.warning((
                f"Correlation save name resolves to a directory. "
                f"Saving under {self.save_prefix}"
                ))

        self.wav_unit = "$\AA$"

        self.rawWav, self.rawSpec, self.rawBpm = self.checkLoad(
            self.fits_list
        )
        self.corrWav, self.corrSpec = self.transform(
            self.rawWav,
            self.rawSpec
        )
        self.spec = np.median(self.corrSpec, axis=1)
        self.normSpec = self.rmvCont(self.spec)

        logging.debug(self.__dict__)
        return
    
    # MARK: Beams property
    @property
    def beams(self) -> str:
        return self._beams
    
    @beams.setter
    def beams(self, mode: str) -> None:
        if mode not in ['O', 'E', 'OE']:
            errMsg = f"Correlation mode '{mode}' not recognized."
            logging.error(errMsg)
            raise ValueError(errMsg)
        
        self._beams = mode

        return
    
    # MARK: Load data
    def load_file(self, filename: Path) -> np.ndarray:
        """
        Loads the data from the given file.

        Parameters
        ----------
        filename : Path
            The path to the file to be loaded.

        Returns
        -------
        np.ndarray
            The data from the file.
        
        """
        with pyfits.open(filename) as hdul:
            spec2D = hdul["SCI"].data
            wav2D = hdul["WAV"].data
            bpm2D = hdul["BPM"].data

        # Return data
        return spec2D, wav2D, bpm2D

    # MARK: Transform spectra
    def transform(self, wav_sol: np.ndarray, spec: np.ndarray) -> np.ndarray:
        """
        Transforms the input wavelength and spectral data based on the given wavelength solution.

        Parameters
        ----------
        wav_sol : np.ndarray
            The wavelength solution.
        spec : np.ndarray
            The spectral data.

        Returns
        -------
        np.ndarray
            The transformed wavelength and spectral data.

        """
        # Create arrays to return
        cw = np.zeros_like(wav_sol)
        cs = np.zeros_like(wav_sol)

        exts = cw.shape[0]
        rows = cw.shape[1]

        for ext in range(exts):
            # Get middle row (to interpolate the rest of the rows to)
            avg_max = [np.where(
                spec[ext][:, col] == spec[ext][:, col].max()
            )[0][0] for col in range(spec[ext].shape[1])]
            avg_max = np.sum(avg_max) // spec[ext].shape[1]

            # Get wavelength values at row with most trace
            wav = wav_sol[ext][avg_max, :]

            # Correct extensions based on wavelength
            # Wavelength ext
            cw[ext][:, :] = wav

            for row in range(rows):
                # Spec extension
                cs[ext][row, :] = np.interp(
                    wav,
                    wav_sol[ext][row, :],
                    spec[ext][row, :]
                )

            if self.cont_plot:
                # Plot results
                fig, ax1 = plt.subplots(figsize=[20, 4])
                ax1.imshow(cs[ext],
                        vmax=cs[0].mean() + 2*cs[0].std(),
                        vmin=cs[0].mean() - 2*cs[0].std(),
                        origin='lower')
                print(f"Average continuum of {'E' if ext else 'O'} at {np.median(np.median(cs[ext], axis=0)):4.3f}")
                ax2 = ax1.twinx()
                ax2.hlines(np.median(np.median(cs[ext], axis=0)), 0, cs[ext].shape[-1], colors='black')
                ax2.plot(cs[ext].mean(axis=0), "k", label=f"mean {'E' if ext else 'O'}")
                ax2.plot(np.median(cs[ext], axis=0), "r", label=f"median {'E' if ext else 'O'}")
                ax2.legend()
                plt.show()

        return cw, cs

    def rmvCont(self):

        return self.spec / self.cont - 1

    def skylines(self,) -> None:
        pass

    # MARK: Process all listed images
    def process(self,) -> None:
        # Load data
        fls = []
        for fl in self.fits_list:
            # spec
            pass

        ## 2. Transform data, skip if filename starts with 't'
        ## 3. Remove continuum
        ## 4. Identify skylines
        ## 5. Fit skylines
        ## 6. Return results
        ## 7. Plot results
        ## 8. Save results

        return


# MARK: Main function
def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])


# Class flow
