"""Module for analyzing the sky lines of a wavelength calibrated image."""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __init__ import __author__, __email__, __version__

import os
import sys
import logging
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import signal

from utils.SharedUtils import find_files

# plt.rcParams['figure.figsize'] = (20, 4)
plt.rcParams['image.origin'] = 'lower'


class Skylines:
    def __init__(
        self,
        data_dir: Path,
        filenames : list[str],
        **kwargs,
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="t", # t[o|e]beam
            ext="fits",
        )
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

    def checkLoad(self, path1 : str) -> np.ndarray:
        # If the path is invalid
        if not os.path.isfile(os.path.expanduser(path1)):
            # Raise File Not Found Error
            raise FileNotFoundError(f"{path1} is invalid")

        # Load data
        with pyfits.open(os.path.expanduser(path1)) as hdul:
            spec2D = hdul["SCI"].data
            wav2D = hdul["WAV"].data
            bpm2D = hdul["BPM"].data

        # Return data
        return spec2D, wav2D, bpm2D

    def transform(wav_sol: np.ndarray, spec: np.ndarray) -> np.ndarray:
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

        return cw, cs

    def rmvCont(self):

        return self.spec / self.cont - 1

    def skylines(self,) -> None:
        pass


def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])
