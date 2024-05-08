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
from scipy import signal, stats, interpolate

from utils.SharedUtils import find_files, continuum
from utils.Constants import SAVE_SKY

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.INFO)
# plt.rcParams['figure.figsize'] = (20, 4)

# MARK: Skylines Class
class Skylines:

    #----------sky0----------

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
    plot : bool, optional
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
    can_plot : bool
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
    
    #----------sky1----------

    # MARK: Skylines init
    def __init__(
        self,
        data_dir: Path,
        filenames : list[str],
        beams: str = "OE",
        split_ccd: bool = False,
        cont_ord: int = 11,
        plot: bool = False,
        transform: bool = True,
        save_prefix: Path | None = None,
        **kwargs,
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="wmxgbp", # t[o|e]beam
            ext="fits",
        )
        self._beams = None
        self.beams = beams

        self.split_ccd = split_ccd
        self.cont_ord = cont_ord
        self.can_plot = plot
        self.must_transform = transform

        self.save_prefix = save_prefix
        # Handle directory save name
        if self.save_prefix and self.save_prefix.is_dir():
            self.save_prefix /= SAVE_SKY
            logging.warning((
                f"Skylines save name resolves to a directory. "
                f"Saving under {self.save_prefix}"
            ))

        self.wav_unit = "$\AA$"

        # self.rawWav, self.rawSpec, self.rawBpm = self.checkLoad(
        #     self.fits_list
        # )
        # self.corrWav, self.corrSpec = self.transform(
        #     self.rawWav,
        #     self.rawSpec
        # )
        # self.spec = np.median(self.corrSpec, axis=1)
        # self.normSpec = self.rmvCont(self.spec)

        logging.debug("__init__ - \n", self.__dict__)
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
    
    # MARK: Load Sky lines
    def load_sky_lines(self, filename: Path | None = None) -> np.ndarray:
        """
        Loads the sky lines from the given file.

        Parameters
        ----------
        filename : Path | None, optional
            The path to the file to be loaded.
            Defaults to loading the skylines from utils/sky.salt

        Returns
        -------
        sky_lines : np.ndarray
            The sky lines from the file.
        
        """
        if not filename:
            filename = Path(__file__).parent.resolve() / 'utils/sky.salt'

        dtype = [('wav', float), ('flux', float)]
        skylines = np.genfromtxt(filename, dtype=dtype, skip_header=3, skip_footer=1)

        if self.can_plot:
            plt.plot(skylines['wav'], skylines['flux'], 'x', label="Model peaks")
            plt.xlabel(f'Wavelength {self.wav_unit}')
            plt.ylabel('Relative intensities')
            plt.title('Known sky lines')
            plt.legend()
            plt.show()

        return skylines

    # MARK: Transform spectra
    @staticmethod
    def transform(wav_sol: np.ndarray, spec: np.ndarray, resPlot: bool = False) -> np.ndarray:
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
        wav, spec : np.ndarray
            The transformed wavelength and spectral data.

        """
        # Create arrays to return
        cw = np.zeros_like(wav_sol)
        cs = np.zeros_like(wav_sol)

        exts = cw.shape[0]
        rows = cw.shape[1]

        for ext in range(exts):
            # Get middle row (to interpolate the rest of the rows to)
            avg_max = [np.where(spec[ext, :, col] == spec[ext, :, col].max())[0][0] for col in range(spec[ext].shape[1])]
            avg_max = np.sum(avg_max) // spec[ext].shape[1]

            # Get wavelength values at row with most trace
            wav = wav_sol[ext, avg_max, :]

            # Correct extensions based on wavelength
            # Wavelength ext
            cw[ext, :, :] = wav

            # Spec ext
            # for row in range(rows):
            #     f_2d = interpolate.interp2d(
            #         wav_sol[ext, row],
            #         np.arange(rows),
            #         spec[ext],
            #     )
            #     cs[ext] = f_2d(cw[ext, row], np.arange(rows))
            for row in range(rows):
                cs[ext][row, :] = np.interp(
                    wav,
                    wav_sol[ext][row, :],
                    spec[ext][row, :]
                )

            # Plot results
            if resPlot:
                fig, ax1 = plt.subplots(figsize=[20, 4])
                ax1.imshow(cs[ext],
                        vmax=cs[0].mean() + 2*cs[0].std(),
                        vmin=cs[0].mean() - 2*cs[0].std()
                )
                print(f"Average continuum of {'E' if ext else 'O'} at {np.median(np.median(cs[ext], axis=0)):4.3f}")
                ax2 = ax1.twinx()
                ax2.hlines(np.median(np.median(cs[ext], axis=0)), 0, cs[ext].shape[-1], colors='black')
                ax2.plot(cs[ext].mean(axis=0), "k", label=f"mean {'E' if ext else 'O'}")
                ax2.plot(np.median(cs[ext], axis=0), "r", label=f"median {'E' if ext else 'O'}")
                ax2.legend()
                plt.show()

        return cw, cs

    # MARK: Remove continuum
    def remove_cont(self, spec: np.ndarray, wav: np.ndarray) -> np.ndarray:
        ctm = continuum(wav, spec, deg=self.cont_ord, plot=self.can_plot)

        return self.spec / ctm - 1

    # MARK: Skylines
    def skylines(self, filename) -> None:
        # raise error if arc image
        with pyfits.open(filename) as hdul:
            if hdul[0].header['OBSTYPE'] == 'ARC':
                logging.warning(f"ARC images, {filename}, contain no sky lines. File skipped.")
                return

            # Load data
            spec2D = hdul["SCI"].data
            wav2D = hdul["WAV"].data
            bpm2D = hdul["BPM"].data

        spec2D *= ~bpm2D

        logging.debug(f"skylines - {filename.name} - spec: {spec2D.shape}")

        # Mask trace
        # TODO@JustinotherGitter: Add trace masking if median is insufficient

        # Save initial feature widths
        # Mean to not filter out skewed features
        spec1D_init = np.mean(spec2D, axis=1)
        # Median to sort for most common wavelength
        wav1D_init = np.median(wav2D, axis=1)

        peaks = [[] for _ in range(2)]
        properties = [[] for _ in range(2)]

        if self.can_plot: fig, axs = plt.subplots(2, 1)
        for ext in range(2):
            peaks[ext], properties[ext] = signal.find_peaks(
                spec1D_init[ext],
                prominence=0.5 * np.std(spec1D_init[ext]),
                width=[1, 1000],
                rel_height=0.3
            )
            peak_width = [properties[ext]['widths'], properties[ext]['width_heights'], properties[ext]['left_ips'], properties[ext]['right_ips']]

            logging.debug(f"skylines - initial features {'E' if ext else 'O'}: {len(peaks[ext])}")
            logging.debug(f"skylines - props: {properties[ext]}")

            if self.can_plot:
                axs[ext].plot(spec1D_init[ext], label=f"{'O' if ext else 'E'} initial")
                axs[ext].vlines(peaks[ext], 0, np.mean(spec1D_init[ext]), color='r', label='Feature positions')
                axs[ext].hlines(*peak_width[1:], color='g', label='Initial widths')
                axs[ext].errorbar(
                    peaks[ext],
                    properties[ext]['prominences'],
                    xerr=np.array([
                        peaks[ext] - properties[ext]['left_ips'],
                        properties[ext]['right_ips'] - peaks[ext]
                        ]),
                    fmt='x',
                    label='Prominences'
                )
                
        if self.can_plot:
            for ax in axs: ax.legend()
            plt.show()

        # Transform data, skip if filename starts with 't'
        if self.must_transform or not filename.name.startswith('t'):
            wav2D, spec2D = self.transform(wav2D, spec2D, self.can_plot)

        # Convert to 1D spectra
        wav1D = np.median(wav2D, axis=1)
        spec1D = np.median(spec2D, axis=1)

        # Remove continuum
        if self.cont_ord > 0:
            for ext in range(2):
                # spec1D_init[ext] = self.remove_cont(spec1D_init[ext], wav1D[ext])
                # spec1D[ext] = self.remove_cont(spec1D[ext], wav1D[ext])
                pass

        if self.can_plot:
            # Plot transformed & normalized feature widths
            plt.plot(spec1D_init[0], label="O, initial")
            plt.plot(spec1D_init[1], label="E, initial")
            plt.plot(spec1D[0], label="O spec")
            plt.plot(spec1D[1], label="E spec")
            plt.legend()
            plt.show()

        # Find observed skylines
        # skyline_cols, properties = find_peaks(sky_norm, prominence=0.2) # 1000
        # prominence is basically the ylimit above which peaks should be found
        skyline_cols, skyline_wavs, properties = [], [], []
        for ext in range(2):
            cols, prop = signal.find_peaks(wav1D[ext], spec1D[ext], prominence=1)
            skyline_cols.append(cols)
            properties.append(prop)
            skyline_wavs.append(wav1D[ext, cols]) # col_to_wav(coeff2, skyline_cols)

            plt.plot(wav1D[ext, cols], cols, label='Observed')
            plt.xlabel('Wavelength ($\AA$)')
            plt.ylabel('Counts')
            plt.legend()
            plt.show()

        for ext in range(2):
            for i in range(len(skyline_cols[ext])):
                logging.debug(("skylines - found features:\n"
                    f"{properties[ext]['right_bases'][i] - properties[ext]['left_bases'][i]:4d}\t"
                    f"{skyline_cols[ext][i]:4d} - {skyline_wavs[ext][i]:.1f} - {properties[ext]['prominences'][i]:.2f}"
                ))


        # Return results
        return
    
    def plot(self, ) -> None:
        plt.style.use(Path(__file__).parent.resolve() / 'utils/STOPS.mplstyle')

        # Find deviation of observed skylines from known skylines

        # Find feature / initial feature widths

        # Load known skylines
        skylines = self.load_sky_lines()

        # Save results
        raise NotImplementedError
        
        return
    
    def show_frame(self, frame: np.ndarray, title: str = None, label: str = None, std: int = 1) -> None:
        if not self.can_plot:
            return
        
        fig, axs = plt.subplots(2, 1)
        axs[0].set_title(title)
        axs[0].imshow(
            frame[0],
            label=f"O beam - {label}",
            vmin=frame[0].mean() - std * frame[0].std(),
            vmax=frame[0].mean() + std * frame[0].std()
        )
        axs[1].imshow(
            frame[1],
            label=f"E beam - {label}",
            vmin=frame[1].mean() - std * frame[0].std(),
            vmax=frame[1].mean() + std * frame[0].std()
        )
        for ax in axs: ax.legend()
        plt.show()

        return
            

    # MARK: Process all listed images
    def process(self,) -> None:
        if self.beams == 'OE':
            for fl in self.fits_list:
                logging.info(f"'OE' skylines of {fl}.")
                self.skylines(fl)
                self.plot()

        if self.beams in ['O', 'E']:
            for fl in self.fits_list:
                logging.info(f"{self.beams} skylines of {fl}.")
                self.skylines(fl)
                self.plot()

        return


# MARK: Main function
def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])
