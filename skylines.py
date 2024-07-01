"""
Module for analyzing the sky lines of a wavelength calibrated image.
"""

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
from utils.Constants import SAVE_SKY, FIND_PEAK_PARAMS, ARC_FILE

# print(
#  [logging.getLogger(name) for name in logging.root.manager.loggerDict]
# )
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.INFO)
pil_logger = logging.getLogger('PIL')
pil_logger.setLevel(logging.INFO)
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
        Transforms the input wavelength and spectral data based on
        the given wavelength solution.
    rmvCont(self) -> np.ndarray:
        Removes the continuum from the spectrum.
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
        self.fits_list, self.arc_list = find_files(
            data_dir=self.data_dir,
            filenames=filenames,
            prefix="wmxgbp", # t[o|e]beam
            ext="fits",
            sep_arc=True,
        )
        self._beams = None
        self.beams = beams
        self.ccds = 1
        if split_ccd:
            # See cross_correlate for initial implementation
            self.ccds = 3
        
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

        self.max_difference = 5

        self.wav_unit = "$\AA$"

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
    
    # MARK: Find Peaks
    def find_peaks(
        self,
        spec: np.ndarray,
        axis: int | None = None,
        min_height: float = 0.5,
        **kwargs,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Finds the peaks in the given spectral data.

        Parameters
        ----------
        spec : np.ndarray
            The spectral data.
        bpm : np.ndarray
            The bad pixel mask.
        min_height : float, optional
            The minimum height of the peaks, by default 0.5.
        rel_height : float, optional
            The relative height of the peaks, by default 0.05.

        Returns
        -------
        peaks, properties : tuple[np.ndarray, dict]
            The peaks and their properties.

        """
        peaks = []
        props = []

        for ext in range(len(self.beams)):
            if axis is None:
                row_mean = spec[ext]
            else:
                row_mean = np.mean(spec[ext], axis=axis)

            peak, property = signal.find_peaks(
                row_mean,
                prominence=min_height * np.max(row_mean),
                width=0,
                **kwargs,
            )
            peaks.append(peak)
            props.append(property)

        if self.can_plot:
            fig, axs = plt.subplots(2, 1)
            for ext in range(len(self.beams)):
                axs[ext].plot(
                    row_mean,
                    label=f"{'E' if ext else 'O'}"
                )
                axs[ext].plot(
                    peak,
                    row_mean[peak],
                    "x",
                    label=f"{'E' if ext else 'O'} peaks"
                )
                axs[ext].legend()
                plt.show()

        logging.debug(f"find_peaks - peaks: {[len(i) for i in peaks]}")
        logging.debug(f"find_peaks - props: {[key for key in props[0].keys()]}")

        return peaks, props
    
    # MARK: Min. of Diff. Matrix
    @staticmethod
    def min_diff_matrix(
        A: np.ndarray,
        B: np.ndarray,
        max_diff: int = 100
    ) -> np.ndarray:
        """
        Find the minimum difference between the elements of two arrays.

        Parameters
        ----------
        A : np.ndarray
            The first 1d array.
        B : np.ndarray
            The second 1d array.
        max_diff : int, optional
            The maximum difference allowed, by default -1.
        
        Returns
        -------
        A : np.ndarray (len(A))
            The elements of the first array.
        min_vals : np.ndarray (len(A))
            The minimum difference between the elements of the two arrays.
        min_idxs : np.ndarray (len(A))
            The indices of the minimum difference between
            the elements of the two arrays.
        
        """
        # Compute the difference matrix using transpose
        diff = np.abs(A - B[:, np.newaxis])
        
        # Find the minimum value in each row (A) of `diff`
        min_vals = np.min(diff, axis=0)
        min_idxs = np.argmin(diff, axis=0)
        # TODO: Recalculate min_val after
        # selecting best min_val and removing the corresponding row/column

        logging.debug(f"min_diff_matrix - min_vals: {np.round(min_vals, 2)}")
        logging.debug(f"min_diff_matrix - min_idxs: {min_idxs}")

        max_mask = min_vals <= max_diff
        
        return A[max_mask], min_vals[max_mask], min_idxs[max_mask]

    # MARK: Load File Data
    def load_file_data(
        self,
        filename: Path
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Loads the data from the given file.

        Parameters
        ----------
        filename : Path
            The path to the file to be loaded.

        Returns
        -------
        spec, wav, bpm : tuple[np.ndarray, np.ndarray, np.ndarray]
            The wavelength, spectral, and bad pixel mask data.

        """
        # Load data from self.beams extension
        with pyfits.open(filename) as hdul:
            exts = [0, 1] if len(self.beams) == 2 else 0 + self.beams == 'E'
            spec2D = np.atleast_3d(hdul["SCI"].data[exts])
            wav2D = np.atleast_3d(hdul["WAV"].data[exts])
            bpm2D = np.atleast_3d(hdul["BPM"].data[exts].astype(bool))

            logging.info(
                f"load_file_data - {filename.name} - shape: {spec2D.shape}"
            )

            return spec2D, wav2D, bpm2D

    # MARK: Load Sky or Arc Lines
    def load_lines(
            self,
            filename: Path | None = None,
            dtype: list[tuple] = [('wav', float), ('flux', float)],
            skip_header: int = 3,
            skip_footer: int = 1
        ) -> np.ndarray:
        """
        Loads the sky or arc lines from the given file.

        Parameters
        ----------
        filename : Path | None, optional
            The path to the file to be loaded.
            Defaults to loading the skylines from `utils/sky.salt`

        Returns
        -------
        sky_lines : np.ndarray['wav', 'flux']
            The sky lines from the file.
        
        """
        usecols = None
        if filename:
            filename = Path(__file__).parent.resolve() / filename
            usecols = (0, 1)
        else:
            filename = Path(__file__).parent.resolve() / 'utils/sky.salt'

        lines = np.genfromtxt(
            filename,
            dtype=dtype,
            skip_header=skip_header,
            skip_footer=skip_footer,
            usecols=usecols
        )

        logging.debug(
            f"load_lines - {filename.name} - shape: {lines.shape}"
        )

        return lines

    # MARK: Mask Traces
    def mask_traces(
            self,
            spec: np.ndarray,
            bpm: np.ndarray,
            max_traces: int = 1,
            tr_pad: int = 5,
            bg_margin: int = 10,
            lr_margins: list[int] = [10, 10],
            h_min: float = 0.5,
            h_rel: float = 1 - 0.05,
        ) -> np.ndarray:
        """
        Masks the traces in the bad pixel mask.

        Parameters
        ----------
        spec : np.ndarray
            The spectral data.
        bpm : np.ndarray
            The bad pixel mask.

        Returns
        -------
        bpm : np.ndarray
            The updated bad pixel mask.

        """
        # Base mask
        bpm[:, :bg_margin] = True
        bpm[:, -bg_margin:] = True
        bpm[:, :, :lr_margins[0]] = True
        bpm[:, :, -lr_margins[1]:] = True

        # Get the traces
        traces, tr_props = self.find_peaks(
            spec,
            axis=1,
            min_height=h_min,
            rel_height=h_rel
        )

        for ext in range(len(self.beams)):
            # Mask the traces
            for i in range(len(traces[ext][:max_traces])):
                lb = max(
                    0,
                    int(tr_props[ext]['left_ips'][i]) - tr_pad
                )
                ub = min(
                    spec.shape[-1],
                    int(tr_props[ext]['right_ips'][i]) + tr_pad
                )
                bpm[ext, lb : ub] = True
                # TODO: Relocate targets after initial masking

        logging.info(f"mask_traces - {min(max_traces, len(traces))} of {len(traces)} traces masked.")

        return bpm

    # MARK: Transform Spectra
    def transform(
        self,
        spec: np.ndarray,
        wav_sol: np.ndarray,
        row_max: int | None = None,
        resPlot: bool = False,
    ) -> np.ndarray:
        """
        Transforms the input wavelength and spectral data
        based on the given wavelength solution.

        Parameters
        ----------
        spec : np.ndarray
            The spectral data.
        wav_sol : np.ndarray
            The wavelength solution.
        resPlot : bool, optional
            Flag indicating whether to plot the results, default False.

        Returns
        -------
        spec, wav : np.ndarray
            The transformed wavelength and spectral data.

        """
        # Create arrays to return
        cs = np.zeros_like(spec)
        cw = np.zeros_like(wav_sol.mean(axis=1))

        for ext in range(len(self.beams)):

            if row_max:
                avg_max = row_max
            else:
                # Get middle row (to interpolate the rest of the rows to)
                avg_max = np.mean(spec[ext], axis=1).argmax()

            # Correct extensions based on wavelength
            # Get wavelength values at row with most trace
            cw[ext] = wav_sol[ext, avg_max]

            # Spec ext
            for row in range(cs.shape[1]):
                cs[ext, row] = np.interp(
                    cw[ext],
                    wav_sol[ext, row],
                    spec[ext, row]
                )
                # f_2d = interpolate.interp2d(
                #     wav_sol[ext, row],
                #     np.arange(rows),
                #     spec[ext],
                # )
                # cs[ext] = f_2d(cw[ext, row], np.arange(rows))

        # Plot results
        if resPlot:
            fig, axs = plt.subplots(2, 1, figsize=[20, 4])
            for ext in range(len(self.beams)):
                axs[ext].imshow(
                    cs[ext],
                    vmax=cs[ext].mean() + 2*cs[ext].std(),
                    vmin=cs[ext].mean() - 2*cs[ext].std()
                )

                logging.debug(f"{'E' if ext else 'O'} Average continuum = {np.median(np.median(cs[ext], axis=0)):4.3f}")

                axx = axs[ext].twinx()
                axx.hlines(
                    np.median(np.median(cs[ext], axis=0)),
                    0,
                    cs[ext].shape[-1],
                    colors='black'
                )
                axx.plot(
                    cs[ext].mean(axis=0),
                    "k",
                    label=f"mean {'E' if ext else 'O'}"
                )
                axx.plot(
                    np.median(cs[ext], axis=0),
                    "r",
                    label=f"median {'E' if ext else 'O'}"
                )
                axx.legend()
            plt.show()

        logging.info(f"transform - {cs.shape} transformed.")

        return cs, cw
    
    # MARK: Plot
    def plot(
            self,
            spectra,
            wavelengths,
            peaks,
            properties,
            arc: bool = False,
        ) -> None:
        plt.style.use(Path(__file__).parent.resolve() / 'utils/STOPS.mplstyle')
        plt.rcParams['figure.subplot.hspace'] *= len(self.beams)
        
        def norm(x):
            return (x - np.min(x)) / (np.max(x) - np.min(x))

        # Load known lines
        if arc:
            lines = self.load_lines(filename=f'utils/RSS_arc_files/{ARC_FILE}')
        else:
            lines = self.load_lines()

        lines = lines[
            (lines['wav'] > wavelengths[1][0][0].min()) &
            (lines['wav'] < wavelengths[1][0][0].max())
        ]


        # Create plot for results
        fig, axs = plt.subplots(2, self.ccds, sharex='col', sharey='row')

        # Convert axs to a 2D array if ccd count is 1
        if self.ccds == 1:
            axs = np.swapaxes(np.atleast_2d(axs), 0, 1)

        for fl in range(len(self.arc_list if arc else self.fits_list)):

            # set color cycle
            color=next(axs[0, 0]._get_lines.prop_cycler)['color']

            for ext in range(len(self.beams)):
                    
                for ccd in range(self.ccds):

                    # spectrum (transformed)
                    ccdrange = spectra[1][fl][ext].shape[-1] // self.ccds
                    axs[0, ccd].plot(
                        wavelengths[1][fl][ext][
                            ccdrange*ccd:ccdrange*(ccd+1)
                        ],
                        norm(spectra[1][fl][ext][
                            ccdrange*ccd:ccdrange*(ccd+1)
                        ]) * 100 + 10 * ext + 30 * fl,
                        color=color,
                        linestyle='dashed' if ext else 'solid',
                        label = f"${{{self.beams[ext]}}}_{{{fl + 1}}}^{{+ {10*ext + 30*fl}}}$" if ccd == 0 else None,
                    )

                    # deviation
                    sky_wavs, dev, peak_idx = self.min_diff_matrix(
                        lines['wav'],
                        wavelengths[1][fl][ext][peaks[1][fl][ext]],
                        max_diff=self.max_difference,
                    )

                    # width/initial width
                    width = properties[1][fl][ext]['widths'][peak_idx]
                    width_i = np.zeros_like(width)

                    sky_i, i_dev, i_idx = self.min_diff_matrix(
                        lines['wav'],
                        wavelengths[0][fl][ext][peaks[0][fl][ext]],
                        max_diff=self.max_difference,
                    )

                    width_i = np.array([
                        properties[0][fl][ext]['widths'][
                            np.where(wav == sky_i)[0][0]
                        ]
                        if wav in sky_i else 1000
                        for wav in sky_wavs
                    ])
                    width_ratio = (width / width_i) - 1
                    width_ratio[width_ratio < 0] = 0

                    ylolims = width_ratio > self.max_difference
                    width_ratio[
                        width_ratio > self.max_difference
                    ] = self.max_difference // 2

                    ok = np.where(
                        (sky_wavs >= wavelengths[1][fl][ext].data[ccdrange*ccd]) &
                        (sky_wavs <= wavelengths[1][fl][ext].data[ccdrange*(ccd+1)])
                    )
                    axs[1, ccd].errorbar(
                        sky_wavs[ok],
                        dev.data[ok],
                        yerr=(width_ratio[ok] * 0, width_ratio[ok]),
                        lolims=ylolims[ok],
                        fmt="." if ext else "x",
                        alpha=0.8,
                        color=color,
                        # markeredgecolor='white',
                        # markeredgewidth=0.5,
                        # label=f"${self.beams[ext]}_{{{fl + 1}}}$",
                    )

                    logging.debug(f"plot - RMS: {np.sqrt(np.mean(dev**2)):2.3f}")

        for ccd in range(self.ccds):
            # spectrum
            ok = np.where(
                (lines['wav'] >= wavelengths[1][fl][0].data[ccdrange*ccd]) &
                (lines['wav'] <= wavelengths[1][fl][0].data[ccdrange*(ccd+1)])
            )
            axs[0, ccd].plot(
                lines['wav'][ok],
                lines['flux'][ok] * 0,
                'x',
                color='C4',
                label="\\textsc{salt}\nModel" if ccd == 0 else None,
            )
            for x in lines['wav'][ok]: axs[0, ccd].axvline(x, ls='dashed', c='0.7')

        axs[0, 0].set_ylabel("Rel. Intensity\n($\%$)")
        axs[1, 0].set_ylabel(
            "Closest Peak\n($|\lambda_{salt} - \lambda_{obs.}|$)"
        )
        # for ax in axs[:, 0]:
        #     ax.legend(loc='upper left', ncols=(fl + 1) * (ext + 1) + 1)
        leg = fig.legend(
            loc='center',
            ncol=min(8, len(spectra[0]) + 1),
            columnspacing=0.5,
            bbox_to_anchor=(
                np.mean((
                    plt.rcParams['figure.subplot.left'],
                    plt.rcParams['figure.subplot.right']
                )),
                np.mean((
                    plt.rcParams['figure.subplot.bottom'],
                    plt.rcParams['figure.subplot.top']
                ))
            ),
        )
        leg.set_draggable(True)
        for ax in axs[1, :]:
            ax.grid(axis='y')
        
        # fig.add_subplot(111, frameon=False)
        # # hide tick and tick label of the big axis
        # plt.tick_params(
        # labelcolor='none',
        #   which='both',
        #   top=False,
        #   bottom=False,
        #   left=False,
        #   right=False
        # )
        axs[-1, 0 if self.ccds == 1 else 1].set_xlabel(
            f"Wavelength ({self.wav_unit})"
        )

        # plt.tight_layout()
        
        plt.show()

        # Save results
        if self.save_prefix:
            fig.savefig(fname=self.save_prefix)
        
        return   

    # MARK: Process all listed images
    def process(self, arc: bool=False) -> None:
        files = self.fits_list
        if arc:
            files = self.arc_list
        
        logging.info(f"Processing '{self.beams}' lines.")

        spectra = [[], []]
        wavs = [[], []]
        peaks = [[], []]
        peak_props = [[], []]

        for fl in files:
            # Load data
            spec2d, wav2d, bpm2d = self.load_file_data(fl)

            # Mask traces in BPM
            bpm2d = self.mask_traces(
                spec2d,
                bpm2d,
                max_traces=0,
                bg_margin=15,
                h_min=0.05
            )
            m_spec2d = np.ma.masked_array(spec2d, mask=bpm2d) # spec2d
            m_wav2d = np.ma.masked_array(wav2d, mask=bpm2d) # wav2d

            # Initial spectra
            spec_i = np.mean(m_spec2d, axis=-2)
            wav_i = np.mean(m_wav2d, axis=-2)

            # Transform data
            t_spec2d, t_wav = self.transform(
                m_spec2d,
                m_wav2d,
                resPlot=self.can_plot
            )

            # Final spectra
            spec_f = np.mean(t_spec2d, axis=-2)
            wav_f = t_wav

            # Find peaks
            peaks_i, props_i = self.find_peaks(
                spec_i,
                **FIND_PEAK_PARAMS
            )
            peaks_f, props_f = self.find_peaks(
                spec_f,
                **FIND_PEAK_PARAMS
            )

            spectra[0].append([*spec_i])
            spectra[1].append([*spec_f])
            wavs[0].append([*wav_i])
            wavs[1].append([*wav_f])
            peaks[0].append([*peaks_i])
            peaks[1].append([*peaks_f])
            peak_props[0].append([*props_i])
            peak_props[1].append([*props_f])

        # Plot results
        self.plot(spectra, wavs, peaks, peak_props, arc=arc)

        if arc:
            return
        elif self.arc_list:
            self.process(arc=True)

        return


# MARK: Main function
def main(argv) -> None:
    return

if __name__ == "__main__":
    main(sys.argv[1:])
