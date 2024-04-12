#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__email__ = "justin.jb78+Masters@gmail.com"

# MARK: Imports

import os
import sys
import logging
import re

import numpy as np
from numpy.polynomial.chebyshev import chebgrid2d as chebgrid2d
from numpy.polynomial.legendre import leggrid2d as leggrid2d
from astropy.io import fits as pyfits
from utils.specpolpy3 import read_wollaston, split_sci

# from lacosmic import lacosmic
from ccdproc import cosmicray_lacosmic as lacosmic


from utils.SharedUtils import get_files, get_arc

# MARK: Constants

DATADIR = os.path.expanduser("~/polsalt-beta/polsalt/data/")
SAVE_PREFIX = {"beam": ["obeam", "ebeam"], "arc": ["oarc", "earc"]}

# CR Cleaning parameters (Deprecated with ccdproc implementation)
CR_CONTRAST = 2
CR_THRESHOLD = 4
CR_NEIGHBOUR_THRESHOLD = 4
# Gain and Readnoise still valid for ccdproc, sourced from
# https://pysalt.salt.ac.za/proposal_calls/current/ProposalCall.html
GAIN = 1
READNOISE = 3.3


class Join:
    # MARK: Join Docstring

    """
    Join class allows for the seperate call of joining the wavelength calibrated O & E beam FITS files

    Parameters
    ----------
    path_name : str
        The path to the data (wmxgbp*.fits files) to be joined
    split_row : int, optional
        The row that the data was split along.
        (The default is 517, the middle row of the CCD's)
    no_arc : bool, optional
        Decides whether the arc frames should be recombined.
        (The default is True, since polsalt only uses the arc frames until spectral extraction)
    save_pref : list of str, optional
        The prefix that the O & E beams are saved as.
        (The default is ["obeam", "ebeam"], which is what split defaults to)

    Returns
    -------
    joined_FITS : list of FITS
        A list of FITS files that were joined and can be returned to polsalt.

    Raises
    ------
    # TODO@JustinotherGitter : Complete docs for which errors are raised and when
    """

    # MARK: Join init

    def __init__(
        self,
        data_dir: str,
        database: str = "database",
        fits_list: list[str] = None,
        solutions_list: list[str] = None,
        split_row: int = 517,
        no_arc: bool = True,
        save_prefix=None,
        verbose: int = 30,
        **kwargs,
    ) -> None:
        self.data_dir = data_dir
        self.database = database
        self.fits_list = get_files(
            data_dir=self.data_dir,
            filenames=fits_list,
            prefix="mxgbp",
            extention="fits",
        )
        self.fc_files, self.custom = self.get_solutions(solutions_list)
        self.split_row = split_row
        self.save_prefix = SAVE_PREFIX
        if type(save_prefix) == dict:
            self.save_prefix = save_prefix

        self.no_arc = no_arc
        self.arc = get_arc(self.fits_list)

        self.verbose = verbose < 30
        return

    def get_solutions(
        self,
        wavlist: list | None,
        prefix: str = "fc"
    ) -> tuple[list[str], bool]:
        # MARK: Find 2D WAV Functions

        # No custom solutions
        if not wavlist:
            # Handle finding solutions
            ws = []
            for fl in os.listdir(
                os.path.join(self.data_dir, self.database)
            ):
                if os.path.isfile(
                    os.path.join(self.data_dir, self.database, fl)
                ) and (prefix == fl[0:2]):
                    ws.append(fl)

            if len(ws) != 2:
                # Handle incorrect number of solutions found
                msg = (
                    f"Incorrect amount of wavelength solutions "
                    f"({len(ws)} fc... files) found in the solution "
                    f"dir.: {os.path.join(self.data_dir, self.database)}"
                )
                logging.error(msg)
                raise FileNotFoundError(msg)

            return (sorted(ws, reverse=True), False)

        # Custom solution
        if len(wavlist) >= 2:
            if len(wavlist) > 2:
                logging.warn(f" Too many solutions: {wavlist}")
                wavlist = wavlist[:2]

            for fl in wavlist:
                if not os.path.isfile(os.path.join(self.data_dir, fl)):
                    msg = (
                        f"{fl} not found in the "
                        f"data directory {self.data_dir}"
                    )
                    logging.error(msg)
                    raise FileNotFoundError(msg)

            return (sorted(wavlist, reverse=True), True) 

    def parse_solution(
        self,
        fc_file,
        xshape: int,
        yshape: int
    ) -> tuple[dict[str, int], np.ndarray]:
        # MARK: Parse 2D WAV Function

        fit_params = {}
        coeff = []

        if self.custom:
            # Load coefficients
            coeff = np.loadtxt(fc_file)

            fit_params["xorder"] = coeff[0].astype(int)
            fit_params["yorder"] = coeff[1].astype(int)
            coeff = coeff[2:]

            f_type = 3
            if "cheb" in str(fc_file): f_type = 1
            elif "leg" in str(fc_file): f_type = 2
            fit_params["function"] = f_type

            fit_params["xmin"], fit_params["xmax"] = 1, xshape
            fit_params["ymin"], fit_params["ymax"] = 1, yshape

        else:
            # Parse IRAF fc database files
            file_contents = []
            with open(self.database + "/" + fc_file) as fcfile:
                for i in fcfile:
                    file_contents.append(re.sub(r"[\n\t\s]*", "", i))

            if file_contents[9] != "1.":  # xterms - Cross-term type
                msg=(
                    "Cross-term not recognised (always 1 for "
                    "FITCOORDS), redo FITCOORDS or change manually."
                )
                raise Exception(msg)

            fit_params["function"] = int(file_contents[6][:-1])

            fit_params["xorder"] = int(file_contents[7][:-1])
            fit_params["yorder"] = int(file_contents[8][:-1])

            fit_params["xmin"] = int(file_contents[10][:-1])
            fit_params["xmax"] = xshape
            # int(file_contents[11][:-1])# stretch fit over x
            fit_params["ymin"] = int(file_contents[12][:-1])
            fit_params["ymax"] = yshape
            # int(file_contents[13][:-1])# stretch fit over y

            coeff = np.array(file_contents[14:], dtype=float)

        coeff = np.reshape(
            coeff,
            (fit_params["xorder"], fit_params["yorder"])
        )

        return (fit_params, coeff)

    def join_file(self, file: os.PathLike) -> None:
        # MARK: Join Files

        # Create empty wavelength appended hdu list
        whdu = pyfits.HDUList()
        primary_ext = ""

        # Handle prefix and names
        pref = "arc" if file == self.arc else "beam"
        o_file = self.save_prefix[pref][0] + file.name[-9:]
        e_file = self.save_prefix[pref][1] + file.name[-9:]

        # Open file
        with pyfits.open(file) as hdu:
            # Check if file has been cropped
            cropsize = self.check_crop(hdu, o_file, e_file)

            y_shape = int(hdu["SCI"].data.shape[0] / 2) - cropsize
            x_shape = hdu["SCI"].data.shape[1]

            # No differences in "PRIMARY" extention header
            primary_ext = hdu["PRIMARY"]
            whdu.append(primary_ext)

            for ext in ["SCI", "VAR", "BPM"]:
                whdu.append(pyfits.ImageHDU(name=ext))
                whdu[ext].header = hdu[ext].header.copy()
                whdu[ext].header["CTYPE3"] = "O,E"

                # Create empty extentions with correct order and format
                if ext == "BPM":
                    whdu[ext].data = np.zeros(
                        (2, y_shape, x_shape),
                        dtype="uint8"
                    )
                    whdu[ext].header["BITPIX"] = "-uint8"
                else:
                    whdu[ext].data = np.zeros(
                        (2, y_shape, x_shape),
                        dtype=">f4"
                    )
                    whdu[ext].header["BITPIX"] = "-32"

                # Fill in empty extentions
                if cropsize:
                    temp_split = split_sci(
                        hdu,
                        self.split_row,
                        ext=ext
                    )[ext].data
                    whdu[ext].data[0] = temp_split[0, cropsize:]
                    whdu[ext].data[1] = temp_split[1, 0:-cropsize]

                else:
                    whdu[ext].data = split_sci(
                        hdu,
                        self.split_row,
                        ext=ext
                    )[ext].data

        # End of hdu calls, close hdu

        # MARK: Wavelength Extension

        # See:
        # https://iraf.net/irafdocs/formats/fitcoords.php,
        # https://numpy.org/doc/stable/reference/generated/numpy.polynomial.chebyshev.chebgrid2d.html
        # https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggrid2d.html

        whdu.append(pyfits.ImageHDU(name="WAV"))
        wav_header = whdu["SCI"].header.copy()
        wav_header["EXTNAME"] = "WAV"
        wav_header["CTYPE3"] = "O,E"
        whdu["WAV"].header = wav_header

        whdu["WAV"].data = np.zeros(
            whdu["SCI"].data.shape,
            dtype=">f4"
        )

        for num, fname in enumerate(self.fc_files):
            pars, chebvals = self.parse_solution(
                fname,
                x_shape,
                y_shape
            )

            if pars["function"] == 1:  # Function type (1 = chebyshev)

                # Set wavelength extention values to function
                whdu["WAV"].data[num] = chebgrid2d(
                    x=np.linspace(-1, 1, pars["ymax"]),
                    y=np.linspace(-1, 1, pars["xmax"]),
                    c=chebvals,
                )

            elif pars["function"] == 2:  # Function type (2 = legendre)

                # Set wavelength extention values to function
                whdu["WAV"].data[num] = leggrid2d(
                    x=np.linspace(-1, 1, pars["ymax"]),
                    y=np.linspace(-1, 1, pars["xmax"]),
                    c=chebvals,
                )

            else:
                # TODO@JustinotherGitter: Handle other functions?
                msg = (
                    "Function type not recognised, please wavelength "
                    "calibrate using either chebychev or legendre."
                )
                raise Exception(msg)

            # MARK: Cosmic Ray Cleaning

            # Constants set work well with base lacosmic.
            # ccdproc lacosmic needs testing as parameter names differ
            # ccdproc lacosmic uses 'gain', defaults to 1
            # https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/08-03-Cosmic-ray-removal.html

            whdu["SCI"].data[num] = lacosmic(
                whdu["SCI"].data[num],
                # contrast=CR_CONTRAST,
                # threshold=CR_THRESHOLD,
                # neighbor_threshold=CR_NEIGHBOUR_THRESHOLD,
                # effective_gain=GAIN,
                # background=None,
                readnoise=READNOISE,
                verbose=self.verbose,
            )[0]

        # MARK: WAV masking

        # Left & Right Crop
        whdu["WAV"].data[whdu["WAV"].data[:] < 3_000] = 0.0
        whdu["WAV"].data[whdu["WAV"].data[:] >= 10_000] = 0.0

        # Top & Bottom Crop (shift\tilt)
        rpix_oc, cols, rbin, lam_c = read_wollaston(
            whdu,
            DATADIR + "wollaston.txt"
        )

        drow_oc = (rpix_oc - rpix_oc[:, int(cols / 2)][:, None]) / rbin

        ## Cropping as suggested
        for c, col in enumerate(drow_oc[0]):
            if np.isnan(col):
                continue

            if int(col) < 0:
                whdu["WAV"].data[0, int(col) :, c] = 0.0
            elif int(col) > cropsize:
                whdu["WAV"].data[0, 0 : int(col) - cropsize, c] = 0.0

        for c, col in enumerate(drow_oc[1]):
            if np.isnan(col):
                continue

            if int(col) > 0:
                whdu["WAV"].data[1, 0 : int(col), c] = 0.0
            elif (int(col) < 0) & (abs(int(col)) > cropsize):
                whdu["WAV"].data[1, int(col) + cropsize :, c] = 0.0

        # MARK: BPM masking

        whdu["BPM"].data[0] = np.where(
            whdu["WAV"].data[0] == 0,
            1,
            whdu["BPM"].data[0]
        )
        whdu["BPM"].data[1] = np.where(
            whdu["WAV"].data[1] == 0,
            1,
            whdu["BPM"].data[1]
        )

        whdu.writeto(f"w{os.path.basename(file)}", overwrite="True")

        return

    def check_crop(self, hdu, o_file, e_file) -> int:
        # MARK: Check Crop
        
        cropsize = 0
        o_y = 0
        e_y = 0

        with pyfits.open(o_file) as o:
            o_y = o[0].data.shape[0]

        with pyfits.open(e_file) as e:
            e_y = e[0].data.shape[0]

        if hdu["SCI"].data.shape[0] != (o_y + e_y):
            # Get crop size, assuming crop same on both sides
            cropsize = int((hdu["SCI"].data.shape[0] - o_y - e_y) / 2)

        return cropsize

    def process(self) -> None:
        # MARK: Process all Listed Images

        for target in self.fits_list:
            logging.debug(f"Processing {target}")
            self.join_file(target)

        return


def main(argv) -> None:
    return


if __name__ == "__main__":
    main(sys.argv[1:])
