#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module for joining the split FITS files with an external wavelength solution."""

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

# from lacosmic import lacosmic # ccdproc ~6x faster
from ccdproc import cosmicray_lacosmic as lacosmic

from utils.specpolpy3 import read_wollaston, split_sci
from utils.SharedUtils import get_files, get_arc
from utils.Constants import DATADIR, SAVE_PREFIX, CR_PARAMS


# MARK: Join Docstring
class Join:
    """
    The `Join` class allows for the joining of previously 
    split files and the appending of an external wavelength 
    solution to the `polsalt` FITS file format.

    Parameters
    ----------
    data_dir : str
        The path to the data to be joined
    database : str, optional
        The name of the `IRAF` database folder.
        (The default is "database")
    fits_list : list[str], optional
        A list of pre-reduced `polsalt` FITS files to be joined within `data_dir`.
        (The default is ``None``, `Join` will search for `mxgbp*.fits` files)
    solutions_list: list[str], optional
        A list of solution filenames from which the wavelength solution is created.
        (The default is ``None``, `Join` will search for `fc*` files within the `database` directory)
    split_row : int, optional
        The row along which the data of each extension in the FITS file was split.
        Necessary when Joining cropped files.
        (The default is 517, the SALT RSS CCD's middle row)
    save_prefix : dict[str, list[str]], optional
        The prefix with which the previously split `O`- & `E`-beams were saved.
        Used for detecting if cropping was applied during the splitting procedure.
        (The default is SAVE_PREFIX (See Notes))
    verbose : int, optional
        The level of verbosity to use for the Cosmic ray rejection
        (The default is 30, I.E. logging.INFO)
    
    Attributes
    ----------
    fc_files : list[str]
        Valid solutions found from `solutions_list`.
    custom : bool
        Internal flag for whether `solutions_list` uses the `IRAF` or a custom format.
        See Notes for custom solution formatting.
        (Default (inherited from `solutions_list`) is False)
    arc : str
        Deprecated. Name of arc FITS file within `data_dir`.
    data_dir
    database
    fits_list
    split_row
    save_prefix


    Methods
    -------
    get_solutions(wavlist: list | None, prefix: str = "fc")
        -> (fc_files, custom): tuple[list[str], bool]
        Parse `solutions_list` and return valid solution files and if they are non-`IRAF` solutions.
    parse_solution(fc_file: str, xshape: int, yshape: int)
        -> tuple[dict[str, int], np.ndarray]
        Loads the wavelength solution file and parses keywords necessary for creating the wavelength extension.
    join_file(file: os.PathLike)
        -> None
        Joins the files, 
        attaches the wavelength solutions, 
        performs cosmic ray cleaning, 
        masks the extension, 
        and checks cropping performed in `Split`.
        Writes the FITS file in a `polsalt` valid format.
    check_crop(hdu: pyfits.HDUList, o_file: str, e_file: str)
        -> int
        Opens the split `O`- and `E`-beam FITS files and returns the amount of cropping that was performed.
    process()
        -> None
        Calls `join_file` on each file in `fits_list` for automation.

    
    Other Parameters
    ----------------
    no_arc : bool, optional
        Deprecated. Decides whether the arc frames should be processed.
        (The default is False, `polsalt` has no use for the arc after wavelength calibrations)
    **kwargs : dict
        keyword arguments. Allows for passing unpacked dictionary to the class constructor.
    
    Notes
    -----
    Constants set are:
        DATADIR
        SAVE_PREFIX
        CR_PARAMS

    Custom wavelength solutions must be formatted as:
        `x`,
        `y`,
        *coefficients...
    where the solutions are of order (`x` by `y`) and contain x*y coefficients.
    The name of the custom wavelength solution file must contain either "cheb" or "leg"
    for Chebychev or Legendre wavelength solutions, respectively.

    Cosmic ray rejection is performed using lacosmic [1]_ implemented in ccdproc via astroscrappy [2]_.

    References
    ----------
    .. [1] van Dokkum 2001, PASP, 113, 789, 1420 (article : http://adsabs.harvard.edu/abs/2001PASP..113.1420V)
    .. [2] https://zenodo.org/records/1482019
    
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

    # MARK: Find 2D WAV Functions
    def get_solutions(
        self,
        wavlist: list[str] | None,
        prefix: str = "fc"
    ) -> tuple[list[str], bool]:
        """
        Get the list of wavelength solution files.

        Parameters
        ----------
        wavlist : list[str] | None
            A list of custom wavelength solutions files.
            If ``None``, `Join` will search for wavelength solutions in the `database` directory.
        prefix : str, optional
            The prefix of the wavelength solution files.
            (Defaults to "fc")

        Returns
        -------
        tuple[list[str], bool]
            A tuple containing the list of wavelength solutions files and 
            a boolean indicating whether custom solutions were provided.

        """
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

    # MARK: Parse 2D WAV Function
    def parse_solution(
        self,
        fc_file: str,
        xshape: int,
        yshape: int
    ) -> tuple[dict[str, int], np.ndarray]:
        """
        Parse the 2D wavelength solution function from `fc_file`.

        Parameters
        ----------
        fc_file : str
            The filename of the wavelength solutions file.
        xshape : int
            The x-order of the 2D solution.
        yshape : int
            The y-order of the 2D solution.

        Returns
        -------
        tuple[dict[str, int], np.ndarray]
            A tuple containing a dictionary of the parameters of the solution function 
            and the function coefficients.

        """
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

    # MARK: Join Files
    def join_file(self, file: os.PathLike) -> None:
        """
        Join the `O`- and `E`-beams, attach the wavelength solutions, 
        perform cosmic ray cleaning, mask the extensions, 
        and checks cropping performed by `Split`.
        Write the FITS file in a `polsalt` valid format.

        Parameters
        ----------
        file : os.PathLike
            The path of the FITS file to be joined.

        See Also
        --------
        IRAF - `fitcoords` task
            https://iraf.net/irafdocs/formats/fitcoords.php,
        numpy.polynomial.chebyshev.chebgrid2d
            https://numpy.org/doc/stable/reference/generated/numpy.polynomial.chebyshev.chebgrid2d.html
        numpy.polynomial.legendre.leggrid2d
            https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggrid2d.html

        """
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

        # MARK: Join (Wav. Ext.)
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
                msg = (
                    "Function type not recognised, please wavelength "
                    "calibrate using either chebychev or legendre."
                )
                raise Exception(msg)

            # MARK: Cosmic Ray Cleaning
            # See utils.Constants for `CR_PARAMS` discussion
            whdu["SCI"].data[num] = lacosmic(
                whdu["SCI"].data[num],
                # contrast=CR_PARAMS['CR_CONTRAST'],
                # threshold=CR_PARAMS['CR_THRESHOLD'],
                # neighbor_threshold=CR_PARAMS['CR_NEIGHBOUR_THRESHOLD'],
                # effective_gain=CR_PARAMS['GAIN'],
                # background=CR_PARAMS['BACKGROUND'],
                readnoise=CR_PARAMS['READNOISE'],
                gain=CR_PARAMS['GAIN'],
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

    # MARK: Check Crop
    def check_crop(
        self,
        hdu: pyfits.HDUList,
        o_file: str,
        e_file: str
    ) -> int:
        """
        Check if cropping is necessary when joining `O`- and `E`-beams.

        Parameters
        ----------
        hdu : astropy.io.fits.HDUList
            The HDUList to check for cropping.
        o_file : str
            The name of the previously split `O`-beam FITS file.
        e_file : str
            The name of the previously split `E`-beam FITS file.

        Returns
        -------
        int
            The number of rows which were cropped by `Split`.

        """
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

    # MARK: Process all Listed Images
    def process(self) -> None:
        """Process all FITS images stored in the `fits_list` attribute"""
        for target in self.fits_list:
            logging.debug(f"Processing {target}")
            self.join_file(target)

        return


def main(argv) -> None:
    """Main function."""
    
    return


if __name__ == "__main__":
    main(sys.argv[1:])
