"""Module for splitting ``polsalt`` FITS files."""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __init__ import __author__, __email__, __version__

# MARK: Imports
import os
import sys
import logging
from copy import deepcopy
from pathlib import Path

import numpy as np
from astropy.io import fits as pyfits

from utils.SharedUtils import find_files, find_arc
from utils.Constants import SAVE_PREFIX, CROP_DEFAULT, SPLIT_ROW


# MARK: Split Class
class Split:
    """
    The `Split` class allows for the splitting of `polsalt` FITS files 
    based on the polarization beam. The FITS files must have basic 
    `polsalt` pre-reductions already applied (`mxgbp...` FITS files).

    Parameters
    ----------
    data_dir : str
        The path to the data to be split
    fits_list : list[str], optional
        A list of pre-reduced `polsalt` FITS files to be split within `data_dir`.
        (The default is None, `Split` will search for `mxgbp*.fits` files)
    split_row : int, optional
        The row along which to split the data of each extension in the FITS file.
        (The default is SPLIT_ROW (See Notes), the SALT RSS CCD's middle row)
    no_arc : bool, optional
        Decides whether the arc frames should be recombined.
        (The default is False, `polsalt` has no use for the arc after wavelength calibrations)
    save_prefix : dict[str, list[str]], optional
        The prefix with which to save the  O & E beams.
        Setting `save_prefix` = ``None`` does not save the split O & E beams.
        (The default is SAVE_PREFIX (See Notes))

    Attributes
    ----------
    arc : str
        Name of arc FITS file within `data_dir`.
        `arc` = `""` if `no_arc` or not detected in `data_dir`.
    o_files, e_files : list[str]
        A list of the `O`- and `E`-beam FITS file names.
        The first entry is the arc file if `arc` defined.
    data_dir
    fits_list
    split_row
    save_prefix

    Methods
    -------
    split_file(file: os.PathLike)
        -> tuple[astropy.io.fits.HDUList]
        Handles creation and saving the separated FITS files 
    split_ext(hdulist: astropy.io.fits.HDUList, ext: str = 'SCI')
        -> astropy.io.fits.HDUList
        Splits the data in the `ext` extension along the `split_row`
    crop_file(hdulist: astropy.io.fits.HDUList, crop: int = CROP_DEFAULT (See Notes))
        -> tuple[numpy.ndarray]
        Crops the data along the edge of the frame, that is,
        `O`-beam cropped as [crop:], and 
        `E`-beam cropped as [:-crop].
    update_beam_lists(o_name: str, e_name: str)
        -> None
        Updates `o_files` and `e_files`.
    save_beam_lists(file_suffix: str = 'frames')
        -> None
        Creates (Overwrites if exists) and writes the `o_files` and `e_files` to files named 
        `o_{file_suffix}` and `e_{file_suffix}`, respectively.
    process()
        -> None
        Calls `split_file` and `save_beam_lists` on each file in `fits_list` for automation.
        
    Other Parameters
    ----------------
    **kwargs : dict
        keyword arguments. Allows for passing unpacked dictionary to the class constructor.

    Notes
    -----
    Constants Imported (See utils.Constants):
        SAVE_PREFIX
        CROP_DEFAULT
        SPLIT_ROW
    
    """
    # MARK: Split init
    def __init__(
        self,
        data_dir: Path,
        fits_list: list[str] = None,
        split_row: int = SPLIT_ROW,
        no_arc: bool = False,
        save_prefix: Path | None = None,
        **kwargs
    ) -> None:
        self.data_dir = data_dir
        self.fits_list = find_files(
            data_dir=data_dir,
            filenames=fits_list,
            prefix="mxgbp",
            ext="fits"
        )
        self.split_row = split_row
        self.save_prefix = SAVE_PREFIX
        if type(save_prefix) == dict:
            self.save_prefix = save_prefix

        self.arc = "" if no_arc else find_arc(self.fits_list)
        self.o_files = []
        self.e_files = []

        logging.debug(self.__dict__)
        return

    # MARK: Split Files
    def split_file(
        self,
        file: os.PathLike
    ) -> tuple[pyfits.HDUList]:
        """
        Split the single FITS file into separated `O`- and `E`- FITS files.

        Parameters
        ----------
        file : os.PathLike
            The name of the FITS file to be split.

        Returns
        -------
        tuple[astropy.io.fits.HDUList]
            Tuple containing the split O and E beam HDULists.
        
        """
        # Create empty HDUList
        O_beam = pyfits.HDUList()
        E_beam = pyfits.HDUList()

        # Open file and split O & E beams
        with pyfits.open(file) as hdul:
            O_beam.append(hdul["PRIMARY"].copy())
            E_beam.append(hdul["PRIMARY"].copy())

            # Split specific extention
            raw_split = self.split_ext(hdul, "SCI")

            # O_beam[0].data = raw_split['SCI'].data[1]
            # E_beam[0].data = raw_split['SCI'].data[0]
            O_beam[0].data, E_beam[0].data = self.crop_file(raw_split)

            # Handle prefix and names
            pref = "arc" if file == self.arc else "beam"
            o_name = self.save_prefix[pref][0] + file.name[-9:]
            e_name = self.save_prefix[pref][1] + file.name[-9:]

            # Add split data to O & E beam lists
            self.update_beam_lists(o_name, e_name, pref == "arc")

            # Handle don't save case
            if self.save_prefix == None:
                return O_beam, E_beam

            # Handle save case
            O_beam.writeto(o_name, overwrite=True)
            E_beam.writeto(e_name, overwrite=True)

            return O_beam, E_beam

    # MARK: Split extensions
    def split_ext(
        self,
        hdulist: pyfits.HDUList,
        ext: str = "SCI"
    ) -> pyfits.HDUList:
        """
        Split the data of the specified extension of `hdulist` into its `O`- and `E`- beams.

        Parameters
        ----------
        hdulist : astropy.io.fits.HDUList
            The FITS HDUList to be split.
        ext : str, optional
            The name of the extension to be split.
            (Defaults to 'SCI')

        Returns
        -------
        astropy.io.fits.HDUList
            The HDUList with the split applied.
        
        """
        hdu = deepcopy(hdulist)
        rows, cols = hdu[ext].data.shape

        # if odd number of rows, strip off the last one
        rows = int(rows / 2) * 2

        # how far split is from center of detector
        offset = int(self.split_row - rows / 2)

        # split arc into o/e images
        ind_rc = np.indices((rows, cols))[0]
        padbins = (ind_rc < offset) | (ind_rc > rows + offset)

        # Roll split_row to be centre row
        image_rc = np.roll(hdu[ext].data[:rows, :], -offset, axis=0)
        image_rc[padbins] = 0.0

        # Split columns equally
        hdu[ext].data = image_rc.reshape((2, int(rows / 2), cols))

        return hdu

    # MARK: Crop files
    def crop_file(
        self,
        hdulist: pyfits.HDUList,
        crop: int = CROP_DEFAULT
    ) -> tuple[np.ndarray]:
        """
        Crop the data with respect to the `O`/`E` beam.

        Parameters
        ----------
        hdulist : astropy.io.fits.HDUList
            The HDUList containing the data to be cropped.
        crop : int, optional
            The number of rows to be cropped from the bottom and top 
            of the `O` and `E` beam, respectively.
            (Defaults to 40)

        Returns
        -------
        tuple[numpy.ndarray]
            Tuple containing the cropped O and E beam data arrays.

        """
        o_data = hdulist["SCI"].data[1, 0:-crop]
        e_data = hdulist["SCI"].data[0, crop:]

        return o_data, e_data

    # MARK: Update beam lists
    def update_beam_lists(
        self,
        o_name,
        e_name,
        arc: bool = True
    ) -> None:
        """
        Update the `o_files` and `e_files` attributes.

        Parameters
        ----------
        o_name : str
            The filename of the O beam.
        e_name : str
            The filename of the E beam.
        arc : bool, optional
            Indicates whether the first entry should be the arc frame.
            (Defaults to True)

        Returns
        -------
        None

        """
        if arc:
            self.o_files.insert(0, o_name)
            self.e_files.insert(0, e_name)
        else:
            self.o_files.append(o_name)
            self.e_files.append(e_name)

        return

    # MARK: Save beam lists
    def save_beam_lists(self, file_suffix: str = 'frames') -> None:
        with open(f"o_{file_suffix}", "w+") as f_o, \
             open(f"e_{file_suffix}", "w+") as f_e:
            for i, j in zip(self.o_files, self.e_files):
                f_o.write(i + "\n")
                f_e.write(j + "\n")

        return

    # MARK: Process all Listed Images
    def process(self) -> None:
        """
        Process all FITS images stored in the `fits_list` attribute
        
        Returns
        -------
        None
        
        """
        for target in self.fits_list:
            logging.debug(f"Processing {target}")
            self.split_file(target)

        self.save_beam_lists()

        return

# MARK: Main function
def main(argv) -> None:
    """Main function."""
    
    return


if __name__ == "__main__":
    main(sys.argv[1:])
