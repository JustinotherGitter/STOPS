#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__email__ = "justin.jb78+Masters@gmail.com"

import os
import sys
import logging
from typing import List
from copy import deepcopy

import numpy as np
from astropy.io import fits as pyfits

from utils.SharedUtils import get_files, get_arc

class Split:
    """
        Split class allows for the seperate call of splitting the polsalt FITS files
        (after basic reductions)

        Parameters
        ----------
        data_dir : str
            The path to the data (mxgbp*.fits files) to be split
        fits_list : list of type FITS, optional
            A list of pre-reduced polsalt FITS files to be split.
            (The default is None, meaning Split will search for files to split in the data directory)
        split_row : int, optional
            The row that the data will be split along.
            (The default is 517, the middle row of the CCD's)
        no_arc : bool, optional
            Decides whether the arc frames should be recombined.
            (The default is False, since polsalt only uses the arc frames until spectral extraction)
        save_prefix : dict of lists of str, optional
            The prefix that the O & E beams are saved as.
            Setting save_prefix = None does not save the split O & E beams
            (The default is {'beam': ["obeam", "ebeam"], 'arc': ["oarc", "earc"]})
        
        Returns
        -------
        split_FITS : list of sets of O & E beam FITS
            A list of sets of FITS files that were split and can be returned to IRAF or python.

        Raises
        ------
            # TODO@JustinotherGitter : Complete docs for which errors are raised and when
    """
    def __init__(self,
                data_dir : str,
                fits_list : list = None, # TODO@JustinotherGitter: Add Lists inner type
                split_row : int = 517,
                no_arc : bool = False,
                save_prefix = None,
                **kwargs
                ) -> None:
        self.data_dir = data_dir
        self.fits_list = get_files(data_dir=data_dir, filenames=fits_list, prefix="mxgbp", extention="fits")
        self.split_row = split_row # TODO@JustinotherGitter: Check valid split and set default to rows // 2 instead of 517
        self.save_prefix = {'beam': ["obeam", "ebeam"], 'arc': ["oarc", "earc"]}
        if type(save_prefix) == dict:
            self.save_prefix = save_prefix # TODO@JustinotherGitter: Check valid list

        self.arc = get_arc(self.fits_list, no_arc)
        self.o_files = []
        self.e_files = []
        return


    def split_file(self, file: str) -> None: # TODO@JustinotherGitter: replace typing return from None to correct type
        # Create empty HDUList
        O_beam = pyfits.HDUList()
        E_beam = pyfits.HDUList()

        # Open file and split O & E beams
        with pyfits.open(file) as hdul:
            O_beam.append(hdul['PRIMARY'].copy())
            E_beam.append(hdul['PRIMARY'].copy())

            # Split specific extention
            raw_split = self.split_ext(hdul, 'SCI')

            # O_beam[0].data = raw_split['SCI'].data[1]
            # E_beam[0].data = raw_split['SCI'].data[0]
            O_beam[0].data, E_beam[0].data = self.crop_file(raw_split)

            # Handle prefix and names
            pref = 'arc' if file == self.arc else 'beam'
            o_name = self.save_prefix[pref][0] + file[-9:]
            e_name = self.save_prefix[pref][1] + file[-9:]
            
            # Add split data to O & E beam lists
            self.update_beam_lists(o_name, e_name, pref == 'arc')

            # Handle don't save case
            if self.save_prefix == None:
                return O_beam, E_beam

            # Handle save case
            O_beam.writeto(o_name, overwrite=True)
            E_beam.writeto(e_name, overwrite=True)

            return O_beam, E_beam


    def split_ext(self, hdulist, ext="SCI") -> np.ndarray: # TODO@JustinotherGitter: Check return is array
        hdu = deepcopy(hdulist)
        rows, cols = hdu[ext].data.shape

        # if odd number of rows, strip off the last one
        rows = int(rows/2) * 2

        # how far split is from center of detector
        offset = int(self.split_row - rows/2)

        # split arc into o/e images
        padbins = (np.indices((rows,cols))[0]<offset) | (np.indices((rows,cols))[0]>rows+offset)

        # Roll split_row to be centre row
        image_rc = np.roll(hdu[ext].data[:rows,:],-offset,axis=0)
        image_rc[padbins] = 0.
        
        # Split columns equally
        hdu[ext].data = image_rc.reshape((2, int(rows/2), cols))

        return hdu


    def crop_file(self, hdulist, crop: int=40) -> None: # TODO@JustinotherGitter: Return type and handle default crop better
        o_data = hdulist['SCI'].data[1, 0:-crop]
        e_data = hdulist['SCI'].data[0, crop:]
        
        return o_data, e_data
    
    
    def update_beam_lists(self, o_name, e_name, arc=True) -> None: # TODO@JustinotherGitter: Add return?
        if arc:
            self.o_files.insert(0, o_name)
            self.e_files.insert(0, e_name)
        else:
            self.o_files.append(o_name)
            self.e_files.append(e_name)
            
        return
    

    def save_beam_lists(self) -> None:
        with open("o_frames", "w+") as f_o:
            for i in self.o_files:
                f_o.write(i + "\n")
                
        with open("e_frames", "w+") as f_e:
            for i in self.e_files:
                f_e.write(i + "\n")

        return


    def process(self) -> None:
        logging.debug(f"Processing the following files: {self.fits_list}")
        for target in self.fits_list:
            self.split_file(target)
        
        self.save_beam_lists()
        return


def main(argv) -> None: # TODO@JustinotherGitter: Handle split.py called directly
    return

if __name__ == "__main__":
    main(sys.argv[1:])
