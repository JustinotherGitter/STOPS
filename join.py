#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__version__ = "20.09.2021"
__email__ = "justin.jb78@gmail.com"

import os
import sys
from typing import List
from copy import deepcopy
import re
import numpy as np
from numpy.polynomial.chebyshev import chebgrid2d as chebgrid2d
import astropy.io as pyfits
from specpolpy3 import read_wollaston
import lacosmic

DATADIR = os.path.expanduser("~/polsalt-beta/polsalt/data/")

class Join:
    """
        Join class allows for the seperate call of joining the wavelength calibrated O & E beam FITS files

        Parameters
        ----------
        path_name : str
            The path to the data (wmxgbp*.fits files) to be joined
        split_row : int, optional
            The row that the data was split along.
            (The default is 517, the middle row of the CCD's)
        verbose : bool, optional
            Decides whether the output should be recorded to the terminal window.
            (The default is False, no output written to the terminal window)
        no_arc : bool, optional
            Decides whether the arc frames should be recombined.
            (The default is False, since polsalt only uses the arc frames until spectral extraction)
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
    def __init__(self,
                data_dir : str,
                database : str = None,
                fits_list : List[str] = None,
                solutions_list : List[str] = None,
                split_row : int = 517,
                verbose : bool = False,
                no_arc : bool = False,
                save_prefix : dict[str, List[str]] = {'beam': ["obeam", "ebeam"], 'arc': ["oarc", "earc"]}
                ) -> None:
        self.data_dir = data_dir # TODO@JustinotherGitter: Check valid path
        self.database = "database" if database == None else database
        self.fits_list = self.get_files(fits_list)
        self.fc_files = self.get_solutions(solutions_list)
        self.split_row = split_row # TODO@JustinotherGitter: Check valid split and set default to rows / 2 instead of 517
        self.verbose = verbose
        self.save_prefix = save_prefix # TODO@JustinotherGitter: Check valid list

        self.no_arc = no_arc
        self.arc = self.get_arc(no_arc)
        return


    def get_files(self, flist: List, prefix: str="m", extention: str="fits") -> List: # TODO@JustinotherGitter: Add Lists inner type
        # Handle recieving list of files
        if flist != None:
            for fl in flist:
                if os.path.isfile(os.path.join(self.data_dir, fl)):
                    continue
                else:
                    raise FileNotFoundError(f"{fl} not found in the data directory {self.data_dir}")
            return flist

        # Handle finding valid files
        flist = []
        for fl in os.listdir(self.data_dir):
            if os.path.isfile(os.path.join(self.data_dir, fl)) and (prefix == fl[0]) and (extention == fl.split(".")[-1]):
                flist.append(fl)
        return flist

    
    def get_solutions(self, wavlist: List, prefix: str="fc"):
        # Handle recieving list of solutions
        if wavlist != None:
            for fl in wavlist:
                if os.path.isfile(os.path.join(self.data_dir, self.database, fl)):
                    continue
                else:
                    raise FileNotFoundError(f"{fl} not found in the data directory {self.data_dir}")
            return wavlist

        # Handle finding solution
        ws = [] # TODO@JustinotherGitter: Check order of solutions
        for fl in os.listdir(os.path.join(self.data_dir, self.database)):
            if os.path.isfile(os.path.join(self.data_dir, self.database, fl)) and (prefix == fl[0:2]):
                ws.append(fl)

        if len(ws) > 0:
            return ws
        else:
            # Handle no solutions found
            raise FileNotFoundError(f"No wavelength solution (fc...) found in the solution directory {os.path.join(self.data_dir, self.database)}")

    
    def get_arc(self, exclude_arc: bool) -> str:

        # Handle finding of arc
        for fl in self.fits_list:
            with pyfits.open(fl) as hdu:
                if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                    return fl

        # Handle arc not found
        raise FileNotFoundError(f"No arc found in the data directory {self.data_dir}")


    def join_file(self, file: str):
        # Handle exclusion of arc
        if self.no_arc:
            return
        
        # Handle prefix and names
        pref = 'arc' if file == self.arc else 'beam'
        o_file = self.save_prefix[pref][0] + file[-9:]
        #e_file = self.save_prefix[pref][1] + file[-9:]

        # Open file
        with pyfits.open(file) as hdu:
            # Check if file has been cropped
            cropsize = False
            with pyfits.open(o_file) as o:
                if hdu["SCI"].data.shape[0] / 2 != o[0].data.shape[0]:
                    cropsize = int(hdu["SCI"].data.shape[0] / 2 - o[0].data.shape[0])
            
            y_shape = int(hdu["SCI"].data.shape[0] / 2) - cropsize
            x_shape = hdu["SCI"].data.shape[1]

            # Create wavelength appended hdu
            whdu = pyfits.HDUList()
            #No differences in "PRIMARY" extention header
            whdu.append(hdu["PRIMARY"])

            for ext in ["SCI", "VAR", "BPM"]:
                whdu.append(pyfits.ImageHDU(name=ext))
                whdu[ext].header = deepcopy(hdu[ext].header)
                whdu[ext].header["CTYPE3"] = "O,E"
                
                if ext == "BPM":
                    whdu[ext].data = np.zeros((2, y_shape, x_shape), dtype='uint8')
                    whdu[ext].header["BITPIX"] = "-uint8"
                else:
                    whdu[ext].data = np.zeros((2, y_shape, x_shape), dtype='>f4')
                    whdu[ext].header["BITPIX"] = "-32"

        whdu.append(pyfits.ImageHDU(name="WAV"))
        wav_header = deepcopy(whdu["SCI"].header)
        wav_header["EXTNAME"] = "WAV"
        wav_header["CTYPE3"] = "O,E"
        whdu["WAV"].header = wav_header
        
        whdu["WAV"].data = np.zeros(whdu["SCI"].data.shape, dtype='>f4')

        for num, fname in enumerate(self.fc_files):

            chebvals = []
            with open(self.database + "/" + fname) as file: # TODO@JustinotherGitter: Check order of fc files here
                for i in file:
                    # TODO: check regex substitution correct
                    chebvals.append(re.sub(r"[\n\t\s]*", "", i))

                    # def stripchars(text):
                    #     if text[0] in ["\t", "\n"]:
                    #         text = text[1:]
                    #     if text[0] in ["\t", "\n"]:
                    #         text = text[1:]
                    #     if text[-1] in ["\t", "\n"]:
                    #         text = text[:-1]
                    #     return text

                    # chebvals.append(stripchars(i))

            if chebvals[6] == "1.": #function - Function type (1=chebyshev, 2=legendre)
                x_ord = int(chebvals[7][:-1]) #xorder - X "order" (highest power of x)
                y_ord = int(chebvals[8][:-1]) #yorder - Y "order" (highest power of y)
                if chebvals[9] == "1.": #xterms - Cross-term type (always 1 for FITCOORDS)
                    xmin = int(float(chebvals[10][:-1])) #xmin - Minimum x over which the fit is defined
                    xmax = int(float(chebvals[11][:-1])) #xmax - Maximum x over which the fit is defined
                    ymin = int(float(chebvals[12][:-1])) #ymin - Minimum y over which the fit is defined
                    ymax = int(float(chebvals[13][:-1])) #ymax - Maximum y over which the fit is defined
                    
                    if ymax != y_shape: # TODO: Fix temporary stretching
                        ymax = y_shape
                        
                c_vals = np.array(chebvals[14:], dtype=float)
                c_vals = np.reshape(c_vals, (x_ord, y_ord))
            
            
                # Set wavelength extention values to function
                whdu["WAV"].data[num] = chebgrid2d(x=np.linspace(-1, 1, ymax),
                                                    y=np.linspace(-1, 1, xmax),
                                                    c=c_vals)

            elif chebvals[6] == "2.":
                # TODO@JustinotherGitter: Handle legendre
                raise NotImplementedError("Legendre functions not yet supported")

            else:
                #TODO@JustinotherGitter: Handle other functions?
                raise Exception("Function type not recognised, please wavelength calibrate using chebychev.")
            
            # Cosmic Ray Cleaning
            # TODO@JustinotherGitter: Hard parameters set. Not good for universal application
            whdu["SCI"].data[num] = lacosmic.lacosmic(whdu["SCI"].data[num], 2, 4, 4, effective_gain=1, readnoise=3.3)[0]
    
        # WAV mask (Left & Right Crop)
        whdu["WAV"].data[whdu["WAV"].data[:] <  3_000] = 0.0
        whdu["WAV"].data[whdu["WAV"].data[:] >= 10_000] = 0.0
            
        # Correct WAV mask shift (Top & Bottom Crop)
        rpix_oc, cols, rbin, lam_c = read_wollaston(whdu, DATADIR + 'wollaston.txt')

        drow_oc = (rpix_oc-rpix_oc[:,int(cols/2)][:,None])/rbin

        ## Cropping as suggested
        for c, col in enumerate(drow_oc[0]):
            if not np.isnan(col):
                if int(col) < 0:
                    whdu["WAV"].data[0, int(col):, c] = 0.0
                elif int(col) > cropsize:
                    whdu["WAV"].data[0, 0:int(col) - cropsize, c] = 0.0
                    
        for c, col in enumerate(drow_oc[1]):
            if not np.isnan(col):
                if int(col) > 0:
                    whdu["WAV"].data[1, 0:int(col), c] = 0.0
                elif (int(col) < 0) & (abs(int(col)) > cropsize):
                    whdu["WAV"].data[1, int(col) + cropsize:, c] = 0.0

        # Mask BPM same as WAV
        whdu["BPM"].data[0] = np.where(whdu["WAV"].data[0] == 0, 1, whdu["BPM"].data[0])
        whdu["BPM"].data[1] = np.where(whdu["WAV"].data[1] == 0, 1, whdu["BPM"].data[1])
        

        whdu.writeto("w" + file, overwrite="True")
    
    def process(self) -> None:
        for target in self.fits_list:
            self.join_file(target)
        
        # TODO@JustinotherGitter: Save here
        return


def main(argv): # TODO@JustinotherGitter: Handle Join.py called directly
    return

if __name__ == "__main__":
    main(sys.argv[1:])
