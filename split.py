#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__version__ = "20.09.2021"
__email__ = "justin.jb78@gmail.com"

import os
import sys
from typing import List
from copy import deepcopy
import numpy as np
from astropy.io import fits as pyfits

class Split:
    """
        Split class allows for the seperate call of splitting the polsalt FITS files
        (after basic reductions)

        Parameters
        ----------
        data_dir : str
            The path to the data (wmxgbp*.fits files) to be joined
        fits_list : list of type FITS, optional
            A list of pre-reduced polsalt FITS files to be split.
            (The default is None, meaning Split will search for files to split in the data directory)
        split_row : int, optional
            The row that the data was split along.
            (The default is 517, the middle row of the CCD's)
        verbose : bool, optional
            Decides whether the output should be recorded to the terminal window.
            (The default is False, no output written to the terminal window)
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
                fits_list : List = None, # TODO@JustinotherGitter: Add Lists inner type
                split_row : int = 517,
                verbose : bool = False,
                no_arc : bool = False,
                save_prefix : dict[str, List[str]] = {'beam': ["obeam", "ebeam"], 'arc': ["oarc", "earc"]}
                ) -> None:
        self.data_dir = data_dir # TODO@JustinotherGitter: Check valid path
        self.fits_list = self.get_files(fits_list)
        self.split_row = split_row # TODO@JustinotherGitter: Check valid split and set default to rows / 2 instead of 517
        self.verbose = verbose
        self.save_prefix = save_prefix # TODO@JustinotherGitter: Check valid list

        self.arc = self.get_arc(no_arc)
        self.o_files = []
        self.e_files = []
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

    
    def get_arc(self, exclude_arc: bool) -> str:
        # Handle exclusion of arc
        if exclude_arc:
            return ''

        # Handle inclusion of arc
        for fl in self.fits_list:
            with pyfits.open(fl) as hdu:
                if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                    return fl

        # Handle arc not found
        raise FileNotFoundError(f"No arc found in the data directory {self.data_dir}")

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
            e_name = self.save_prefix[pref][0] + file[-9:]
            
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


    def crop_file(hdulist, crop: int=40) -> None: # TODO@JustinotherGitter: Fix return type and Handle default crop better
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


def main(): # TODO@JustinotherGitter: Handle Split.py called directly
    return

if __name__ == "__main__":
    main()

########################
# TODO@JustinotherGitter: Finish refactoring/rewriting and remove below
# not to be included in any distributions
def split(pathname, splitrow=517, display=False, no_arc=False, save_pref=["obeam", "ebeam"]):
    if display: print(f"Path to data given as: {pathname}")
    
    # Get ARC and DATA FITS files from path given
    arcfile = ""
    infilelist = []
    o_files = []
    e_files = []
    
    # Get all valid FITS files
    for fl in os.listdir(pathname):
        if os.path.isfile(os.path.join(pathname, fl)) and ("m" == fl[0]) and ("fits" == fl.split(".")[-1]):
            #print(fl)
            infilelist.append(fl)
    
    # TODO: Allow specification of ARC file in case multiple ARCS
    # Get arcfile from valid FITS files
    for fl in infilelist:
        with pyfits.open(fl) as hdu:
            if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                arcfile = fl
    
    if (arcfile == "") & (no_arc == False):
        cont = input("No arc found. Continue without it? [Y/n]\t")
        if cont.lower() in ("y", "yes", ""):
            no_arc = True
        elif cont.lower() in ("n", "no"):
            print("Quitting wav_replacement.py")
            sys.exit()
            
    if no_arc:
        arcfile = ""
    else:
        if display: print("Arcfile:\t\t\t" + arcfile)
        
        arcO, arcE = outfiles(arcfile, splitrow, save_O="oarc" + arcfile[-9:], save_E="earc" + arcfile[-9:], display=display)
        
        o_files.append("oarc" + arcfile[-9:])
        e_files.append("earc" + arcfile[-9:])
        if display: print(f"{arcfile} split into {o_files[-1]} & {e_files[-1]}")
        
    # Load existing FITS file and create copied single-extention FITS file
    if display: print("Pipeline m*.fits files found:\t" + str(infilelist))
    
    tar_list = []
    for i in infilelist:
        # Create blank FITS file for wavelength calibration and transformation
        if i != arcfile:
            tar_list.append(outfiles(i, splitrow, save_O=save_pref[0] + i[-9:], save_E=save_pref[1] + i[-9:], display=display))
            
            o_files.append(save_pref[0] + i[-9:])
            e_files.append(save_pref[1] + i[-9:])
        
        if display: print(f"{i} split into {o_files[-1]} & {e_files[-1]}")
    
    
    ####################################
            
    ## Crop fits files for reidentify to work
    # TODO: Find crop size dynamically / user given
    cut_size = 40
    
    
    # Obeam
    for i in o_files:
        with pyfits.open(i) as hdu:
            #print(hdu.info())
            #fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=[18, 8])
            #ax1.imshow(hdu[0].data, vmax=np.std(hdu[0].data), origin="lower")
            #ax2.imshow(hdu[0].data[0:-cut_size], vmax=np.std(hdu[0].data), origin="lower")
    
            hdu[0].data = hdu[0].data[0:-cut_size]
            hdu.writeto(i, overwrite=True)
    
    #Ebeam
    for i in e_files:
        with pyfits.open(i) as hdu:
            #print(hdu.info())
            #fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=[18, 8])
            #ax1.imshow(hdu[0].data, vmax=np.std(hdu[0].data), origin="lower")
            #ax2.imshow(hdu[0].data[0:-cut_size], vmax=np.std(hdu[0].data), origin="lower")
    
            hdu[0].data = hdu[0].data[cut_size:]
            hdu.writeto(i, overwrite=True)

    # Create text files listing o and e beams
    with open("o_frames", "w+") as f_o:
        for i in o_files:
            f_o.write(i + "\n")
            
    with open("e_frames", "w+") as f_e:
        for i in e_files:
            f_e.write(i + "\n")
    
    if display: print("Wrote o and e frames to o_frames and e_frames, respectfully.")
    print("Finished with no errors raised.")