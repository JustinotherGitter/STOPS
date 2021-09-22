#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__version__ = "20.09.2021"
__email__ = "justin.jb78@gmail.com"

import sys
import os

from typing import List

import astropy.io as pyfits

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

        # Handle inclusion of solution
        ws = []
        for fl in os.listdir(os.path.join(self.data_dir, self.database)):
            if os.path.isfile(os.path.join(self.data_dir, self.database, fl)) and ("fc" == fl[0:2]):
                ws.append(fl)

        if len(ws) > 0:
            return ws
        else:
            # Handle no solutions found
            raise FileNotFoundError(f"No wavelength solution (fc...) found in the solution directory {os.path.join(self.data_dir, self.database)}")

    
    def get_arc(self, exclude_arc: bool) -> str:
        # Handle exclusion of arc
        if exclude_arc:
            # TODO@JustinotherGitter: Check whether to remove arc here or not.
            return ''

        # Handle inclusion of arc
        for fl in self.fits_list:
            with pyfits.open(fl) as hdu:
                if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                    return fl

        # Handle arc not found
        raise FileNotFoundError(f"No arc found in the data directory {self.data_dir}")


















def main(argv): # TODO@JustinotherGitter: Handle Split.py called directly
    return

if __name__ == "__main__":
    main(sys.argv[1:])
