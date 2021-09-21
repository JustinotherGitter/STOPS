#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__version__ = "20.09.2021"
__email__ = "justin.jb78@gmail.com"

import sys
import os

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
    def __init__():
        pass


def main(argv): # TODO@JustinotherGitter: Handle Split.py called directly
    return

if __name__ == "__main__":
    main(sys.argv[1:])
