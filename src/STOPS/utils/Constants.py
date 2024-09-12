"""File tracking constants for consistency across modules."""

import os


"""The CCD row along which the $O$- and $E$-beams are split."""
SPLIT_ROW = 517

# Note order of save prefixes: O first, E second. Always.
"""The prefix with which the split $O$- and $E$-beam science files are saved."""
PREFIX = ["beamo", "beame"]

"""The prefix with which the split $O$- and $E$-beam FITS files are saved."""
SAVE_PREFIX = {"beam": PREFIX, "arc": ["arco", "arce"]}

"""The prefix with which the correlation plots are saved."""
SAVE_CORR = "corr"

"""The prefix with which the skyline plots are saved."""
SAVE_SKY = "sky"


"""The `polsalt` directory containing the calibration data for the Wollaston prism."""
DATADIR = os.path.expanduser("~/polsalt-beta/polsalt/data/")

"""The default directory to search for data files. DEPRECATED."""
DEFAULT_DIR = os.getcwd()


"""The cropping applied during `split`, based on best cropping of `PG0300` with `Ar` arc lamp."""
CROP_DEFAULT = 40


"""The default parameters for the CLI parser."""
PARSE = {
    "VERBOSE": 0,  # 0, 1, 2 see ParserUtils.parse_loglevel
    "DATA_DIR": DEFAULT_DIR,
    "CONT_ORD": 11,
    "OFFSET": 0,
    "BEAMS": "OE",  # 'O', 'E', 'OE'
}

# ccdproc lacosmic needs testing as parameter names differ
# https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/08-03-Cosmic-ray-removal.html
"""
The default parameters for cosmic ray cleaning.
Note deprecations arising from ccdproc implementation of the lacosmic algorithm.
Gain and readnoise sourced from https://pysalt.salt.ac.za/proposal_calls/current/ProposalCall.html
"""
CR_PARAMS = {
    "READNOISE": 3.3,
    "GAIN": 1,
    "CR_CONTRAST": 2,  # Deprecated
    "CR_THRESHOLD": 4,  # Deprecated
    "CR_NEIGHBOUR_THRESHOLD": 4,  # Deprecated
    "BACKGROUND": None,
}

"""The default parameters for the `skylines` `find_peaks` function."""
FIND_PEAK_PARAMS = {
    'rel_height': 0.5,
    'min_height': 0.01,
    'distance': 2,
}

"""The default `skylines` Arc lamp."""
ARC_FILE = 'Argon_lores.txt' # 'NeAr.txt' | 'ThAr.txt' | 'Xe.txt'


"""The vertical offset of the spectra in the `correlate` plot."""
OFFSET = 0.3
