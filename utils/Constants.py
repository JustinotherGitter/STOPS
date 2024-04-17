#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

"""File to keep constants consistent across modules."""

# Directory to `polsalt` calibration data,
# specifically the wollaston correction calibration data
DATADIR = os.path.expanduser("~/polsalt-beta/polsalt/data/")

SPLIT_ROW = 517 # Half of the CCD pixel `Height`
CROP_DEFAULT = 40 # Best crop from usage, PG0300 with Ar arc lamp

# Naming convention for FITS prefixes, Note order: O first, E second. Always.
PREFIX = ["obeam", "ebeam"]
SAVE_PREFIX = {"beam": PREFIX, "arc": ["oarc", "earc"]}


# CR Cleaning parameters
    # Gain and Readnoise still valid for ccdproc, sourced from
    # https://pysalt.salt.ac.za/proposal_calls/current/ProposalCall.html
    # Constants (Deprecated) work well with base lacosmic.
    # ccdproc lacosmic needs testing as parameter names differ
    # https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/08-03-Cosmic-ray-removal.html
CR_PARAMS = {
    'READNOISE': 3.3,
    'GAIN': 1,
    'CR_CONTRAST': 2, # Deprecated
    'CR_THRESHOLD': 4, # Deprecated
    'CR_NEIGHBOUR_THRESHOLD': 4, # Deprecated
    'BACKGROUND': None,
}
