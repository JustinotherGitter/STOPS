#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "2024.02.04"
__author__ = "Justin Cooper"
__email__ = "justin.jb78@gmail.com"

# General Imports
import os
import sys
import argparse
import logging

import split
import join
import cross_correlate

from utils import ParserUtils as pu


PROG = (
    "Supplementary TOols for Polsalt Spectropolarimetry (STOPS)"
    + " OR "
    + "Spectropolarimetric Wavelength calibration Alternative for Polsalt (SWAP)"
)
DESCRIPTION = """
Supplementary tools created for SALT's POLSALT pipeline, allowing for wavelength calibrations with IRAF.
Additional tools provide support for cross correlating complementary polarimetric beams.
Scripts created for and as part of Master thesis (2024).

DOI: 10.22323/1.401.0056
"""

SPLITROW = 517
PREFIX = ["obeam", "ebeam"]

# TODO@JustinotherGitter: Add plot return options?
# TODO@JustinotherGitter: Add type of file check (I.E. FITS)
# TODO@JustinotherGitter: Implement saving log to --log (logging in all files)


# Universal parser
parser = argparse.ArgumentParser(
    description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "-V",
    "--version",
    action="version",
    version=f"%(prog)s as of {__version__}",
)
parser.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="Enables verbose mode. Use -v or -vv for greater verbosity levels.",
)
parser.add_argument(
    "-l",
    "--log",
    action="store",
    type=pu.parse_logfile,
    help="Filename to save logging to. Defaults to None.",
)


# Parent parser template
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument(
    "data_dir",
    action="store",
    nargs="?",
    default=os.getcwd(),
    type=pu.parse_path,
    help="Path to directory containing data. Defaults to the CWD.",
)
parent_parser.add_argument(
    "-n",
    "--no_arc",
    action="store_true",
    help="Flag to exclude arc files from processing.",
)
parent_parser.add_argument(
    "-s",
    "--split_row",
    default=SPLITROW,
    type=int,
    help=f"Row along which to split the O and E beams at. Defaults to {SPLITROW}.",
)
parent_parser.add_argument(
    "-p",
    "--save_prefix",
    nargs=2,
    default=PREFIX,
    help=f"Prefix with which to save the O and E beams. Defaults to {PREFIX}.",
)


# Create subparser modes
subparsers = parser.add_subparsers(
    dest="mode", help="Operational mode of supplementary tools"
)


# Split subparser mode
split_parser = subparsers.add_parser(
    "split", aliases=["s"], help="Split mode", parents=[parent_parser]
)
split_parser.set_defaults(mode="split", func=split.Split)


# Join subparser mode
join_parser = subparsers.add_parser(
    "join", aliases=["j"], help="Join mode", parents=[parent_parser]
)
join_parser.set_defaults(mode="join", func=join.Join)


# Correlate subparser mode
corr_parser = subparsers.add_parser(
    "correlate", aliases=["x"], help="Cross correlation mode"
)
corr_parser.add_argument(
    "filenames",
    action="store",
    nargs="+",
    type=pu.parse_file,
    help="Filenames to be compared. Provide only one filename for O/E beam comparisons.",
)
corr_parser.add_argument(
    "-ccd",
    "--split_ccd",
    action="store_false",
    help="Flag to NOT split CCD's. Recommended to leave off unless data cleaned of chip gaps.",
)
corr_parser.add_argument(
    "-c",
    "--continuum_order",
    type=int,
    default=11,
    help="Order of continuum to remove from spectra. Higher orders recommended to remove most variation, leaving only significant features.",
)
corr_parser.add_argument(
    "-p",
    "--continuum_plot",
    action="store_true",
    help="Flag to plot fitting of continuum. Used to confirm only notable features left in spectrum.",
)
corr_parser.add_argument(
    "-o",
    "--offset",
    type=int,
    default=0,
    help="Offset introduction when correcting for known offset in spectra or for testing purposes. (For testing, not used in regular operation.) ",
)
corr_parser.add_argument(
    "-s",
    "--save_name",
    help="Name to save plot with. If left undefined,d plot will not be saved.",
)
corr_parser.set_defaults(mode="correlate", func=cross_correlate.CrossCorrelate)


# Parse mode and arguments + and keyword clean up
args = parser.parse_args()
args.verbose = pu.parse_loglevel(args.verbose)

# Begin logging
logfile = pu.parse_logfile(args.log)
# logFormatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
logging.basicConfig(filename=logfile, level=args.verbose)
# format=logFormatter
# handlers=[logging.StreamHandler(sys.stdout)]

# Run mode using arguments
logging.debug(f"Argparse namespace: {args}")
logging.info(f"Mode:{args.mode}")
args.func(args).process()

# Confirm all processes completed
logging.info("All done! Come again!\n")
