#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "2024.03.01"
__author__ = "Justin Cooper"
__email__ = "justin.jb78+Masters@gmail.com"

# General Imports
import os
import argparse
import logging

import split
import join
import cross_correlate
import skylines

from utils import ParserUtils as pu


PROG = "STOPS"
DESCRIPTION = """
Supplementary TOols for Polsalt Spectropolarimetry (STOPS) is a
collection of supplementary tools created for SALT's POLSALT pipeline,
allowing for wavelength calibrations with IRAF. The tools provide
support for slitting and joining polsalt formatted data as well as
cross correlating complementary polarimetric beams.

DOI: 10.22323/1.401.0056
"""
# Scripts created for and as part of Master thesis (2024).

SPLITROW = 517
PREFIX = ["obeam", "ebeam"]

# TODO@JustinotherGitter: Add plot return options?
# TODO@JustinotherGitter: Add type of file check (I.E. FITS)
# TODO@JustinotherGitter: Implement saving log to --log (logging in all files)


# Universal parser
parser = argparse.ArgumentParser(
    prog=PROG,
    description=DESCRIPTION,
    formatter_class=argparse.RawDescriptionHelpFormatter
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
parser.add_argument(
    "data_dir",
    action="store",
    nargs="?",
    default=os.getcwd(),
    type=pu.parse_path,
    help="Path to directory containing data. Defaults to the CWD.",
)


# Split / Join specific arguments
split_join_args = argparse.ArgumentParser(add_help=False)
split_join_args.add_argument(
    "-n",
    "--no_arc",
    action="store_true",
    help="Flag to exclude arc files from processing.",
)
split_join_args.add_argument(
    "-s",
    "--split_row",
    default=SPLITROW,
    type=int,
    help=f"Row along which to split the O and E beams at. Defaults to {SPLITROW}.",
)
split_join_args.add_argument(
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
    "split", aliases=["s"], help="Split mode", parents=[split_join_args]
)
# 'children' split args here
# Change defaults here
split_parser.set_defaults(mode="split", func=split.Split)


# Join subparser mode
join_parser = subparsers.add_parser(
    "join", aliases=["j"], help="Join mode", parents=[split_join_args]
)
# 'children' join args here
# Change defaults here
join_parser.set_defaults(mode="join", func=join.Join)


# Correlate subparser mode
corr_parser = subparsers.add_parser(
    "correlate", aliases=["x"], help="Cross correlation mode"
)
# 'children' correlate args here
corr_parser.add_argument(
    "filenames",
    action="store",
    # dest="fits_list",
    nargs="+",
    type=pu.parse_corr_file,
    help="Filenames to be compared. Provide only one filename for O/E beam comparisons. Include relative path to working directory.",
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
    dest="cont_ord",
    help="Order of continuum to remove from spectra. Higher orders recommended to remove most variation, leaving only significant features.",
)
corr_parser.add_argument(
    "-p",
    "--continuum_plot",
    action="store_true",
    dest="cont_plot",
    help="Flag to plot fitting of continuum. Used to confirm only notable features left in spectrum.",
)
corr_parser.add_argument(
    "-o",
    "--offset",
    type=int,
    default=0,
    help="Offset introduction when correcting for known offset in spectra or for testing purposes. (For testing, not used during regular operation.)",
)
corr_parser.add_argument(
    "-s",
    "--save_name",
    help="Name to save plot with. If left undefined, plot will not be saved.",
)
# Change defaults here
corr_parser.set_defaults(mode="correlate", func=cross_correlate.CrossCorrelate)


# Skyline subparser mode
sky_parser = subparsers.add_parser(
    "skylines", aliases=["sky"], help="Sky line check mode"
)
# 'children' skyline args here
sky_parser.add_argument(
    "filenames",
    action="store",
    nargs="+",
    type=pu.parse_file,
    help="File name(s) of FITS files to be checked using SALT sky atlas. At least one filename expected.",
)
sky_parser.add_argument(
    "-s",
    "--save_prefix",
    action="store",
    nargs="?",
    # default=None,
    const="sky",
    help="Prefix used when saving plot. Excluding option does not save, flag usage of option uses 'sky' default prefix, and a provided prefix overwrites default prefix.",
)
sky_parser.add_argument(
    "-b",
    "--beams",
    choices=["o", "e", "oe"],
    type=str.lower,
    default="both",
    help="Beam(s) for skyline checking. Defaults to both beams overplotted, but may be given either 'o', 'e', or 'oe' for separating and excluding beam plots.",
)
sky_parser.add_argument(
    "-t",
    "--transform",
    action="store_true",
    help="Flag to NOT transform image. Defaults to true. Recommended to only use when input image(s) already transformed.",
)
sky_parser.set_defaults(mode="skyline", func=skylines.Skylines)


# Parse 'mode' and arguments + any keyword clean up
args = parser.parse_args()
args.verbose = pu.parse_loglevel(args.verbose)
args.log = args.data_dir / args.log
if "filenames" in args:
    args.filenames = pu.flatten(args.filenames)

# Begin logging
logging.basicConfig(
    filename=args.log,
    format='%(asctime)s - %(module)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=args.verbose
)


# Run mode using arguments
logging.debug(f"Argparse namespace: {args}")
logging.info(f"Mode:{args.mode}")
# args.func(**vars(args)).process()


# Confirm all processes completed and exit without error
logging.info("All done! Come again!\n")
