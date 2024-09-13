#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The Command Line Interface (CLI) for STOPS, handling argument parsing,
class instantiation, processing, and logging.

Usage
-----
`python -m STOPS`

Suffix the above command with `-h` for help on usage.

Description
-----------
This module is the entry point for the STOPS package. It provides
the CLI argument parser for the supplementary tools provided by STOPS.
The parser is built using the `argparse` module and is designed to be
user-friendly.

See Also
--------
IRAF:
    For more information on the IRAF package, see the IRAF website:
    https://iraf-community.github.io/

POLSALT:
    For more information on the POLSALT pipeline, see the POLSALT website:
    https://github.com/saltastro/polsalt

PYSALT:
    For more information on the PYSALT package, see the PYSALT website:
    https://astronomers.salt.ac.za/software/pysalt-documentation/

PYRAF:
    For more information on the PYRAF package, see the PYRAF website:
    https://pyraf.readthedocs.io/en/latest/

"""

# MARK: Imports
import sys
import argparse
import logging
from pathlib import Path

from STOPS import __version__
from STOPS import Split, Join, CrossCorrelate, Skylines
from STOPS.utils import ParserUtils as Parser
from STOPS.utils.Constants import SPLIT_ROW, PREFIX, PARSE, SAVE_CORR, SAVE_SKY


# MARK: Constants
PROG = "STOPS"
DESCRIPTION = """
Supplementary TOols for Polsalt Spectropolarimetry (STOPS) is a
collection of supplementary tools created for SALT's POLSALT pipeline,
allowing for wavelength calibrations with IRAF. The tools provide
support for splitting and joining polsalt formatted data as well as
cross correlating complementary polarimetric beams.

Pre-release Reference:
    DOI: 10.22323/1.401.0056
"""


# MARK: Universal Parser
parser = argparse.ArgumentParser(
    prog=PROG,
    description=DESCRIPTION,
    formatter_class=argparse.RawDescriptionHelpFormatter,
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
    default=PARSE['VERBOSE'],
    help=(
        "Counter flag which enables and increases verbosity. "
        "Use -v or -vv for greater verbosity levels."
    ),
)
parser.add_argument(
    "-l",
    "--log",
    action="store",
    type=Parser.parse_logfile,
    help=(
        "Filename of the logging file. "
        "File is created if it does not exist. Defaults to None."
    ),
)
parser.add_argument(
    "data_dir",
    action="store",
    nargs="?",
    default=PARSE['DATA_DIR'],
    type=Parser.parse_path,
    help=(
        "Path of the directory which contains the working data. "
        f"Defaults to the cwd -> `{PARSE['DATA_DIR']}` (I.E. '.')."
    ),
)


# MARK: Split\Join Parent
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
    default=SPLIT_ROW,
    type=int,
    help=(
        "Row along which the O and E beams are split. "
        f"Defaults to polsalt's default -> {SPLIT_ROW}."
    ),
)
split_join_args.add_argument(
    "-p",
    "--save_prefix",
    nargs=2,
    default=PREFIX,
    help=(
        "Prefix appended to the filenames, "
        "with which the O and E beams are saved. "
        f"Defaults to {PREFIX}."
    ),
)


# MARK: Corr.\Sky. Parent
corr_sky_args = argparse.ArgumentParser(add_help=False)
corr_sky_args.add_argument(
    "filenames",
    action="store",
    nargs="+",
    type=Parser.parse_file,
    help=(
        "File name(s) of FITS file(s) to be processed."
        "A minimum of one filename is required."
    ),
)
corr_sky_args.add_argument(
    "-b",
    "--beams",
    choices=["O", "E", "OE"],
    type=str.upper,
    default=PARSE['BEAMS'],
    help=(
        "Beams to process. "
        f"Defaults to {PARSE['BEAMS']}, but "
        "may be given 'O', 'E', or 'OE' to "
        "determine which beams are processed."
    ),
)
corr_sky_args.add_argument(
    "-ccd",
    "--split_ccd",
    action="store_false",
    help=(
        "Flag to NOT split CCD's. "
        "Recommended to leave off unless the chip gaps "
        "have been removed from the data."
    ),
)
corr_sky_args.add_argument(
    "-c",
    "--continuum_order",
    type=int,
    default=PARSE['CONT_ORD'],
    dest="cont_ord",
    help=(
        "Order of continuum to remove from spectra. "
        "Higher orders recommended to remove most variation, "
        "leaving only significant features."
    ),
)
corr_sky_args.add_argument(
    "-p",
    "--plot",
    action="store_true",
    help="Flag for additional plot outputs.",
)


# MARK: Create subparser modes
subparsers = parser.add_subparsers(
    dest="mode",
    help="Operational mode of supplementary tools",
)


# MARK: Split Subparser
split_parser = subparsers.add_parser(
    "split",
    aliases=["s"],
    help="Split mode",
    parents=[split_join_args],
)
# 'children' split args here
# Change `split` defaults here
split_parser.set_defaults(
    mode="split",
    func=Split,
)


# MARK: Join Subparser
join_parser = subparsers.add_parser(
    "join",
    aliases=["j"],
    help="Join mode",
    parents=[split_join_args],
)
# 'children' join args here
join_parser.add_argument(
    "-c",
    "--coefficients",
    dest="solutions_list",
    nargs='*',
    type=Parser.parse_file,
    help=(
        "Custom coefficients to use instead of the `IRAF` fitcoords "
        "database. Use as either '-c <o_solution> <e_solution>' or "
        "a regex descriptor '-c <*solution*extention>'."
    ),
)
# Change `join` defaults here
join_parser.set_defaults(
    mode="join",
    func=Join,
)


# MARK: Correlate Subparser
corr_parser = subparsers.add_parser(
    "correlate",
    aliases=["x"],
    help="Cross correlation mode",
    parents=[corr_sky_args],
)
# 'children' correlate args here
corr_parser.add_argument(
    "-o",
    "--offset",
    type=int,
    default=PARSE['OFFSET'],
    help=(
        "Introduces an offset when correcting for "
        "known offset in spectra or for testing purposes. "
        f"Defaults to {PARSE['OFFSET']}. "
        "(For testing, not used during regular operation.)"
    ),
)
corr_parser.add_argument(
    "-s",
    "--save_prefix",
    action="store",
    nargs="?",
    type=lambda path: Path(path).expanduser().resolve(),
    const=SAVE_CORR,
    help=(
        "Prefix used when saving plot. "
        "Excluding flag does not save output plot, "
        f"flag usage of option uses default prefix {SAVE_CORR}, "
        "and a provided prefix overwrites default prefix."
    ),
)
# Change `correlate` defaults here
corr_parser.set_defaults(
    mode="correlate",
    func=CrossCorrelate,
)


# MARK: Skyline Subparser
sky_parser = subparsers.add_parser(
    "skylines",
    aliases=["sky"],
    help="Sky line check mode",
    parents=[corr_sky_args],
)
# 'children' skyline args here
sky_parser.add_argument(
    "-t",
    "--transform",
    action="store_false",
    help=(
        "Flag to force transform images. "
        "Recommended to use only when input image(s) "
        "are prefixed 't' but are not yet transformed."
    ),
)
sky_parser.add_argument(
    "-s",
    "--save_prefix",
    action="store",
    nargs="?",
    type=lambda path: Path(path).expanduser().resolve(),
    const=SAVE_SKY,
    help=(
        "Prefix used when saving plot. "
        "Excluding flag does not save output plot, "
        f"flag usage of option uses default prefix {SAVE_SKY}, "
        "and a provided prefix overwrites default prefix."
    ),
)
# Change `skylines` defaults here
sky_parser.set_defaults(
    mode="skyline",
    func=Skylines,
)


# MARK: Keyword Clean Up
args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

args.verbose = Parser.parse_loglevel(args.verbose)

if 'log' in args and args.log not in ["", None]:
    args.log = args.data_dir / args.log

if "filenames" in args:
    args.filenames = Parser.flatten(args.filenames)

if "solutions_list" in args and isinstance(args.solutions_list, list):
    args.solutions_list = Parser.flatten(args.solutions_list)

# MARK: Begin logging
logging.basicConfig(
    filename=args.log,
    format="%(asctime)s - %(module)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=args.verbose,
)

# MARK: Call Relevant Class(Args)
logging.debug(f"Argparse namespace: {args}")
logging.info(f"Mode:{args.mode}")
args.func(**vars(args)).process()


# Confirm all processes completed and exit without error
logging.info("All done! Come again!\n")
