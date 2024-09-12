"""Utility functions for the CLI Parser"""

# MARK: Imports
import os
import logging
from pathlib import Path


# MARK: Parse Path
def parse_path(path: str) -> Path:
    """
    Parse a given string path into a valid Path object.

    Parameters
    ----------
    path : str
        A relative or absolute path to a directory.

    Returns
    -------
    Path
        A valid Path object to a directory.

    Raises
    ------
    NotADirectoryError
        Raised if the path is not a valid directory.
    """
    # Expand and resolve the path
    path = Path(path).expanduser().resolve()

    # Check if the path is a directory
    if path.is_dir():
        os.chdir(path)
        return path

    # Raise Directory not found error
    msg = f" The directory path `{path}` is not valid."
    raise NotADirectoryError(msg)


# MARK: Parse File
def parse_file(filename: str) -> Path | list[Path]:
    """
    Parse a given, filename or regex, string into a valid Path object.

    Parameters
    ----------
    filename : str
        A relative or absolute path, terminating in a file or a regex pattern.

    Returns
    -------
    Path | list[Path]
        A valid Path object to a file, or a list of Path objects matching a regex pattern.

    Raises
    ------
    IsADirectoryError
        Raised if the filename is a directory.
    FileNotFoundError
        Raised if the filename is not found.
    """
    # Expand and resolve the filename
    filename = Path(filename).expanduser().resolve()

    # Check if the file is a directory
    if filename.is_dir():
        errMsg = f"Filename, {filename}, is a directory."
        raise IsADirectoryError(errMsg)

    # Check if the file exists
    if filename.is_file():
        return filename

    # Check if the file is a regex pattern
    filelist = Path(filename).parent.glob(Path(filename).name)
    filelist = sorted(filelist)

    # Check if the filelist is empty
    if not filelist:
        errMsg = f"No file(s), {filename}, found."
        raise FileNotFoundError(errMsg)

    return filelist


# MARK: Parse Logging Level
def parse_loglevel(loglevel: int) -> int:
    """
    Parse the logging level into a valid `logging` level.

    Parameters
    ----------
    loglevel : int
        The logging level to parse, as a count of `v`'s.

    Returns
    -------
    int
        A log level parsable by the `logging` module.
    """
    # Set order of desired -v -> -vvv logging
    # https://docs.python.org/3/library/logging.html#logging-levels
    loglist = [logging.WARNING, logging.INFO, logging.DEBUG]

    return loglist[max(0, min(loglevel, len(loglist) - 1))]


# MARK: Parse Logging File
# lambda name: '' if name is None else '.'.join([*[i.upper() for i in name.split('.')[:-1 if name.endswith('log') else len(name.split('.'))]], 'log'])
def parse_logfile(logfile: str) -> str:
    """
    Parse the logfile into a valid `logging` file.

    Parameters
    ----------
    logfile : str
        A string representing the logfile to parse, located in the data directory.

    Returns
    -------
    str
        A valid logfile parsable by the `logging` module, or an empty string.
    """
    # Handle no logfile
    if logfile in [None, ""]:
        return ""

    # Return valid 'LOGFILE.log'
    # logfile = logfile.upper()
    if logfile.endswith(".LOG") or logfile.endswith(".log"):
        logfile = logfile[:-4]
    logfile += ".log"

    return logfile


# MARK: Flatten
def flatten(filelist: list[list | Path]) -> list[Path]:
    """
    A utility function to flatten a list of lists into a single list.

    Parameters
    ----------
    filelist : list
        A list of lists to flatten, irrespective of inner type.

    Returns
    -------
    list[Path]
        A list of Path objects, sorted and unique.
    """
    flatlist = []

    # Flatten the list
    for item in filelist:
        if type(item) != list:
            flatlist.append(item)
            continue

        for subitem in item:
            flatlist.append(subitem)

    # Sort and remove duplicates
    flatlist = sorted(set(flatlist))

    return flatlist
