"""Utility functions for the Parser"""

# MARK: Imports
import os
import logging
from pathlib import Path


# MARK: Parse Path
def parse_path(path: str) -> Path:
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
def parse_file(filename: str) -> Path:
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
    # Set order of desired -v -> -vvv logging
    loglist = [logging.WARNING, logging.INFO, logging.DEBUG]

    return loglist[max(0, min(loglevel, len(loglist) - 1))]


# MARK: Parse Logging File
# lambda name: '' if name is None else '.'.join([*[i.upper() for i in name.split('.')[:-1 if name.endswith('log') else len(name.split('.'))]], 'log'])
def parse_logfile(logfile: str) -> str | None:
    # Handle no logfile
    if logfile in [None, ""]:
        return ""

    # Return valid 'LOGFILE.log'
    logfile = logfile.upper()
    if logfile.endswith(".LOG"):
        logfile = logfile[:-4]
    logfile += ".log"

    return logfile


# MARK: Flatten
def flatten(filelist: list) -> list[Path]:
    flatlist = []

    for item in filelist:
        if type(item) != list:
            flatlist.append(item)

        for subitem in item:
            flatlist.append(subitem)

    flatlist = sorted(set(flatlist))

    return flatlist
