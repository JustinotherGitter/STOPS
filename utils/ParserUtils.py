"""Utility functions for the Parser"""

# MARK: Imports
import os
import logging
from pathlib import Path


# MARK: Parse Path
def parse_path(path: str) -> Path:
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        os.chdir(path)
        return path

    # Raise Directory not found error
    msg = f" The directory path `{path}` is not valid."
    raise NotADirectoryError(msg)


# MARK: Parse File
def parse_file(filename: str) -> Path:
    filename = Path(filename).expanduser().resolve()
    if filename.is_file():
        return filename

    # Raise File not found error
    msg = f" The filename `{filename}` is not valid."
    raise FileNotFoundError(msg)


# MARK: Parse Correlation File
def parse_corr_file(filename: str) -> Path:
    filelike = Path(filename)
    filelist = Path(filelike.parent).glob(filelike.name)
    filelist = sorted(filelist)

    if not filelist:
        errMsg = f"No file, {filename}, found."
        raise FileNotFoundError(errMsg)

    return filelist


# MARK: Parse Surface Function File
def parse_coeff_file(filename: str) -> list[Path]:
    try:
        return [parse_file(filename)]

    except:
        return list(sorted(Path().glob(filename)))


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
        fname = fname[:-4]
    fname += ".log"

    return fname


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
