# MARK: Imports

import os
from pathlib import Path
import logging


def parse_path(path: str) -> os.PathLike:
    # MARK: Parse Path

    pathlike = ""

    # Parse 'quick' path
    if path in ["", "."]:
        return Path(os.getcwd())
    
    # Parse 'relative' path
    npath = Path(path).expanduser().resolve()
    pathlike = npath if npath.is_dir() else ""

    # Validate path
    if pathlike != "":
        os.chdir(pathlike)
        return pathlike

    # Raise Directory not found error
    msg = (
        f"The given directory, {path}, does not resolve to a valid "
        "directory."
    )
    logging.error(msg)
    raise FileNotFoundError(msg)


def parse_file(file: str) -> os.PathLike:
    # MARK: Parse File

    filelike = ""

    # Parse 'relative' file path
    nfile = Path(file).expanduser().resolve()
    filelike = nfile if nfile.is_file() else ""

    # Validate file path
    if filelike != "":
        return filelike
    
    # Raise File not found error
    msg = (
        f"The given filename, {file}, does not resolve to a valid "
        "file."
    )
    logging.error(msg)
    raise FileNotFoundError(msg)


def parse_corr_file(filename: str) -> os.PathLike:
    # MARK: Parse Correlation File

    filelike = Path(filename)
    filelist = Path(filelike.parent).glob(filelike.name)
    filelist = list(sorted(filelist))

    if len(filelist) == 0:
        raise FileNotFoundError(f"No file, {filename}, found.")

    return filelist


def parse_coeff_file(file:str) -> os.PathLike | list[os.PathLike]:
    # MARK: Parse Surface Function File

    try:
        return parse_file(file)
    except:
        return list(sorted(Path().glob(file)))


def parse_loglevel(loglevel: int) -> int:
    # MARK: Parse Logging Level

    # Set order of desired -v -> -vvv logging
    loglist = [logging.WARNING, logging.INFO, logging.DEBUG]
    # Return desired level
    return loglist[max(0, min(loglevel, 2))]


def parse_logfile(logfile: str) -> str | None:
    # MARK: Parse Logging File

    # Handle no logfile
    if logfile in [None, ""]:
        return ""
    
    # Return valid 'LOGFILE.log'
    fname = logfile.upper()
    if fname[-4:] == ".LOG":
        fname = fname[:-4]
    fname += ".log"

    return fname


def flatten(filelist: list) -> list[os.PathLike]:
    # MARK: Flatten
    flatlist = []

    for item in filelist:
        if type(item) != list:
            flatlist.append(item)

        for subitem in item:
            flatlist.append(subitem)

    flatlist = sorted(set(flatlist))

    return flatlist
