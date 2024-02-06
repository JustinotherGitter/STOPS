# Parser helper functions for convenience
import os
from pathlib import Path
import logging

# TODO@JustinotherGitter: Implement logfile name validation and relative pathing


def parse_path(path: str) -> os.PathLike:
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
    raise FileNotFoundError(f"The given directory, {path}, does not resolve to a valid directory.")


def parse_file(file: str) -> os.PathLike:
    filelike = ""

    # Parse 'relative' file path
    nfile = Path(file).expanduser().resolve()
    filelike = nfile if nfile.is_file() else ""

    # Validate file path
    if filelike != "":
        return filelike
    
    # Raise File not found error
    raise FileNotFoundError(f"The given filename, {file}, does not resolve to a valid file.")


def parse_loglevel(loglevel: int) -> int:
    loglist = [logging.WARNING, logging.INFO, logging.DEBUG]
    return loglist[max(0, min(loglevel, 2))]

def parse_logfile(logfile: str) -> str | None:
    if logfile in [None, ""]:
        return None
    
    return logfile.upper() + ".log"
