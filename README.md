# [STOPS] (*Supplementary TOols for POLSALT Spectropolarimetry*)

[STOPS] is a Python3-based pipeline designed to provide additional functionality of the [POLSALT] pipeline. This pipeline provides a command-line interface (CLI) that can call four main classes, each implementing a unique function: `Split`, `Join`, `Skylines`, and `CrossCorrelate`.

## Installation

To install the STOPS pipeline, clone the repository and install the dependencies:

```bash
git clone https://github.com/JustinotherGitter/STOPS.git
cd STOPS
# Additionally create a STOPS environment
pip install -r requirements.txt
```

## Classes

* [**Split**](/split.py): Separate the pre-processed [POLSALT] files by the perpendicular beams.
* [**Join**](/join.py): Combine external wavelength files and with the perpendicular beams, as expected by [POLSALT]'s `spectral extraction`
* [**Skylines**](/skylines.py): Automatically identify and calculate the difference between [known](http://pysalt.salt.ac.za/lineatlas/sky_strengths.txt "SALT identified sky lines") and observed sky lines (or [arc lines](https://astronomers.salt.ac.za/data/salt-longslit-line-atlas/ "SALT arc lines available for calibration")) in spectropolarimetric data.
* [**CrossCorrelate**](/cross_correlate.py): Perform cross-correlation analysis on spectropolarimetric data, for each file comparing the perpendicular polarization beams, or across multiple files comparing a singular polarization beam.

## Procedure
A simplistic workflow is provided below, for further information and implementation please see the [STOPS writeup](https://github.com/JustinotherGitter/Masters-Thesis/Thesis.pdf "Justin Cooper - Master Thesis").

1. Pre-reductions may be performed using [POLSALT], or downloaded directly.
1. The $O$- & $E$-beams are split into to separate FITS files using `split`, into the format required by [IRAF].
1. Wavelength calibrations are performed via [IRAF]
    * Alternative external tools, such as Python may also be used for wavelength calibrations, but must be formatted into the format specified in the [Wavelength section](#wavelength-calibration).
1. The $O$- & $E$-beams, along with their respective wavelength solutions, are recombined into a single file using `join`, with the respective header and extensions updated for [POLSALT], and cosmic-ray cleaning performed via the `lacosmic` algorithm implemented through `ccdproc`.
1. The wavelength calibrations for the $O$- & $E$-beams may be compared using `skylines`, highlighting variations between the individual wavelength solutions.
1. The extracted spectra may be compared using `correlate`, allowing the correlation between the perpendicular polarization beams within a file to be correlated, or for a polarization beam across multiple files.
1. Spectral extraction and polarization calculations are performed via [POLSALT].
1. Flux calibration may be performed using the astropy and scipy packages, assuming a standard is available for the observational setup.

## Wavelength Calibration


## CLI Usage

The [STOPS] pipeline can be controlled via a command-line interface. Below are the details of the available commands and their options.

### General Options

- `-V`, `--version`: Show the version of [STOPS].
- `-v`, `--verbose`: Enable and increase verbosity. Use `-v` or `-vv` for greater verbosity levels.
- `-l`, `--log`: Specify the filename of the logging file. The file is created if it does not exist.

### Commands

#### Split

```bash
[py-path]python3 [rel-path]STOPS [General options] [data-dir] split [Options] [data-files]
```

**Options:**
- `-n`, `--no_arc`: Exclude arc files from processing.
- `-s`, `--split_row`: Row along which the O and E beams are split. Defaults to the pipeline's default.
- `-p`, `--save_prefix`: Prefix appended to the filenames for saving the O and E beams. Defaults to the pipeline's default.

#### Join

Combine multiple data sets into a single cohesive unit.

```bash
<py-path>python3 <rel-path>STOPS join [OPTIONS] [data_dir]
```

**Options:**
- `-n`, `--no_arc`: Exclude arc files from processing.
- `-s`, `--split_row`: Row along which the O and E beams are split. Defaults to the pipeline's default.
- `-p`, `--save_prefix`: Prefix appended to the filenames for saving the O and E beams. Defaults to the pipeline's default.
- `-c`, `--coefficients`: Custom coefficients to use instead of the `IRAF` fitcoords database. Use as `-c <o_solution> <e_solution>` or a regex descriptor `-c <*solution*extension>`.

#### Skylines

Identify and process sky lines in spectropolarimetric data.

```bash
<py-path>python3 <rel-path>STOPS skylines [OPTIONS] filenames
```

**Options:**
- `-b`, `--beams`: Beams to process. Defaults to `OE`, but may be given `O`, `E`, or `OE`.
- `-ccd`, `--split_ccd`: Flag to NOT split CCD's.
- `-c`, `--continuum_order`: Order of continuum to remove from spectra.
- `-p`, `--plot`: Flag for additional plot outputs.
- `-s`, `--save_prefix`: Prefix used when saving plot.
- `-t`, `--transform`: Force transform images.

#### CrossCorrelate

Perform cross-correlation analysis on spectropolarimetric data.

```bash
<py-path>python3 <rel-path>STOPS crosscorrelate [OPTIONS] filenames
```

**Options:**
- `-b`, `--beams`: Beams to process. Defaults to `OE`, but may be given `O`, `E`, or `OE`.
- `-ccd`, `--split_ccd`: Flag to NOT split CCD's.
- `-c`, `--continuum_order`: Order of continuum to remove from spectra.
- `-p`, `--plot`: Flag for additional plot outputs.
- `-s`, `--save_prefix`: Prefix used when saving plot.
- `-o`, `--offset`: Introduce an offset when correcting for known offset in spectra or for testing purposes.

### Help Command

For detailed information on each command and its options, use the help command:

```bash
<py-path>python3 <rel-path>STOPS <command> --help
```

Replace `<command>` with `split`, `join`, `skylines`, or `crosscorrelate` to get detailed information about the available options for each function.


[IRAF]: <https://iraf-community.github.io/> (IRAF GitHub Repository)
[POLSALT]: <https://github.com/saltastro/polsalt> (POLSALT GitHub Repository)
[STOPS]: <https://github.com/JustinotherGitter/STOPS> (STOPS GitHub Repository)