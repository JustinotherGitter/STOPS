<!-- markdownlint-disable MD033 -->
# *<ins>S</ins>upplementary <ins>To</ins>ols for <span style="font-variant:small-caps;">[<ins>p</ins>olsalt](https://github.com/saltastro/polsalt "POLSALT GitHub Repository")</span> <ins>S</ins>pectropolarimetry* ([<span style="font-variant:small-caps;">stops</span>])

[<span style="font-variant:small-caps;">stops</span>] is a [Python]3-based pipeline designed to provide additional functionality to the [<span style="font-variant:small-caps;">polsalt</span>] pipeline.
This pipeline provides a command-line interface (CLI) that can call four main classes, each implementing a unique function: `Split`, `Join`, `Skylines`, and `CrossCorrelate`.
For an in-depth discussion and usage, see the [<span style="font-variant:small-caps;">stops</span>] [writeup].

<!-- MARK: Install -->
## Installation

To install the [<span style="font-variant:small-caps;">stops</span>] pipeline, simply download the wheel file and run:

```console
cd ~/Downloads
# Activate any conda environment
pip install STOPS
```

The repository may also be cloned and the [dependencies](/requirements.txt "Requirements file") installed:

```console
git clone https://github.com/JustinotherGitter/STOPS.git # or manually download the repository
cd STOPS
# Additionally create and activate a python3 `STOPS` venv
pip install -r requirements.txt
```

**Note:** [<span style="font-variant:small-caps;">stops</span>] also requires $\LaTeX$ to be installed as well as to have [<span style="font-variant:small-caps;">polsalt</span>] installed in the Home, `~/`, directory.

<!-- MARK: Classes -->
## Classes

* [**Split**](/src/STOPS/split.py): Separate the pre-processed [<span style="font-variant:small-caps;">polsalt</span>] [FITS] files by their perpendicular polarization beams into two [<span style="font-variant:small-caps;">iraf</span>] parsable [FITS] files.
* [**Join**](/src/STOPS/join.py): Combine external wavelength calibration solutions with the perpendicular beams, as expected by [<span style="font-variant:small-caps;">polsalt</span>]'s `spectral extraction`
* [**Skylines**](/src/STOPS/skylines.py): Automatically identify and calculate the difference between [known](http://pysalt.salt.ac.za/lineatlas/sky_strengths.txt "SALT identified sky lines") and observed sky lines (or [arc lines](https://astronomers.salt.ac.za/data/salt-longslit-line-atlas/ "SALT arc lines available for calibration")) in wavelength calibrated spectropolarimetric data.
* [**CrossCorrelate**](/src/STOPS/cross_correlate.py): Perform cross-correlation analysis on the spectropolarimetric data (possible after [<span style="font-variant:small-caps;">polsalt</span>]'s `spectral extraction`), either for each file comparing the perpendicular polarization beams, or across multiple files comparing a singular polarization beam.

<!-- MARK: Procedure -->
## Procedure

A simplistic workflow is provided below, for further information and implementation please see the [write-up] of [<span style="font-variant:small-caps;">stops</span>].

1. Pre-reductions may be performed on the raw data using [<span style="font-variant:small-caps;">polsalt</span>], or downloaded directly alongside the raw data.
2. The $O$- & $E$-beams are split into two separate [FITS] files using `split`, into the format required by [<span style="font-variant:small-caps;">iraf</span>].
3. Wavelength calibrations are performed via [<span style="font-variant:small-caps;">iraf</span>], replacing the [<span style="font-variant:small-caps;">polsalt</span>] `wavelength calibration`.
    * Alternate external tools, such as [Python] may also be used for wavelength calibrations.
    * Formatting of non-[<span style="font-variant:small-caps;">iraf</span>] wavelength solutions must be formatted as described in the [Wavelength section](#wavelength-calibration).
4. The $O$- & $E$-beams, along with their respective wavelength solutions, are recombined into a single file using `join`, and cosmic-ray cleaning performed via the `lacosmic` algorithm implemented through `ccdproc`.
5. The wavelength calibrations for the $O$- & $E$-beams may be compared using `skylines`, highlighting variations between the individual wavelength solutions.
6. Spectral extraction is performed via [<span style="font-variant:small-caps;">polsalt</span>].
7. The extracted spectra may be compared using `correlate`, allowing the correlation between the perpendicular polarization beams within a file to be correlated, or for a polarization beam across multiple files.
8. If either `skylines` or `correlate` show poor wavelength calibrations, the wavelength calibration procedure may be repeated.
9. The files generated, excluding the wavelength solution, for wavelength calibrations may be moved or deleted.
10. Polarization calculations are performed via [<span style="font-variant:small-caps;">polsalt</span>].
11. Flux calibration may be performed using the `astropy` and `scipy` packages, assuming a standard is available for the observational setup.

<!-- MARK: Wav Cal -->
## Wavelength Calibration

Wavelength calibrations are ideally intended for [<span style="font-variant:small-caps;">iraf</span>].
The [<span style="font-variant:small-caps;">iraf</span>] wavelength calibration procedure uses `identify`, `reidentify`, `fitcoords`, and optionally `transform`.
Any preferred parameters may be used during calibration with only the polynomial type of the resultant wavelength solution (as produced by `fitcoords`) being limited to either `Chebyshev` or `Legendre` polynomials.

Alternate wavelength solutions may be used but must be:

* Saved in the working directory,
* Have separate files for the $O$- and $E$-beams,
* Have a name containing:
  * The polynomial type (e.g. `cheb` for Chebyshev, or `leg` for Legendre), and
  * The polarization beam contained within (e.g. `O` or `E`),
* Contain on:
  * line 1 → the $x$-order of the 2D wavelength solution,
  * line 2 → the $y$-order of the 2D wavelength solution,
  * lines 3+ → the ($x * y$) solution coefficients, in order, separated by line, and
  * all lines → no delimiters.

<br>
e.g.

`cheb_params_e.txt`

```text cheb_params_e.txt
5
5
7419.096745914063
1510.03933621895
-21.10886852752348
-2.079553916887794
0.06772631420528228
0.7720164913117386
'...'
```

<!-- MARK: CLI Usage -->
## CLI Usage

The [<span style="font-variant:small-caps;">stops</span>] pipeline is most generally controlled via a CLI.
The basic format of commands follows:

```console
<py_dir>python<3> <STOPS_dir>STOPS (General Options) [data_dir] MODE (Options) [File names]
```

where:

* `<>` parameters are optional depending on the system setup, e.g. if [Python] or [<span style="font-variant:small-caps;">stops</span>] has been added to `$PATH`, etc. (for simplicity, these parameters will be left out of the usage examples below),
* `MODE` refers to the operational mode of [<span style="font-variant:small-caps;">stops</span>], as listed in [Modes](#modes),
* `()` parameters are optional, and
* `[]` parameters are compulsory (unless otherwise stated).

Below are the details of the [<span style="font-variant:small-caps;">stops</span>] options, the available Modes, and their respective sub-options.

<!-- MARK: STOPS Opts. -->
### [<span style="font-variant:small-caps;">stops</span>] Options

**Optional:**

* `-V`, `--version`: Show the version of [<span style="font-variant:small-caps;">stops</span>].
* `-v`, `--verbose`: Enable and increase verbosity.
  Use `-v` or `-vv` for greater verbosity levels.
* `-l`, `--log`: Specify the filename of the logging file, which is created if it does not exist.

**Compulsory:**

* `data_dir` : The Path (absolute or relative) to the directory containing the Working data.
  `.` may be used to indicate the current directory that the CLI is running in.

### Help Commands

For detailed information on [<span style="font-variant:small-caps;">stops</span>], use the help command:

```console
python STOPS . --help
```

Help for each mode is also accessible, and may be viewed using the help command:

```console
python STOPS . MODE --help
```

where `MODE` may be replaced with [`split`](#split), [`join`](#join), [`skylines`](#skylines), or [`crosscorrelate`](#correlate).

### Modes

---

<!-- MARK: STOPS Split -->
#### <ins>s</ins>plit

```console
python STOPS . split (Options) [mxgbp*.fits]
```

##### <ins>s</ins>plit Options

**Optional:**

* `-n`, `--no_arc`: Exclude arc files from processing.
* `-s`, `--split_row`: Row along which the $O$- & $E$-beams are split.
  Defaults to the [<span style="font-variant:small-caps;">polsalt</span>]'s default.
* `-p`, `--save_prefix`: Prefix appended to the filenames for saving the $O$- & $E$-beams.
  Defaults to the [<span style="font-variant:small-caps;">polsalt</span>]'s default.

**Compulsory:**

* Filenames to be split.
  May be excluded if the [RegEx] pattern in the example matches the desired files.

---

<!-- MARK: STOPS Join -->
#### <ins>j</ins>oin

```console
python STOPS . join (Options) []
```

##### <ins>j</ins>oin Options

**Optional:**

* `-n`, `--no_arc`: Exclude arc files from processing.
* `-s`, `--split_row`: Row along which the $O$- & $E$-beams are split.
  Defaults to the pipeline's default.
* `-p`, `--save_prefix`: Prefix appended to the filenames for saving the $O$- & $E$-beams.
  Defaults to the pipeline's default.
* `-c`, `--coefficients`: Custom coefficients to use instead of the `IRAF` fitcoords database.
  Use as `-c <o_solution> <e_solution>` or a [RegEx] descriptor `-c <*solution*extension>`.

**Compulsory:**

* May be excluded as `join` will identify all split files and wavelength solution database entries to recombine.

---

<!-- MARK: STOPS Sky. -->
#### <ins>sky</ins>lines

```console
python STOPS . skylines (Options) [Filenames]
```

##### <ins>sky</ins>lines Options

**Optional:**

* `-b`, `--beams`: Beams to process.
  Defaults to `OE`, but may be given `O`, `E`, or `OE`.
* `-ccd`, `--split_ccd`: Flag to NOT split CCD's.
* `-c`, `--continuum_order`: Order of continuum to remove from spectra.
* `-p`, `--plot`: Flag for additional debug plot outputs.
* `-s`, `--save_prefix`: Prefix used when saving plot.
* `-t`, `--transform`: Force transform images, for [<span style="font-variant:small-caps;">iraf</span>] `transform` [FITS] file debugging.

**Compulsory:**

* Filenames to be considered for `skyline` comparisons.
  May either be:
  * the `wmxgbp*.fits` [RegEx] pattern for [<span style="font-variant:small-caps;">polsalt</span>] formatted, wavelength calibrated, [FITS] files, or
  * the `tbeam*.fits` [RegEx] pattern for [<span style="font-variant:small-caps;">iraf</span>] formatted, `transform` output, [FITS] files.

---

<!-- MARK: STOPS Corr. -->
#### <ins>correlate</ins>

```console
python STOPS . correlate (Options) [ecwmxgbp*.fits]
```

##### <ins>correlate</ins> Options

**Optional:**

* `-b`, `--beams`: Beams to process.
  Defaults to `OE`, but may be given `O`, `E`, or `OE`.
* `-ccd`, `--split_ccd`: Flag to NOT split CCD's.
* `-c`, `--continuum_order`: Order of continuum to remove from spectra.
* `-p`, `--plot`: Flag for additional debug plot outputs.
* `-s`, `--save_prefix`: Prefix used when saving plot.
* `-o`, `--offset`: Introduce an offset when correcting for known offset in spectra or for testing purposes.
  Deprecated keyword.

**Compulsory:**

* Filenames to be considered for `correlate` cross-correlation.
  May be excluded if the [RegEx] pattern in the example matches the desired files.

<br>

<!-- MARK: Contribute -->
## Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes.
The following styles are broadly implemented, and serve as a collection of references used when creating the [<span style="font-variant:small-caps;">stops</span>] pipeline:

* General [Python] [project structure](https://docs.python-guide.org/writing/structure/ "Structuring a Python project") applies
* Docstrings follow the NumPy [documentation style](https://numpydoc.readthedocs.io/en/latest/format.html "NumPy style guide")
  * Possible `sphinx` docs?
* Classes and methods implement [typing](https://docs.python.org/3/library/typing.html "typing in Python") for type hinting
* [Logging](https://realpython.com/python-logging/ "Logging in Python") is implemented
* [Tests](https://docs.python-guide.org/writing/tests/ "Testing in Python") are planned but not implemented
  * `pytest` or `unittest` (`pytest-cov` | `tox`)?
* [requirements.txt](./requirements.txt "Python requirements") was generated using [`pipreqs`](https://pypi.org/project/pipreqs/ "`pipreqs` PyPI page")
  * Found using:

    ```console
    cd /media/justin/Transcend/STOPS
    pipreqs
    ```
  
* The minimum required [Python] version was found using [`vermin`](https://pypi.org/project/vermin/ "`vermin` PyPI page")
  * Found using:

    ```console
    cd /media/justin/Transcend/STOPS
    vermin --eval-annotations --feature union-types ./src/STOPS
    ```
  
* [<span style="font-variant:small-caps;">stops</span>] management
  * Building generates `tar.gz` and `whl` files, overwrites if version not updated.
    Built using:

    ```console
    cd /media/justin/Transcend/STOPS
    python3 -m build
    ```

  * Installed using `pip` in an editable format (*exclude* `-e` for `stable` release):

    ```console
    cd /media/justin/Transcend/STOPS
    pip install -e ../STOPS
    ```

<!-- MARK: See Also -->
## See Also

* Package creation:
  * <https://packaging.python.org/en/latest/guides/writing-pyproject-toml/>
  * <https://packaging.python.org/en/latest/tutorials/packaging-projects/>
  * <https://packaging.python.org/en/latest/tutorials/installing-packages/>

* Package resources
  * <https://docs.python.org/3.11/library/importlib.resources.html>
  * <https://importlib-resources.readthedocs.io/en/latest/using.html>
  * <https://setuptools.pypa.io/en/latest/userguide/datafiles.html>

* [SALT] Line Atlases
  * <http://pysalt.salt.ac.za/lineatlas/>

<!-- MARK: Liscensing -->
## License

This project is licensed under the BSD 3-Clause License.
See the [LICENSE](/LICENSE "STOPS License") for further details.

<!-- MARK: Links/Refs -->
[writeup]: <https://github.com/JustinotherGitter/Masters-Thesis/Thesis.pdf> (Justin Cooper - Master Thesis)

[<span style="font-variant:small-caps;">iraf</span>]: <https://iraf-community.github.io/> (IRAF GitHub Repository)
[<span style="font-variant:small-caps;">polsalt</span>]: <https://github.com/saltastro/polsalt> (POLSALT GitHub Repository)
[<span style="font-variant:small-caps;">stops</span>]: <https://github.com/JustinotherGitter/STOPS> (STOPS GitHub Repository)

[SALT]: <https://astronomers.salt.ac.za/> (SALT for astronomers Homepage)
[FITS]: <https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf> (FITS file standard)
[RegEx]: <https://regexr.com/> (Basic RegEx test environment)
[Python]: <https://www.python.org/> (The Python website)
