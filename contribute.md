# \_\_main\_\_.py
  * Add type of file check (I.E. FITS)
  * Add plot return options? [open | open and save to file | file]

```mermaid
  graph LR;
    A[Universal Parser] -->|Split\Join Parent Args| B[Split\Join Subparsers]
    B -->|Split mode| C[Split Subparser]
    B -->|Join mode| D[Join Subparser]
    A -->|Cross correlation mode| E[Correlate Subparser]
    A -->|Sky line check mode| F[Skyline Subparser]
    D -->|Custom coefficients, etc.| G[Join Subparser Children]
    E -->|Filenames, options, etc.| H[Correlate Subparser Children]
    F -->|Filenames, options, etc.| I[Skyline Subparser Children]
    C -->|Split mode args| J[Split Subparser Children]
    J -->|Defaults| K[Split Defaults]
    D -->|Join mode args| L[Join Subparser Children]
    L -->|Defaults| M[Join Defaults]
    E -->|Cross correlation mode args| N[Correlate Subparser Children]
    N -->|Defaults| O[Correlate Defaults]
    F -->|Sky line check mode args| P[Skyline Subparser Children]
    P -->|Defaults| Q[Skyline Defaults]
    B -.-> A
```


# split.py
  * self.split_row → Check valid split row at assignment
  * self.save_prefix → Check valid dict at assignment

```mermaid
  classDiagram
    class Split {
      str Data_dir
      list[str] fits_list  = None
      int split_row = 517
      bool no_arc = False
      dict save_prefix = None
      
      + split_file(file: os.PathLike) tuple[pyfits.HDUList]
      + split_ext(hdulist, ext: str = "SCI") pyfits.HDUList
      + crop_file(hdulist, crop: int = 40) tuple[np.ndarray]
      + update_beam_lists(o_name, e_name, arc: bool = True) None
      + save_beam_lists() None
      + process() None
    }

```


# join.py
  * self.split_row → Check valid split row at assignment
  * self.save_prefix → Check valid dict at assignment

```mermaid
  classDiagram
    class Join {
      str data_dir
      str database = "database"
      list[str] fits_list = None
      list[str] | None fc_files = None
      bool custom = `fc_files != None`
      int split_row = 517
      dict save_prefix = None
      bool no_arc = True
      str arc = None
      int verbose = `False`

      + get_solutions(wavlist: list | None, prefix: str = "fc") tuple[list[str], bool]
      + parse_solution(fc_file, x_shape, y_shape) tuple[dict[str, int], np.ndarray]
      + join_file(file: os.PathLike) None
      + check_crop(hdu, o_file, e_file) int
      + process() None
    }
```


# correlate.py
  * Update correlate (fits_list ← in1/in2), (cont_ord ← cont)
  * LoadFile: Recheck return of o and e beams.
  * Complete correlate


# skylines.py
  * Complete skylines


# Using general python project structure
  * Handle direct calls in main()?
  * Create Makefile to clear up installation
  * Create license (CC BY-SA 4.0 international(?))
    * https://github.com/santisoler/cc-licenses
    * https://creativecommons.org/licenses/by-sa/4.0/deed.en
    * check NC, etc.

  * Refer to:
      * https://docs.python-guide.org/writing/structure/


# Using NumPy docstrings / type hinting
  * Make sure implemented throughout
  * Make sure type hinting using correct types

  * Refer to:
      * https://numpydoc.readthedocs.io/en/latest/format.html
      * https://docs.python.org/3/library/typing.html


# Implement logging
  * Make sure implemented throughout
  * Add verbosity logging statements

  * Refer to:
      * https://realpython.com/python-logging/


# Testing
  * Add tests for \_\_main\_\_
  * Add tests for each module
  * Add tests for utils

  * Refer to:
      * https://docs.python-guide.org/writing/tests/