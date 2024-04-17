# Look at all TODO comments.
* Move all TODO comments here
* Handle direct call in main()
* Add verbosity logging statements (check log to logfile)
* Complete `NumPy styled` docstrings

## __main__.py
  * Add plot return options? [open | open and save to file | file]
  * Add type of file check (I.E. FITS)

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

## split.py
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

## join.py
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

## correlate.py
  *

## skylines.py
  *


# Using general python project structure
  * make sure implemented throughout
  * create dependencies.txt
  * create Makefile to clear up installation
  * update docs to reflect updated structure

  * refer to:
      * https://docs.python-guide.org/writing/structure/

# Using NumPy docstring formatting as well as type hinting
  * make sure implemented throughout
  * make sure type hinting using correct types (I.E. find type for fits files and implement)

  * refer to:
      * https://numpydoc.readthedocs.io/en/latest/format.html
      * https://docs.python.org/3/library/typing.html

# Implement logging
  * make sure implemented throughout
  * decide on log structure (I.E. [Datetime] Function called and returned {output} / raised FileNotFound error)

  * refer to:
      * https://realpython.com/python-logging/
  