# Masters
### Substitution of the wavelength calibration of the SALT Spectropolarimetric Pipeline using python 3 and IRAF.

## Installation:
wav_replacement.py and cross_corr.py are added to the main polsalt-beta directory and the specpolpy3.py is added to the polsalt sub-directory.

## Procedure
Basic reduction steps using this workflow are as follows:

  * Pre-reductions are performed using polsalt
  * The O & E beams are split into to separate FITS files and cropped to remove any data-less rows, into the format required by IRAF
  * Wavelength calibrations are performed in IRAF
  * The wavelength solution for the O & E beams may be compared, highlighting any variations in the individual wavelength solutions
  * The O & E beam FITS files are recombined into a single file, the header and extensions are updated for polsalt, and cosmic-ray cleaning is performed with the lacosmic package
  * Spectral extraction and polarization calculations are performed with polsalt
  * Flux calibration may be performed using the astropy and scipy packages, assuming a standard is available for the observations.
