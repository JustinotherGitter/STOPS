"""Functions modified from POLSALT, updated for Python 3 compatibility."""

# MARK: Imports
import os
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift

from STOPS.utils.Constants import DATADIR


# MARK: DATED LINE
def datedline(filename, date):
    line_l = [ll for ll in open(filename) if ll[8:10] == "_v"]
    datever_l = [line_l[l].split()[0] for l in range(len(line_l))]

    line = ""
    for (l, datever) in enumerate(datever_l):
        if int(date) < int(datever[:8]):
            continue

        for (v, vdatever) in enumerate(datever_l[l:]):
            if int(vdatever[:8]) > int(datever[:8]):
                continue

            datever = datever_l[l + v]
            line = line_l[l + v]

    return line


# MARK: RSS TR ALIGN
def rssdtralign(datobs, trkrho):
    # optic axis is center of imaging mask.  In columns, same as longslit position
    rc0_pd = np.loadtxt(DATADIR + "RSSimgalign.txt", usecols=(1, 2))
    flex_p = np.array([np.sin(np.radians(trkrho)),
                      np.cos(np.radians(trkrho)) - 1.0])
    rcflex_d = (rc0_pd[0:2] * flex_p[:, None]).sum(axis=0)

    row0, col0, C0 = np.array(
        datedline(DATADIR + "RSSimgalign.txt", datobs).split()[1:]
    ).astype(float)
    row0, col0 = np.array([row0, col0]) + rcflex_d  # sign changed 20190610

    return row0, col0, C0


# MARK: RSSMODELWAVE
def rssmodelwave(grating, grang, artic, trkrho, cbin, cols, datobs):
    row0, col0 = rssdtralign(datobs, trkrho)[:2]
    spec_dp = np.array(
        datedline(DATADIR + "RSSspecalign.txt", datobs).split()[1:]
    ).astype(float)
    Grat0, Home0, ArtErr, T2Con, T3Con = spec_dp[:5]
    FCampoly = spec_dp[5:]

    grname = np.loadtxt(
        DATADIR + "gratings.txt",
        dtype=str,
        usecols=(0,),
        skiprows=2
    )
    grlmm, grgam0 = np.loadtxt(
        DATADIR + "gratings.txt",
        usecols=(1, 2),
        skiprows=2,
        unpack=True
    )
    grnum = np.where(grname == grating)[0][0]
    lmm = grlmm[grnum]
    alpha_r = np.radians(grang + Grat0)
    beta0_r = np.radians(artic * (1 + ArtErr) + Home0) \
        - 0.015 * col0 / FCampoly[0] - alpha_r
    gam0_r = np.radians(grgam0[grnum])
    lam0 = 1e7 * np.cos(gam0_r) * (np.sin(alpha_r) + np.sin(beta0_r)) / lmm

    # image center (unbinned pixels) for wavelength calibration model
    modelcenter = 3162.
    ww = lam0 / 1000.0 - 4.0
    fcam = np.polyval(FCampoly[::-1], ww)
    disp = (1e7 * np.cos(gam0_r) * np.cos(beta0_r) / lmm) / (fcam / 0.015)
    dfcam = (modelcenter / 1000.0) * disp \
        * np.polyval([FCampoly[5 - x] * (5 - x) for x in range(5)], ww)

    T2 = -0.25 * (1e7 * np.cos(gam0_r) * np.sin(beta0_r) / lmm) \
        / (fcam / 47.43)**2 + T2Con * disp * dfcam
    T3 = (-1.0 / 24.0) * modelcenter * disp / (fcam / 47.43)**2 + T3Con * disp
    T0 = lam0 + T2
    T1 = modelcenter * disp + 3 * T3
    X = (np.array(range(cols)) - cols / 2) * cbin / modelcenter
    lam_X = T0 + T1 * X + T2 * (2 * X**2 - 1) + T3 * (4 * X**3 - 3 * X)

    return lam_X


# MARK: READ WOLLASTON
def read_wollaston(hdu, wollaston_file):
    # set up data
    data = hdu['SCI'].data[0]
    rows, cols = data.shape
    grating = hdu[0].header['GRATING'].strip()
    grang = hdu[0].header['GR-ANGLE']
    artic = hdu[0].header['CAMANG']
    trkrho = hdu[0].header['TRKRHO']
    date = hdu[0].header['DATE-OBS'].replace('-', '')
    cbin, rbin = [int(x) for x in hdu[0].header['CCDSUM'].split(" ")]

    # load data from wollaston file
    lam_m = np.loadtxt(wollaston_file, dtype=float, usecols=(0,))
    rpix_om = np.loadtxt(
        wollaston_file,
        dtype=float,
        unpack=True,
        usecols=(1, 2)
    )
    lam_c = rssmodelwave(grating, grang, artic, trkrho, cbin, cols, date)

    return interp1d(lam_m, rpix_om, kind='cubic', bounds_error=False)(lam_c), cols, rbin, lam_c


# MARK: CORRECT WOLLASTON
def correct_wollaston(data, drow_shift):
    rows, cols = data.shape
    sdata = np.zeros(data.shape, dtype='float32')
    for c in range(cols):
        shift(data[:, c], drow_shift[c], sdata[:, c], order=1)
    return sdata


############################################################################################

# MARK: SPLIT_SCI
def split_sci(hdulist, splitrow, ext="SCI"):
    """
    split pipeline output FITS files at a certain splitrow
    roll over the array so that the split is in the middle of the FITS frame
    """

    hdu = deepcopy(hdulist)
    rows, cols = hdu[ext].data.shape

    # if odd number of rows, strip off the last one
    rows = int(rows / 2) * 2

    # how far split is from center of detector
    offset = int(splitrow - rows / 2)

    # split arc into o/e images
    padbins = (np.indices((rows, cols))[0] < offset) | (
        np.indices((rows, cols))[0] > rows+offset)

    image_rc = np.roll(hdu[ext].data[:rows, :], -offset, axis=0)
    image_rc[padbins] = 0.0

    # print(padbins)
    # print(hdulist['SCI'].data[:int(rows / 2)] == image_rc)

    hdu[ext].data = image_rc.reshape((2, int(rows/2), cols))

    return hdu


# MARK: OUTFILES
def outfiles(fname, splitrow, save_O=None, save_E=None, display=False):
    """
    save split FITS files
    """

    O_beam = pyfits.HDUList()
    E_beam = pyfits.HDUList()

    with pyfits.open(fname) as hdul:
        O_beam.append(hdul['PRIMARY'].copy())
        E_beam.append(hdul['PRIMARY'].copy())

        temp = split_sci(hdul, splitrow)

        O_beam[0].data = temp['SCI'].data[1]
        E_beam[0].data = temp['SCI'].data[0]

        arc = False

        if (hdul['PRIMARY'].header['OBJECT'] == 'ARC') & (display == True):
            print("Arc lamp: ", hdul['PRIMARY'].header['LAMPID'])
            print("Grating : ", hdul['PRIMARY'].header['GRATING'])
            print("GR angle: ", hdul['PRIMARY'].header['GR-ANGLE'])
            print("AR angle: ", hdul['PRIMARY'].header['AR-ANGLE'])
            # print("CR VAL  : ", hdul['PRIMARY'].header['CRVAL'])
            # print("C DELT  : ", hdul['PRIMARY'].header['CDELT'])

            arc = True  # skips output for arc

        if display & (not arc):
            hdul.info()
            O_beam.info()
            E_beam.info()

        if save_O and save_E:
            O_beam.writeto(save_O, overwrite=True)
            E_beam.writeto(save_E, overwrite=True)

    return O_beam, E_beam
