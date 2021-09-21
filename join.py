#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Justin Cooper"
__version__ = "20.09.2021"
__email__ = "justin.jb78@gmail.com"

class Join:
    """
    Join class allows for the seperate call of joining the wavelength calibrated O & E beam FITS files

    Parameters
    ----------
    path_name : str
        The path to the data (wmxgbp*.fits files) to be joined
    split_row : int, optional
        The row that the data was split along.
        (The default is 517, the middle row of the CCD's)
    verbose : bool, optional
        Decides whether the output should be recorded to the terminal window.
        (The default is False, no output written to the terminal window)
    no_arc : bool, optional
        Decides whether the arc frames should be recombined.
        (The default is False, since polsalt only uses the arc frames until spectral extraction)
    save_pref : list of str, optional
        The prefix that the O & E beams are saved as.
        (The default is ["obeam", "ebeam"], which is what split defaults to)
    
    Returns
    -------
    joined_FITS : list of FITS
        A list of FITS files that were joined and can be returned to polsalt.

    Raises
    ------
    # TODO@JustinotherGitter : Complete docs for which errors are raised and when
    """
    def __init__():


########################
# TODO@JustinotherGitter: Finish refactoring/rewriting and remove below
# not to be included in any distributions
def join(pathname, splitrow=517, display=False, no_arc=False, save_pref=["obeam", "ebeam"]):
    DATADIR = os.path.expanduser("~/polsalt-beta/polsalt/data/")
    db = "database"
    
    print(f"Path to data given as: {pathname}")
    
    # Get pipeline FITS files from path given
    infilelist = []
    arcfile = ""
    fc_files = []
    
    # Get valid FITS files
    for fl in os.listdir(pathname):
        if os.path.isfile(os.path.join(pathname, fl)) and ("m" == fl[0]) and ("fits" == fl.split(".")[-1]):
            #print(fl)
            infilelist.append(fl)
            
    # TODO: Allow specification of ARC file in case multiple ARCS
    # Get arcfile from valid FITS files
    for fl in infilelist:
        with pyfits.open(fl) as hdu:
            if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                if no_arc:
                    arcfile = fl
                    infilelist.remove(fl)
                else:
                    arcfile = fl
                    
    for fl in os.listdir(os.path.join(pathname, db)):
        if os.path.isfile(os.path.join(pathname, db, fl)) and ("fc" == fl[0:2]):
            fc_files.append(fl)
            if display: print(f"Found fitcoords solution {fl}")
    
    if display:
        print("Arcfile:\t\t\t" + arcfile)
        print("Pipeline m*.fits files found:\t" + str(infilelist))
        if no_arc: print("Removed arcfile from list to join")


    for fits in infilelist:
        # TODO: Make sure non-t files have wav solution
        if fits != arcfile:
            o_file = save_pref[0] + fits[-9:]
            e_file = save_pref[1] + fits[-9:]
        else:
            o_file = "oarc" + fits[-9:]
            e_file = "earc" + fits[-9:]
        
        with pyfits.open(fits) as hdu:
            with pyfits.open(o_file) as o:
                with pyfits.open(e_file) as e:
                    print(o[0].data.shape, e[0].data.shape)
                    
                    # Check if files have been cropped
                    if hdu["SCI"].data.shape[0] / 2 != o[0].data.shape[0]:
                        # get crop size
                        cropsize = int(hdu["SCI"].data.shape[0] / 2 - o[0].data.shape[0])
                        
                        print(f"The original shape was {hdu['SCI'].data.shape} and was croppped by {cropsize} to {o[0].data.shape}")
                        
                    else:
                        # haven't been cropped
                        cropsize = False
                    
                    y_shape = int(hdu["SCI"].data.shape[0] / 2) - cropsize
                    x_shape = hdu["SCI"].data.shape[1]
                    
                    whdu = pyfits.HDUList()
                    # No differences in "PRIMARY" extention header
                    whdu.append(hdu["PRIMARY"])
                    
                    
                    for ext in ["SCI", "VAR", "BPM"]:
                        whdu.append(pyfits.ImageHDU(name=ext))
                        whdu[ext].header = deepcopy(hdu[ext].header)
                        whdu[ext].header["CTYPE3"] = "O,E"
                        
                        if ext == "BPM":
                            whdu[ext].data = np.zeros((2, y_shape, x_shape), dtype='uint8')
                            whdu[ext].header["BITPIX"] = "-uint8"
                        else:
                            whdu[ext].data = np.zeros((2, y_shape, x_shape), dtype='>f4')
                            whdu[ext].header["BITPIX"] = "-32"
    
                        
                        if cropsize:
                            temp_split = split_sci(hdu, splitrow, ext=ext)[ext].data
                            whdu[ext].data[0] = temp_split[0, cropsize:]
                            whdu[ext].data[1] = temp_split[1, 0:-cropsize]
                            
                        else:
                            whdu[ext].data = split_sci(hdu, splitrow, ext=ext)[ext].data
                    
                    
                    whdu.append(pyfits.ImageHDU(name="WAV"))
                    wav_header = deepcopy(whdu["SCI"].header)
                    wav_header["EXTNAME"] = "WAV"
                    wav_header["CTYPE3"] = "O,E"
                    whdu["WAV"].header = wav_header
                    
                    whdu["WAV"].data = np.zeros(whdu["SCI"].data.shape, dtype='>f4')
    
                    for num, fname in enumerate(fc_files):
                        print(f"Using {fname} as solution")
    
                        chebvals = []
                        with open(db + "/" + fname) as file:
                            for i in file:
                                # TODO: replace with strip
                                def stripchars(text):
                                    if text[0] in ["\t", "\n"]:
                                        text = text[1:]
                                    if text[0] in ["\t", "\n"]:
                                        text = text[1:]
                                    if text[-1] in ["\t", "\n"]:
                                        text = text[:-1]
                                    return text
    
                                chebvals.append(stripchars(i))
    
    
    
    
                        if chebvals[6] == "1.": #function - Function type (1=chebyshev, 2=legendre)
                            x_ord = int(chebvals[7][:-1]) #xorder - X "order" (highest power of x)
                            y_ord = int(chebvals[8][:-1]) #yorder - Y "order" (highest power of y)
                            if chebvals[9] == "1.": #xterms - Cross-term type (always 1 for FITCOORDS)
                                xmin = int(float(chebvals[10][:-1])) #xmin - Minimum x over which the fit is defined
                                xmax = int(float(chebvals[11][:-1])) #xmax - Maximum x over which the fit is defined
                                ymin = int(float(chebvals[12][:-1])) #ymin - Minimum y over which the fit is defined
                                ymax = int(float(chebvals[13][:-1])) #ymax - Maximum y over which the fit is defined
                                
                                if ymax != y_shape: # TODO: Fix temporary stretching
                                    print(f"xmin {xmin}, xmax {xmax}, ymin {ymin}, ymax {ymax}")
                                    ymax = y_shape
                                    print(f"new:\txmin {xmin}, xmax {xmax}, ymin {ymin}, ymax {ymax}")
                                    
                            c_vals = np.array(chebvals[14:], dtype=float)
                            c_vals = np.reshape(c_vals, (x_ord, y_ord))
                        
                        
                            # Set wavelength extention values to function
                            whdu["WAV"].data[num] = chebgrid2d(x=np.linspace(-1, 1, ymax),
                                                               y=np.linspace(-1, 1, xmax),
                                                               c=c_vals)
    
                        elif chebvals[6] == "2.":
                            # TODO: Handle legendre
                            pass
    
                        else:
                            #TODO: Handle other functions???
                            pass
                        
                        # Cosmic Ray Cleaning ## TODO: Check parameters
                        whdu["SCI"].data[num] = lacosmic.lacosmic(whdu["SCI"].data[num], 2, 4, 4, effective_gain=1, readnoise=3.3)[0]
                        
                        
                    # WAV mask (Left & Right Crop)
                    whdu["WAV"].data[whdu["WAV"].data[:] <  3_000] = 0.0
                    whdu["WAV"].data[whdu["WAV"].data[:] >= 10_000] = 0.0
                        
                    # Correct WAV mask shift (Top & Bottom Crop)
                    rpix_oc, cols, rbin, lam_c = read_wollaston(whdu, DATADIR + 'wollaston.txt')
    
                    drow_oc = (rpix_oc-rpix_oc[:,int(cols/2)][:,None])/rbin
    
                    ## Cropping as suggested
                    for c, col in enumerate(drow_oc[0]):
                        if not np.isnan(col):
                            if int(col) < 0:
                                whdu["WAV"].data[0, int(col):, c] = 0.0
                            elif int(col) > cropsize:
                                whdu["WAV"].data[0, 0:int(col) - cropsize, c] = 0.0
                                
                    for c, col in enumerate(drow_oc[1]):
                        if not np.isnan(col):
                            if int(col) > 0:
                                whdu["WAV"].data[1, 0:int(col), c] = 0.0
                            elif (int(col) < 0) & (abs(int(col)) > cropsize):
                                whdu["WAV"].data[1, int(col) + cropsize:, c] = 0.0
    
                    # Mask BPM same as WAV
                    whdu["BPM"].data[0] = np.where(whdu["WAV"].data[0] == 0, 1, whdu["BPM"].data[0])
                    whdu["BPM"].data[1] = np.where(whdu["WAV"].data[1] == 0, 1, whdu["BPM"].data[1])
                    
    
                    whdu.writeto("w" + fits, overwrite="True")
                    print(f"{fits} joined and saved to {'w' + fits}\n")
                    
    print("\nCOMPLETED WITH NO ERRORS\n")