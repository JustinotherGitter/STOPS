# Code replaced from main py file. Kept for comparison while rewriting but to be deleted as still in main branch

#region Split PIPELINE files for IRAF 
def split(pathname, splitrow=517, display=False, no_arc=False, save_pref=["obeam", "ebeam"]):
    if display: print(f"Path to data given as: {pathname}")
    
    # Get ARC and DATA FITS files from path given
    arcfile = ""
    infilelist = []
    o_files = []
    e_files = []
    
    # Get all valid FITS files
    for fl in os.listdir(pathname):
        if os.path.isfile(os.path.join(pathname, fl)) and ("m" == fl[0]) and ("fits" == fl.split(".")[-1]):
            #print(fl)
            infilelist.append(fl)
    
    # TODO: Allow specification of ARC file in case multiple ARCS
    # Get arcfile from valid FITS files
    for fl in infilelist:
        with pyfits.open(fl) as hdu:
            if hdu['PRIMARY'].header['OBJECT'] == 'ARC':
                arcfile = fl
    
    if (arcfile == "") & (no_arc == False):
        cont = input("No arc found. Continue without it? [Y/n]\t")
        if cont.lower() in ("y", "yes", ""):
            no_arc = True
        elif cont.lower() in ("n", "no"):
            print("Quitting wav_replacement.py")
            sys.exit()
            
    if no_arc:
        arcfile = ""
    else:
        if display: print("Arcfile:\t\t\t" + arcfile)
        
        arcO, arcE = outfiles(arcfile, splitrow, save_O="oarc" + arcfile[-9:], save_E="earc" + arcfile[-9:], display=display)
        
        o_files.append("oarc" + arcfile[-9:])
        e_files.append("earc" + arcfile[-9:])
        if display: print(f"{arcfile} split into {o_files[-1]} & {e_files[-1]}")
        
    # Load existing FITS file and create copied single-extention FITS file
    if display: print("Pipeline m*.fits files found:\t" + str(infilelist))
    
    tar_list = []
    for i in infilelist:
        # Create blank FITS file for wavelength calibration and transformation
        if i != arcfile:
            tar_list.append(outfiles(i, splitrow, save_O=save_pref[0] + i[-9:], save_E=save_pref[1] + i[-9:], display=display))
            
            o_files.append(save_pref[0] + i[-9:])
            e_files.append(save_pref[1] + i[-9:])
        
        if display: print(f"{i} split into {o_files[-1]} & {e_files[-1]}")
    
    
    ####################################
            
    ## Crop fits files for reidentify to work
    # TODO: Find crop size dynamically / user given
    cut_size = 40
    
    
    # Obeam
    for i in o_files:
        with pyfits.open(i) as hdu:
            #print(hdu.info())
            #fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=[18, 8])
            #ax1.imshow(hdu[0].data, vmax=np.std(hdu[0].data), origin="lower")
            #ax2.imshow(hdu[0].data[0:-cut_size], vmax=np.std(hdu[0].data), origin="lower")
    
            hdu[0].data = hdu[0].data[0:-cut_size]
            hdu.writeto(i, overwrite=True)
    
    #Ebeam
    for i in e_files:
        with pyfits.open(i) as hdu:
            #print(hdu.info())
            #fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=[18, 8])
            #ax1.imshow(hdu[0].data, vmax=np.std(hdu[0].data), origin="lower")
            #ax2.imshow(hdu[0].data[0:-cut_size], vmax=np.std(hdu[0].data), origin="lower")
    
            hdu[0].data = hdu[0].data[cut_size:]
            hdu.writeto(i, overwrite=True)

    # Create text files listing o and e beams
    with open("o_frames", "w+") as f_o:
        for i in o_files:
            f_o.write(i + "\n")
            
    with open("e_frames", "w+") as f_e:
        for i in e_files:
            f_e.write(i + "\n")
    
    if display: print("Wrote o and e frames to o_frames and e_frames, respectfully.")
    print("Finished with no errors raised.")

#end region

#region Join IRAF files for PIPELINE
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

#end region