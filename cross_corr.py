#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "17.05.2022"

import sys
import os
import getopt
import itertools as iters

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from numpy.polynomial import chebyshev
from scipy import signal

class CrossCorrelate:
    """
        Cross correlate allows for comparing the extensions of multiple
        FITS files, or comparing the O and E beams of a single FITS file.

        Parameters
        ----------


        Returns
        -------


        Raises
        ------

        
    """
    def __init__(self,
                in1 : str,
                in2 : str = None,
                split_ccd : bool = True,
                cont : int = 11,
                offset : int = 0,
                save_name : str = None) -> None:
        
        self.invert = False
        self.wavUnits = "Å"
        self.wav1, self.spec1, self.bpm1 = self.checkLoad(in1)
        self.wav2, self.spec2, self.bpm2 = self.checkLoad(in2, in1)

        self.exts = self.spec1.shape[0]
        self.ccds = 1
        # Bounds shape [extensions, ccds, lower / upper bound]
        self.bounds1 = np.array([[[0, self.spec1[0].shape[-1]]], [[0, self.spec1[1].shape[-1]]]], dtype=int) 
        self.bounds2 = np.array([[[0, self.spec2[0].shape[-1]]], [[0, self.spec2[1].shape[-1]]]], dtype=int)
        if split_ccd:
            self.splitCCD()
        
        self.cont = cont
        if cont > 0:
            self.rmvCont()

        # Add an offset to the spectra to test cross correlation
        self.spec1 = np.insert(self.spec1, [0] * offset, self.spec1[:, :offset], axis=-1)[:, :self.spec1.shape[-1]]
        
        self.corrdb = []
        self.lagsdb = []
        self.correlate()

        self.save_name = save_name
        self.checkPlot()

        return

    def checkLoad(self, path1 : str, path2 : str = None) -> np.ndarray:
        # If the first path is invalid
        if not os.path.isfile(path1):
            # And the second path is not defined, raise an error
            if path2 == None: raise FileNotFoundError(f"{path1} is invalid")

            # Use the second path but swap the O and E beams
            path1 = path2
            self.invert = True

        # Load data
        with pyfits.open(path1) as hdu:
            spec = hdu['SCI'].data.sum(axis=1)
            wav = np.arange(spec.shape[-1]) * hdu["SCI"].header["CDELT1"] + hdu["SCI"].header["CRVAL1"]
            bpm = hdu['BPM'].data.sum(axis=1)

            if "Angstroms" not in hdu["SCI"].header["CTYPE1"]: self.wavUnits = hdu["SCI"].header["CTYPE1"]

        # Return data and implement swap if necessary
        return (wav, spec[::-1], bpm[::-1]) if self.invert else (wav, spec, bpm)

    def splitCCD(self) -> None:
        # Assumed BPM has a value of 2 near the center of each CCD (i.e. sum(bpm == 2) = count(ccd))
        self.ccds = sum(self.bpm1[0] == 2)
        # update bounds to reflect ccds
        self.bounds1 = np.zeros([self.exts, self.ccds, 2], dtype=int)
        self.bounds2 = np.zeros([self.exts, self.ccds, 2], dtype=int)
        
        # Get lower and upper bound for each ccd, save to bounds
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            mid1 = np.where(self.bpm1[ext] == 2)[0][ccd]
            mid2 = np.where(self.bpm2[ext] == 2)[0][ccd]

            # Lower bound, min non-zero
            lowb1 = max(mid1 - self.bpm1.shape[-1] // (self.ccds * 2), 0)
            uppb1 = min(mid1 + self.bpm1.shape[-1] // (self.ccds * 2), self.bpm1.shape[-1])

            # Upper bound, max bpm length
            lowb2 = max(mid2 - self.bpm2.shape[-1] // (self.ccds * 2), 0)
            uppb2 = min(mid2 + self.bpm2.shape[-1] // (self.ccds * 2), self.bpm2.shape[-1])

            self.bounds1[ext, ccd] = (lowb1, uppb1)
            self.bounds2[ext, ccd] = (lowb2, uppb2)
    
    def rmvCont(self) -> None:
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Get the range for current extension, ccd combination
            ccdBound1 = range(*self.bounds1[ext][ccd])
            ccdBound2 = range(*self.bounds2[ext][ccd])

            # Mask out the bad pixels for fitting continua
            okwav1 = np.where(self.bpm1[ext][ccdBound1] != 1)
            okwav2 = np.where(self.bpm2[ext][ccdBound2] != 1)
            
            # Define continua
            ctm1 = continuum(self.wav1[ccdBound1][okwav1], self.spec1[ext][ccdBound1][okwav1], deg=self.cont)
            ctm2 = continuum(self.wav2[ccdBound2][okwav2], self.spec2[ext][ccdBound2][okwav2], deg=self.cont)
            
            # Normalise spectra
            self.spec1[ext][ccdBound1] /= chebyshev.chebval(self.wav1[ccdBound1], ctm1)
            self.spec1[ext][ccdBound1] -= 1

            self.spec2[ext][ccdBound2] /= chebyshev.chebval(self.wav2[ccdBound2], ctm2)
            self.spec2[ext][ccdBound2] -= 1

            return
    
    def correlate(self) -> None:
        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Get the range for current extension, ccd combination
            ccdBound1 = range(*self.bounds1[ext][ccd])
            ccdBound2 = range(*self.bounds2[ext][ccd])

            # Add rows/cols for correlation and lags data
            if len(self.corrdb) <= ext:
                self.corrdb.append([])
                self.lagsdb.append([])
            if len(self.corrdb[ext]) <= ccd:
                self.corrdb[ext].append([])
                self.lagsdb[ext].append([])
            
            # Invert BPM (and account for 2 in BPM) to zero bad pixels
            sig1 = self.spec1[ext][ccdBound1] * abs(self.bpm1[ext][ccdBound1] * -1 + 1)
            sig2 = self.spec2[ext][ccdBound2] * abs(self.bpm2[ext][ccdBound2] * -1 + 1)

            # Finally(!!!) cross correlate signals
            corr = signal.correlate(sig1, sig2)
            corr /= np.max(corr) # Scales array so that the maximum correlation is at 1
            lags = signal.correlation_lags(sig1.shape[-1], sig2.shape[-1])
            
            self.corrdb[ext][ccd] = corr
            self.lagsdb[ext][ccd] = lags
        
        return
    
    def checkPlot(self, default_name : str = "OEcorr.pdf") -> None:
        # Plot
        fig = plt.figure()
        for ccd in range(self.ccds):
            # Add cross correlation to plots
            ax = fig.add_subplot(self.exts + 1, self.ccds, ccd + 1)
            for ext in range(self.exts):
                ax.plot(self.lagsdb[ext][ccd], self.corrdb[ext][ccd], label=f"max lag @ {self.lagsdb[ext][ccd][self.corrdb[ext][ccd].argmax()]}")


        for ext, ccd in iters.product(range(self.exts), range(self.ccds)):
            # Add wav, spec to plots
            ax = fig.add_subplot(self.exts + 1, self.ccds, (ext + 1) * self.exts + (ccd + 1))

            ax.plot(self.wav1[ext][ccd], self.spec1[ext][ccd], label="sig1")
            ax.plot(self.wav2[ext][ccd], self.spec2[ext][ccd], label="sig2")

            if ext == self.exts - 1: ax.set_xlabel(f"Wavelength ({self.wavUnits})")
            if ccd == 0: ax.set_ylabel("Normalised Intensity (Counts)")
        
        plt.legend()
        plt.show()


        # Handle do not save
        if self.save_name == None:
            return
        
        # Handle lazy save_name
        if self.save_name == ".":
            self.save_name = os.getcwd()

        # Handle save name directory, use a default name (overwrite with warning)
        if self.save_name[-1] == "/" or os.path.isdir(self.save_name):
            self.save_name += default_name
            print(f"Save name is a directory. Saving cross correlation results as {default_name}")

        # Check save location valid
        save_dir = os.path.expanduser('/'.join(self.save_name.split('/')[:-1]))
        if not os.path.isdir(save_dir):
            raise FileNotFoundError(f"The path ({save_dir}) does not exist")

        # Save
        if self.save_name != None:
            fig.savefig(fname=self.save_name)
        
        return

    # def self.continuum(self) -> np.array:
    #     pass



###############################################################

debug = True

# Functions
def continuum(w,c,deg=11,std=1.6,steps=5,pos=False, plot=False):
    
    if plot:
        plt.plot(w,c)
    
    nw = w.copy()
    nc = c.copy()
    
    for i in range(steps):
        
        p = chebyshev.chebfit(nw,nc,deg=deg)
        ch = chebyshev.chebval(nw,p)
        diff = nc - ch
        sigma = np.std(diff)
        
        if(pos):
            ok = np.where( np.abs(diff) > std*sigma)
        else:
            ok = np.where( diff > -1*std*sigma)
        
        nw = nw[ok]
        nc = nc[ok]
        
        if plot:
            plt.plot(w,chebyshev.chebval(w,p),label=str(i))
            
    if plot:
        plt.legend()
        plt.show()

        plt.plot(w,c/chebyshev.chebval(w,p))
        plt.show()
        
    return p


def load_data(a_dir=os.path.expanduser("~/polsalt-beta/masters_pol/sci/"),
              b_dir=os.path.expanduser("~/polsalt-beta/masters_pol/sci_comp/"),
              name_pref="e",
              llim=4000,
              ulim=9000):
    
    sci_list = {"ewhdu": [], "comp": []}
    
    #IRAF
    a_infilelist = []
    for fl in os.listdir(a_dir):
        if os.path.isfile(os.path.join(a_dir, fl)) and (name_pref == fl[0]) and ("fits" == fl.split(".")[-1]):
            a_infilelist.append(fl)
    
    
    for n, fname in enumerate(a_infilelist):#["ecwmxgbpP201912210018.fits", "ecwmxgbpP201912210019.fits", "ecwmxgbpP201912210020.fits", "ecwmxgbpP201912210021.fits"]):
        with pyfits.open(a_dir + fname) as hdu:
            sci_list["ewhdu"].append([0, 0, 0])
            for ext in [0, 1]:
                cs = hdu["SCI"].data[ext][0]
                ws = np.arange(cs.size)  * hdu["SCI"].header["CDELT1"] + hdu["SCI"].header["CRVAL1"]
    
                #limit wavelength range
                ok = np.where( (llim < ws) & (ws < ulim) )[0]
                ws = ws[ok]
                cs = cs[ok]
    
                #cut out chipgaps
                ok = np.where( (4500 > ws) | (ws > 4660) )[0]
                ws = ws[ok]
                cs = cs[ok]
    
                ok = np.where( (7750 > ws) | (ws > 7895) )[0]
                ws = ws[ok]
                cs = cs[ok]
                
                sci_list["ewhdu"][n][ext] = cs
                sci_list["ewhdu"][n][2]   = ws
                
    
    #PURE PIPELINE
    b_infilelist = []
    for fl in os.listdir(b_dir):
        if os.path.isfile(os.path.join(b_dir, fl)) and (name_pref == fl[0]) and ("fits" == fl.split(".")[-1]):
            b_infilelist.append(fl)
            
    
    for n, fname in enumerate(b_infilelist):#["ecwmxgbpP201912210018.fits", "ecwmxgbpP201912210019.fits", "ecwmxgbpP201912210020.fits", "ecwmxgbpP201912210021.fits"]):
        with pyfits.open(b_dir + fname) as hdu:
            sci_list["comp"].append([0, 0, 0])
            for ext in [0, 1]:
                cs = hdu["SCI"].data[ext][0]
                ws = np.arange(cs.size)  * hdu["SCI"].header["CDELT1"] + hdu["SCI"].header["CRVAL1"]
    
                #limit wavelength range
                ok = np.where( (llim < ws) & (ws < ulim) )[0]
                ws = ws[ok]
                cs = cs[ok]
    
                #cut out chipgaps
                ok = np.where( (4500 > ws) | (ws > 4660) )[0]
                ws = ws[ok]
                cs = cs[ok]
    
                ok = np.where( (7750 > ws) | (ws > 7895) )[0]
                ws = ws[ok]
                cs = cs[ok]
                
                sci_list["comp"][n][ext] = cs
                sci_list["comp"][n][2]   = ws

    return sci_list

# Compare E extension to O extension for Iraf and then for Pure wavelength.
def cross_corr_oe(sci_list, sub_cont, split_ccds, save_plots, save_path, offset=0):
    if not split_ccds:
        for d in ["ewhdu", "comp"]:
        
            fig, axs = plt.subplots(5, 1, figsize=(16, 12))
            fig.suptitle("IRAF WAV" if d == "ewhdu" else "Pure pipeline")
            fig.tight_layout(pad=1.5, h_pad=3)
            
            axs[4].set_title(f'{"IRAF WAV" if d == "ewhdu" else "Pure pipeline"} Cross-correlated signal')
            axs[4].set_xlabel('Lag')
            axs[4].margins(0, 0.1)
            
            for f in range(4):
                #Compare E against O
                ext1 = sci_list[d][f][0]
                ext2 = sci_list[d][f][1]
                w    = sci_list[d][f][2]
                
                #Subtract continuum
                if sub_cont:
                    p1    = continuum(w, ext1)
                    ext1  = ext1 / chebyshev.chebval(w, p1)
                    ext1 -= 1
                    
                    p2    = continuum(w, ext2)
                    ext2  = ext2 / chebyshev.chebval(w, p2)
                    ext2 -= 1
                    
                if offset != 0:
                    ext1  = np.insert(ext1, 0, ext1[0: offset])[0:w.size]
        
                corr = signal.correlate(ext2, ext1)
                lags = signal.correlation_lags(len(ext1), len(ext2))
                corr /= np.max(corr)
        
                axs[f].plot(w, ext1, label="E-beam")
                axs[f].plot(w, ext2, label="O-beam")
                
                axs[f].set_title(f"...{f + 18}")
                axs[f].legend()
        
                axs[4].plot(lags, corr, alpha=0.5, label=f"max lag @ {lags[np.where(corr == 1.0)][0]}")
                axs[4].axvline(lags[np.where(corr == 1.0)], c="orange", alpha=0.5)
                axs[4].legend()
                
                
                axs[f].margins(0, 0.1)
            
            if save_plots:
                file_name = ("bg_sub_" if sub_cont else "") + ("iraf_" if d == "ewhdu" else "salt_")
                fig.savefig(fname=f"cross_corr_results/{file_name}cross_correlation.pdf")
                
    elif split_ccds:
        for d in ["ewhdu", "comp"]:
            
            fig, axs = plt.subplots(5, 3, figsize=(16, 12))
            fig.suptitle("IRAF WAV" if d == "ewhdu" else "Pure pipeline")
            fig.tight_layout(pad=1.5, h_pad=3)
            
            for f in range(4):
                for ccd_num in range(3):
                    
                    #Compare E against O
                    ext1 = sci_list[d][f][0]
                    ext2 = sci_list[d][f][1]
                    w    = sci_list[d][f][2]
                    
                    if ccd_num == 0:
                        ccd = np.where(w < 4500)[0]
                    elif ccd_num == 1:
                        ccd = np.where((w > 4500) & (w < 7800))[0]
                    else:
                        ccd = np.where(w > 7800)[0]
                        
                    ext1 = ext1[ccd]
                    ext2 = ext2[ccd]
                    w    =    w[ccd]
        
                    #Subtract continuum
                    if sub_cont:
                        p1    = continuum(w, ext1)
                        ext1  = ext1 / chebyshev.chebval(w, p1)
                        ext1 -= 1
        
                        p2    = continuum(w, ext2)
                        ext2  = ext2 / chebyshev.chebval(w, p2)
                        ext2 -= 1
                        
                    if offset != 0:
                        ext1  = np.insert(ext1, 0, ext1[0: offset])[0:w.size]
        
                    corr = signal.correlate(ext2, ext1)
                    lags = signal.correlation_lags(len(ext1), len(ext2))
                    corr /= np.max(corr)
        
                    axs[f, ccd_num].plot(w, ext1, label="E-beam")
                    axs[f, ccd_num].plot(w, ext2, label="O-beam")
        
                    axs[f, ccd_num].set_title(f"...{f + 18}")
                    axs[f, ccd_num].legend()
        
                    axs[4, ccd_num].plot(lags, corr, alpha=0.5, label=f"max lag @ {lags[np.where(corr == 1.0)][0]}")#label=f"{round(corr[np.where(lags == 0.0)][0], 5)} @lag=0")
                    axs[4, ccd_num].axvline(lags[np.where(corr == 1.0)], c="orange", alpha=0.5)
                    axs[4, ccd_num].legend()
        
        
                    axs[f, ccd_num].margins(0, 0.1)
            
            if save_plots:
                file_name = "ccd_" + ("bg_sub_" if sub_cont else "") + ("iraf_" if d == "ewhdu" else "salt_")
                fig.savefig(fname=f"cross_corr_results/{file_name}cross_correlation.pdf")

# Compare E extensions for Iraf vs Pure wavelength and then O extensions.
def cross_corr_ab(sci_list, sub_cont, split_ccds, save_plots, save_path, offset=0):
    if not split_ccds:
        for f in range(4):
            fig, axs = plt.subplots(3, 1, figsize=(16, 8))
            fig.suptitle(f"IRAF vs Pure pipeline ...{f + 18}")
            fig.tight_layout(pad=1.5, h_pad=3)
            
            axs[2].set_title('Cross-correlated signal')
            axs[2].set_xlabel('Lag')
            axs[2].margins(0, 0.1)
            
            for i in [0, 1]:
                # Compare comp vs ewhdu
                ext1 = sci_list["ewhdu"][f][i]
                ext2 = sci_list["comp"][f][i]
                w1   = sci_list["ewhdu"][f][2]
                w2   = sci_list["comp"][f][2]
                
                #Subtract continuum
                if sub_cont:
                    p1    = continuum(w1, ext1)
                    ext1  = ext1 / chebyshev.chebval(w1, p1)
                    ext1 -= 1
                    
                    p2    = continuum(w2, ext2)
                    ext2  = ext2 / chebyshev.chebval(w2, p2)
                    ext2 -= 1
                    
                if offset != 0:
                    ext1  = np.insert(ext1, 0, ext1[0: offset])[0:w1.size]
        
                corr = signal.correlate(ext2, ext1)
                lags = signal.correlation_lags(len(ext1), len(ext2))
                corr /= np.max(corr)
        
                axs[i].plot(w1, ext1, label="IRAF")
                axs[i].plot(w2, ext2, label="COMP")
                axs[i].set_title(f"ext. {i}")
                axs[i].legend()
        
                axs[2].plot(lags, corr, label=f"ext {i}, max lag @ {lags[np.where(corr == 1.0)][0]}")# \n( {round(corr[np.where(lags == 0.0)][0], 5)} )")
                axs[2].axvline(lags[np.where(corr == 1.0)], c="orange")
                
                axs[2].legend()
                
                axs[i].margins(0, 0.1)
        
            if save_plots:
                file_name = ("bg_sub_" if sub_cont else "") + str(f + 18)
                fig.savefig(fname=f"cross_corr_results/{file_name}_cross_correlation.pdf")
        
    elif split_ccds:
        for f in range(4):
            fig, axs = plt.subplots(3, 3, figsize=(16, 8))
            fig.suptitle(f"IRAF vs Pure pipeline ...{f + 18}")
            fig.tight_layout(pad=1.5, h_pad=3)
            
            for i in [0, 1]:
                for ccd_num in [0, 1, 2]:
                    
                    # Compare comp vs ewhdu
                    ext1 = sci_list["ewhdu"][f][i]
                    ext2 = sci_list["comp"][f][i]
                    w1   = sci_list["ewhdu"][f][2]
                    w2   = sci_list["comp"][f][2]
                    
                    if ccd_num == 0:
                        ccd1 = np.where(w1 < 4500)[0]
                        ccd2 = np.where(w2 < 4500)[0]
                    elif ccd_num == 1:
                        ccd1 = np.where((w1 > 4500) & (w1 < 7800))[0]
                        ccd2 = np.where((w2 > 4500) & (w2 < 7800))[0]
                    else:
                        ccd1 = np.where(w1 > 7800)[0]
                        ccd2 = np.where(w2 > 7800)[0]
                        
                    ext1 = ext1[ccd1]
                    ext2 = ext2[ccd2]
                    w1 = w1[ccd1]
                    w2 = w2[ccd2]
        
                    #Subtract continuum
                    if sub_cont:
                        p1    = continuum(w1, ext1)
                        ext1  = ext1 / chebyshev.chebval(w1, p1)
                        ext1 -= 1
        
                        p2    = continuum(w2, ext2)
                        ext2  = ext2 / chebyshev.chebval(w2, p2)
                        ext2 -= 1
                        
                    if offset != 0:
                        ext1  = np.insert(ext1, 0, ext1[0: offset])[0:w1.size]
        
                    corr = signal.correlate(ext2, ext1)
                    lags = signal.correlation_lags(len(ext1), len(ext2))
                    corr /= np.max(corr)
        
                    axs[i, ccd_num].plot(w1, ext1, label="IRAF")
                    axs[i, ccd_num].plot(w2, ext2, label="COMP")
                    axs[i, ccd_num].set_title(f"ext. {i}")
                    axs[i, ccd_num].legend()
        
                    axs[2, ccd_num].plot(lags, corr, label=f"ext {i}, max lag @ {lags[np.where(corr == 1.0)][0]}")# \n( {round(corr[np.where(lags == 0.0)][0], 5)} )")
                    axs[2, ccd_num].axvline(lags[np.where(corr == 1.0)], c="orange")
        
                    axs[2, ccd_num].legend()
        
                    axs[i, ccd_num].margins(0, 0.1)
        
            if save_plots:
                file_name = "ccd_" + ("bg_sub_" if sub_cont else "") + str(f + 18)
                fig.savefig(fname=f"cross_corr_results/{file_name}_cross_correlation.pdf")


# Parse args
def main(argv):
    a_path    = ""
    b_path    = ""
    prefix    = ""
    mode      = ""
    cont      = False
    split     = False
    save_plot = False
    save_path = ""
    offset    = 0
    
    try:
        opts, args = getopt.getopt(argv, "hcda:b:p:m:s:o:", ["help",
                                                           "cont_subtract",
                                                           "distinct_ccds",
                                                           "a_path=", "b_path=", "prefix=",
                                                           "mode=",
                                                           "save_plots=",
                                                           "offset="])
        
    except getopt.GetoptError:
        print('cross_corr.py -h\n')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("USAGE:\tcross_corr.py -a <path-to-IRAF-data> -b <path-to-PURE-data> -m <oe | ab | both>")
            print("Optional:\n-p <prefix-of-data> -c (subtract continuum) -d (split ccd cross correlations) -o <offset to set if cross_corr working>")
            print("-s <save directory [a | b | directory]>")
         
            sys.exit()

        elif opt in ("-a", "--a_path"):
            if arg == ".":
                a_path = os.getcwd()
            
            elif os.path.isdir(os.path.expanduser(arg)):
                a_path = os.path.expanduser(arg)
                
            if debug: print(f"set a to {a_path}")
     
        elif opt in ("-b", "--b_path"):
            if arg == ".":
                b_path = os.getcwd()
            
            elif os.path.isdir(os.path.expanduser(arg)):
                b_path = os.path.expanduser(arg)
                
            if debug: print(f"set b to {b_path}")
            
                
        elif opt in ("-p", "--prefix"):
            prefix = arg
            
            if debug: print(f"set prefix to {prefix}")
         
        elif opt in ("-m", "--mode"):
            if arg == "both":
                mode = "both"
            elif arg == "oe":
                mode = "oe"
            elif arg == "ab":
                mode = "ab"
                
            if debug: print(f"set mode to {mode}")
          
        elif opt in ("-c", "--cont_subtract"):
            cont = True
            
            if debug: print(f"set cont to {cont}")
          
        elif opt in ("-d", "--distinct_ccds"):
            split = True
            
            if debug: print(f"set split to {split}")
          
        elif opt in ("-s", "--save_plots"):
            save_plot = True
          
            if arg == ".":
                save_path = os.getcwd()
            
            elif arg in ("a", "b"):
                save_path = arg
            
            elif os.path.isdir(os.path.expanduser(arg)):
                save_path = os.path.expanduser(arg)
                
            if debug: print(f"set save_path to {save_path}")
            
        elif opt in ("-o", "--offset"):
            offset = int(arg)
            
            if debug: print(f"set offset to {offset}")
    
    if a_path == "":
        a_path = os.getcwd()
    if b_path == "":
        b_path = os.getcwd()
        
    if save_path == "a":
        save_path = a_path
    elif save_path == "b":
        save_path = b_path
    elif save_path == "":
        save_path = os.getcwd()
        
    try:
        os.mkdir(os.path.join(save_path, "cross_corr_results"))
        
    except FileExistsError:
        print("File exists so older figures may be deleted")
    
    if prefix == "":
        sci_list = load_data(a_path, b_path)
    else:
        sci_list = load_data(a_path, b_path, prefix)
        
    if mode == "both":
        print("Running both")
        cross_corr_oe(sci_list, cont, split, save_plot, save_path, offset)
        cross_corr_ab(sci_list, cont, split, save_plot, save_path, offset)
        
    elif mode == "oe":
        print("Running O vs E")
        cross_corr_oe(sci_list, cont, split, save_plot, save_path, offset)
        
    elif mode == "ab":
        print("Running A vs B")
        cross_corr_ab(sci_list, cont, split, save_plot, save_path, offset)


if __name__ == "__main__":
    main(sys.argv[1:])