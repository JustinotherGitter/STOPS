#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 13:53:44 2021

@author: justin
"""
# General Imports
import sys
import os
import getopt
from copy import deepcopy

from astropy.io import fits as pyfits
import numpy as np
from numpy.polynomial.chebyshev import chebgrid2d as chebgrid2d
import lacosmic

import split
import join
# Polsalt Imports, PYTHON 2 edited to PYTHON 3
from polsalt.specpolpy3 import read_wollaston
from polsalt.specpolpy3 import outfiles, split_sci


def main(argv):
    pathname = ""
    verbose = False
    mode = ""
    no_arc = False
    split_row = 517
    prefix = ["obeam", "ebeam"]
    
    try:
        # TODO: Add plot option
        opts, args = getopt.getopt(argv, "hvm:d:ns:p:", ["help",
                                                       "verbose",
                                                       "mode=",
                                                       "dir=",
                                                       "no_arc",
                                                       "split_row=",
                                                       "prefix="])
        
    except getopt.GetoptError:
        print('wav_replacement.py -m < split / join > -d <data-dir [.]>')
        sys.exit(2)
        
    for opt, arg in opts:
      if opt in ("-h", "--help"):
         print("USAGE:")
         print('wav_replacement.py -m <split / join> -p <path-to-data> -v -n -s\n')
         print('Directory can be given as "." for current terminal directory\n')
         print('Use -v/--verbose for verbose output, including plots')
         print('Use -n/--no_arc for no arc input')
         print('Use -s/--split_row <integer> to set the split row')
         print('Use -p/--prefix <[o_name,ename]> to set the target o/e beam fits prefixes')
         sys.exit()
         
      elif opt in ("-v", "--verbose"):
          verbose = True
         
      elif opt in ("-m", "--mode"):
         if arg.lower() == "split":
             mode = "split"
             
         elif arg.lower() == "join":
             mode = "join"
             
      elif opt in ("-d", "--dir"):
         if arg == ".":
            pathname = os.getcwd()
            
         elif os.path.isdir(os.path.expanduser(arg)):
            pathname = os.path.expanduser(arg)
            os.chdir(pathname)
            
         else:
            raise FileNotFoundError(f"The given directory, {arg}, does not exist.")
            
      elif opt in ("-n", "--no_arc"):
         no_arc = True
         
      elif opt in ("-s", "--splitrow"):
          # TODO@JustinotherGitter: Pass in splitrow automatically
         try:
             split_row = int(arg)
         except ValueError:
             print("Split row usage: -s|-split_row < integer value >")
             
      elif opt in ("-p", "--prefix"):
         prefix = list(arg.strip("[]()").split(",").strip())
    
    if pathname == "":
        pathname = os.getcwd()
        
    if mode == "split":
        print("Running split")
        split.main(pathname, split_row, verbose, no_arc, prefix)
        
    elif mode == "join":
        print("Running join")
        join.main(pathname, split_row, verbose, no_arc, prefix)

if __name__ == "__main__":
    main(sys.argv[1:])