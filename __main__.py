#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Oct 18 2021

@author: Justin Cooper
"""

# General Imports
import os
import sys
import getopt

import split, join, cross_corr


argv = sys.argv[1:]
pathname = ""
verbose = False
mode = ""
no_arc = False
split_row = 517
prefix = ["obeam","ebeam"]

parent_folder = os.path.split(os.path.dirname(os.path.realpath(__file__)))[-1]

try:
    # TODO@JustinotherGitter: Add plot option and update docs
    opts, args = getopt.getopt(argv, "hm:d:vns:p:", ["help",
                                                    "mode=",
                                                    "dir=",
                                                    "verbose",
                                                    "no_arc",
                                                    "split_row=",
                                                    "prefix="])
    
except getopt.GetoptError:
    print(f"python .../{parent_folder} -m < split / join > -d <data-dir [.]>")
    print(f"python .../{parent_folder} -h")
    sys.exit(2)
    
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print("USAGE:")
        print(f"python .../{parent_folder} -m <split / join> -d <directory of data> -v -n -s <integer> -p <o_name,e_name>\n")

        print("Use -m/--mode to set whether splitting or joining FITS files")
        print("Use -d/--dir to give data directory, can be given as <.> for current terminal directory\n")

        print("OPTIONAL")
        print("Use -v/--verbose for verbose output, including plots")
        print("Use -n/--no_arc for no arc input")
        print("Use -s/--split_row <integer> to manually set the split row")
        print("Use -p/--prefix <[o_name,ename]> to set the target o/e beam fits prefixes")
        sys.exit()

    elif opt in ("-m", "--mode"):
        if arg.lower() in ("split", "s"):
            mode = "split"
            
        elif arg.lower() in ("join", "j"):
            mode = "join"
        
        elif arg.lower() in ("cross_corr", "c"):
            mode = "xcorr"
            
    elif opt in ("-d", "--dir"):
        if arg == ".":
            pathname = os.getcwd()
        
        elif os.path.isdir(os.path.expanduser(arg)):
            pathname = os.path.expanduser(arg)
            os.chdir(pathname)
        
        else:
            raise FileNotFoundError(f"The given directory, {arg}, does not exist.")

    elif opt in ("-v", "--verbose"):
        verbose = True
        
    elif opt in ("-n", "--no_arc"):
        no_arc = True
        
    elif opt in ("-s", "--split_row"):
        # TODO@JustinotherGitter: Find splitrow automatically
        try:
            split_row = int(arg)
        except ValueError:
            print("Split row usage: -s|--split_row <positive integer value>")
            
    elif opt in ("-p", "--prefix"):
        prefix = list(arg.strip("[]()").split(","))
        if len(prefix) < 2:
            raise IndexError(f"prefix usage: -p|--prefix <[o_name,e_name]>")

if pathname == "":
    pathname = os.getcwd()
    
if mode == "split":
    print("Running split")
    split.Split(data_dir=pathname, split_row=split_row, verbose=verbose, no_arc=no_arc, save_prefix=prefix).process()
    
elif mode == "join":
    print("Running join")
    join.Join(data_dir=pathname, split_row=split_row, verbose=verbose, no_arc=no_arc, save_prefix=prefix).process()

elif mode == "xcorr":
    print("Running cross correlation")
    cross_corr.CrossCorrelate().process() # TODO@JustinotherGitter: Add relevant args to Correlate object creation
    