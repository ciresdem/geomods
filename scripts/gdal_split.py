#!/usr/bin/env python
### gdal_split.py
##
## Copyright (c) 2018 - 2020 Matthew Love <matthew.love@colorado.edu>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## Split a gdal grid by value.
##
### Code:

import os
import sys
import numpy as np
from gdalconst import *
from osgeo import osr
from osgeo import gdal

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

import geomods.gdalfun

_version = "0.0.8"

_license = """
version %s
    """ %(_version)

_usage = """
gdal_split.py: Split the topo from a gdal file (>0)

usage: gdal_split.py [ -s [ args ] ] [ file ]

Options:
  file\t\tThe input DEM file-name

  -split\tValue to split [0]
  -help\t\tPrint the usage text
  -version\tPrint the version information

Example:
gdal_split.py input.tif
gdal_split.py input.tif -split -1

gdal_split.py v.%s 
""" %(_version)

if __name__ == '__main__':
    
    elev = None
    split_value = 0
    outFormat = 'GTiff'
    
    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '-s' or arg == '-split' or arg == '--split':
            split_value = float(sys.argv[i+1])
            i = i + 1

        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '-version' or arg == '--version':
            print('gdal_split.py v.%s' %(_version))
            print(_license)
            sys.exit(1)

        elif elev is None:
            elev = arg

        else:
            print(_usage)
            sys.exit(1)

        i = i + 1

    if elev is None:
        print(_usage)
        sys.exit(1)

    if not os.path.exists(elev):
        print("Error: %s is not a valid file" %(elev))
    else: geomods.gdalfun.split(elev)

### End
