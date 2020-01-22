#!/usr/bin/env python
### gdal_crop.py
##
## Copyright (c) 2018 Matthew Love <matthew.love@colorado.edu>
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
## apply a mask to a gdal file
##
### Code:

import os
import sys
import numpy as np
from osgeo import gdal

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

_version = "0.0.3"

_license = """
version %s
    """ %(_version)

_usage = """
gdal_mask.py: apply a mask to a gdal dataset

usage: gdal_mask.py [ -mask [ args ] ] [ in-file ] [ out-file ]

Options:
  in-file\tThe input DEM file-name
  out-file\tThe output DEM file-name

  -mask\t\tThe mask grid file-name

  -help\t\tPrint the usage text
  -version\tPrint the version information

Example:
gdal_mask.py input.tif

gdal_mask.py v.%s 
""" %(_version)

#--
def writeArray(outArray, outfile, outdriver, ycount, xcount, gt, ndata, overwrite, verbose):
    # Create the output GDAL Raster
    if outdriver == "AAIGrid":
        driver = gdal.GetDriverByName("MEM")
    else:
        driver = gdal.GetDriverByName(outdriver)

    if os.path.exists(outfile):
        if not overwrite:
            sys.exit("Error - %s already exists!  Use -O or --overwrite to force an \
                    overwrite of the output file." %(outfile))
        else:
            driver.Delete(outfile)
            if verbose:
                print "Warning - Overwriting file %s " %(outfile)

    dst_ds = driver.Create(outfile, xcount, ycount, 1, gdal.GDT_Float64)

    if dst_ds is None:
        sys.exit("failed to open output file...%s" %(outfile))

    dst_ds.SetGeoTransform(gt)
    dst_band = dst_ds.GetRasterBand(1)

    dst_band.SetNoDataValue( ndata )

    # Write Numpy Array to the GDAL Raster Band1
    dst_band.WriteArray(outArray)

    if outdriver == "AAIGrid":
        driver = gdal.GetDriverByName(outdriver)
        dst_ds_aii = driver.CreateCopy(outfile, dst_ds)

    dst_ds = None
    dst_ds_aii = None

#--

def maskArray(srcfile, srcfile2):

    ds = gdal.Open(srcfile2)
    gt2 = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    ndata = band.GetNoDataValue()
    dsArray2 = ds.ReadAsArray()

    ##maskArray = dsArray2 == 1
    maskArray = dsArray2 != ndata
    dsArray2 = None

    ds = gdal.Open(srcfile)
    gt = ds.GetGeoTransform()
    dsArray = ds.ReadAsArray()

    np.place(dsArray, maskArray, ndata)
    
    maskArray = None
    ds = None

    return dsArray, gt, ndata

if __name__ == "__main__":

    inmask=None
    ingrd=None
    outgrd=None
    overwrite = True
    verbose=False
    out_nodata=None
    outdriver = "GTiff"
    
    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-mask':
            inmsk = argv[i+1]
            i = i + 1

        elif arg == '-d_nodata':
            out_nodata = float(argv[i+1])
            i = i + 1
            
        elif arg == '-verbose':
            verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(1)

        elif ingrd is None:
            ingrd = arg

        elif outgrd is None:
            outgrd = arg

        else:
            print(_usage)
            sys.exit(1)

        i = i + 1

    if ingrd is None or outgrd is None:
        print(_usage)
        sys.exit(0)

    progress( 0.0 )
    outArray, gt, ndata = maskArray(ingrd, inmsk)
    progress( 0.5 )
    writeArray(outArray, outgrd, outdriver, outArray.shape[0], outArray.shape[1], gt, ndata, overwrite, verbose)
    progress( 1.0 )

### END
