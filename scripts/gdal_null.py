#!/usr/bin/env python
### gdal_null.py
##
## Copyright (c) 2018, 2019 Matthew Love <matthew.love@colorado.edu>
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
## Create a nodata grid
##
### Code:

import sys
import os
import numpy as np
from osgeo import gdal

try:
    progress = gdal.TermProgress_nocb
except:
    progress = gdal.TermProgress

_version = "0.1.6"

_usage = """
usage: gdal_null.py [-region xmin xmax ymin ymax] [-cell_size value]
                    [-t_nodata value] [-d_format grid-format] [-overwrite]
                    [-copy grid] [-verbose] output_grid

gdal_null v.%s 
""" %(_version)

def verbosePrint(xcount, ycount, extent, cellsize, outf):
    print '--'
    print 'xcols: ',xcount,' ycols: ',ycount
    print 'xmax: ',extent[1],'xmin: ',extent[0]
    print 'ymax: ',extent[3],'ymin: ',extent[2]
    print 'cellsize: ',cellsize
    print 'output grid format: ',outf
    print '--'

def createNullCopy(srcfile, outfile, nodata, outformat, verbose, overwrite):
    ds = gdal.Open(srcfile)
    gt = ds.GetGeoTransform()
    #print ds.RasterXSize,ds.RasterYSize
    
    progress( 0.0 )
    
    (ycount,xcount) = ds.RasterYSize,ds.RasterXSize

    if nodata is None:
        nodata = ds.GetRasterBand(1).GetNoDataValue()
    if nodata is None:
        nodata = -9999

    #if verbose:
    #    print nodata

    ds = None
    dsArray = np.zeros([ycount,xcount])

    progress( 0.2 )

    dsArray[:] = float(nodata)

    progress( 0.4 )

    # Create the output GDAL Raster
    driver = gdal.GetDriverByName("GTiff")

    if os.path.exists(outfile):
        if not overwrite:
            progress( 1.0 )
            sys.exit("gdal_null: Error - %s already exists!  Use -overwrite to force an \
overwrite of the output file." %(outfile))
        else:
            driver.Delete(outfile)
            if verbose:
                print >> sys.stderr, "gdal_null: Warning - Overwriting file %s " %(outfile)

    dst_ds = driver.Create(outfile, xcount, ycount, 1, gdal.GDT_Float32)

    if dst_ds is None:
        progress( 1.0 )
        sys.exit("gdal_null: failed to open output file...%s" %(outfile))

    dst_ds.SetGeoTransform(gt)
    dst_band = dst_ds.GetRasterBand(1)

    progress( 0.75 )

    dst_band.SetNoDataValue(nodata)

    # Write Numpy Array to the GDAL Raster Band
    dst_band.WriteArray(dsArray)
    dst_ds = None
    progress( 1.0 )
    
def createGrid(outfile, extent, cellsize, nodata, outformat, verbose, overwrite):

    # if verbose:
    #     print extent

    # Set the rows and columns for the output grid
    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize/cellsize)+1
    ycount = int(ysize/cellsize)+1

    if verbose:
        verbosePrint(xcount, ycount, extent, cellsize, outformat)

    # Create the output GDAL Raster
    if outformat == "AAIGrid":
        driver = gdal.GetDriverByName("MEM")
    else:
        driver = gdal.GetDriverByName(outformat)

    if os.path.exists(outfile):
        if not overwrite:
            sys.exit("gdal_null: Error - %s already exists! Use -overwrite to force an \
overwrite of the output file." %(outfile))
        else:
            driver.Delete(outfile)
            if verbose:
                print >> sys.stderr, "gdal_null: Warning - Overwriting file %s " %(outfile)

    dst_ds = driver.Create(outfile, xcount, ycount, 1, gdal.GDT_Float32)
    if dst_ds is None:
        sys.exit("gdal_null: failed to open output file...%s" %(outfile))

    gt = (extent[0],cellsize,0,extent[3],0,(cellsize*-1.))
    dst_ds.SetGeoTransform(gt)
    dst_band = dst_ds.GetRasterBand(1)

    if nodata is None:
        nodata = 0

    dst_band.SetNoDataValue(nodata)

    # Create a Numpy Array filled with 0's
    nullArray = np.zeros( (ycount, xcount) )

    nullArray[nullArray==0]=nodata

    # Write Numpy Array to the GDAL Raster Band
    dst_band.WriteArray(nullArray)

    if outformat == "AAIGrid":
        driver = gdal.GetDriverByName(outformat)
        dst_ds_aai = driver.CreateCopy(outfile, dst_ds)
    
    # Clear the GDAL Raster(s)x
    dst_ds = None
    dst_ds_aai = None

    if verbose:
        print '\ngdal_null: Grid created.'

if __name__ == '__main__':

    extent = None
    cellsize = 1
    d_format = "GTiff"
    cpgrd = None
    overwrite = False
    verbose = False
    nodata = None
    output = None

    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-region' or arg == '-r' or arg == '--region':
            extent = (float(argv[i+1]),float(argv[i+2]),
                      float(argv[i+3]),float(argv[i+4]))
            i = i + 4

        elif arg == '-cell_size' or arg == '-s' or arg == '--cell_size':
            cellsize = float(argv[i+1])
            i = i + 1

        elif arg == '-d_format' or arg == 'd' or arg == '--d_format':
            d_format = str(argv[i+1])
            i = i + 1

        elif arg == '-t_nodata' or arg == '-t' or arg == '--t_nodata':
            nodata = float(argv[i+1])
            i = i + 1

        elif arg == '-copy' or arg == '-c' or arg == '--copy':
            cpgrd = argv[i+1]
            i = i + 1
            
        elif arg == '-overwrite' or arg == '--overwrite':
            overwrite = True
            
        elif arg == '-verbose' or arg == '--verbose':
            verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(1)

        elif output is None:
            output = arg

        else:
            print(_usage)
            sys.exit(1)

        i = i + 1

    if output is None:
        print(_usage)
        sys.exit(0)
    if extent == None:
        extent = '1'

    #Run the program given the user input
    if cpgrd is not None:
        createNullCopy(cpgrd, output, nodata, d_format, verbose, overwrite)
    else:
        createGrid(output, extent, cellsize, nodata, d_format, verbose, overwrite)

### END
