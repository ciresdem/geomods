#!/usr/bin/env python
### gdal_crop.py
##
## Copyright (c) 2018 - 2020 CIRES Coastal DEM Team
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
## crop a gdal grid by the nodata value
##
### Code:

import os
import sys

import gdal
import geomods

_version = '0.0.4'

_usage = '''gdal_crop.py ({}): crop a gdal grid by the nodata value

usage: gdal_crop.py [ file ]

 Options:
  file\t\tThe input DEM file-name

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_crop.py input.tif

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

if __name__ == '__main__':    
    elev = None
    split_value = 0

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '-version' or arg == '--version':
            print('gdal_crop.py, version {}\n{}'.format(_version, geomods._license))
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
        print('Error: %s is not a valid file' %(elev))
    else:
        output_name=elev[:-4] + '_crop.tif'

        out_array, out_config = geomods.gdalfun.gdal_crop(elev)
        outsize = out_array.shape

        #Export Tif
        ods = gdal.GetDriverByName('GTiff').Create(output_name, outsize[1], outsize[0], 1, out_config['dt'])
        ods.SetGeoTransform(out_config['geoT'])
        ods.SetProjection(out_config['proj'])
        ods.GetRasterBand(1).SetNoDataValue(out_config['ndv'])
        ods.GetRasterBand(1).WriteArray(out_array)

        ods = None
### End
