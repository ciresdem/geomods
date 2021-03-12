#!/bin/env python
### gdal_gquery.py
##
## Copyright (c) 2018 - 2021 CIRES Coastal DEM Team
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
## transform z values in a gdal grid using an associated transformation grid
##
### Code:

import os
import sys
import gdal

## ==============================================
## import geomods
## ==============================================
from geomods import utils
from geomods import regions
from geomods import gdalfun

_version = '0.0.2'
_usage = '''gdal_gquery.py ({}): transform band data

usage: gdal_gquery.py [ src_gdal src_trans [ OPTIONS ] ]

  src_gdal\tThe input raster file-name
  src_trans\tThe input transformation raster file-name

 Options:

  -o, --output\tGenerate new grid as output.
  -p, --op\ttransformation operation [d/m/...]
  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_gquery.py input.tif vA2vB.tif -o input_vB.tif -p 'd'

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == '__main__':    
    i_ds = None
    d_ds = None
    t_ds = None
    d_form = 'xyd'
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-o' or arg == '--output':
            d_ds = sys.argv[i+1]
            i = i + 1
        elif arg == '-p' or arg == '--op':
            d_form = 'xy' + sys.argv[i+1]
            i = i + 1
        elif i_ds is None:
            i_ds = arg
        elif t_ds is None:
            t_ds = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(1)
        i = i + 1
        
    if i_ds is None or not os.path.exists(i_ds):
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter a valid input raster file')
        sys.exit(1)

    if t_ds is None or not os.path.exists(t_ds):
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must enter a valid input raster file')
        sys.exit(1)

    if d_ds is None:
        d_ds = i_ds.split('.')[0] + '_trans.' + i_ds.split('.')[-1]

    ds = gdal.Open(i_ds)
    ds_config = gdalfun.gdal_gather_infos(ds)
    ds_region = regions.gt2region(ds_config)
    ds_inc = ds_config['geoT'][1]

    gp = gdalfun.gdal_parse(ds)
    gq = gdalfun.gdal_yield_query(gp, t_ds, d_form)
    out, status = gdalfun.gdal_xyz2gdal(gq, d_ds, ds_region, ds_inc, mode='m')
    
    ds = None
## End
