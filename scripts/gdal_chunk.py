#!/usr/bin/env python
### gdal_chunk.py
##
## Copyright (c) 2019 - 2020 CIRES Coastal DEM Team
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
## chunk a gdal grid
##
### Code:

import os
import sys

import geomods

_version = '0.0.2'

_usage = '''gdal_chunk.py ({}): chunk a gdal grid

usage: gdal_chunk.py [ file ]

  file\t\tThe input grid file-name

 Options:
  -c, --chunk\tThe chunk size
      output will be in chunks of `chunk X chunk` cells.

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % gdal_chunk.py input.tif --chunk 1000

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

if __name__ == '__main__':
    src_fn = None
    chunk_value = 1000

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '-c' or arg == '-chunk' or arg == '--chunk':
            try:
                chunk_value = int(sys.argv[i + 1])
            except: pass
            i = i + 1

        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '-version' or arg == '--version':
            print('gdal_chunk.py, version {}\n{}'.format(_version, geomods._license))
            sys.exit(1)

        elif src_fn is None:
            src_fn = arg

        else:
            print(_usage)
            sys.exit(1)

        i = i + 1

    if src_fn is None:
        print(_usage)
        sys.exit(1)

    if not os.path.exists(src_fn):
        print('Error: {} is not a valid file'.format(src_fn))
    else: geomods.gdalfun.chunks(src_fn, chunk_value)

### End
