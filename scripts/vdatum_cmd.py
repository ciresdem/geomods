#!/usr/bin/env python
### vdatum_cmd.py
##
## Copyright (c) 2019 - 2020 Matthew Love <matthew.love@colorado.edu>
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

_usage = '''vdatum_cmd.py ({}): run NOAAs vdatum

usage: vdatum [ file ]

  file\t\tThe input file to translate

 Options:

  -i the input vertical datum (default mllw:m:height) while
  -o specifies the output vertical datum (defafult navd88:m:height).

  -r specifies the input horizontal datum (default nad83:geo:deg) while
  -z specifies the output horizontal datum (default nad83:geo:deg).

  -g specifies the vdatum region (default 3)

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % vdatum_cmd.py elev_data.xyz -i mllw -o navd88 -g 1
 % vdatum_cmd.py elev_data.xyz -i lmsl:ft:sounding -o navd88:m:height -g 1
 % vdatum_cmd.py elev_data.xyz -i navd88:m:height -o navd88:m:height -r nad83:utm:ft:10 -z nad83:geo:deg

See: vdatum.noaa.gov/docs/userguide_cmd.html for specifics on cli parameters.
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

def main():
    src_fn = None

    ihorz = 'NAD83_2011'
    ohorz = 'NAD83_2011'
    ivert = 'navd88:m:height'
    overt = 'mhw:m:height'
    region = '3'

    i = 1

    argv = sys.argv
    while i < len(argv):
        arg = argv[i]

        if arg == '-i' or arg == '--ivert':
            try:
                ivert = argv[i + 1]
            except: pass
            i = i + 1

        if arg == '-o' or arg == '--overt':
            try:
                overt = argv[i + 1]
            except: pass
            i = i + 1

        if arg == '-r' or arg == '--ihorz':
            try:
                ihorz = argv[i + 1]
            except: pass
            i = i + 1

        if arg == '-z' or arg == '--ohorz':
            try:
                ohorz = argv[i + 1]
            except: pass
            i = i + 1

        if arg == '-e' or arg == '--region':
            try:
                region = argv[i + 1]
            except: pass
            i = i + 1

        elif arg == '-help' or arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '-version' or arg == '--version':
            print('vdatum_cmd.py, version {}\nVDatum, version {}\n{}'.format(_version, geomods.utils.vdatum()._version, geomods._license))
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
    else: 
        #vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:space,0,1,2:{}:result region:{}'.format(ihorz, ivert, ohorz, overt, src_fn, region)
        #geomods.utils.run_vdatum(vdc)
        geomods.utils.vdatum().run_vdatum(src_fn)

if __name__ == '__main__':
    main()

### End
