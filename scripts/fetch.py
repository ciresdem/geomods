#!/usr/bin/env python
### fetch.py
##
## Copyright (c) 2012 - 2020 Matthew Love <matthew.love@colorado.edu>
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
### Code:

import os
import sys
import time
import threading

import osgeo.ogr as ogr
import geomods

_version = '0.1.0'

fetch_infos = { 
    'dc':[lambda x, c: geomods.fetches.dc(x, callback = c), 'digital coast'],
    'nos':[lambda x, c: geomods.fetches.nos(x, callback = c), 'noaa nos bathymetry'],
    'mb':[lambda x, c: geomods.fetches.mb(x, callback = c), 'noaa multibeam'],
    'gmrt':[lambda x, c: geomods.fetches.gmrt(x, callback = c), 'gmrt'],
    'srtm':[lambda x, c: geomods.fetches.srtm_cgiar(x, callback = c), 'srtm from cgiar'],
    'charts':[lambda x, c: geomods.fetches.charts(x, callback = c), 'noaa nautical charts'],
    'tnm':[lambda x, c: geomods.fetches.tnm(x, callback = c), 'the national map'],
    'ngs':[lambda x, c: geomods.fetches.ngs(x, callback = c), 'ngs monuments'],
    'usace':[lambda x, c: geomods.fetches.usace(x, callback = c), 'usace bathymetry'] }    

def fetch_desc(x):
    fd = []
    for key in x: 
        fd.append('%-10s\t\t%s' %(key, x[key][1]))
    return fd

_license = '''
fetch.py version %s
Copyright (c) 2012 - 2020 Matthew Love <matthew.love@colorado.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
''' %(_version)

_usage = '''fetch.py (%s): Fetch geographic elevation data.

usage: fetch.py [ -fhlpRuv [ args ] ] module(s) ...

Modules:
  %s

Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -l, --list-only\tOnly fetch a list of surveys in the given region.
  -f, --filter\t\tSQL style attribute filter; use ? to see available field names.
  -p, --process\t\tProcess fetched data to ASCII XYZ format.

  --update\t\tUpdate the stored list of surveys.
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Examples:
 %% fetch.py nos charts -R -90.75/-88.1/28.7/31.25 -f "Date > 2000"
 %% fetch.py --update

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>''' %(_version, '\n  '.join(fetch_desc(fetch_infos)))


## =============================================================================
##
## Run fetch from command-line
##
## =============================================================================
def main():
    extent = None
    poly = False
    want_list = False
    want_update = False
    want_proc = False
    f = []
    i = 1
    fetch_class = []
    these_regions = []
    status = 0

    ## process Command-Line
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '--region' or arg == '-R':
            extent = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            iregion = str(arg[2:])

        elif arg == '--list-only' or arg == '-l':
            want_list = True

        elif arg == '--filter' or arg == '-f':
            f.append(sys.argv[i + 1])
            i = i + 1

        elif arg == '--process' or arg == '-p':
            want_proc = True

        elif arg == '--update' or arg == '-u':
            want_update = True

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print(_license)
            sys.exit(1)

        else:
            fetch_class.append(sys.argv[i])

        i = i + 1

    ## if no fetch source is specified, fetch from them all.
    if len(fetch_class) == 0:
        fetch_class = fetch_infos.keys()

    ## print the available fields for filtering
    if len(f) > 0 and f[0] == '?':
        for fc in fetch_class:
            fl = fetch_infos[fc][0](None, None)
            if fl._ref_vector is not None:
                print('%s: %s' %(fc, geomods.gdalfun.ogr_get_fields(fl._ref_vector)))
        sys.exit(1)

    ## update the reference vector
    if want_update: 
        for fc in fetch_class:
            pb = geomods.clis.prog_bar('updating %s reference vector' %(fc))
            fl = fetch_infos[fc][0](None, pb)
            try:
                fl_update = threading.Thread(target = fl._update)
                fl_update.start()

                while True:
                    time.sleep(5)
                    pb._update()
                    if not fl_update.is_alive():
                        break
            except: status = -1

            pb._end(status)
        sys.exit(status)

    if extent is None:
        print(_usage)
        sys.exit(1)

    ## process input region(s)
    tw = geomods.clis.prog_bar("processing region(s)")
    try: 
        these_regions = [geomods.regions.region(extent)]
    except:
        if os.path.exists(extent):
            _poly = ogr.Open(extent)
            _player = _poly.GetLayer(0)
            for pf in _player:
                _pgeom = pf.GetGeometryRef()
                these_regions.append(geomods.regions.region("/".join(map(str, _pgeom.GetEnvelope()))))

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
            status = -1

    tw.pm = 'processing %s region(s)' %(len(these_regions))
    tw._end(status)

    ## fetch some data
    for this_region in these_regions:
        for fc in fetch_class:
            pb = geomods.clis.prog_bar('initializing %s' %(fc))
            fl = fetch_infos[fc][0](this_region, pb)
            pb._update()

            #try: 
            fl.search_gmt(f)
            #except: pass
            pb._end(fl._status)

            pb = geomods.clis.prog_bar('fetching %s data' %(fc))

            if want_list: 
                fl.print_results()
            else:
                fl_fetch = threading.Thread(target = fl.fetch_results)
                fl_fetch.start()

                while True:
                    time.sleep(5)
                    pb._update()
                    if not fl_fetch.is_alive():
                        break

            pb._end(fl._status)

            if want_proc:
                pb = geomods.clis.prog_bar('processing %s data' %(fc))
                fl_proc = threading.Thread(target = fl.proc_results)

                fl_proc.start()

                while True:
                    time.sleep(5)
                    pb._update()
                    if not fl_proc.is_alive():
                        break

                pb._end(fl._status)

if __name__ == '__main__':
    main()

### End
