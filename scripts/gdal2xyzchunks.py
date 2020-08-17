#!/usr/bin/env python


### Code:
import os
import sys
from geomods import waffles

_version = '0.0.2'
_usage = '''gdal2xyzchunks.py ({}): chunk a gdal grid and output xyz data

usage: gdal2xyzchunks.py [ FILE ] [ OPTIONS ]

  file\t\tThe input grid file-name

 Options:
  -c, --chunk\tThe chunk size
      output will be in chunks of `chunk X chunk` cells.
  -E, --inc\tThe increment to resample bag to in degrees.
  -V, --vdatum\tThe output vertical datum input:output string
      e.g. -V 'mllw:navd88'
  -P, --epsg\tThe output EPSG
  -D, --datalist\tThe chunked xyz datalist

  --help\tPrint the usage text
  --version\tPrint the version information

 Examples:
 % bag2xyzchunks.py input_mllw.bag --chunk 1000 --increment 0.00003064

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_version)

if __name__ == '__main__':
    src_fn = None
    chunk_value = 1000
    inc = None
    vdatum = None
    epsg = None
    datalist = 'gdal_chunks.datalist'
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-C' or arg == '-c' or arg == '-chunk' or arg == '--chunk':
            try:
                chunk_value = int(sys.argv[i + 1])
            except: pass
            i = i + 1
        elif arg == '-E' or arg == '-e' or arg == '-inc' or arg == '--inc':
            try:
                inc = float(sys.argv[i + 1])
            except: pass
            i = i + 1
        elif arg == '-V' or arg == '-v' or arg == '-vdatum' or arg == '--vdatum':
            vdatum = sys.argv[i + 1]
            i = i + 1
        elif arg == '-D' or arg == '-d' or arg == '-datalist' or arg == '--datalist':
            datalist = sys.argv[i + 1]
            i = i + 1
        elif arg == '-P' or arg == '-p' or arg == '-epsg' or arg == '--epsg':
            try:
                epsg = int(sys.argv[i + 1])
            except: pass
            i = i + 1
        elif arg == '-help' or arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(-1)
        elif arg == '-version' or arg == '--version':
            sys.stdout.write('{}\n'.format(_version))
            sys.exit(-1)
        elif src_fn is None: src_fn = arg
        else:
            sys.stderr.write(_usage)
            sys.exit(-1)
        i = i + 1

    if src_fn is None:
        waffles.echo_error_msg('you must enter an input file')
        sys.stderr.write(_usage)
        sys.exit(1)

    if not os.path.exists(src_fn):
        waffles.echo_error_msg('{} is not valid'.format(src_fn))
    else: waffles.gdal2xyz_chunks(src_fn, chunk_value, inc, epsg, vdatum, datalist, True)

### End
