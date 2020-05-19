### configs.py
##
## Copyright (c) 2010 - 2020 CIRES Coastal DEM Team
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
import ConfigParser

import utils

## =============================================================================
##
## config-file - configs.py
##
## The waffles config file holds system information,
## such as host system, python version, external programs and their paths, etc.
##
## This needs to be updated for python3 and to work better in general!
##
## =============================================================================

CONFIG_FILE = os.path.expanduser('~/geomods.ini')
_waff_co = ConfigParser.ConfigParser()

def check_config(recheck = False, verbose = False):

    waffles_co = {}
    vdatum_path = None
    gmt_vers = None
    mbgrid_vers = None
    gdal_vers = None
    bounds_vers = None
    
    _waff_co.read(CONFIG_FILE)
    
    try:
        vdatum_path = _waff_co.get('VDATUM', 'jar')
        vdatum_path = None if vdatum_path == 'None' else vdatum_path
    except: _waff_co.add_section('VDATUM')

    if vdatum_path is None or recheck:
        vd = vdatum(verbose=verbose)
        vdatum_path = vd.vdatum_path
        vdatum_path = None if vdatum_path == 'None' else vdatum_path
        _waff_co.set('VDATUM', 'jar', vdatum_path)
        
    waffles_co['VDATUM'] = vdatum_path
        
    try:
        gmt_vers = _waff_co.get('GMT', 'vers')
        gmt_vers = None if gmt_vers == 'None' else gmt_vers
    except: _waff_co.add_section('GMT')

    if gmt_vers is None or recheck:
        gmt_vers = cmd_check('gmt', 'gmt --version')
        gmt_vers = None if gmt_vers == 'None' else gmt_vers
        _waff_co.set('GMT', 'vers', gmt_vers)
        
    waffles_co['GMT'] = gmt_vers
        
    try:
        gdal_vers = _waff_co.get('GDAL', 'vers')
        gdal = None if gdal_vers == 'None' else gdal_vers
    except: _waff_co.add_section('GDAL')

    if gdal_vers is None or recheck:
        gdal_vers = cmd_check('gdal-config', 'gdal-config --version')
        gdal = None if gdal_vers == 'None' else gdal_vers
        waff_co.set('GDAL', 'vers', gdal_vers)
        
    waffles_co['GDAL'] = gdal_vers
    
    try:
        mbgrid_vers = _waff_co.get('MBGRID', 'vers')
        mbgrid_vers = None if mbgrid_vers == 'None' else mbgrid_vers
    except: _waff_co.add_section('MBGRID')

    if mbgrid_vers is None or recheck:
        mbgrid_vers = cmd_check('mbgrid', 'mbgrid -version | grep Version')
        mbgrid_vers = None if mbgrid_vers == 'None' else mbgrid_vers
        _waff_co.set('MBGRID', 'vers', mbgrid_vers)
        
    waffles_co['MBGRID'] = mbgrid_vers
    
    try:
        bounds_vers = _waff_co.get('BOUNDS', 'vers')
        bounds_vers = None if bounds_vers == 'None' else bounds_vers
    except: _waff_co.add_section('BOUNDS')

    if bounds_vers is None or recheck:
        bounds_vers = cmd_check('bounds', 'bounds --version')
        bounds_vers = None if bounds_vers == 'None' else bounds_vers
        _waff_co.set('BOUNDS', 'vers', bounds_vers)
        
    waffles_co['BOUNDS'] = bounds_vers
    
    with open(CONFIG_FILE, 'w') as conf:
        _waff_co.write(conf)

    return(waffles_co)

## geomods-config console script

gc_version = '0.0.1'
gc_usage = '''{} ({}): geomods configuration

usage: {} [ -hv [ args ] ] ...

Options:
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            gc_version, 
            os.path.basename(sys.argv[0]))


def geomods_config():

    want_verbose = False
    
    argv = sys.argv
        
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================

    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '--help' or arg == '-h':
            print(gc_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), gc_version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True
            
        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        i = i + 1

    gc = check_config(True, want_verbose)
    utils.echo_msg(gc)

    import pkg_resources
    pkg = pkg_resources.get_distribution('geomods')
    utils.echo_msg('Geomods Console Scripts:')
    for src in pkg._get_metadata('entry_points.txt'):
        utils.echo_msg(src)
    for src in pkg._get_metadata('SOURCES.txt'):
        if 'scripts' in src:
            utils.echo_msg(src)
### End
