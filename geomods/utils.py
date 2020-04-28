### utils.py
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
import glob
import time
import subprocess
import threading

import ConfigParser

_version = '0.1.6'

_license = '''geomods, version {}
Copyright (c) 2010 - 2020 CIRES Coastal DEM Team

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

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

## =============================================================================
##
## Config File
##
## The config file holds system info, including available software and versions...
##
## =============================================================================


CONFIG_FILE = os.path.expanduser('~/geomods.ini')

def init_config():
    co = ConfigParser.ConfigParser()
    cf = os.path.expanduser(CONFIG_FILE)

    co['PATHINFO'] = {
        'vdatum': None
    }

    co['VERSINFO'] = {
        'vdatum': None,
        'gdal': None,
        'gmt': None,
        'mbgrid': None,
        'lastools': None,
    }
    
    with open(cf, 'w') as conf:
        co.write(conf)

def config_set_gmt():
    co = ConfigParser.ConfigParser()
    co.read(CONFIG_FILE)
    
    try:
        gmt_vers = co.get('VERSINFO', 'gmt')
    except:
        gmt_vers = _cmd_check('gmt', 'gmt --version')
        co.set('VERSINFO', 'gmt', gmt_vers)

    with open(CONFIG_FILE, 'w') as conf:
        co.write(conf)

def config_get_vers(cmd = None):
    if cmd is not None:
        co = ConfigParser.ConfigParser()
        co.read(CONFIG_FILE)

        try:
            cmd_vers = co.get(cmd, 'vers')
        except: cmd_vers = None
    else: cmd_vers = None
    
    return(cmd_vers)    
        
def check_config(recheck = False, verbose = False):

    geomods_co = {}
    vdatum_path = None
    gmt_vers = None
    mbgrid_vers = None
    gdal_vers = None
    bounds_vers = None
    
    co = ConfigParser.ConfigParser()
    co.read(CONFIG_FILE)
    
    ## VDATUM
    import vdatum
    try:
        vdatum_path = co.get('VDATUM', 'jar')
        vdatum_path = None if vdatum_path == 'None' else vdatum_path
    except: co.add_section('VDATUM')

    if vdatum_path is None or recheck:
        vd = vdatum.vdatum(verbose=verbose)
        vdatum_path = vd.vdatum_path
        vdatum_path = None if vdatum_path == 'None' else vdatum_path
        co.set('VDATUM', 'jar', vdatum_path)
        
    geomods_co['VDATUM'] = vdatum_path
        
    ## GMT
    try:
        gmt_vers = co.get('GMT', 'vers')
        gmt_vers = None if gmt_vers == 'None' else gmt_vers
    except: co.add_section('GMT')

    if gmt_vers is None or recheck:
        gmt_vers = _cmd_check('gmt', 'gmt --version')
        gmt_vers = None if gmt_vers == 'None' else gmt_vers
        co.set('GMT', 'vers', gmt_vers)
        
    geomods_co['GMT'] = gmt_vers
        
    ## GDAL
    try:
        gdal_vers = co.get('GDAL', 'vers')
        gdal = None if gdal_vers == 'None' else gdal_vers
    except: co.add_section('GDAL')

    if gdal_vers is None or recheck:
        gdal_vers = _cmd_check('gdal-config', 'gdal-config --version')
        gdal = None if gdal_vers == 'None' else gdal_vers
        co.set('GDAL', 'vers', gdal_vers)
        
    geomods_co['GDAL'] = gdal_vers
    
    ## mbgrid
    try:
        mbgrid_vers = co.get('MBGRID', 'vers')
        mbgrid_vers = None if mbgrid_vers == 'None' else mbgrid_vers
    except: co.add_section('MBGRID')

    if mbgrid_vers is None or recheck:
        mbgrid_vers = _cmd_check('mbgrid', 'mbgrid -version | grep Version')
        mbgrid_vers = None if mbgrid_vers == 'None' else mbgrid_vers
        co.set('MBGRID', 'vers', mbgrid_vers)
        
    geomods_co['MBGRID'] = mbgrid_vers
    
    ## bounds
    try:
        bounds_vers = co.get('BOUNDS', 'vers')
        bounds_vers = None if bounds_vers == 'None' else bounds_vers
    except: co.add_section('BOUNDS')

    if bounds_vers is None or recheck:
        bounds_vers = _cmd_check('bounds', 'bounds --version')
        bounds_vers = None if bounds_vers == 'None' else bounds_vers
        co.set('BOUNDS', 'vers', bounds_vers)
        
    geomods_co['BOUNDS'] = bounds_vers
    
    with open(CONFIG_FILE, 'w') as conf:
        co.write(conf)

    return(geomods_co)

        
## =============================================================================
##
## Command execution, OS functions, et cetra
##
## OS System commands and general OS functions and CLI tools.
## run a command with a progress bar with 'run_cmd'
## check if a command exists on the system with 'cmd_exists'
## Et Cetra...
##
## =============================================================================

def remove_glob(glob_str):
    '''glob `glob_str` and os.remove results'''

    globs = glob.glob(glob_str)
    if len(globs) > 0:
        for g in globs:
            try:
                os.remove(g)
            except: pass

cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

def run_cmd(cmd, verbose = False, prog = True, data_fun = None):
    '''Run a command with or without a progress bar while passing data'''

    if prog: pb = _progress('running cmd: \033[1m{}\033[m...'.format(cmd[:84]))
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    
    p = subprocess.Popen(cmd, shell = True, stdin = pipe_stdin, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)    

    if data_fun is not None:
        if prog: pb_stdin = _progress('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
        if prog: pb_stdin.end(0, 'piped data to cmd subprocess.')
    
    while p.poll() is None:
        if verbose:
            rl = p.stderr.readline()
            sys.stderr.write('\x1b[2K\r')
            sys.stderr.write(rl)
    if verbose: sys.stderr.write(p.stderr.read())

    out = p.stdout.read()
    p.stderr.close()
    p.stdout.close()

    if prog: pb.end(p.returncode, 'ran cmd: \033[1m{}\033[m.'.format(cmd[:84]))

    return out, p.returncode

def _cmd_check(cmd_str, cmd_vers_str):

    cmd_vers = None
    pb = _progress('checking for {}...'.format(cmd_str))
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str), prog = False)
        cmd_vers = cmd_vers.split()[-1].rstrip()
    else:
        cmd_vers = None
        status = -1
    pb.end(status, 'found {} version {}'.format(cmd_str, cmd_vers))

    return(cmd_vers)

## =============================================================================
##
## Progress indicator, messaging, etc.
##
## =============================================================================
    
def _error_msg(msg):
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('geomods: error, {}\n'.format(msg))

def _msg(msg):
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('geomods: {}\n'.format(msg))

class _progress:
    '''geomods minimal progress indicator'''

    def __init__(self, message = None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.opm = message
        #self.pm = self.opm

        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one

        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
            #sys.stderr.flush()
        
    def _switch_way(self):
        if self.spin_way == self.add_one:
            self.spin_way = self.sub_one
        else: self.spin_way = self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def err_msg(self, msg):
        self._clear_stderr()
        sys.stderr.write('geomods: error, {}\n'.format(msg))

    def msg(self, msg):
        self._clear_stderr()
        sys.stderr.write('geomods: {}\n'.format(msg))

    def update(self):
        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw+1))
        self._clear_stderr()

        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], self.opm))

        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one

        self.count = self.spin_way(self.count)

    def end(self, status, end_msg = None):
        self._clear_stderr()
        if end_msg is None:
            end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))

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
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), gc_version, _license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True
            
        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        i = i + 1

    gc = check_config(True, want_verbose)
    _msg(gc)

    import pkg_resources
    pkg = pkg_resources.get_distribution('geomods')
    _msg('Geomods Console Scripts:')
    for src in pkg._get_metadata('entry_points.txt'):
        _msg(src)
    for src in pkg._get_metadata('SOURCES.txt'):
        if 'scripts' in src:
            _msg(src)
    

    
        
### End
