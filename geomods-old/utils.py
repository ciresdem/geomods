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
### Commentary:
##
## General geomods utility functions
##
### Code:

import os
import sys
import subprocess

_version = '0.3.0'

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

## ==============================================
## General utility functions - utils.py
## ==============================================

def con_dec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''
    if x is None: return(x)
    return("%." + str(dec) + "f" % x)

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)'''
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_year():
    '''return the current year'''
    import datetime
    return(datetime.datetime.now().strftime('%Y'))

def remove_glob(glob_str):
    '''glob `glob_str` and os.remove results'''
    import glob
    globs = glob.glob(glob_str)
    if len(globs) > 0:
        for g in globs:
            try:
                os.remove(g)
                return(0)
            except: return(None)
    else: return(None)

cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

def run_cmd(cmd, data_fun = None, verbose = False):
    '''Run a command with or without a progress bar while passing data'''
    if verbose: echo_msg('running cmd: \033[1m{}\033[m...'.format(cmd.rstrip()))    
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    p = subprocess.Popen(cmd, shell = True, stdin = pipe_stdin, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)    

    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
    
    while p.poll() is None:
        if verbose:
            rl = p.stderr.readline()
            sys.stderr.write('\x1b[2K\r')
            sys.stderr.write(rl.decode('utf-8'))
    if verbose: sys.stderr.write(p.stderr.read().decode('utf-8'))

    out = p.stdout.read()
    p.stderr.close()
    p.stdout.close()
    if verbose: echo_msg('ran cmd: \033[1m{}\033[m and returned {}.'.format(cmd.rstrip(), p.returncode))
    return(out, p.returncode)

def cmd_check(cmd_str, cmd_vers_str):
    '''check system for availability of 'cmd_str' and return it's version'''
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str))
        return(cmd_vers.split()[-1].rstrip())
    else: return(None)

def config_check(chk_vdatum = False, verbose = False):
    '''check for needed geomods external software'''
    _waff_co = {}
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform
    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers
    if host_os == 'win32': ae = '.exe'
    else: ae = ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal_grid --version')
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version')
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version | grep Version')
    _waff_co['BOUNDS'] = cmd_check('bounds{}'.format(ae), 'bounds --version')
    
    return(_waff_co)
    
## ==============================================
## stderr messaging
## ==============================================

def echo_error_msg(msg, prefix = 'geomods'):
    '''echo error msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: error, {}\n'.format(prefix, msg))

def echo_msg(msg, prefix = 'geomods'):
    '''echo msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: {}\n'.format(prefix, msg))
    
### End
