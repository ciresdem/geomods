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

_version = '0.1.5'

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

def init_config():
    co = ConfigParser.ConfigParser()
    cf = os.path.expanduser('~/geomods.ini')

    co['PATHINFO'] = {
        'vdatum': None
    }
    
    with open(cf, 'w') as conf:
        co.write(conf)

## =============================================================================
##
## Command execution, et cetra
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

def run_cmd_with_input(cmd, data_fun, verbose = False, prog = True):
    '''Run a command with or without a progress bar while passing data'''

    if prog:
        pb = _progress('running cmd: \033[1m{}\033[m...'.format(cmd[:44]))

    p = subprocess.Popen(cmd, shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)    
    
    c = lambda x: map(p.stdin.write, x)
    t = threading.Thread(target = data_fun, args = (p.stdin,))
    t.start()

    if prog:
        while True:
            time.sleep(3)
            pb.update()
            if not t.is_alive():
                break

    out, err = p.communicate()

    if verbose:
        _progress()._clear_stderr()
        sys.stdout.write(out)
        sys.stderr.write(err)
    
    if prog:
        pb.opm = 'ran cmd: \033[1m{}\033[m.'.format(cmd[:44])
        pb.end(p.returncode)

    return out, p.returncode

def run_cmd(cmd, verbose = False, prog = True):
    '''Run a command with or without a progress bar.'''

    if prog:
        pb = _progress('running cmd: \033[1m{}\033[m...'.format(cmd[:64]))
    
    p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)

    if prog:
        while True:
            time.sleep(3)
            pb.update()
            if p.poll() is not None:
                break

    out, err = p.communicate()

    if verbose:
        _progress()._clear_stderr()
        sys.stdout.write(out)
        sys.stderr.write(err)

    if prog:
        pb.opm = 'ran cmd: \033[1m{}\033[m.'.format(cmd[:64])
        pb.end(p.returncode)

    return out, p.returncode

def _cmd_check(cmd_str, cmd_vers_str):

    cmd_vers = None
    pb = _progress('checking for {}...'.format(cmd_str))
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str), prog = False)
    else: status = -1
    pb.opm = 'found {} version {}'.format(cmd_str, cmd_vers.split()[-1].rstrip())
    pb.end(status)

    return(cmd_vers)

def _module_check(mod_name, check_func, vers_func):

    pb = _progress('checking for {}...'.format(mod_name))
    if check_func():
        status = 0
        vers = vers_func
    else: status = -1
    pb.end(status)

    return(vers)
    
    # pb = _progress('checking system status...')
    # platform = sys.platform
    # pb.opm = 'platform is {}'.format(platform)
    # pb.end(status)

    # pb.opm = 'checking for VDatum'
    # if vdatum().vdatum_path is not None:
    #     status = 0
    #     vdatum_vers = vdatum()._version 
    # else: status = -1
    # pb.end(status)

class _progress:
    '''geomods minimal progress indicator'''

    def __init__(self, message=''):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw

        self.opm = message 
        self.opl = len(self.opm)
        self.pm = self.opm

        self._clear_stderr()

        sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw-1), self.opm))
        sys.stderr.flush()

        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one

    def _switch_way(self):
        if self.spin_way == self.add_one:
            self.spin_way = self.sub_one
        else: self.spin_way = self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def err_msg(self, msg):
        self._clear_stderr()
        sys.stderr.write('{}\n'.format(msg))

    def update(self):
        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw+1))
        self._clear_stderr()

        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], self.pm))
        sys.stderr.flush()

        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one

        self.count = self.spin_way(self.count)

    def end(self, status):
        self._clear_stderr()

        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', self.opm))
        else:
            sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', self.opm))

        sys.stderr.flush()
    
## =============================================================================
##
## VDatum class for working with NOAA's VDatum program
## Currently only compatible with VDatum > 4.0
##
## =============================================================================

class vdatum:
    '''vdatum object to communicate with NOAA's VDatum'''

    def __init__(self, vdatum_path = None, verbose = False):
        self.verbose = verbose
        self.status = 0
        if vdatum_path is None:

            co = ConfigParser.ConfigParser()
            try:
                co.read(os.path.expanduser('~/geomods.ini'))
                self.vdatum_path = co.get('VDATUM', 'jar')
            except:

                self.vdatum_path = self._find_vdatum()[0]
                if self.status == 0:
                    co.add_section('VDATUM')
                    co.set('VDATUM', 'jar', self.vdatum_path)
            
                    with open(os.path.expanduser('~/geomods.ini'), 'w') as conf:
                        co.write(conf)
                else: self.vdatum_path = None
        else: self.vdatum_path = vdatum_path

        self._version = None
        self._get_version()
        
        self.ivert = 'navd88:m:height'
        self.overt = 'mhw:m:height'
        self.ihorz = 'NAD83_2011'
        self.ohorz = 'NAD83_2011'

        self.region = '3'
        self.ft = 'txt'
        self.fd = 'space'

        self.ds_dir = 'result'

    def _get_version(self):
        if self.vdatum_path is not None:
            out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, '-'), prog = False)
            for i in out.split('\n'):
                if '- v' in i.strip():
                    self._version = i.strip().split('v')[-1]
                    break

    def _find_vdatum(self):
        '''Find the VDatum executable on the system and 
        return a list of found vdatum.jar paths'''

        results = []
        if self.verbose:
            pb = _progress('checking for vdatum.')
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.join(root, 'vdatum.jar'))
        if len(results) <= 0:
            self.status = -1
        if self.verbose:
            pb.opm = '{}..{}'.format(pb.opm, results[0])
            pb.end(self.status)

        return results

    def run_vdatum(self, src_fn):
        '''Run vdatum on src_fn which is an XYZ file'''
        
        if self.verbose:
            pb = 'transforming data with VDatum'
        else: pb = None
        
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}\
        '.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.fd, src_fn, self.ds_dir, self.region)
        out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, self.verbose)
        #out, status = run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, True)

        return status


### End
