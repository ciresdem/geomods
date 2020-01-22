### utils.py
##
## Copyright (c) 2012 - 2020 CIRES Coastal DEM Team
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
import subprocess
import threading

import ConfigParser

_version = '0.1.0'

_license = '''geomods, version {}
Copyright (c) 2012 - 2020 CIRES Coastal DEM Team

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
## OS System commands and checks.
## run a command with a progress bar with 'run_cmd'
## check if a command exists on the system with 'cmd_exists'
##
## =============================================================================

cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

def run_cmd_with_input(cmd, data_fun, verbose = False, prog = None):
    '''Run a command with or without a progress bar while passing data'''

    want_poll = lambda a, b: a or b is not None

    if prog is not None:
        prog = _progress('{}'.format(prog))

    p = subprocess.Popen(cmd, shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)    
    
    c = lambda x: map(p.stdin.write, x)
    t = threading.Thread(target = data_fun, args = (p.stdin,))
    t.start()

    if prog is not None:
        while True:
            time.sleep(4)
            prog.update()

            if not t.is_alive():
                break

    out, err = p.communicate()

    if verbose:
        _progress()._clear_stderr()
        sys.stdout.write(out)
        sys.stderr.write(err)
    
    if prog is not None:
        prog.end(p.returncode)

    return out, p.returncode

def run_cmd(cmd, verbose = False, prog = None):
    '''Run a command with or without a progress bar.'''

    want_poll = lambda a, b: a or b is not None

    if prog is not None:
        prog = _progress('{}'.format(prog))
    
    p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)

    if prog is not None:
        while True:
            time.sleep(3)
            prog.update()
            if p.poll() is not None:
                break

    # if verbose:
    #     _progress()._clear_stderr()
    #     for line in iter(p.stderr.readline, ''):
    #         sys.stderr.write(line)
    #         #l = p.stderr.readline()
    #         #if l: sys.stderr.write(l)

    out, err = p.communicate()

    if verbose:
        _progress()._clear_stderr()
        sys.stdout.write(out)
        sys.stderr.write(err)

    # if verbose:
    #     _progress()._clear_stderr()
    #     sys.stderr.write(err)
    #     p.stderr.close()
    if prog is not None:
        prog.end(p.returncode)

    return out, p.returncode

## =============================================================================
##
## VDatum class for working with NOAA's VDatum program
## Currently only compatible with VDatum > 4.0
##
## =============================================================================

class vdatum:
    def __init__(self, vdatum_path = None, verbose = False):
        '''vdatum object to communicate with NOAA's VDatum'''

        self.verbose = verbose
        
        if vdatum_path is None:

            co = ConfigParser.ConfigParser()
            try:
                co.read(os.path.expanduser('~/geomods.ini'))
                self.vdatum_path = co.get('VDATUM', 'jar')
            except:

                self.vdatum_path = self._find_vdatum()[0]
                co.add_section('VDATUM')
                co.set('VDATUM', 'jar', self.vdatum_path)
            
                with open(os.path.expanduser('~/geomods.ini'), 'w') as conf:
                    co.write(conf)
        else: self.vdatum_path = vdatum_path
        
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
        #if len(self.vdatum_paths) > 0:
        if self.vdatum_path is not None:
            self._version = None
            out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, '-'))
            for i in out.split('\n'):
                if '- v' in i.strip():
                    #print i.strip().split('v')[-1]
                    self._version = i.strip().split('v')[-1]
                    break

    def _find_vdatum(self):
        '''Find the VDatum executable on the system and 
        return a list of found vdatum.jar paths'''

        results = []
        status = 0
        if self.verbose:
            pb = _progress('checking for vdatum.')
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.join(root, 'vdatum.jar'))
        if len(results) <= 0:
            status = -1
        if self.verbose:
            pb.opm = '{}..{}'.format(pb.opm, results[0])
            pb.end(status)

        return results

    def run_vdatum(self, src_fn):
        '''Run vdatum on src_fn which is an XYZ file'''
        
        if self.verbose:
            pb = 'transforming data with VDatum'
        else: pb = None
        
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}'.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.fd, src_fn, self.ds_dir, self.region)
        #out, status = run_cmd('java -jar {} {}'.format(self.vdatum_paths[0], vdc), False, 'transforming data with vdatum')
        out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, vdc), False, pb)

        return status

## =============================================================================
##
## Progress Bar
##
## with 'prog_message' print a simple progress bar and message.
## use the 'prog_bar.pm' variable to update the message while running.
##
## =============================================================================

class _progress:
    def __init__(self, message=''):
        self.tw = 5
        self.count = 0
        self.pc = self.count % self.tw

        self.opm = message 
        self.opl = len(self.opm)
        self.pm = self.opm

        self._clear_stderr()

        sys.stderr.write('\r[{}] {:40}\r'.format(" " * (self.tw-1), self.opm))
        sys.stderr.flush()

        self.spinner = ['*   ', '**  ', '*** ', ' ***', '  **', '   *']
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one

    def _switch_way(self):
        if self.spin_way == self.add_one:
            self.spin_way = self.sub_one
        else: self.spin_way = self.add_one

    def _terminal_size(self):
        return os.popen('stty size', 'r').read().split()

    def _clear_stderr(self, slen = 79):
        #sys.stderr.write('\r{}\r'.format(' ' * int(self._terminal_size()[1])))
        sys.stderr.write('\x1b[2K')
        sys.stderr.flush()

    def update(self):
        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw+1))
        self._clear_stderr()

        sys.stderr.write('\r[\033[32m{:4}\033[m] {:40}\r'.format(self.spinner[self.sc], self.pm))
        sys.stderr.flush()

        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one

        self.count = self.spin_way(self.count)

    def end(self, status):
        self._clear_stderr()

        if status != 0:
            sys.stderr.write('\r[\033[31mfail\033[m] {:40}\n'.format(self.opm))
        else:
            sys.stderr.write('\r[\033[32m{:^4}\033[m] {:40}\n'.format('ok', self.opm))

        sys.stderr.flush()

### End
