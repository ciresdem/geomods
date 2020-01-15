### utils.py
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
import subprocess

## =============================================================================
##
## Command execution, et cetra
#
## OS System commands and checks.
## run a command with a progress bar with 'run_cmd'
## check if a command exists on the system with 'cmd_exists'
##
## =============================================================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

def run_cmd(cmd, verbose = False, prog = None):
    '''Run a command with or without a progress bar.'''

    want_poll = lambda a, b: a or b is not None

    if prog is not None:
        prog = _progress('geomods: `{}`'.format(prog))
    
    p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)
    
    if want_poll(verbose, prog):
        while p.poll() is None:
            if verbose:
                l = p.stderr.readline()
                _progress()._clear_stderr()
                if l: sys.stderr.write('{}\n'.format(l.strip()))
            elif prog is not None:
                prog.update()
                time.sleep(3)

        if verbose:
            l = p.stderr.read()
            _progress()._clear_stderr()
            if l: sys.stderr.write('{}\n'.format(l.strip()))

    out, err = p.communicate()

    if verbose:
        p.stderr.close()

    if prog is not None:
        prog.end(p.returncode)

    return out, p.returncode

class vdatum:
    def __init__(self, vdatum_path = None):
        if vdatum_path is None:
            self.vdatum_paths = self._find_vdatum()
        else: self.vdatum_paths = [vdatum_path]

        self.ivert = 'navd88:m:height'
        self.overt = 'mhw:m:height'
        self.ihorz = 'NAD83_2011'
        self.ohorz = 'NAD83_2011'

        self.region = '3'
        self.ft = 'txt'
        self.fd = 'space'

        self.ds_dir = 'result'

    def _find_vdatum(self):
        results = []
        status = 0
        pb = _progress('geomods: attempting to locate a vdatum')
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.join(root, 'vdatum.jar'))
        if len(results) <= 0:
            status = -1
        pb.opm = '{}: {}'.format(pb.opm, results[0])
        pb.end(status)
        return results

    def run_vdatum(self, src_fn):
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}'.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.fd, src_fn, self.ds_dir, self.region)
        out, status = run_cmd('java -jar {} {}'.format(self.vdatum_paths[0], vdc), False, 'transforming data with vdatum')
        return status

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

    def _clear_stderr(self):
        sys.stderr.write('\r{}\r'.format(' ' * int(self._terminal_size()[1])))
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
