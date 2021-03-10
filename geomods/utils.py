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
### Code:

import os
import sys
import time
import subprocess
import glob
import math
import zipfile
import gzip
import datetime
import numpy as np

## ==============================================
## General utility functions - utils.py
## ==============================================
_version = '0.4.0'

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    returns a str representation of float(inc)'''
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_date():
    '''return the current date'''
    return(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))

def this_year():
    '''return the current year'''
    return(datetime.datetime.now().strftime('%Y'))

def rm_f(f_str):
    '''os.remove f_str, pass if error'''
    try:
        if os.path.exists(f_str):
            os.remove(f_str)
    except: pass
    return(0)
        
def remove_glob(*args):
    '''glob `glob_str` and os.remove results, pass if error'''
    for glob_str in args:
        try:
            globs = glob.glob(glob_str)
            for g in globs:
                if os.path.isdir(g):
                    remove_globs('{}/*'.format(g))
                    os.removedirs(g)
                else: os.remove(g)
        except: pass
    return(0)

def base_name(instr, extension):
    '''Return the basename of a file-string'''
    return(instr[:-len(extension)])

def args2dict(args, dict_args = {}):
    '''convert list of arg strings to dict.
    args are a list of ['key=val'] pairs

    returns a dictionary of the key/values'''
    for arg in args:
        p_arg = arg.split('=')
        dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else True if p_arg[1].lower() == 'true' else None if p_arg[1].lower() == 'none' else p_arg[1]
    return(dict_args)

def int_or(val, or_val = None):
    '''returns val as int otherwise returns or_val'''
    try:
        return(int(val))
    except: return(or_val)

def euc_dst(pnt0, pnt1):
    '''return the distance between pnt0 and pnt1,
    using the euclidean formula.
    `pnts` are geographic and result is in meters.'''
    rad_m = 637100
    distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(pnt0, pnt1)]))
    return(rad_m * distance)
    
def hav_dst(pnt0, pnt1):
    '''return the distance between pnt0 and pnt1,
    using the haversine formula.
    `pnts` are geographic and result is in meters.'''
    x0 = float(pnt0[0])
    y0 = float(pnt0[1])
    x1 = float(pnt1[0])
    y1 = float(pnt1[1])
    rad_m = 637100
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    a = math.sin(dx / 2) * math.sin(dx / 2) + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy / 2) * math.sin(dy / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return(rad_m * c)

def _clean_zips(zip_files):
    '''remove all files\directories in `zip_files`'''
    for i in zip_files:
        if os.path.isfile(i):
            os.remove(i)
            zip_files = [x for x in zip_files if x != i]
    if len(zip_files) > 0:
        for i in zip_files:
            if os.path.isdir(i):
                try:
                    os.removedirs(i)
                except: pass
    return(0)

def unzip(zip_file):
    '''unzip (extract) `zip_file`

    return a list of extracted file names.'''
    zip_ref = zipfile.ZipFile(zip_file)
    zip_files = zip_ref.namelist()
    zip_ref.extractall()
    zip_ref.close()
    return(zip_files)

def gunzip(gz_file):
    '''gunzip `gz_file`

    return the extracted file name.'''
    if os.path.exists(gz_file):
        gz_split = gz_file.split('.')[:-1]
        guz_file = '{}.{}'.format(gz_split[0], gz_split[1])
        with gzip.open(gz_file, 'rb') as in_gz, \
             open(guz_file, 'wb') as f:
            while True:
                block = in_gz.read(65536)
                if not block:
                    break
                else: f.write(block)
    else:
        echo_error_msg('{} does not exist'.format(gz_file))
        guz_file = None
    return(guz_file)

def p_unzip(src_file, exts = None):
    '''unzip/gunzip src_file based on `exts`'''
    src_procs = []
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for ext in exts:
                for zf in zfs:
                    if ext == zf.split('.')[-1]:
                        #if ext in zf:
                        src_procs.append(os.path.basename(zf))
                        with open(os.path.basename(zf), 'wb') as f:
                            f.write(z.read(zf))
    elif src_file.split('.')[-1] == 'gz':
        tmp_proc = gunzip(src_file)
        if tmp_proc is not None:
            for ext in exts:
                if ext == tmp_proc.split('.')[-1]:
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                else: remove_glob(tmp_proc)
    else:
        for ext in exts:
            if ext == src_file.split('.')[-1]:
                src_procs.append(src_file)
                break
        
    return(src_procs)
    
def procs_unzip(src_file, exts):
    '''unzip/gunzip src_file based on `exts`
    
    return the file associated with `exts`'''
    zips = []
    src_proc = None
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for ext in exts:
                for zf in zfs:
                    if ext in zf:
                        #src_proc = os.path.join(os.path.dirname(src_file), zf)
                        src_proc = os.path.basename(zf)
                        with open(src_proc, 'wb') as f:
                            f.write(z.read(zf))
                        break
    elif src_file.split('.')[-1] == 'gz':
        tmp_proc = gunzip(src_file)
        if tmp_proc is not None:
            for ext in exts:
                if ext in tmp_proc:
                    src_proc = os.path.basename(tmp_proc)
                    os.rename(tmp_proc, src_proc)
                    break
    else:
        for ext in exts:
            if ext in src_file:
                src_proc = src_file
                break
    return([src_proc, zips])

def err_fit_plot(xdata, ydata, out, fitfunc, dst_name = 'unc', xa = 'distance'):
    '''plot a best fit plot'''
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        plt.plot(xdata, ydata, 'o')
        plt.plot(xdata, fitfunc(out, xdata), '-')
        plt.xlabel(xa)
        plt.ylabel('error (m)')
        out_png = '{}_bf.png'.format(dst_name)
        plt.savefig(out_png)
        plt.close()
        
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

def err_scatter_plot(error_arr, dist_arr, dst_name = 'unc', xa = 'distance'):
    '''plot a scatter plot'''
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        plt.scatter(dist_arr, error_arr)
        #plt.title('Scatter')
        plt.xlabel(xa)
        plt.ylabel('error (m)')
        out_png = '{}_scatter.png'.format(dst_name)
        plt.savefig(out_png)
        plt.close()
        
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

def err2coeff(err_arr, coeff_guess = [0, 0.1, 0.2], dst_name = 'unc', xa = 'distance'):
    '''calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist`'''
    from scipy import optimize
    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    max_int_dist = np.max(distance)
    nbins = 10
    n, _ = np.histogram(distance, bins = nbins)
    while 0 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins = nbins)
    serror, _ = np.histogram(distance, bins = nbins, weights = error)
    serror2, _ = np.histogram(distance, bins = nbins, weights = error**2)
    mean = serror / n
    std = np.sqrt(serror2 / n - mean * mean)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)
    #print(xdata)
    xdata[xdata - 0 < 0.0001] = 0.0001
    #print(xdata)
    fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    #f = lambda x, a, b, c: a/1+b*x**c
    #g = lambda x, a, b, c: b/a*x**c+1/a
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)
    #if out[2]<0.001: out[2]=0.001
    #try:
    #popt, pcov = optimize.curve_fit(f, xdata, ydata, p0=optimize.curve_fit(g, xdata, 1/ydata)[0])
    #except: popt = coeff_guess
    #print(popt, np.sqrt(np.diag(pcov)))
    
    try:
        err_fit_plot(xdata, ydata, out, fitfunc, dst_name, xa)
        err_scatter_plot(error, distance, dst_name, xa)
    except: echo_error_msg('unable to generate error plots, please check configs.')
    return(out)

## ==============================================
## system cmd verification and configs.
##
## run_cmd - run a system command as a subprocess
## and return the return code and output
##
## yield_cmd - run a system command as a subprocess
## and yield the output
## ==============================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

def run_cmd(cmd, data_fun = None, verbose = False):
    '''Run a system command while optionally passing data.
    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    returns [command-output, command-return-code]'''
    if verbose:
        #echo_msg('running cmd: {}...'.format(cmd.rstrip()))
        _prog = _progress('running cmd: {}...'.format(cmd.rstrip()[:20]))
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    if verbose:
        p = subprocess.Popen(cmd, shell = True, stdin = pipe_stdin, stdout = subprocess.PIPE, close_fds = True)
    else: p = subprocess.Popen(cmd, shell = True, stdin = pipe_stdin, stdout = subprocess.PIPE, stderr = subprocess.PIPE, close_fds = True)

    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
    
    while p.poll() is None:
        if verbose:
            #rl = p.stderr.readline()
            time.sleep(2)
            #sys.stderr.write('\x1b[2K\r')
            #sys.stderr.write(rl.decode('utf-8'))
            _prog.update()
    #if verbose: sys.stderr.write(p.stderr.read().decode('utf-8'))

    out = p.stdout.read()
    if not verbose: p.stderr.close()
    p.stdout.close()
    if verbose: _prog.end(p.returncode, 'ran cmd: {}... and returned {}.'.format(cmd.rstrip()[:20], p.returncode))
    return(out, p.returncode)

def yield_cmd(cmd, data_fun = None, verbose = False):
    '''Run a system command while optionally passing data.
    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    returns [command-output, command-return-code]'''
    
    if verbose: echo_msg('running cmd: {}...'.format(cmd.rstrip()))    
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    p = subprocess.Popen(cmd, shell = True, stdin = pipe_stdin, stdout = subprocess.PIPE, close_fds = True)
    
    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()

    while p.poll() is None:
        #while True:
        line = p.stdout.readline().decode('utf-8')
        if not line: break
        else: yield(line)
    line = p.stdout.read().decode('utf-8')
    #yield(line)
    p.stdout.close()
    if verbose: echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))

def cmd_check(cmd_str, cmd_vers_str):
    '''check system for availability of 'cmd_str' 

    returns the commands version or None'''
    
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str))
        return(cmd_vers.rstrip())
    else: return(None)

def config_check(chk_vdatum = False, verbose = False):
    '''check for needed waffles external software.
    waffles external software: gdal, gmt, mbsystem
    also checks python version and host OS and 
    records waffles version

    returns a dictionary of gathered results.'''
    
    _waff_co = {}
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform
    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers[0]
    ae = '.exe' if host_os == 'win32' else ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal_grid --version').decode()
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version').decode()
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version 2>&1 | grep Version').decode()
    _waff_co['LASTOOLS'] = cmd_check('las2txt{}'.format(ae), 'las2txt -version 2>&1 | awk \'{print $5}\'').decode()
    _waff_co['GEOMODS'] = str(_version)
    return(_waff_co)
    
## ==============================================
## stderr messaging
## ==============================================
def echo_warning_msg2(msg, prefix = 'waffles'):
    '''echo warning msg to stderr using `prefix`
    >> echo_warning_msg2('message', 'test')
    test: warning, message'''
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1mwarining\033[m, {}\n'.format(prefix, msg))

def echo_error_msg2(msg, prefix = 'waffles'):
    '''echo error msg to stderr using `prefix`
    >> echo_error_msg2('message', 'test')
    test: error, message'''
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles', nl = True):
    '''echo `msg` to stderr using `prefix`
    >> echo_msg2('message', 'test')
    test: message'''
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
    sys.stderr.flush()

## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg_inline = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), nl = False)

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix = os.path.basename(sys.argv[0]))

class _progress:
    '''geomods minimal progress indicator'''

    def __init__(self, message = None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.opm = message
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
        
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def update_perc(self, p, msg = None):
        if len(p) == 2:
            self._clear_stderr()
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {:40}\r'.format(self.perc(p), msg if msg is not None else self.opm))
        else: self.update()
        
    def update(self, msg = None):

        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], msg if msg is not None else self.opm))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)

    def end(self, status, end_msg = None):
        self._clear_stderr()
        if end_msg is None: end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))

### End
