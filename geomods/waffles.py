### waffles.py
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
## WAFFLES - Generate Digital Elevation Models and derivatives using a variety of algorithms, etc.
##
## GDAL and gdal-python are required to run waffles.
##
## Recommended external software for full functionality:
## - GMT (FOSS)
## - MBSystem (FOSS)
## - VDatum (US/NOAA)
## - LASTools (Non-Free) - data processing
##
## see/set `_waffles_grid_info` dictionary to run a grid.
##
## Current DEM modules:
## surface (GMT), triangulate (GMT/GDAL), nearneighbor (GMT/GDAL), mbgrid (MBSYSTEM), num (waffles), average (GDAL)
##
## optionally, clip, filter, buffer the resulting DEM.
##
## find data to grid with GEOMODS' fetch.py
##
## DATALIST and REGION functions - datalists.py
##
## MBSystem/Waffles style datalists.
## Recurse through a datalist file and process the results.
##
## a datalist '*.datalist' file should be formatted as in MBSystem:
## ~path ~format ~weight ~metadata,list ~etc
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
## a format of 200 represents GDAL data
## a format of 300 represents LAS/LAZ data <not implemented>
## a format of 400 represents a FETCHES module - e.g. `nos:datatype=bag`
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
## or via `waffles -M datalists:infos=True`
##
## GDAL/LIDAR/FETCHES data don't need inf files.
##
## if 'region' is specified, will only process data that falls within
## the given region
##
### TODO:
## Add remove/replace module
## Add source uncertainty to uncertainty module
## Add LAS/LAZ support to datalits
## -W weight-range for datalist processing
## -B for 'breakline' (densify line/exract nodes/add to datalist)
##
### Code:
import sys
import os
import io
import time
import glob
import math
import copy
import shutil
import subprocess

try:
    import Queue as queue
except: import queue as queue
import threading
## ==============================================
## import gdal, etc.
## ==============================================
import numpy as np
import json
import gdal
import ogr
import osr

#from geomods import fetches
import fetches

## ==============================================
## General utility functions - utils.py
## ==============================================
_version = '0.6.0'

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    returns a str representation of float(inc)'''
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_date():
    '''return the current date'''
    
    import datetime
    return(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))

def this_year():
    '''return the current year'''
    
    import datetime
    return(datetime.datetime.now().strftime('%Y'))

def rm_f(f_str):
    '''os.remove f_str, pass if error'''
    
    try:
        if os.path.exists(f_str):
            os.remove(f_str)
    except: pass
    return(0)
        
def remove_glob(glob_str):
    '''glob `glob_str` and os.remove results, pass if error'''

    try:
        globs = glob.glob(glob_str)
    except: globs = None
    if globs is None: return(0)
    for g in globs:
        try:
            os.remove(g)
        except: pass
    return(0)

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

def path_exists_or_url(src_str):
    if os.path.exists(src_str): return(True)
    if src_str[:4] == 'http': return(True)
    if src_str.split(':')[0] in _known_datalist_fmts[400]: return(True)
    echo_error_msg('invalid datafile/datalist: {}'.format(src_str))
    return(False)

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
    
    import zipfile
    zip_ref = zipfile.ZipFile(zip_file)
    zip_files = zip_ref.namelist()
    zip_ref.extractall()
    zip_ref.close()
    return(zip_files)

def gunzip(gz_file):
    '''gunzip `gz_file`

    return the extracted file name.'''
    
    import gzip
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

def procs_unzip(src_file, exts):
    '''unzip/gunzip src_file based on `exts`
    
    return the file associated with `exts`'''

    zips = []
    src_proc = None
    if src_file.split('.')[-1] == 'zip':
        zips = unzip(src_file)
        for ext in exts:
            for zf in zips:
                if ext in zf:
                    src_proc = zf
                    break
                #else: remove_glob(zf)
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
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')
    
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel(xa)
    plt.ylabel('error (m)')
    out_png = '{}_bf.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

def err_scatter_plot(error_arr, dist_arr, dst_name = 'unc', xa = 'distance'):
    '''plot a scatter plot'''
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

    plt.scatter(dist_arr, error_arr)
    #plt.title('Scatter')
    plt.xlabel(xa)
    plt.ylabel('error (m)')
    out_png = '{}_scatter.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

def err2coeff(err_arr, coeff_guess = [0, 0.1, 0.2], dst_name = 'unc', xa = 'distance'):
    '''calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist`'''

    from scipy import optimize

    #np.seterr(divide='ignore', invalid='ignore')
    
    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    max_int_dist = np.max(distance)
    nbins = 10
    n, _ = np.histogram(distance, bins = nbins)
    #print(n)
    # want at least 2 values in each bin?
    while 0 in n:
        nbins -= 1
        #print(nbins)
        n, _ = np.histogram(distance, bins = nbins)
    #print(0 in n)
    serror, _ = np.histogram(distance, bins = nbins, weights = error)
    serror2, _ = np.histogram(distance, bins = nbins, weights = error**2)
    mean = serror / n
    std = np.sqrt(serror2 / n - mean * mean)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)
    fitfunc = lambda p, x: p[0] + p[1] * (abs(x) ** abs(p[2]))
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)
    err_fit_plot(xdata, ydata, out, fitfunc, dst_name, xa)
    err_scatter_plot(error, distance, dst_name, xa)
    return(out)

## ==============================================
## system cmd verification and configs.
## ==============================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

def run_cmd(cmd, data_fun = None, verbose = False):
    '''Run a system command while optionally passing data.
    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    returns [command-output, command-return-code]'''
    
    if verbose: echo_msg('running cmd: {}...'.format(cmd.rstrip()))    
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
    if verbose: echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))
    return(out, p.returncode)

def yield_cmd(cmd, data_fun = None, verbose = False):
    '''Run a system command while optionally passing data.
    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    returns [command-output, command-return-code]'''
    
    #if verbose: echo_msg('running cmd: {}...'.format(cmd.rstrip()))    
    #if data_fun is not None:
    #    pipe_stdin = subprocess.PIPE
    #else: pipe_stdin = None
    p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, close_fds = True)

    #if data_fun is not None:
    #    if verbose: echo_msg('piping data to cmd subprocess...')
    #    data_fun(p.stdin)
    #    p.stdin.close()

    while True:
        line = p.stdout.readline().decode('utf-8')
        if not line: break
        else: yield(line)
    p.stdout.close()
    #if verbose: echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))

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
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal_grid --version')
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version')
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version | grep Version')
    _waff_co['WAFFLES'] = str(_version)
    return(_waff_co)
    
## ==============================================
## stderr messaging
## ==============================================
def echo_error_msg2(msg, prefix = 'waffles'):
    '''echo error msg to stderr using `prefix`
    >> echo_error_msg2('message', 'test')
    test: error, message'''
    
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: error, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles', nl = True):
    '''echo `msg` to stderr using `prefix`
    >> echo_msg2('message', 'test')
    test: message'''
    
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))

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
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
        
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

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
        if end_msg is None: end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))

## ==============================================
## regions - regions are a bounding box list:
## [w, e, s, n]
## -- regions.py --
## ==============================================
def region_valid_p(region):
    '''return True if `region` [xmin, xmax, ymin, ymax] appears to be valid'''
    
    if region is not None:
        if region[0] <= region[1] and region[2] <= region[3]: return(True)
        else: return(False)
    else: return(False)

def region_center(region):
    '''find the center point [xc, yc] of the `region` [xmin, xmax, ymin, ymax]

    returns the center point [xc, yc]'''
    
    xc = region[0] + (region[1] - region[0] / 2)
    yc = region[2] + (region[3] - region[2] / 2)
    return([xc, yc])

def region_pct(region, pctv):
    '''calculate a percentage buffer for the `region` [xmin, xmax, ymin, ymax]

    returns the pctv buffer val of the region'''
    
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv = 0, pct = False):
    '''return the region buffered by buffer-value `bv`
    if `pct` is True, attain the buffer-value via: region_pct(region, bv)

    returns the buffered region [xmin, xmax, ymin, ymax]'''
    
    if pct: bv = region_pct(region, bv)
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def regions_reduce(region_a, region_b):
    '''combine two regions and find their minimum combined region.
    if the regions don't overlap, will return an invalid region.
    check the result with region_valid_p()
    
    return the minimum region [xmin, xmax, ymin, ymax] when combining `region_a` and `region_b`'''
    
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] > region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] < region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] > region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] < region_b[3] else region_b[3]
    return(region_c)

def regions_merge(region_a, region_b):
    '''combine two regions and find their maximum combined region.

    returns maximum region [xmin, xmax, ymin, ymax] when combining `region_a` `and region_b`'''
    
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] < region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] > region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] < region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] > region_b[3] else region_b[3]
    if len(region_a) > 4 and len(region_b) > 4:
        region_c.append(region_a[4] if region_a[4] < region_b[4] else region_b[4])
        region_c.append(region_a[5] if region_a[5] > region_b[5] else region_b[5])
    return(region_c)

def regions_intersect_p(region_a, region_b):
    '''check if two regions intersect.
    region_valid_p(regions_reduce(region_a, region_b))

    return True if `region_a` and `region_b` intersect else False'''
    
    if region_a is not None and region_b is not None:
        return(region_valid_p(regions_reduce(region_a, region_b)))
    else: return(False)
    
def regions_intersect_ogr_p(region_a, region_b):
    '''check if two regions intersect.
    region_a_ogr_geom.Intersects(region_b_ogr_geom)
    
    return True if `region_a` and `region_b` intersect else False.'''
    
    if region_a is not None and region_b is not None:
        geom_a = gdal_region2geom(region_a)
        geom_b = gdal_region2geom(region_b)
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(True)

def region_format(region, t = 'gmt'):
    '''format region to string, defined by `t`
    t = 'str': xmin/xmax/ymin/ymax
    t = 'gmt': -Rxmin/xmax/ymin/ymax
    t = 'bbox': xmin,ymin,xmax,ymax
    t = 'te': xmin ymin xmax ymax
    t = 'ul_lr': xmin ymax xmax ymin
    t = 'fn': ymax_xmin

    returns the formatted region as str'''

    if t == 'str': return('/'.join([str(x) for x in region[:4]]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region[:4]]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'te': return(' '.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'ul_lr': return(' '.join([str(region[0]), str(region[3]), str(region[1]), str(region[2])]))
    elif t == 'fn':
        ns = 's' if region[3] < 0 else 'n'
        ew = 'e' if region[0] > 0 else 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100)) % 100, 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100)) % 100))
    elif t == 'inf': return(' '.join([str(x) for x in region]))

def region_chunk(region, inc, n_chunk = 10):
    '''chunk the region [xmin, xmax, ymin, ymax] into 
    n_chunk by n_chunk cell regions, given inc.

    returns a list of chunked regions.'''
    
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk
    o_chunks = []
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    
    while True:
        y_chunk = n_chunk
        while True:
            this_x_origin = x_chunk - n_chunk
            this_y_origin = y_chunk - n_chunk
            this_x_size = x_chunk - this_x_origin
            this_y_size = y_chunk - this_y_origin
            
            geo_x_o = region[0] + this_x_origin * inc
            geo_x_t = geo_x_o + this_x_size * inc
            geo_y_o = region[2] + this_y_origin * inc
            geo_y_t = geo_y_o + this_y_size * inc

            if geo_y_t > region[3]: geo_y_t = region[3]
            if geo_y_o < region[2]: geo_y_o = region[2]
            if geo_x_t > region[1]: geo_x_t = region[1]
            if geo_x_o < region[0]: geo_x_o = region[0]
            o_chunks.append([geo_x_o, geo_x_t, geo_y_o, geo_y_t])
        
            if y_chunk < ycount:
                y_chunk += n_chunk
                i_chunk += 1
            else: break
        if x_chunk < xcount:
            x_chunk += n_chunk
            x_i_chunk += 1
        else: break
    return(o_chunks)

def regions_sort(trainers):
    '''sort regions by distance; regions is a list of regions [xmin, xmax, ymin, ymax].

    returns the sorted region-list'''
    
    train_sorted = []
    for z, train in enumerate(trainers):
        train_d = []
        np.random.shuffle(train)
        while True:
            if len(train) == 0: break
            this_center = region_center(train[0][0])
            train_d.append(train[0])
            train = train[1:]
            if len(train) == 0: break
            dsts = [hav_dst(this_center, region_center(x[0])) for x in train]
            min_dst = np.percentile(dsts, 50)
            d_t = lambda t: hav_dst(this_center, region_center(t[0])) > min_dst
            np.random.shuffle(train)
            train.sort(reverse=True, key=d_t)
        #echo_msg(' '.join([region_format(x[0], 'gmt') for x in train_d[:25]]))
        train_sorted.append(train_d)
    return(train_sorted)

def region_warp(region, s_warp = 4326, t_warp = 4326):
    '''warp region from source `s_warp` to target `t_warp`, using EPSG keys

    returns the warped region'''
    
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(int(s_warp))

    if t_warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(t_warp))
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)        
        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[0], region[2]))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[1], region[3]))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(region)

def region2ogr(region, dst_ogr, append = False):
    '''convert a region string to an OGR vector'''

    dst_wkt = gdal_region2wkt(region)
    
    driver = ogr.GetDriverByName('ESRI Shapefile')

    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)
        
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPolygon)

    dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
    dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(dst_wkt))
    dst_feat.SetField('id', 1)
    dst_lyr.CreateFeature(dst_feat)
    dst_feat = None

    dst_ds = None

def z_region_valid_p(z_region):
    '''return True if z_region appears to be valid'''
    
    if len(z_region) < 2: return(False)
    if z_region[0] > z_region[1]: return(False)
    return(True)

def z_region_pass(region, upper_limit = None, lower_limit = None):
    '''return True if extended region [xmin, xmax, ymin, ymax, zmin, zmax] is 
    within upper and lower z limits'''
    
    if region is not None:
        z_region = region[4:]
        if z_region is not None and len(z_region) >= 2:
            if upper_limit is not None:
                if z_region[0] >= upper_limit:
                    return(False)
            if lower_limit is not None:
                if z_region[1] <= lower_limit:
                    return(False)
    return(True)

def z_pass(z, upper_limit = None, lower_limit = None):
    '''return True if z-value is between upper and lower z limits'''

    if z is None: return(True)
    if upper_limit is not None:
        if z >= upper_limit:
            return(False)
    if lower_limit is not None:
        if z <= lower_limit:
            return(False)
    return(True)

## =============================================================================
##
## VDatum - vdatumfun.py
## wrapper functions for NOAA's VDatum
##
## Currently only compatible with VDatum >= 4.0
##
## TODO: add all vdatum cli options
## =============================================================================
_vd_config = {
    'jar': None,
    'ivert': 'navd88:m:height',
    'overt': 'mhw:m:height',
    'ihorz': 'NAD83_2011',
    'ohorz': 'NAD83_2011',
    'region': '3',
    'fmt': 'txt',
    'xyzl': '0,1,2',
    'skip': '0',
    'delim': 'space',
    'result_dir': 'result',
    'verbose': False,
}

def vdatum_locate_jar():
    '''Find the VDatum executable on the local system.

    returns a list of found vdatum.jar system paths'''
    
    results = []
    for root, dirs, files in os.walk('/'):
        if 'vdatum.jar' in files:
            results.append(os.path.abspath(os.path.join(root, 'vdatum.jar')))
            break
    if len(results) == 0:
        return(None)
    else: return(results)

def vdatum_get_version(vd_config = _vd_config):
    '''run vdatum and attempt to get it's version
    
    return the vdatum version or None'''
    
    if vd_config['jar'] is None:
        vd_config['jar'] = vdatum_locate_jar()
    if vd_config['jar'] is not None:
        out, status = run_cmd('java -jar {} {}'.format(vd_config['jar'], '-'), verbose = self.verbose)
        for i in out.decode('utf-8').split('\n'):
            if '- v' in i.strip():
                return(i.strip().split('v')[-1])
    return(None)

def vdatum_xyz(xyz, vd_config = _vd_config):
    '''run vdatum on an xyz list [x, y, z]

    returns the transformed xyz list'''
    
    if vd_config['jar'] is None: vd_config['jar'] = vdatum_locate_jar()[0]
    if vd_config['jar'] is not None:
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -pt:{},{},{} region:{}\
        '.format(vd_config['ihorz'], vd_config['ivert'], vd_config['ohorz'], vd_config['overt'], \
                 xyz[0], xyz[1], xyz[2], vd_config['region'])
        out, status = run_cmd('java -Djava.awt.headless=false -jar {} {}'.format(vd_config['jar'], vdc), verbose = False)
        for i in out.split('\n'):
            if 'Height/Z' in i:
                z = float(i.split()[2])
                break
        return([xyz[0], xyz[1], z])
    else: return(xyz)

def vdatum_clean_result(result_f = 'result'):
    '''clean the vdatum 'result' folder'''
    
    remove_glob('{}/*'.format(result_f))
    try:
        os.removedirs(result_f)
    except: pass
    
def run_vdatum(src_fn, vd_config = _vd_config):
    '''run vdatum on src_fn which is an XYZ file
    use vd_config to set vdatum parameters.

    returns [command-output, command-return-code]'''
    
    if vd_config['jar'] is None: vd_config['jar'] = vdatum_locate_jar()[0]
    if vd_config['jar'] is not None:
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},{},skip{}:{}:{} region:{}\
        '.format(vd_config['ihorz'], vd_config['ivert'], vd_config['ohorz'], vd_config['overt'], \
                 vd_config['delim'], vd_config['xyzl'], vd_config['skip'], src_fn, vd_config['result_dir'], vd_config['region'])
        #return(run_cmd('java -jar {} {}'.format(vd_config['jar'], vdc), verbose = True))
        return(run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(vd_config['jar'], vdc), verbose = vd_config['verbose']))
    else: return([], -1)
    
## ==============================================
## GMT Wrapper Functions - gmtfun.py
## wrapper functions to GMT system commands
##
## GMT must be installed on the system to run these
## functions and commands.
## ==============================================
def gmt_inf(src_xyz):
    '''generate an info (.inf) file from a src_xyz file using GMT.

    returns [cmd-output, cmd-return-code]'''
    
    return(run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), verbose = False))

def gmt_grd_inf(src_grd):
    '''generate an info (.inf) file from a src_gdal file using GMT.

    returns [cmd-output, cmd-return-code]'''
    
    return(run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_grd, src_grd), verbose = False))

def gmt_inc2inc(inc_str):
    '''convert a GMT-style `inc_str` (6s) to geographic units
    c/s - arc-seconds
    m - arc-minutes

    return float increment value.'''
    
    if inc_str is None or inc_str.lower() == 'none': return(None)
    units = inc_str[-1]
    if units == 'c': inc = float(inc_str[:-1]) / 3600
    elif units == 's': inc = float(inc_str[:-1]) / 3600
    elif units == 'm': inc = float(inc_str[:-1]) / 360
    else:
        try:
            inc = float(inc_str)
        except ValueError as e:
            echo_error_msg('could not parse increment {}, {}'.format(inc_str, e))
            return(None)
    return(inc)

def gmt_grd2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''convert the grd file to tif using GMT

    returns the gdal file name or None'''
    
    dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
    grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
    '.format(src_grd, dst_gdal, dst_fmt))
    out, status = run_cmd(grd2gdal_cmd, verbose = verbose)
    if status == 0:
        return(dst_gdal)
    else: return(None)

def gmt_grdinfo(src_grd, verbose = False):
    '''gather infos about src_grd using GMT grdinfo.

    return an info list of `src_grd`'''
    
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
    out, status = run_cmd(grdinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)

def gmt_gmtinfo(src_xyz, verbose = False):
    '''gather infos about src_xyz using GMT gmtinfo

    return an info list of `src_xyz`'''
    
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    gmtinfo_cmd = ('gmt gmtinfo {} -C'.format(src_xyz))
    out, status = run_cmd(gmtinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)
        
def gmt_select_split(o_xyz, sub_region, sub_bn, verbose = False):
    '''split an xyz file into an inner and outer region.

    returns [inner_region, outer_region]'''
    
    out_inner = None
    out_outer = None
    gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, region_format(sub_region, 'gmt'), sub_bn)
    out, status = run_cmd(gmt_s_inner, verbose = verbose)
    if status == 0: out_inner = '{}_inner.xyz'.format(sub_bn)
    gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, region_format(sub_region, 'gmt'), sub_bn)
    out, status = run_cmd(gmt_s_outer, verbose = verbose)
    if status == 0:  out_outer = '{}_outer.xyz'.format(sub_bn)
    return([out_inner, out_outer])
        
def gmt_grdcut(src_grd, src_region, dst_grd, verbose = False):
    '''cut `src_grd` to `src_region` using GMT grdcut
    
    returns [cmd-output, cmd-return-code]'''
    
    cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
    return(run_cmd(cut_cmd1, verbose = True))

def gmt_grdfilter(src_grd, dst_grd, dist = '3s', verbose = False):
    '''filter `src_grd` using GMT grdfilter

    returns [cmd-output, cmd-return-code]'''
    
    ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1'.format(src_grd, dst_grd, src_grd, dist))
    return(run_cmd(ft_cmd1, verbose = verbose))

def gmt_nan2zero(src_grd, node = 'pixel', verbose = False):
    '''convert nan values in `src_grd` to zero

    returns status code (0 == success) '''
    
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = tmp.tif=gd+n-9999:GTiff'.format(src_grd))
    out, status = run_cmd(num_msk_cmd, verbose = True)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_grdcut(src_grd, region, verbose = False):
    '''cut a grid to region using GMT grdcut

    return status code (0 == success)'''
    
    cut_cmd = ('gmt grdcut -V {} -Gtmp.grd {}\
    '.format(src_grd, region_format(region, 'gmt')))
    out, status = run_cmd(cut_cmd, verbose = True)
    if status == 0:
        remove_glob(src_grd)
        os.rename('tmp.grd', '{}'.format(src_grd))
    return(status)

def gmt_slope(src_dem, dst_slp, verbose = False):
    '''generate a Slope grid from a DEM with GMT

    return status code (0 == success)'''
    
    o_b_name = '{}'.format(src_dem.split('.')[0])
    slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
    '.format(src_dem, o_b_name, src_dem))
    out, status = run_cmd(slope_cmd0, verbose = verbose)
    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}=gd+n-9999:GTiff\
        '.format(o_b_name, dst_slp))
        out, status = run_cmd(slope_cmd1, verbose = verbose)
    remove_glob('{}_pslp.grd'.format(o_b_name))
    return(status)

def gmt_num_msk(num_grd, dst_msk, verbose = False):
    '''generate a num-msk from a NUM grid using GMT grdmath

    returns [cmd-output, cmd-return-code]'''
    
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    return(run_cmd(num_msk_cmd, verbose = verbose))

def gmt_sample_gnr(src_grd, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node
    grid registration.

    returns status code (0 == success)'''
    
    out, status = run_cmd('gmt grdsample -T {} -Gtmp.tif=gd+n-9999:GTiff'.format(src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_sample_inc(src_grd, inc = 1, verbose = False):
    '''resamele src_grd to increment `inc` using GMT grdsample

    returns status code (0 == success)'''
    
    out, status = run_cmd('gmt grdsample -I{:.10f} {} -R{} -Gtmp.tif=gd+n-9999:GTiff'.format(inc, src_grd, src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

## ==============================================
## MB-System Wrapper Functions - mbsfun.py
##
## MS-System must be installed on the system to run
## these functions and commands.
## ==============================================
def mb_inf(src_xyz, src_fmt = 168):
    '''generate an info (.inf) file from a src_xyz file using MBSystem.

    return inf_parse(inf_file)'''
    
    run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz))
    return(inf_parse('{}.inf'.format(src_xyz)))

## ==============================================
## GDAL Wrappers and Functions - gdalfun.py
## ==============================================
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.UseExceptions()
_gdal_progress = gdal.TermProgress #crashes on osgeo4w
_gdal_progress_nocb = gdal.TermProgress_nocb
    
def gdal_sr_wkt(epsg, esri = False):
    '''convert an epsg code to wkt

    returns the sr Wkt or None'''
    
    try:
        int(epsg)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)

def gdal_fext(src_drv_name):
    '''find the common file extention given a GDAL driver name
    older versions of gdal can't do this, so fallback to known standards.

    returns list of known file extentions or None'''
    
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER): fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None: return(fexts.split()[0])
        else: return(None)
    except:
        if src_drv_name == 'GTiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        return(fext)

def gdal_prj_file(dst_fn, epsg):
    '''generate a .prj file given an epsg code

    returns 0'''
    
    with open(dst_fn, 'w') as out:
        out.write(gdal_sr_wkt(int(epsg), True))
    return(0)
    
def gdal_set_epsg(src_fn, epsg = 4326):
    '''set the projection of gdal file src_fn to epsg

    returns status-code (0 == success)'''
    
    ds = gdal.Open(src_fn, gdal.GA_Update)
    if ds is not None:
        ds.SetProjection(gdal_sr_wkt(int(epsg)))
        ds = None
        return(0)
    else: return(None)

def gdal_set_nodata(src_fn, nodata = -9999):
    '''set the nodata value of gdal file src_fn

    returns 0'''
    
    ds = gdal.Open(src_fn, gdal.GA_Update)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    ds = None
    return(0)

def gdal_infos(src_fn, scan = False):
    '''scan gdal file src_fn and gather region info.

    returns region dict.'''
    
    if os.path.exists(src_fn):
        ds = gdal.Open(src_fn)
        if ds is not None:
            dsc = gdal_gather_infos(ds, scan = scan)
            ds = None
            return(dsc)
        else: return(None)
    else: return(None)

def gdal_gather_infos(src_ds, scan = False):
    '''gather information from `src_ds` GDAL dataset

    returns gdal_config dict.'''

    src_band = src_ds.GetRasterBand(1)
    ds_config = {
        'nx': src_ds.RasterXSize,
        'ny': src_ds.RasterYSize,
        'nb':src_ds.RasterCount,
        'geoT': src_ds.GetGeoTransform(),
        'proj': src_ds.GetProjectionRef(),
        'dt': src_band.DataType,
        'dtn': gdal.GetDataTypeName(src_band.DataType),
        'ndv': src_band.GetNoDataValue(),
        'fmt': src_ds.GetDriver().ShortName,
        'zr': None,
    }
    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
    if scan:
        try:
            ds_config['zr'] = src_band.ComputeRasterMinMax()
        except: ds_config = None    
    return(ds_config)

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''set a datasource config dictionary

    returns gdal_config dict.'''
    
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def gdal_cpy_infos(src_config):
    '''copy src_config

    returns copied src_config dict.'''
    
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def gdal_write (src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff'):
    '''write src_arr to gdal file dst_gdal using src_config

    returns [output-gdal, status-code]'''
    
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception as e:
            echo_error_msg(e)
            remove_glob(dst_gdal)
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        ds.SetProjection(ds_config['proj'])
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        ds.GetRasterBand(1).WriteArray(src_arr)
        ds = None
        return(dst_gdal, 0)
    else: return(None, -1)

def gdal_null(dst_fn, region, inc, nodata = -9999, outformat = 'GTiff'):
    '''generate a `null` grid with gdal

    returns [output-gdal, status-code]'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    null_array = np.zeros((ycount, xcount))
    null_array[null_array == 0] = nodata
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(4326), gdal.GDT_Float32, -9999, outformat)
    return(gdal_write(null_array, dst_fn, ds_config))
    
def gdal_cut(src_gdal, region, dst_fn):
    '''cut src_fn gdal file to srcwin and output dst_fn gdal file

    returns [output-gdal, status-code]'''
    
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        ds_config = gdal_gather_infos(src_ds)
        srcwin = gdal_srcwin(src_ds, region)
        gt = ds_config['geoT']
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
        ds_config = gdal_set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, gdal_sr_wkt(4326), ds_config['dt'], ds_config['ndv'], ds_config['fmt'])
        src_ds = None
        return(gdal_write(ds_arr, dst_fn, ds_config))
    else: return(None, -1)

def gdal_clip(src_gdal, src_ply = None, invert = False):
    '''clip dem to polygon `src_ply`, optionally invert the clip.

    returns [gdal_raserize-output, gdal_rasterize-return-code]'''

    # clip src_ply to src_gdal extent
    gi = gdal_infos(src_gdal)
    g_region = gdal_gt2region(gi)
    tmp_ply = 'tmp_clp_ply.shp'
    run_cmd('ogr2ogr {} {} -clipsrc {} {} {} {}'.format(tmp_ply, src_ply, g_region[0], g_region[3], g_region[1], g_region[2]), verbose = True)
    if gi is not None and src_ply is not None:
        gr_inv = '-i' if invert else ''
        gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'\
                 .format(gi['ndv'], gr_inv, os.path.basename(tmp_ply).split('.')[0], tmp_ply, src_gdal)
        out, status = run_cmd(gr_cmd, verbose = True)
    else: return(None)
    remove_glob('tmp_clp_ply.*')
    return(out, status)

def gdal_query(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    returns array of values'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()

        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                #print(g)
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    xyzl.append(np.array(outs))
        dsband = ds = None
        out_array = np.array(xyzl)
    return(out_array)

def gdal_yield_query(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()
        dsband = ds = None
        
        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    yield(outs)

def gdal_query2(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()

        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                #print(g)
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    #c = con_dec(math.fabs(d / (g + 0.00000001) * 100), 2)
                    #s = con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    #d = con_dec(d, 4)
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    yield(outs)
        dsband = ds = None
                    
def np_split(src_arr, sv = 0, nd = -9999):
    '''split numpy `src_arr` by `sv` (turn u/l into `nd`)

    returns [upper_array, lower_array]'''
    
    try:
        sv = int(sv)
    except: sv = 0
    u_arr = np.array(src_arr)
    l_arr = np.array(src_arr)
    u_arr[u_arr <= sv] = nd
    l_arr[l_arr >= sv] = nd
    return(u_arr, l_arr)

def gdal_split(src_gdal, split_value = 0):
    '''split raster file `src_gdal`into two files based on z value

    returns [upper_grid-fn, lower_grid-fn]'''
    
    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_u.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_l.tif'.format(os.path.basename(src_gdal)[:-4]))
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = gdal_gather_infos(src_ds)
        dst_config = gdal_cpy_infos(src_config)
        dst_config['fmt'] = 'GTiff'
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny'])
        ua, la = np_split(ds_arr, split_value, src_config['ndv'])
        gdal_write(ua, dst_upper, dst_config)
        gdal_write(la, dst_lower, dst_config)
        ua = la = ds_arr = src_ds = None
        return([dst_upper, dst_lower])
    else: return(None)

def gdal_crop(src_ds):
    '''crop `src_ds` GDAL datasource by it's NoData value. 

    returns [cropped array, cropped_gdal_config].'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_arr = src_ds.GetRasterBand(1).ReadAsArray()

    src_arr[elev_array == ds_config['ndv']] = np.nan
    nans = np.isnan(src_arr)
    nancols = np.all(nans, axis=0)
    nanrows = np.all(nans, axis=1)

    firstcol = nancols.argmin()
    firstrow = nanrows.argmin()        
    lastcol = len(nancols) - nancols[::-1].argmin()
    lastrow = len(nanrows) - nanrows[::-1].argmin()

    dst_arr = src_arr[firstrow:lastrow,firstcol:lastcol]
    src_arr = None

    dst_arr[np.isnan(dst_arr)] = ds_config['nv']
    GeoT = ds_config['geoT']
    dst_x_origin = GeoT[0] + (GeoT[1] * firstcol)
    dst_y_origin = GeoT[3] + (GeoT[5] * firstrow)
    dst_geoT = [dst_x_origin, GeoT[1], 0.0, dst_y_origin, 0.0, GeoT[5]]
    ds_config['geoT'] = dst_geoT
    return(dst_arr, ds_config)
        
def gdal_gdal2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, dst_gdal = None, co = True):
    '''convert the gdal file to gdal using gdal

    return output-gdal-fn'''
    
    if os.path.exists(src_grd):
        if dst_gdal is None:
            dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
        if not co:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {}\
            '.format(src_grd, dst_gdal, dst_fmt))
        else:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {} -co TILED=YES -co COMPRESS=DEFLATE\
            '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(gdal2gdal_cmd, verbose = False)
        if status == 0: return(dst_gdal)
        else: return(None)
    else: return(None)

def gdal_sum(src_gdal):
    '''sum the z vale of src_gdal

    return the sum'''
    
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_array = ds.GetRasterBand(1).ReadAsArray() 
        sums = np.sum(ds_array)
        ds = ds_array = None
        return(sums)
    else: return(None)

def gdal_percentile(src_gdal, perc = 95):
    '''calculate the `perc` percentile of src_fn gdal file.

    return the calculated percentile'''
    
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        x_dim = ds_array.shape[0]
        ds_array_flat = ds_array.flatten()
        p = np.percentile(ds_array_flat, perc)
        percentile = 2 if p < 2 else p
        ds = ds_array = None
        return(percentile)
    else: return(None)

def gdal_mask(src_gdal, dst_gdal, invert = False):
    '''transform src_gdal to a raster mask (1 = data; 0 = nodata)
    if invert is True, 1 = nodata, 0 = data.'''
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_band = ds.GetRasterBand(1)
        ds_array = ds_band.ReadAsArray()
        ds_config = gdal_gather_infos(ds)
        ds_config['dtn'] = 'Int32'
        ndv = ds_band.GetNoDataValue()
        #print(ndv)
        if not invert:
            ds_array[ds_array != ndv] = 1
            ds_array[ds_array == ndv] = 0
        else:
            ds_array[ds_array != ndv] = 0
            ds_array[ds_array == ndv] = 1
            
        return(gdal_write(ds_array, dst_gdal, ds_config))
    else: return(-1, -1)
    
def gdal_mask_analysis(mask = None):
    '''mask is a GDAL mask grid of 0/1

    returns [sum, max, percentile]'''
    
    msk_sum = gdal_sum(mask)
    msk_gc = gdal_infos(mask)
    msk_max = float(msk_gc['nx'] * msk_gc['ny'])
    msk_perc = float((msk_sum / msk_max) * 100.)
    return(msk_sum, msk_max, msk_perc)
    
def gdal_proximity(src_fn, dst_fn):
    '''compute a proximity grid via GDAL

    return 0 if success else None'''
    
    prog_func = None
    src_ds = gdal.Open(src_fn)
    dst_ds = None
    if src_ds is not None:
        src_band = src_ds.GetRasterBand(1)
        ds_config = gdal_gather_infos(src_ds)

        if dst_ds is None:
            drv = gdal.GetDriverByName('GTiff')
            dst_ds = drv.Create(dst_fn, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'], [])
        dst_ds.SetGeoTransform(ds_config['geoT'])
        dst_ds.SetProjection(ds_config['proj'])
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(ds_config['ndv'])
        gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)
        dst_band = src_band = dst_ds = src_ds = None
        return(0)
    else: return(None)
    
def gdal_gt2region(ds_config):
    '''convert a gdal geo-tranform to a region [xmin, xmax, ymin, ymax] via a data-source config dict.

    returns region of gdal data-source'''
    
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + geoT[1] * ds_config['nx'], geoT[3] + geoT[5] * ds_config['ny'], geoT[3]])

def gdal_region2gt(region, inc):
    '''return a count info and a gdal geotransform based on extent and cellsize

    returns a list [xcount, ycount, geot]'''

    dst_gt = (region[0], inc, 0, region[3], 0, (inc * -1.))
    
    this_origin = _geo2pixel(region[0], region[3], dst_gt)
    this_end = _geo2pixel(region[1], region[2], dst_gt)
    #this_size = ((this_end[0] - this_origin[0]) + 1, (this_end[1] - this_origin[1]) + 1)
    this_size = (this_end[0] - this_origin[0], this_end[1] - this_origin[1])
    
    return(this_size[0], this_size[1], dst_gt)

def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class'''
    
    if dst_defn is None: dst_defn = src_layer.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    feats = len(src_layer)
    echo_msg('unioning {} features'.format(feats))
    for n, f in enumerate(src_layer):
        _gdal_progress_nocb((n+1 / feats) * 100)
        if f.GetField(src_field) == 0:
            src_layer.DeleteFeature(f.GetFID())
        elif f.GetField(src_field) == 1:
            f.geometry().CloseRings()
            wkt = f.geometry().ExportToWkt()
            multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
            src_layer.DeleteFeature(f.GetFID())
    #union = multi.UnionCascaded() ## slow on large multi...
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometryDirectly(multi)
    #union = multi = None
    return(out_feat)

def gdal_ogr_regions(src_ds):
    '''return the region(s) of the ogr dataset'''
    
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                these_regions.append(pgeom.GetEnvelope())
        poly = None
    return(these_regions)

def gdal_create_polygon(coords):
    '''convert coords to Wkt

    returns polygon as wkt'''
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[1], coord[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def gdal_region2wkt(region):

    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    return(gdal_create_polygon(eg))
    
def gdal_region2geom(region):
    '''convert an extent [west, east, south, north] to an OGR geometry

    returns ogr geometry'''
    
    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    geom = ogr.CreateGeometryFromWkt(gdal_create_polygon(eg))
    return(geom)

def gdal_getEPSG(src_ds):
    '''returns the EPSG of the given gdal data-source'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_region = gdal_gt2region(ds_config)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)

    return(srs_auth)

def gdal_region(src_ds, warp = None):
    '''return the extent of the src_fn gdal file.
    warp should be an epsg to warp the region to.

    returns the region of the gdal data-source'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_region = gdal_gt2region(ds_config)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)
    
    if srs_auth is None or srs_auth == warp: warp = None

    if warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(warp))
        #dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)

        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(ds_region[0], ds_region[2]))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(ds_region[1], ds_region[3]))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        ds_region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(ds_region)

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''convert a geographic x,y value to a pixel location of geoTransform'''
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = ((geo_x - geoTransform[0]) / geoTransform[1])# + .5
        pixel_y = ((geo_y - geoTransform[3]) / geoTransform[5])# + .5
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
    return(int(pixel_x), int(pixel_y))

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''convert a pixel location to geographic coordinates given geoTransform'''
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    '''apply geotransform to in_x,in_y'''
    
    out_x = geoTransform[0] + (in_x + 0.5) * geoTransform[1] + (in_y + 0.5) * geoTransform[2]
    out_y = geoTransform[3] + (in_x + 0.5) * geoTransform[4] + (in_y + 0.5) * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    '''invert the geotransform'''
    
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001: return
    invDet = 1.0 / det
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet
    return(outGeoTransform)

def gdal_srcwin(src_ds, region):
    '''given a gdal file src_fn and a region [w, e, s, n],
    output the appropriate gdal srcwin.

    returns the gdal srcwin'''
    
    ds_config = gdal_gather_infos(src_ds)
    this_origin = [0 if x < 0 else x for x in _geo2pixel(region[0], region[3], ds_config['geoT'])]
    this_end = [0 if x < 0 else x for x in _geo2pixel(region[1], region[2], ds_config['geoT'])]
    this_size = [0 if x < 0 else x for x in ((this_end[0] - this_origin[0]), (this_end[1] - this_origin[1]))]
    if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
    if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
    #this_size = [0 if x < 0 else x for x in this_size]
    return(this_origin[0], this_origin[1], this_size[0], this_size[1])

def xyz2gdal_ds(src_xyz, dst_ogr):
    '''Make a point vector OGR DataSet Object from src_xyz

    returns the in-memory GDAL data-source'''
    
    ds = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
    layer = ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPoint25D)
    fd = ogr.FieldDefn('long', ogr.OFTReal)
    fd.SetWidth(10)
    fd.SetPrecision(8)
    layer.CreateField(fd)
    fd = ogr.FieldDefn('lat', ogr.OFTReal)
    fd.SetWidth(10)
    fd.SetPrecision(8)
    layer.CreateField(fd)
    fd = ogr.FieldDefn('elev', ogr.OFTReal)
    fd.SetWidth(12)
    fd.SetPrecision(12)
    layer.CreateField(fd)
    f = ogr.Feature(feature_def = layer.GetLayerDefn())
    for this_xyz in src_xyz:
        f.SetField(0, this_xyz[0])
        f.SetField(1, this_xyz[1])
        f.SetField(2, this_xyz[2])
        wkt = 'POINT({:.8f} {:.8f} {:.10f})'.format(this_xyz[0], this_xyz[1], this_xyz[2])
        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)
    return(ds)

def gdal_xyz2gdal(src_xyz, dst_gdal, region, inc, dst_format = 'GTiff', mode = 'n', epsg = 4326, verbose = False):
    '''Create a GDAL supported grid from xyz data 
    `mode` of `n` generates a num grid
    `mode` of `m` generates a mean grid
    `mode` of `k` generates a mask grid
    `mode` of 'w' generates a wet/dry mask grid

    returns output, status'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    if verbose:
        echo_msg('gridding data with mode: {} to {}'.format(mode, dst_gdal))
        echo_msg('grid size: {}/{}'.format(ycount, xcount))
    if mode == 'm' or mode == 'w':
        sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    #else: gdt = gdal.GDT_Int32
    ptArray = np.zeros((ycount, xcount))
    if mode == 'w': ptArray[ptArray == 0] = np.nan
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(epsg), gdt, -9999, dst_format)
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    if mode == 'm' or mode == 'w':
                        sumArray[ypos, xpos] += z
                    if mode == 'n' or mode == 'm':
                        ptArray[ypos, xpos] += 1
                    else: ptArray[ypos, xpos] = 1
                except Exception as e: echo_error_msg(e)
    if mode == 'm' or mode == 'w':
        ptArray[ptArray == 0] = np.nan
        outarray = sumArray / ptArray
        if mode == 'w':
            outarray[outarray >= 0] = 1
            outarray[outarray < 0] = 0
    elif mode == 'n': outarray = ptArray
    else: outarray = ptArray
    outarray[np.isnan(outarray)] = -9999
    return(gdal_write(outarray, dst_gdal, ds_config))

def gdal_xyz_mask(src_xyz, dst_gdal, region, inc, dst_format='GTiff', epsg = 4326):
    '''Create a num grid mask of xyz data. The output grid
    will contain 1 where data exists and 0 where no data exists.

    yields the xyz data'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    ptArray = np.zeros((ycount, xcount))
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(epsg), gdal.GDT_Float32, -9999, 'GTiff')
    for this_xyz in src_xyz:
        yield(this_xyz)
        x = this_xyz[0]
        y = this_xyz[1]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                except: pass
    out, status = gdal_write(ptArray, dst_gdal, ds_config)

def np_gaussian_blur(in_array, size):
    '''blur an array using fftconvolve from scipy.signal
    size is the blurring scale-factor.

    returns the blurred array'''
    
    from scipy.signal import fftconvolve
    from scipy.signal import convolve
    import scipy.fftpack._fftpack as sff
    padded_array = np.pad(in_array, size, 'symmetric')
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(in_array.dtype)
    in_array = None
    #try:
    out_array = fftconvolve(padded_array, g, mode = 'valid')
    #except:
    #print('switching to convolve')
    #out_array = convolve(padded_array, g, mode = 'valid')
    return(out_array)

def gdal_blur(src_gdal, dst_gdal, sf = 1):
    '''gaussian blur on src_gdal using a smooth-factor of `sf`
    runs np_gaussian_blur(ds.Array, sf)'''
    
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
        msk_array = np.array(ds_array)
        msk_array[msk_array != ds_config['ndv']] = 1
        msk_array[msk_array == ds_config['ndv']] = np.nan
        ds_array[ds_array == ds_config['ndv']] = 0
        smooth_array = np_gaussian_blur(ds_array, int(sf))
        smooth_array = smooth_array * msk_array
        mask_array = ds_array = None
        smooth_array[np.isnan(smooth_array)] = ds_config['ndv']
        return(gdal_write(smooth_array, dst_gdal, ds_config))
    else: return([], -1)

def gdal_smooth(src_gdal, dst_gdal, fltr = 10, split_value = None, use_gmt = False):
    '''smooth `src_gdal` using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero)

    return 0 for success or -1 for failure'''
    
    if os.path.exists(src_gdal):
        if split_value is not None:
            dem_u, dem_l = gdal_split(src_gdal, split_value)
        else: dem_l = src_gdal
        if use_gmt:
            out, status = gmt_grdfilter(dem_l, 'tmp_fltr.tif=gd+n-9999:GTiff', dist = fltr, verbose = True)
        else: out, status = gdal_blur(dem_l, 'tmp_fltr.tif', fltr)
        if split_value is not None:
            ds = gdal.Open(src_gdal)
            ds_config = gdal_gather_infos(ds)
            msk_arr = ds.GetRasterBand(1).ReadAsArray()
            msk_arr[msk_arr != ds_config['ndv']] = 1
            msk_arr[msk_arr == ds_config['ndv']] = np.nan
            ds = None
            u_ds = gdal.Open(dem_u)
            if u_ds is not None:
                l_ds = gdal.Open('tmp_fltr.tif')
                if l_ds is not None:
                    u_arr = u_ds.GetRasterBand(1).ReadAsArray()
                    l_arr = l_ds.GetRasterBand(1).ReadAsArray()
                    u_arr[u_arr == ds_config['ndv']] = 0
                    l_arr[l_arr == ds_config['ndv']] = 0
                    ds_arr = (u_arr + l_arr) * msk_arr
                    ds_arr[np.isnan(ds_arr)] = ds_config['ndv']
                    gdal_write(ds_arr, 'merged.tif', ds_config)
                    l_ds = None
                    remove_glob(dem_l)
                u_ds = None
                remove_glob(dem_u)
            os.rename('merged.tif', 'tmp_fltr.tif')
        os.rename('tmp_fltr.tif', dst_gdal)
        return(0)
    else: return(-1)

def gdal_sample_inc(src_grd, inc = 1, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node grid registration.'''
    
    out, status = run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te -R{} -r -Gtmp.tif=gd+n-9999:GTiff'.format(inc, inc, src_grd, src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)
        
def gdal_polygonize(src_gdal, dst_layer, verbose = False):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''
    
    ds = gdal.Open('{}'.format(src_gdal))
    ds_arr = ds.GetRasterBand(1)
    if verbose: echo_msg('polygonizing {}'.format(src_gdal))
    status = gdal.Polygonize(ds_arr, None, dst_layer, 0, callback = _gdal_progress if verbose else None)
    ds = ds_arr = None
    return(0, 0)
    
def gdal_chunks(src_fn, n_chunk = 10):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.

    returns a list of chunked filenames or None'''
    
    band_nums = []
    o_chunks = []
    if band_nums == []: band_nums = [1]
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk

    try:
        src_ds = gdal.Open(src_fn)
    except: src_ds = None
    if src_ds is not None:
        ds_config = gdal_gather_infos(src_ds)
        band = src_ds.GetRasterBand(1)
        gt = ds_config['geoT']

        while True:
            y_chunk = n_chunk
            while True:
                if x_chunk > ds_config['nx']:
                    this_x_chunk = ds_config['nx']
                else: this_x_chunk = x_chunk

                if y_chunk > ds_config['ny']:
                    this_y_chunk = ds_config['ny']
                else: this_y_chunk = y_chunk

                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = this_x_chunk - this_x_origin
                this_y_size = this_y_chunk - this_y_origin

                ## chunk size aligns with grid cellsize
                if this_x_size == 0 or this_y_size == 0: break
                
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)
                dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
                
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                if not np.all(band_data == band_data[0,:]):
                    o_chunk = '{}_chnk{}x{}.tif'.format(os.path.basename(src_fn).split('.')[0], x_i_chunk, i_chunk)
                    dst_fn = os.path.join(os.path.dirname(src_fn), o_chunk)
                    o_chunks.append(dst_fn)

                    dst_config = gdal_cpy_infos(ds_config)
                    dst_config['nx'] = this_x_size
                    dst_config['ny'] = this_y_size
                    dst_config['geoT'] = dst_gt                    
                    gdal_write(band_data, dst_fn, dst_config)

                band_data = None

                if y_chunk > ds_config['ny']:
                    break
                else: 
                    y_chunk += n_chunk
                    i_chunk += 1
            if x_chunk > ds_config['nx']:
                break
            else:
                x_chunk += n_chunk
                x_i_chunk += 1
        src_ds = None
        return(o_chunks)
    else: return(None)

def gdal_slope(src_gdal, dst_gdal, s = 111120):
    '''generate a slope grid with GDAL

    return cmd output and status'''
    
    gds_cmd = 'gdaldem slope {} {} {} -compute_edges'.format(src_gdal, dst_gdal, '' if s is None else '-s {}'.format(s))
    return(run_cmd(gds_cmd))
    
## ==============================================
## inf files (data info) inf.py
## mbsystem/gmt/waffles infos
## ==============================================
def inf_generate(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    
    return(inf_entry([data_path, data_fmt]), True)

def inf_parse(src_inf):
    '''parse an inf file (mbsystem or gmt)
    
    returns region: [xmin, xmax, ymin, ymax, zmin, zmax]'''
    
    minmax = [0, 0, 0, 0, 0, 0]
    with open(src_inf) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                try: 
                    minmax = [float(x) for x in til]
                except:
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            minmax[0] = til[2]
                            minmax[1] = til[5]
                        elif til[1] == 'Latitude:':
                            minmax[2] = til[2]
                            minmax[3] = til[5]
                        # elif til[1] == 'Altitude:':
                        #     minmax[4] = til[2]
                        #     minmax[5] = til[5]
                        elif til[1] == 'Depth:':
                            minmax[4] = float(til[5])*-1
                            minmax[5] = float(til[2])*-1
    return([float(x) for x in minmax])

def inf_entry(src_entry, overwrite = False):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file or the 
    result of `gmt gmtinfo file.xyz -C`, which is a 6 column line 
    with minmax info, etc.

    returns the region of the inf file.'''
    
    minmax = None
    if os.path.exists(src_entry[0]):
        path_i = src_entry[0] + '.inf'
        if not os.path.exists(path_i) or overwrite:
            minmax = _dl_inf_h[src_entry[1]](src_entry)
        else: minmax = inf_parse(path_i)
    if not region_valid_p(minmax): minmax = None
    return(minmax)

## ==============================================
## gdal processing (datalist fmt:200)
## ==============================================
def gdal_parse(src_ds, dump_nodata = False, srcwin = None, mask = None, warp = None, verbose = False, z_region = None, step = 1):
    '''parse the data from gdal dataset src_ds (first band only)
    optionally mask the output with `mask` or transform the coordinates to `warp` (epsg-code)

    yields the parsed xyz data'''

    #if verbose: sys.stderr.write('waffles: parsing gdal file {}...'.format(src_ds.GetDescription()))
    ln = 0
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_gather_infos(src_ds)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)
    
    if srs_auth is None or srs_auth == warp: warp = None

    if warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(warp))
        ## GDAL 3+
        #dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        
    gt = ds_config['geoT']
    msk_band = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
    if band.GetNoDataValue() is not None: nodata.append('{:g}'.format(band.GetNoDataValue()))
    if dump_nodata: nodata = []
    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        if z_region is not None:
            if z_region[0] is not None:
                band_data[band_data < z_region[0]] = -9999
            if z_region[1] is not None:
                band_data[band_data > z_region[1]] = -9999
        if msk_band is not None:
            msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            band_data[msk_data==0]=-9999
        band_data = np.reshape(band_data, (srcwin[2], ))
        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            z = band_data[x_i]
            if '{:g}'.format(z) not in nodata:
                ln += 1
                geo_x,geo_y = _pixel2geo(x, y, gt)
                if warp is not None:
                    point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(geo_x, geo_y))
                    point.Transform(dst_trans)
                    pnt = point.GetPoint()
                    line = [pnt[0], pnt[1], z]
                else: line = [geo_x, geo_y, z]
                yield(line)
    band = None
    src_mask = None
    msk_band = None
    if verbose: echo_msg('parsed {} data records from {}'.format(ln, src_ds.GetDescription()))

def xyz_transform(src_fn, xyz_c = None, inc = None, vdatum = None, region = None, datalist = None, verbose = False):
    '''transform xyz data file src_fn to xyz data'''
    
    src_xyz, src_zips = procs_unzip(src_fn, _known_datalist_fmts[168])
    if xyz_c is None: xyz_c = copy.deepcopy(_xyz_config)
    if not os.path.exists(os.path.join(os.path.dirname(src_fn), 'xyz')): os.mkdir(os.path.join(os.path.dirname(src_fn), 'xyz'))

    if vdatum is not None:
        vds = vdatum.split(',')
        iv = vds[0]
        ov = vds[1]
        vdc = _vd_config
        if vdc['jar'] is None: vdc['jar'] = vdatum_locate_jar()[0]
        vdc['ivert'] = iv
        vdc['overt'] = ov
        vdc['delim'] = 'comma' if xyz_c['delim'] == ',' else 'space' if xyz_c['delim'] == ' ' else xyz_c['delim']
        vdc['skip'] = xyz_c['skip']
        vdc['xyzl'] = ','.join([str(x) for x in [xyz_c['xpos'], xyz_c['ypos'], xyz_c['zpos']]])
        out, status = run_vdatum(src_xyz, vdc)
        src_result = os.path.join('result', os.path.basename(src_xyz))
    else: src_result = src_xyz
    xyz_final = os.path.join('xyz', os.path.basename(src_xyz))

    if datalist is not None:
        datalist = os.path.join(os.path.dirname(src_fn), 'xyz', datalist)
    
    if os.path.exists(src_result):
        with open(src_result, 'r') as in_n, open(xyz_final, 'w') as out_n:
            for xyz in xyz_parse(in_n, xyz_c = xyz_c, region = region, verbose = verbose):
                xyz_line(xyz, out_n)

    if os.path.exists(xyz_final):
        if datalist is not None:
            mb_inf(xyz_final)
            datalist_append_entry([os.path.basename(xyz_final), 168, 1], datalist)
            if verbose: echo_msg('appended xyz file {} to datalist {}'.format(xyz_final, datalist))
                
    if src_xyz != src_fn: remove_glob(src_xyz)
    
def gdal2xyz_chunks(src_fn, chunk_value = 1000, inc  = None, epsg = None, vdatum = None, datalist = None, verbose = False):
    '''chunk and transform gdal file src_fn to xyz file'''
    
    if verbose:
        echo_msg('------------------------------------')
        echo_msg('input gdal:\t\t{}'.format(src_fn))
        echo_msg('chunk size:\t\t{}'.format(chunk_value))
        echo_msg('output epsg:\t\t{}'.format(epsg))
        echo_msg('output increment:\t{}'.format(inc))
        echo_msg('vdatum string:\t{}'.format(vdatum))
        echo_msg('output datalist:\t{}'.format(datalist))
        echo_msg('------------------------------------')

    if not os.path.exists('xyz'): os.mkdir('xyz')

    src_gdal, src_zips = procs_unzip(src_fn, _known_datalist_fmts[200])
    src_c = gdal_infos(src_gdal)
    if verbose:
        echo_msg('{}'.format(src_c))
        echo_msg('chunking grid file {}...'.format(src_gdal))
    chunks = gdal_chunks(src_gdal, chunk_value)
    if verbose: echo_msg('generated {} chunks from {}'.format(len(chunks), src_gdal))
    if src_gdal != src_fn: remove_glob(src_gdal)

    if vdatum is not None:
        vds = vdatum.split(',')
        if len(vds) < 2:
            echo_error_msg('bad vdatum string {}'.format(vdatum))
            vdatum = None
        else:
            iv = vds[0]
            ov = vds[1]
            vdc = _vd_config
            if vdc['jar'] is None: vdc['jar'] = vdatum_locate_jar()[0]
            vdc['ivert'] = iv
            vdc['overt'] = ov
            vdc['delim'] = 'space'
            vdc['skip'] = '0'
            vdc['xyzl'] = '0,1,2'

    if datalist is not None:
        datalist = os.path.join('xyz', datalist)
        
    for i,chunk in enumerate(chunks):
        echo_msg('* processing chunk {} [{}/{}]...'.format(chunk, i+1, len(chunks)))
        xyz_chunk = '{}.xyz'.format(chunk.split('.')[0])
        xyz_chunk_final = os.path.join('xyz', os.path.basename(xyz_chunk))

        if epsg is not None or inc is not None:
            tif_chunk = '{}_warp.tif'.format(chunk.split('.')[0])
            gdw = 'gdalwarp {} -dstnodata -9999 -overwrite {}'.format(chunk, tif_chunk)
            if epsg is not None: gdw += ' -t_srs EPSG:{}'.format(epsg)
            if inc is not None: gdw += ' -tr {} {}'.format(inc, inc)
            out, status = run_cmd(gdw, verbose = verbose)
            remove_glob(chunk)
        else: tif_chunk = chunk

        with open(xyz_chunk, 'w') as xyz_c:
            gdal_dump_entry([tif_chunk, 200, None], dst_port = xyz_c, verbose = verbose)
        remove_glob(tif_chunk)

        if vdatum is not None:
            out, status = run_vdatum(xyz_chunk, vdc)
            remove_glob(xyz_chunk)
            xyz_chunk = os.path.join('result', os.path.basename(xyz_chunk)) 
            os.rename(xyz_chunk, xyz_chunk_final)
            vdatum_clean_result()

            if verbose: echo_msg('transformed {} chunk to {}'.format(iv, ov))
        else: os.rename(xyz_chunk, xyz_chunk_final)

        if datalist is not None:
            mb_inf(xyz_chunk_final)
            datalist_append_entry([os.path.basename(xyz_chunk_final), 168, 1], datalist)
            if verbose: echo_msg('appended xyz chunk {} to datalist {}'.format(xyz_chunk_final, datalist))
    
def gdal_inf(src_ds, warp = None):
    '''generate an info (.inf) file from a src_gdal file using gdal

    returns the region [xmin, xmax, ymin, ymax] of src_ds'''
    
    minmax = gdal_region(src_ds, warp)
    try: zr = src_ds.GetRasterBand(1).ComputeRasterMinMax()
    except: zr = [None, None]
    minmax = minmax + list(zr)
    with open('{}.inf'.format(src_ds.GetDescription()), 'w') as inf:
        echo_msg('generating inf file for {}'.format(src_ds.GetDescription()))
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)
                    
def gdal_inf_entry(entry, warp = None):
    ''' scan a gdal entry and find it's region

    returns the region [xmin, xmax, ymin, ymax] of the gdal entry'''
    
    ds = gdal.Open(entry[0])
    minmax = gdal_inf(ds, warp)
    ds = None
    return(minmax)

def gdal_yield_entry(entry, region = None, verbose = False, epsg = None, z_region = None):
    '''yield the xyz data from the datalist entry.

    yields [x, y, z, <w, ...>]'''
    
    ds = gdal.Open(entry[0])
    if region is not None:
        srcwin = gdal_srcwin(ds, region)
    else: srcwin = None
    for xyz in gdal_parse(ds, dump_nodata = False, srcwin = srcwin, warp = epsg, verbose = verbose, z_region = z_region):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)
    ds = None
    
def gdal_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, epsg = None, z_region = None):
    '''dump the xyz data from the gdal entry to dst_port'''
    
    for xyz in gdal_yield_entry(entry, region, verbose, epsg, z_region):
        xyz_line(xyz, dst_port, True)

## ==============================================
## fetches processing (datalists fmt:400 - 499)
## ==============================================
def fetch_yield_entry(entry = ['nos:datatype=xyz'], region = None, verbose = False):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''
    
    fetch_mod = entry[0].split(':')[0]
    fetch_args = entry[0].split(':')[1:]
    
    fl = fetches.fetch_infos[fetch_mod][0](region_buffer(region, 5, pct = True), [], lambda: False)
    args_d = args2dict(fetch_args, {})
    fl._verbose = verbose

    for xyz in fl._yield_results_to_xyz(**args_d):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dump_entry(entry = ['nos:datatype=nos'], dst_port = sys.stdout, region = None, verbose = False):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    
    for xyz in fetch_yield_entry(entry, region, verbose):
        xyz_line(xyz, dst_port, True)
        
def fetch_module_yield_entry(entry, region = None, verbose = False, module = 'dc'):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''
    
    fl = fetches.fetch_infos[module][0](region_buffer(region, 5, pct = True), [], lambda: False)
    fl._verbose = verbose
    fetch_entry = [entry[0], entry[0].split('/')[-1], module]
    
    for xyz in fl._yield_xyz(fetch_entry):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dc_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, module = 'dc'):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    
    for xyz in fetch_module_yield_entry(entry, region, verbose, module):
        xyz_line(xyz, dst_port, True)        
        
## ==============================================
## xyz processing (datalists fmt:168)
## ==============================================
_xyz_config = {
    'delim': None,
    'xpos': 0,
    'ypos': 1,
    'zpos': 2,
    'skip': 0,
    'name': '<xyz-data-stream>',
    'upper_limit': None,
    'lower_limit': None}

_known_delims = [',', ' ', '\t', '/', ':']

def xyz_line_delim(xyz):
    for delim in _known_delims:
        this_xyz = xyz.split(delim)
        if len(this_xyz) > 1: return(delim)
    return(None)

def xyz_parse_line(xyz, xyz_c = _xyz_config):
    '''parse an xyz line-string, using _xyz_config

    returns [x, y, z]'''
    
    this_line = xyz.strip()
    if xyz_c['delim'] is None:
        xyz_c['delim'] = xyz_line_delim(this_line)
    this_xyz = this_line.split(xyz_c['delim'])
    try:
        o_xyz = [float(this_xyz[xyz_c['xpos']]), float(this_xyz[xyz_c['ypos']]), float(this_xyz[xyz_c['zpos']])]
    except IndexError as e:
        echo_error_msg(e)
        return(None)
    except Exception as e:
        echo_error_msg(e)
        return(None)
    return(o_xyz)
    
def xyz_parse(src_xyz, xyz_c = _xyz_config, region = None, verbose = False):
    '''xyz file parsing generator
    `src_xyz` is a file object or list of xyz data.

    yields each xyz line as a list [x, y, z, ...]'''
    
    ln = 0
    skip = int(xyz_c['skip'])
    xpos = xyz_c['xpos']
    ypos = xyz_c['ypos']
    zpos = xyz_c['zpos']
    pass_d = True
    #if verbose: echo_msg('parsing xyz data from {}...'.format(xyz_c['name']))
    for xyz in src_xyz:
        pass_d = True
        if ln >= skip:
            this_xyz = xyz_parse_line(xyz, xyz_c)
            if this_xyz is not None:
                if region is not None:
                    if not xyz_in_region_p(this_xyz, region): pass_d = False
                if xyz_c['upper_limit'] is not None or xyz_c['lower_limit'] is not None:
                    if not z_pass(this_xyz[2], upper_limit = xyz_c['upper_limit'], lower_limit = xyz_c['lower_limit']): pass_d = False
            else: pass_d = False
            
            if pass_d:
                ln += 1
                yield(this_xyz)
        else: skip -= 1
    if verbose: echo_msg('parsed {} data records from {}'.format(ln, xyz_c['name']))

def xyz2py(src_xyz):
    '''return src_xyz as a python list'''
    
    xyzpy = []
    return([xyzpy.append(xyz) for xyz in xyz_parse(src_xyz)])

def xyz_block(src_xyz, region, inc, dst_xyz = sys.stdout, weights = False, verbose = False):
    '''block the src_xyz data to the mean block value

    yields the xyz value for each block with data'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    #xcount += 1
    #ycount += 1
    sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    ptArray = np.zeros((ycount, xcount))
    if weights: wtArray = np.zeros((ycount, xcount))
    if verbose: echo_msg('blocking data to {}/{} grid'.format(ycount, xcount))
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if weights:
            w = this_xyz[3]
            z = z * w
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    sumArray[ypos, xpos] += z
                    ptArray[ypos, xpos] += 1
                    if weights: wtArray[ypos, xpos] += w
                except: pass
    ptArray[ptArray == 0] = np.nan
    if weights:
        wtArray[wtArray == 0] = 1
        outarray = (sumArray / wtArray) / ptArray
    else: outarray = sumArray / ptArray

    sumArray = ptArray = None
    if weights: wtArray = None

    outarray[np.isnan(outarray)] = -9999
    
    for y in range(0, ycount):
        for x in range(0, xcount):
            geo_x, geo_y = _pixel2geo(x, y, dst_gt)
            z = outarray[y,x]
            if z != -9999:
                yield([geo_x, geo_y, z])
    
def xyz_line(line, dst_port = sys.stdout, encode = False):
    '''write "xyz" `line` to `dst_port`
    `line` should be a list of xyz values [x, y, z, ...].'''
    delim = _xyz_config['delim'] if _xyz_config['delim'] is not None else ' '
    
    l = '{}\n'.format(delim.join([str(x) for x in line]))
    if encode: l = l.encode('utf-8')
    dst_port.write(l)

def xyz_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n], else False.'''
    
    if src_xy[0] < src_region[0]: return(False)
    elif src_xy[0] > src_region[1]: return(False)
    elif src_xy[1] < src_region[2]: return(False)
    elif src_xy[1] > src_region[3]: return(False)
    else: return(True)
    
def xyz_inf(src_xyz):
    '''scan an xyz file and find it's min/max values and
    write an associated inf file for the src_xyz file.

    returns region [xmin, xmax, ymin, ymax] of the src_xyz file.'''
    
    minmax = []
    for i,l in enumerate(xyz_parse(src_xyz)):
        if i == 0:
            minmax = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            try:
                if l[0] < minmax[0]: minmax[0] = l[0]
                elif l[0] > minmax[1]: minmax[1] = l[0]
                if l[1] < minmax[2]: minmax[2] = l[1]
                elif l[1] > minmax[3]: minmax[3] = l[1]
                if l[2] < minmax[4]: minmax[4] = l[2]
                elif l[2] > minmax[5]: minmax[5] = l[2]
            except: pass
    if len(minmax) == 6:
        with open('{}.inf'.format(src_xyz.name), 'w') as inf:
            #echo_msg('generating inf file for {}'.format(src_xyz.name))
            inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
        return(minmax)
    else: return(0,0,0,0,0,0)

def xyz_inf_entry(entry):
    '''find the region of the xyz datalist entry
    
    returns the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry'''
    
    with open(entry[0]) as infile:
        try:
            minmax = mb_inf(infile)
        except: minmax = xyz_inf(infile)
    return(minmax)        

def xyz_yield_entry(entry, region = None, verbose = False, z_region = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''

    xyzc = copy.deepcopy(_xyz_config)
    xyzc['name'] = entry[0]
    if z_region is not None and len(z_region) >= 2:
        xyzc['lower_limit'] = z_region[0]
        xyzc['upper_limit'] = z_region[1]
    
    with open(entry[0]) as infile:
        for line in xyz_parse(infile, xyz_c = xyzc, region = region, verbose = verbose):
            yield(line + [entry[2]] if entry[2] is not None else line)

def gmt_yield_entry(entry, region = None, verbose = False, z_region = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''
    ln = 0
    delim = ' '
    if z_region is not None: z_region = ['-' if x is None else str(x) for x in z_region]
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    for line in yield_cmd('gmt gmtselect {} {} {}\
    '.format(entry[0], '' if region is None else region_format(region, 'gmt'), '' if z_region is None else '-Z{}'.format('/'.join(z_region))),
                          data_fun = None, verbose = False):
        ln += 1
        #if delim is None: delim = xyz_line_delim(line)
        yield(line)
        #xyz = [float(x) for x in line.split(delim)]
        #yield(xyz + [entry[2]] if entry[2] is not None else xyz)
    if verbose: echo_msg('read {} data points from {}'.format(ln, entry[0]))
        
def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, z_region = None):
    '''dump the xyz data from the xyz datalist entry to dst_port'''
    
    for xyz in xyz_yield_entry(entry, region, verbose, z_region):
        xyz_line(xyz, dst_port, True, None)

## ==============================================
## datalists and entries - datalists.py
##
## datalist processing (datalists fmt:-1)
## entry processing fmt:*
## TODO -> pass_h to list for all entry passing
## ==============================================
_known_dl_delims = [' ']
_known_datalist_fmts = {-1: ['datalist', 'mb-1'], 168: ['xyz', 'csv', 'dat', 'ascii'], 200: ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'], 400: ['nos', 'dc', 'gmrt', 'srtm', 'charts', 'mb']}
_known_datalist_fmts_short_desc = lambda: '\n  '.join(['{}\t{}'.format(key, _known_datalist_fmts[key]) for key in _known_datalist_fmts])
_dl_inf_h = {
    -1: lambda e: datalist_inf_entry(e),
    168: lambda e: xyz_inf_entry(e),
    200: lambda e: gdal_inf_entry(e),
    201: lambda e: gdal_inf_entry(e, 4269)
}
_dl_pass_h = [lambda e: path_exists_or_url(e[0])]

def datalist_default_hooks():
    return([lambda e: path_exists_or_url(e[0])])

def datalist_inf(dl, inf_file = True, overwrite = False):
    '''return the region of the datalist and generate
    an associated `.inf` file if `inf_file` is True.'''
    
    out_regions = []
    minmax = None
    dp_h = _dl_pass_h
    #dp_h.append(lambda e: region_valid_p(inf_entry(e, True)) if e[1] == -1 else region_valid_p(inf_entry(e)))

    echo_msg('generating inf for datalist {}'.format(dl))
    
    for entry in datalist(dl, pass_h = dp_h):
        #entry_inf = inf_entry(entry, True) if entry[1] == -1 else inf_entry(entry)
        if entry[1] == 200:
            entry_inf = inf_entry(entry, True)
        else: entry_inf = inf_entry(entry)
        if entry_inf is not None:
            out_regions.append(inf_entry(entry)[:6])
    
    out_regions = [x for x in out_regions if x is not None]
    if len(out_regions) == 0:
        minmax = None
    elif len(out_regions) == 1:
        minmax = out_regions[0]
    else:
        out_region = out_regions[0]
        for x in out_regions[1:]:
            out_region = regions_merge(out_region, x)
        minmax = out_region
    if minmax is not None and inf_file:
        with open('{}.inf'.format(dl), 'w') as inf:
            inf.write('{}\n'.format(region_format(minmax, 'inf')))
    return(minmax)

def datalist_inf_entry(e):
    '''write an inf file for datalist entry e
    
    return the region [xmin, xmax, ymin, ymax]'''
    
    return(datalist_inf(e[0]))

def datalist_append_entry(entry, datalist):
    '''append entry to datalist file `datalist`'''
    
    with open(datalist, 'a') as outfile:
        outfile.write('{}\n'.format(' '.join([str(x) for x in entry])))

def datalist_archive_yield_entry(entry, dirname = 'archive', region = None, inc = 1, weight = None, verbose = None, z_region = None):
    '''archive a datalist entry.
    a datalist entry is [path, format, weight, ...]

    yield the xyz line data'''
    
    if region is None:
        a_name = entry[-1]
    else: a_name = '{}_{}_{}'.format(entry[-1], region_format(region, 'fn'), this_year())
    
    i_dir = os.path.dirname(entry[0])
    i_xyz = os.path.basename(entry[0]).split('.')[0]
    i_xyz = ''.join(x for x in i_xyz if x.isalnum())
    a_dir = os.path.join(dirname, a_name, 'data', entry[-1])
    a_xyz_dir = os.path.join(a_dir, 'xyz')
    a_xyz = os.path.join(a_xyz_dir, i_xyz + '.xyz')
    a_dl = os.path.join(a_xyz_dir, '{}.datalist'.format(entry[-1]))
    
    if not os.path.exists(a_dir): os.makedirs(a_dir)
    if not os.path.exists(a_xyz_dir): os.makedirs(a_xyz_dir)

    with open(a_xyz, 'w') as fob:
        for xyz in datalist_yield_entry(entry, region = region, verbose = verbose, z_region = z_region):
            xyz_line(xyz, fob)
            yield(xyz)
            
    mb_inf(a_xyz)
    datalist_append_entry([i_xyz + '.xyz', 168, entry[2] if entry[2] is not None else 1], a_dl)
        
def datalist_echo(entry):
    '''echo datalist entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))
    datalist(entry[0])

def datafile_echo(entry):
    '''echo datafile entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))

def datalist_major(dls, major = '.mjr.datalist', region = None):
    '''set the major datalist
    `dls` is a list of datalist entries, minimally: ['datafile.xyz']

    returns the major datalist filename'''
    
    with open(major, 'w') as md:        
        for dl in dls:
            entries = datalist2py(dl, region)
            #print(entries)
            #sys.exit()
            for entry in entries:
                md.write('{}\n'.format(' '.join([str(e) for e in entry])))
    if os.stat(major).st_size == 0:
        remove_glob(major)
        echo_error_msg('bad datalist/entry, {}'.format(dls))
        return(None)    
    return(major)

def entry2py(dl_e):
    '''convert a datalist entry to python

    return the entry as a list [fn, fmt, wt, ...]'''
    
    this_entry = dl_e.rstrip().split()
    try:
        entry = [x if n == 0 else int(x) if n < 2 else float(x) if n < 3 else x for n, x in enumerate(this_entry)]
    except Exception as e:
        echo_error_msg('could not parse entry {}'.format(dl_e))
        return(None)
    if len(entry) < 2:
        for key in _known_datalist_fmts.keys():
            se = entry[0].split('.')
            if len(se) == 1:
                see = entry[0].split(':')[0]
            else: see = se[-1]
            if see in _known_datalist_fmts[key]:
                entry.append(int(key))
    if len(entry) < 3: entry.append(1)
    return(entry)

def datalist2py(dl, region = None):
    '''convert a datalist to python data
    
    returns a list of datalist entries.'''
    
    these_entries = []
    #print(this_dir)
    this_entry = entry2py(dl)
    this_dir = os.path.dirname(this_entry[0])
    if this_entry[1] == -1:
        if os.path.exists(this_entry[0]):
            with open(this_entry[0], 'r') as op:
                for this_line in op:
                    if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                        these_entries.append([os.path.join(this_dir, x) if n == 0 else x for n,x in enumerate(entry2py(this_line.rstrip()))])
                        #these_entries.append(entry2py(this_line.rstrip()))
        else: echo_error_msg('could not open datalist/entry {}'.format(this_entry[0]))
    elif this_entry[1] == 400:
        fetch_mod = this_entry[0].split(':')[0]
        fetch_args = this_entry[0].split(':')[1:]
        if fetch_mod in fetches.fetch_infos.keys():
            fl = fetches.fetch_infos[fetch_mod][0](region_buffer(region, 5, pct = True), [], lambda: False)
            args_d = args2dict(fetch_args, {})
            fl._verbose = True

            results = fl.run(**args_d)
            if len(results) > 0:
                with open('{}.datalist'.format(fetch_mod), 'w') as fdl:
                    for r in results:
                        e = [r[0], fl._datalists_code, 1]
                        fdl.write('{} {} {}\n'.format(e[0], e[1], e[2]))
                these_entries.append(['{}.datalist'.format(fetch_mod), -1, 1])
                
    else: these_entries.append(this_entry)
    return(these_entries)
            
def datalist(dl, wt = None, pass_h = _dl_pass_h,
             yield_dl_entry = False, verbose = False):
    '''recurse a datalist/entry
    for entry in datalist(dl): do_something_with entry

    yields entry [path, fmt, wt, ...]'''

    #this_dir = os.path.dirname(dl)
    these_entries = datalist2py(dl)
    if len(these_entries) == 0: these_entries = [entry2py(dl)]
    for this_entry in these_entries:
        if this_entry is not None:
            #this_entry[0] = os.path.join(this_dir, this_entry[0])
            this_entry[2] = wt if wt is None else wt * this_entry[2] 
            this_entry = this_entry[:3] + [' '.join(this_entry[3:]).split(',')] + [os.path.basename(dl).split('.')[0]]
            if not False in [x(this_entry) for x in pass_h]:
                if verbose and this_entry[1] == -1: echo_msg('parsing datalist ({}) {}'.format(this_entry[2], this_entry[0]))
                if this_entry[1] == -1:
                    if yield_dl_entry: yield(this_entry)
                    for entry in datalist(this_entry[0], wt = this_entry[2], pass_h = pass_h, yield_dl_entry = yield_dl_entry, verbose = verbose):
                        yield(entry)
                else: yield(this_entry)

def datalist_yield_entry(this_entry, region = None, verbose = False, z_region = None, w_region = None):
    '''yield the xyz data from the datalist entry [entry/path, entry-format, entry-weight]

    yields xyz line data [x, y, z, ...]'''

    if this_entry[1] == 168:
        for xyz in xyz_yield_entry(this_entry, region = region, verbose = verbose, z_region = z_region):
            #for xyz in gmt_yield_entry(this_entry, region = region, verbose = verbose, z_region = z_region):
            yield(xyz)
    elif this_entry[1] == 200:
        for xyz in gdal_yield_entry(this_entry, region = region, verbose = verbose, z_region = z_region):
            yield(xyz)
    elif this_entry[1] == 201:
        for xyz in gdal_yield_entry(this_entry, verbose = verbose, z_region = z_region, epsg = 4269):
            yield(xyz)
    elif this_entry[1] == 400:
        for xyz in fetch_yield_entry(this_entry, region = region, verbose = verbose):
            yield(xyz)
    elif this_entry[1] == 401:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'nos'):
            yield(xyz)
    elif this_entry[1] == 402:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'dc'):
            yield(xyz)
    elif this_entry[1] == 403:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'charts'):
            yield(xyz)
    elif this_entry[1] == 404:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'srtm'):
            yield(xyz)
    elif this_entry[1] == 406:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'mb'):
            yield(xyz)
    elif this_entry[1] == 408:
        for xyz in fetch_module_yield_entry(this_entry, region, verbose, 'gmrt'):
            yield(xyz)

def datalist_yield_queue(q):
    while True:
        this_entry_info = q.get()
        this_entry = this_entry_info[0]
        this_region = this_entry_info[1]
        this_z_region = this_entry_info[2]
        this_verbose = this_entry_info[3]
        this_archive = this_entry_info[4]
        this_wt = this_entry_info[5]
        echo_msg('parsing {}'.format(this_entry[0]))
        dly = datalist_yield_entry(this_entry, this_region, verbose = this_verbose, z_region = this_z_region)
        if this_archive: dly = datalist_archive_yield_entry(this_entry, dirname = 'archive', region = this_region, weight = this_wt, verbose = this_verbose, z_region = this_z_region)
        for xyz in dly:
            yield(xyz)
            
        q.task_done()    

def datalist_yield_xyz_queue(dl, wt = None, pass_h = _dl_pass_h,
                             region = None, archive = False,
                             mask = False, verbose = False, z_region = None):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''

    q = queue.Queue()

    echo_msg('starting 3 parsing threads')
    for _ in range(3):
        t = threading.Thread(target = datalist_yield_queue, args= (q, ))
        t.daemon = True
        t.start()
        
    for this_entry in datalist(dl, wt = wt, pass_h = pass_h, verbose = verbose):
        q.put([this_entry, region, z_region, verbose, archive, wt])
    q.join()
        
def datalist_yield_xyz(dl, wt = None, pass_h = _dl_pass_h,
                       region = None, archive = False,
                       mask = False, verbose = False, z_region = None):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''
    for this_entry in datalist(dl, wt = wt, pass_h = pass_h, verbose = verbose):
        dly = datalist_yield_entry(this_entry, region, verbose = verbose, z_region = z_region)
        if archive: dly = datalist_archive_yield_entry(this_entry, dirname = 'archive', region = region, weight = wt, verbose = verbose, z_region = z_region)
        for xyz in dly:
            yield(xyz)

def datalist_dump_xyz(dl, wt = None,  pass_h = _dl_pass_h,
                      region = None, archive = False, mask = False,
                      verbose = False, dst_port = sys.stdout, z_region = None):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''

    for xyz in datalist_yield_xyz(dl, wt, pass_h, region, archive, mask, verbose, z_region):
        xyz_line(xyz, dst_port, verbose)
                
## ==============================================
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT, GDAL and VDatum for full functionality
## ==============================================
_waffles_grid_info = {
    'datalist': None,
    'datalists': [],
    'region': None,
    'inc': None,
    'name': 'waffles_dem',
    'node': 'pixel',
    'fmt': 'GTiff',
    'extend': 0,
    'extend_proc': 20,
    'weights': None,
    'z_region': None,
    'w_region': None,
    'fltr': None,
    'sample': None,
    'clip': None,
    'chunk': None,
    'epsg': 4326,
    'mod': 'help',
    'mod_args': (),
    'verbose': False,
    'archive': False,
    'spat': False,
    'mask': False,
    'unc': False,
    'gc': config_check()
}

## ==============================================
## the default waffles config dictionary.
## lambda returns dictionary with default waffles
## ==============================================
#waffles_config = lambda: copy.deepcopy(_waffles_grid_info)
waffles_config_copy = lambda wg: copy.deepcopy(wg)

def waffles_config(datalist = None, datalists = [], region = None, inc = None, name = 'waffles_dem',
                   node = 'pixel', fmt = 'GTiff', extend = 0, extend_proc = 20, weights = None,
                   z_region = None, w_region = None, fltr = None, sample = None, clip = None, chunk = None, epsg = 4326,
                   mod = 'help', mod_args = (), verbose = False, archive = False, spat = False, mask = False,
                   unc = False, gc = None):
    wg = waffles_config_copy(_waffles_grid_info)
    wg['datalist'] = datalist
    wg['datalists'] = datalists
    wg['region'] = region
    wg['inc'] = gmt_inc2inc(str(inc))
    wg['name'] = name
    wg['node'] = node
    wg['fmt'] = 'GTiff'
    wg['extend'] = int_or(extend, 0)
    wg['extend_proc'] = int_or(extend_proc, 20)
    wg['weights'] = weights
    wg['z_region'] = z_region
    wg['w_region'] = w_region
    wg['fltr'] = fltr
    wg['sample'] = gmt_inc2inc(str(sample))
    wg['clip'] = clip
    wg['chunk'] = chunk
    wg['epsg'] = int_or(epsg, 4326)
    wg['mod'] = mod
    wg['mod_args'] = mod_args
    wg['verbose'] = verbose
    wg['archive'] = archive
    wg['spat'] = spat
    wg['mask'] = mask
    wg['unc'] = unc
    wg['gc'] = config_check()

    if wg['datalists'] is None:
        if wg['datalist'] is not None:
            wg['datalists'] = [x[0] for x in datalist2py(wg['datalist'])]
        else: wg['datalists'] = None
    
    if wg['datalist'] is None and len(wg['datalists']) > 0:
        wg['datalist'] = datalist_major(wg['datalists'], region = wg['region'], major = '{}_major.datalist'.format(wg['name']))
        
    #if wg['mod'].lower() != 'vdatum' and wg['mod'].lower() != 'coastline':
    if _waffles_modules[wg['mod']][3]:
        if wg['datalist'] is None:
            echo_error_msg('invalid datalist/s entry')
            return(None)

    if wg['region'] is None or not region_valid_p(wg['region']):
        #if wg['mod'].lower() == 'datalists':
        #    wg['region'] = [-180, 180, -90, 90]
        #else:
        echo_error_msg('invalid region {}'.format(wg['region']))
        return(None)

    if wg['inc'] is None: wg['inc'] = (wg['region'][1] - wg['region'][0]) / 500
    
    return(wg)

## ==============================================
## The waffles modules
## { 'module-name': [module-lambda, module description], ... }
## the module lambda should point to a function that takes at least
## the waffles config as its first option (e.g. def mod(wg))
## ==============================================
_waffles_modules = {
    'surface': [lambda args: waffles_gmt_surface(**args), '''SPLINE DEM via GMT surface
    < surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
     :tension=[0-1] - Spline tension.''', 'raster', True],
    'triangulate': [lambda args: waffles_gmt_triangulate(**args), '''TRIANGULATION DEM via GMT triangulate''', 'raster', True],
    'cudem': [lambda args: waffles_cudem(**args), '''Generate a CUDEM Bathy/Topo DEM <beta>''', 'raster', True],
    'update': [lambda args: waffles_update_dem(**args), '''Update a CUDEM DEM with data from datalist <beta>''', 'raster', True],
    'nearest': [lambda args: waffles_nearneighbor(**args), '''NEAREST NEIGHBOR DEM via GMT or gdal_grid
    < nearest:radius=6s:use_gdal=False >
     :radius=[value] - Nearest Neighbor search radius
     :use_gdal=[True/False] - use gdal grid nearest algorithm''', 'raster', True],
    'num': [lambda args: waffles_num(**args), '''Uninterpolated DEM populated by <mode>.
    < num:mode=n >
     :mode=[key] - specify mode of grid population: k (mask), m (mean) or n (num)''', 'raster', True],
    'vdatum': [lambda args: waffles_vdatum(**args), '''VDATUM transformation grid
    < vdatum:ivert=navd88:overt=mhw:region=3:jar=None >
     :ivert=[vdatum] - Input VDatum vertical datum.
     :overt=[vdatum] - Output VDatum vertical datum.
     :region=[0-10] - VDatum region (3 is CONUS).
     :jar=[/path/to/vdatum.jar] - VDatum jar path - (auto-locates by default)''', 'raster', False],
    'mbgrid': [lambda args: waffles_mbgrid(**args), '''Weighted SPLINE DEM via mbgrid
    < mbgrid:tension=35:dist=10/3:use_datalists=False >
     :tension=[0-100] - Spline tension.
     :dist=[value] - MBgrid -C switch (distance to fill nodata with spline)
     :use_datalists=[True/False] - use waffles built-in datalists''', 'raster', True],
    'invdst': [lambda args: waffles_invdst(**args), '''INVERSE DISTANCE DEM via gdal_grid
    < invdst:power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1 >''', 'raster', True],
    'average': [lambda args: waffles_moving_average(**args), '''Moving AVERAGE DEM via gdal_grid
    < average:radius1=0.01:radius2=0.01 >''', 'raster', True],
    'linear': [lambda args: waffles_linear(**args), '''LINEAR DEM via gdal_grid
    < linear:radius=0.01 >''', 'raster', True],
    'spat-meta': [lambda args: waffles_spatial_metadata(**args), '''generate SPATIAL-METADATA''', 'vector', True],
    #'uncertainty': [lambda args: waffles_interpolation_uncertainty(**args), '''generate DEM UNCERTAINTY
    #< uncertainty:mod=surface:dem=None:msk=None:prox=None:slp=None:sims=2 >''', 'raster', False],
    'help': [lambda args: waffles_help(**args), '''display module info''', None, False],
    'coastline': [lambda args: waffles_coastline(**args), '''generate a coastline (landmask)''', 'vector', False],
    'datalists': [lambda args: waffles_datalists(**args), '''recurse the DATALIST
    < datalists:dump=False:echo=False:infos=False:recurse=True >
     :dump=[True/False] - dump the data from the datalist(s)
     :echo=[True/False] - echo the data entries from the datalist(s)
     :infos=[True/False] - generate inf files for the datalists datalist entries.
     :recurse=[True/False] - recurse the datalist (default = True)''', None, True],
}

## ==============================================
## module descriptors (used in cli help)
## ==============================================
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['{:14}{}\n'.format(key, x[key][1]) for key in x]) + '\n'
_waffles_module_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

## ==============================================
## the grid-node region
## ==============================================
waffles_grid_node_region = lambda wg: region_buffer(wg['region'], wg['inc'] * .5)

## ==============================================
## the "proc-region" region_buffer(wg['region'], (wg['inc'] * 20) + (wg['inc'] * wg['extend']))
## ==============================================
waffles_proc_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * wg['extend_proc']) + (wg['inc'] * wg['extend']))
waffles_coast_region = lambda wg: region_buffer(waffles_proc_region(wg), (wg['inc'] * 200))
waffles_proc_str = lambda wg: region_format(waffles_proc_region(wg), 'gmt')
waffles_proc_bbox = lambda wg: region_format(waffles_proc_region(wg), 'bbox')
waffles_proc_ul_lr = lambda wg: region_format(waffles_proc_region(wg), 'ul_lr')

## ==============================================
## the "dist-region" region_buffer(wg['region'], (wg['inc'] * wg['extend']))
## ==============================================
waffles_dist_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * wg['extend']) if wg['sample'] is None else (wg['sample'] * wg['extend']))
waffles_dist_ul_lr = lambda wg: region_format(waffles_dist_region(wg), 'ul_lr')

## ==============================================
## the datalist dump function, to use in run_cmd()
## ==============================================
waffles_dl_func = lambda wg: lambda p: waffles_dump_datalist(wg, dst_port = p)

## ==============================================
## grid registration string for use in GTM programs
## ==============================================
waffles_gmt_reg_str = lambda wg: '-r' if wg['node'] == 'pixel' else ''

## ==============================================
## the 'long-name' used from prefix
## ==============================================
waffles_append_fn = lambda bn, region, inc: '{}{}_{}_{}v1'.format(bn, inc2str_inc(inc), region_format(region, 'fn'), this_year())

def waffles_wg_valid_p(wg = _waffles_grid_info):
    '''return True if wg_config appears valid'''
    
    if wg['datalists'] is None: return(False)
    if wg['mod'] is None: return(False)
    else: return(True)

def waffles_help(wg = _waffles_grid_info):
    sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
    return(0, 0)

def waffles_dlp_hooks(wg = _waffles_grid_info):
    
    region = waffles_proc_region(wg)
    dlp_hooks = datalist_default_hooks()
    
    if region is not None: dlp_hooks.append(lambda e: regions_intersect_ogr_p(region, inf_entry(e)))
    if wg['z_region'] is not None:
        dlp_hooks.append(lambda e: z_region_pass(inf_entry(e), upper_limit = wg['z_region'][1], lower_limit = wg['z_region'][0]))

    if wg['w_region'] is not None:
        dlp_hooks.append(lambda e: z_pass(e[2], upper_limit = wg['w_region'][1], lower_limit = wg['w_region'][0]))

    return(dlp_hooks)

def waffles_yield_datalist(wg = _waffles_grid_info):
    '''recurse the datalist and do things to it.
    
    yield the xyz line data'''
    
    region = waffles_proc_region(wg)
    dlp_hooks = waffles_dlp_hooks(wg)
    
    dly = datalist_yield_xyz(wg['datalist'], pass_h = dlp_hooks, wt = 1 if wg['weights'] else None, region = region, archive = wg['archive'], verbose = wg['verbose'], z_region = wg['z_region'])
    if wg['mask']: dly = gdal_xyz_mask(dly, '{}_msk.tif'.format(wg['name']), region, wg['inc'], dst_format = wg['fmt'])
    for xyz in dly: yield(xyz)
    
    if wg['archive']:
        a_dl = os.path.join('archive', '{}.datalist'.format(wg['name']))

        for dir_, _, files in os.walk('archive'):
            for f in files:
                if '.datalist' in f:
                    rel_dir = os.path.relpath(dir_, 'archive')
                    rel_file = os.path.join(rel_dir, f)
                    datalist_append_entry([rel_file, -1, 1], a_dl)

def waffles_dump_datalist(wg = _waffles_grid_info, dst_port = sys.stdout):
    '''dump the xyz data from datalist.'''

    for xyz in waffles_yield_datalist(wg):
        xyz_line(xyz, dst_port, True)
        
def waffles_datalist_list(wg = _waffles_grid_info):
    '''list the datalist entries in the given region'''
        
    for this_entry in datalist(wg['datalist'], wt = 1, pass_h = waffles_dlp_hooks(wg)):
        print(' '.join([','.join(x) if i == 3 else os.path.abspath(str(x)) if i == 0 else str(x) for i,x in enumerate(this_entry[:-1])]))
        
def waffles_datalists(wg = _waffles_grid_info, dump = False, echo = False, infos = False, recurse = True):
    '''dump the xyz data from datalist and generate a data mask while doing it.'''

    if echo: waffles_datalist_list(wg)
    if infos: print(datalist_inf(wg['datalist'], inf_file = True))
    if dump:
        recurse = True
        pass_func = lambda xyz: xyz_line(xyz, sys.stdout, True)
    else: pass_func = lambda xyz: None

    if recurse:
        for xyz in waffles_yield_datalist(wg): pass_func(xyz)
    return(0,0)

## ==============================================
## Waffles Spatial Metadata
## Polygonize each datalist entry into vector data source
## ==============================================
def waffles_spatial_metadata(wg):
    '''generate spatial-metadata'''

    dst_vector = '{}_sm.shp'.format(wg['name'])
    dst_layer = '{}_sm'.format(wg['name'])
    v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
    t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
    remove_glob('{}.*'.format(dst_layer))
    gdal_prj_file('{}.prj'.format(dst_layer), wg['epsg'])
    
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer('{}'.format(dst_layer), None, ogr.wkbMultiPolygon)
        [layer.CreateField(ogr.FieldDefn('{}'.format(f), t_fields[i])) for i, f in enumerate(v_fields)]
        [layer.SetFeature(feature) for feature in layer]
    else: layer = None
    defn = layer.GetLayerDefn()

    for this_entry in datalist(wg['datalist'], wt = 1 if wg['weights'] else None, pass_h = waffles_dlp_hooks(wg), yield_dl_entry = True, verbose = wg['verbose']):
        if this_entry[1] == -1 or this_entry[-1] == wg['datalist'].split('.')[0]:
            echo_msg(this_entry[0])

            defn = None if layer is None else layer.GetLayerDefn()            
            twg = waffles_config_copy(wg)
            twg['datalist'] = this_entry[0]
            twg['name'] = '{}_{}_msk'.format(os.path.basename(this_entry[0]).split('.')[0], region_format(twg['region'], 'fn'))
            if twg['inc'] < gmt_inc2inc('.3333333s'):
                twg['inc'] = gmt_inc2inc('.3333333s')
                twg['extend'] = twg['extend'] / 3
            twg = waffles_config(**twg)

            if len(this_entry[3]) == 8:
                o_v_fields = this_entry[3]
            else: o_v_fields = [twg['name'], 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

            out, status = waffles_num(twg, mode='k')
            waffles_gdal_md(twg)
            ng = '{}.tif'.format(twg['name'])
            if gdal_infos(ng, True)['zr'][1] == 1:
                tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(twg['name']))
                if tmp_ds is not None:
                    tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(twg['name']), None, ogr.wkbMultiPolygon)
                    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                    gdal_polygonize(ng, tmp_layer, verbose = twg['verbose'])

                    if len(tmp_layer) > 1:
                        if defn is None: defn = tmp_layer.GetLayerDefn()
                        out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                        [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(v_fields)]
                        layer.CreateFeature(out_feat)
                tmp_ds = None
                remove_glob('{}_poly.*'.format(twg['name']))
            remove_glob(ng)
    ds = None
    return(0, 0)
    
## ==============================================
## Waffles CUDEM generation module
## ==============================================
def waffles_cudem(wg = _waffles_grid_info, coastline = None):
    '''generate bathy/topo DEM'''
    
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the CUDEM module')
        return(None, -1)

    ## ==============================================
    ## generate the bathy-surface
    ## using 'surface' with upper_limit of -0.1
    ## at 1 arc-second spacing
    ## ==============================================
    b_wg = waffles_config_copy(wg)
    ul = -0.1
    if coastline is not None:

        # fetch gmrt to sample coastline
        #gmrt_res = fetches.gmrt(extent = waffles_coast_region(b_wg)).run()
        #print(gmrt_res)
        #gmrt = 'gmrt_{}.tif'.format(region_format(waffles_coast_region(b_wg), 'fn'))
        #status = fetches.fetch_file(gmrt_res[0][0], gmrt, verbose = True)
        #print(status)
        #echo_msg('----------------------------')
        #echo_msg('generating coastline surface')
        #echo_msg('----------------------------')
        #c_wg = waffles_config_copy(wg)
        #c_wg['inc'] = gmt_inc2inc('3s')
        #c_wg['name'] = 'coast_{}'.format(wg['name'])
        #c_wg['mod'] = 'surface'
        #c_wg['mod_args'] = ('tension=1',)
        #c_wg = waffles_config(**c_wg)
        #coast_surf = waffles_run(c_wg)
        out, status = run_cmd('coastline2xyz.sh -I {} -O {}_coast.xyz -Z {} -W {} -E {} -S {} -N {}'\
                              .format(coastline, b_wg['name'], 0, waffles_coast_region(b_wg)[0], waffles_coast_region(b_wg)[1], waffles_coast_region(b_wg)[2], waffles_coast_region(b_wg)[3]), verbose = True)
        #remove_glob(coast_surf)
        coast_xyz = '{}_coast.xyz'.format(b_wg['name'])
        coast_region = xyz_inf_entry([coast_xyz])
        print(coast_region)
        #ul = coast_region[5]

    echo_msg('----------------------------')
    echo_msg('generating bathy surface')
    echo_msg('----------------------------')
    b_wg['z_region'] = [None, ul + .5]
    b_wg['name'] = 'bathy_{}'.format(wg['name'])
    b_wg['spat'] = False
    b_wg['fltr'] = None
    b_wg['datalist'] = None
    b_wg['datalists'].append('{} 168 .1'.format(coast_xyz))
    b_wg['mod'] = 'surface'
    b_wg['mod_args'] = ('upper_limit={}'.format(ul),)
    #b_wg['mod'] = 'triangulate'
    b_wg['sample'] = wg['inc']
    b_wg['inc'] = gmt_inc2inc('1s')
    b_wg['name'] = 'bathy_{}'.format(wg['name'])
    b_wg['clip'] = '{}:invert=True'.format(coastline)
    b_wg['extend_proc'] = 40
    b_wg['mask'] = False
    
    b_wg = waffles_config(**b_wg)
    
    bathy_surf = waffles_run(b_wg)
    remove_glob(coast_xyz)

    echo_msg('----------------------------')
    echo_msg('generating integrated bathy-topo surface')
    echo_msg('----------------------------')
    ## ==============================================
    ## append the bathy-surface to the datalist and
    ## generate final DEM using 'surface'
    ## ==============================================
    wg['datalist'] = None
    wg['datalists'].append('{} 200 .5'.format(bathy_surf))
    wg['w_region'] = [.4, None]
    wg = waffles_config(**wg)

    print(wg)
    
    return(waffles_gmt_surface(wg))

## ==============================================
## Waffles Update module <beta>
## ==============================================
def waffles_update_dem(wg = _waffles_grid_info, dem = None):
    ## grid datalist with nearneighbor
    nn_wg = waffles_config_copy(wg)
    nn_wg['mod'] = 'surface'
    nn_wg['mod_args'] = ()
    nn_wg['name'] = 'surf_{}'.format(wg['name'])
    nn_wg['spat'] = True
    nn_wg['mask'] = False
    nn_wg['fltr'] = None
    nn_wg['sample'] = wg['inc']
    nn_wg['inc'] = gmt_inc2inc('1s')
    nn_wg['clip'] = '{}_sm.shp:invert=True'.format(wg['name'])

    nn_dem = waffles_run(nn_wg)
    
    ## mask nearneighbor grid
    dst_msk = 'msk_{}.tif'.format(nn_wg['name'])
    gmt_num_msk(nn_dem, '{}=gd:GTiff'.format(dst_msk), verbose = wg['verbose'])
    gdal_set_nodata('msk_{}.tif'.format(nn_wg['name']), -9999)
    msk_inf = gdal_infos(dst_msk, scan = False)
    print(msk_inf)
    
    dem_inf = gdal_infos(dem, scan = False)
    print(dem_inf)
    nd = dem_inf['ndv']
    
    #dst_msk1 = 'msk_{}.tif'.format(nn_wg['name'])
    #out, status = run_cmd('gmt grdmath {} {} MUL = {}=gd:GTiff'.format(dst_msk, nd, dst_msk1), verbose = True)
    #gdal_set_nodata('msk_{}.tif'.format(nn_wg['name']), -9999)
        
    ### clip DEM to nearneighbor mask
    #msk_dem = 'msk_{}'.format(dem)
    ##out, status = run_cmd('gmt grdmath -N {} {} SUM = {}=gd:GTiff'.format(dem, dst_msk1, msk_dem), verbose = True)
    #out, status = run_cmd('gdal_calc.py -A {} -B {} --calc A+B --outfile {}'.format(dem, dst_msk1, msk_dem), verbose = True)
    
    ## add clipped DEM to datalist
    f_wg = waffles_config_copy(wg)
    f_wg['mod'] = 'surface'
    f_wg['mod_args'] = ()
    f_wg['datalist'] = None
    #f_wg['datalists'] =[]
    f_wg['datalists'].append('{} 200 10'.format(nn_dem))
    f_wg['datalists'].append('{} 200 .25'.format(dem)) #msk_dem))
    
    f_wg = waffles_config(**f_wg)
    
    ## grid datalist with surface
    #dem_update = waffles_run(f_wg)
    dem_update = waffles_gmt_surface(f_wg)
    return(0, 0)
    
## ==============================================
## Waffles MBGrid module
## ==============================================
def waffles_mbgrid(wg = _waffles_grid_info, dist = '10/3', tension = 35, use_datalists = False):
    '''Generate a DEM with MBSystem's mbgrid program.
    if `use_datalists` is True, will parse the datalist through
    waffles instead of mbsystem.'''
    
    if wg['gc']['MBGRID'] is None:
        echo_error_msg('MBSystem must be installed to use the MBGRID module')
        return(None, -1)
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the MBGRID module')
        return(None, -1)

    if use_datalists:
        #datalist_archive(wg, arch_dir = '.mb_tmp_datalist', verbose = True)
        archive = wg['archive']
        wg['archive'] = True
        for xyz in waffles_yield_datalist(wg): pass
        wg['datalist'] = datalist_major(['archive/{}.datalist'.format(wg['name'])])
        wg['archive'] = archive

    region = waffles_proc_region(wg)
    region = region_buffer(region, wg['inc'] * -.5)
    xsize, ysize, gt = gdal_region2gt(waffles_proc_region(wg), wg['inc'])
    
    if len(dist.split('/')) == 1: dist = dist + '/2'
    mbgrid_cmd = ('mbgrid -I{} {} -D{}/{} -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} \
    '.format(wg['datalist'], region_format(region, 'gmt'), xsize, ysize, wg['name'], dist, tension, '-M' if wg['mask'] else ''))
    #for out in yield_cmd(mbgrid_cmd, verbose = wg['verbose']): sys.stderr.write('{}'.format(out))
    out, status = run_cmd(mbgrid_cmd, verbose = wg['verbose'])

    remove_glob('*.cmd')
    remove_glob('*.mb-1')
    gmt_grd2gdal('{}.grd'.format(wg['name']))
    remove_glob('{}.grd'.format(wg['name']))
    if use_datalists and not wg['archive']: shutil.rmtree('archive')
    if wg['mask']:
        remove_glob('*_sd.grd')
        num_grd = '{}_num.grd'.format(wg['name'])
        dst_msk = '{}_msk.tif=gd+n-9999:GTiff'.format(wg['name'])
        out, status = gmt_num_msk(num_grd, dst_msk, verbose = wg['verbose'])
        remove_glob(num_grd)
    if not use_datalists:
        if wg['spat'] or wg['archive']:
            for xyz in waffles_yield_datalist(wg): pass
    return(0, 0)

## ==============================================
## Waffles GMT surface module
## ==============================================
def waffles_gmt_surface(wg = _waffles_grid_info, tension = .35, relaxation = 1.2,
                        lower_limit = 'd', upper_limit = 'd'):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the SURFACE module')
        return(None, -1)

    dem_surf_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} -r\
'.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), \
             wg['inc'], wg['name'], tension, relaxation, lower_limit, upper_limit))
    return(run_cmd(dem_surf_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

## ==============================================
## Waffles GMT triangulate module
## ==============================================
def waffles_gmt_triangulate(wg = _waffles_grid_info):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the TRIANGULATE module')
        return(None, -1)
    dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt triangulate {} -I{:.10f} -V -G{}.tif=gd+n-9999:GTiff -r\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), wg['inc'], wg['name']))
    return(run_cmd(dem_tri_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

## ==============================================
## Waffles nearest neighbor module
## GMT if available else GDAL
## ==============================================
def waffles_nearneighbor(wg = _waffles_grid_info, radius = None, use_gdal = False):
    '''genearte a DEM with GMT nearneighbor or gdal_grid nearest'''
    
    radius = wg['inc'] * 3 if radius is None else gmt_inc2inc(radius)
    if wg['gc']['GMT'] is not None and not use_gdal:
        dem_nn_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), \
                 wg['inc'], radius, wg['name']))
        return(run_cmd(dem_nn_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))
    else: return(waffles_gdal_grid(wg, 'nearest:radius1={}:radius2={}:nodata=-9999'.format(radius, radius)))

## ==============================================
## Waffles 'NUM grid' module
## ==============================================
def waffles_num(wg = _waffles_grid_info, mode = 'n'):
    '''Generate an uninterpolated num grid.
    mode of `k` generates a mask grid
    mode of `m` generates a mean grid
    mode of `n` generates a num grid
    mode of `A<mode> ` generates grid using GMT xyz2grd'''

    if mode[0] == 'A':
        if wg['gc']['GMT'] is None:
            echo_error_msg('GMT must be installed to use the SURFACE module')
            return(None, -1)

        dem_xyz2grd_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt xyz2grd -{} -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', mode, waffles_proc_str(wg), wg['inc'], wg['name']))
        return(run_cmd(dem_xyz2grd_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))
    else:
        dly = waffles_yield_datalist(wg)
        if wg['weights']: dly = xyz_block(dly, waffles_proc_region(wg), wg['inc'], weights = True)
        return(gdal_xyz2gdal(dly, '{}.tif'.format(wg['name']), waffles_proc_region(wg), wg['inc'], dst_format = wg['fmt'], mode = mode, verbose = wg['verbose']))

## ==============================================
## Waffles GDAL_GRID module
## ==============================================
def waffles_gdal_grid(wg = _waffles_grid_info, alg_str = 'linear:radius=1'):
    '''run gdal grid using alg_str
    parse the data through xyz_block to get weighted mean before
    building the GDAL dataset to pass into gdal_grid'''

    region = waffles_proc_region(wg)
    dly = xyz_block(waffles_yield_datalist(wg), region, wg['inc'], weights = False if wg['weights'] is None else True)
    ds = xyz2gdal_ds(dly, '{}'.format(wg['name']))
    if ds.GetLayer().GetFeatureCount() == 0: return(-1,-1)
    xcount, ycount, dst_gt = gdal_region2gt(region, wg['inc'])
    gd_opts = gdal.GridOptions(outputType = gdal.GDT_Float32, noData = -9999, format = 'GTiff', \
                               width = xcount, height = ycount, algorithm = alg_str, callback = _gdal_progress if wg['verbose'] else None, \
                               outputBounds = [region[0], region[3], region[1], region[2]])
    gdal.Grid('{}.tif'.format(wg['name']), ds, options = gd_opts)
    ds = None
    gdal_set_nodata('{}.tif'.format(wg['name']), -9999)
    return(0, 0)

## ==============================================
## Waffles GDAL_GRID invdist module
## ==============================================
def waffles_invdst(wg = _waffles_grid_info, power = 2.0, smoothing = 0.0,
                   radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999):
    '''Generate an inverse distance grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmt_inc2inc(radius2)
    gg_mod = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
                             .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_moving_average(wg = _waffles_grid_info, radius1 = None, radius2 = None,
                           angle = 0.0, min_points = 0, nodata = -9999):
    '''generate a moving average grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmt_inc2inc(radius2)
    gg_mod = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'.format(radius1, radius2, angle, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_linear(wg = _waffles_grid_info, radius = None, nodata = -9999):
    '''generate a moving average grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius is None else gmt_inc2inc(radius)
    gg_mod = 'linear:radius={}:nodata={}'.format(radius, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles VDATUM 'conversion grid' module
## ==============================================
def waffles_vdatum(wg = _waffles_grid_info, ivert = 'navd88', overt = 'mhw',
                   region = '3', jar = None):
    '''generate a 'conversion-grid' with vdatum.
    output will be the differences (surfaced) between 
    `ivert` and `overt` for the region'''
    
    vc = _vd_config
    if jar is None:
        vc['jar'] = vdatum_locate_jar()[0]
    else: vc['jar'] = jar
    vc['ivert'] = ivert
    vc['overt'] = overt
    vc['region'] = region

    gdal_null('empty.tif', waffles_proc_region(wg), 0.00083333, nodata = 0)
    with open('empty.xyz', 'w') as mt_xyz:
        gdal_dump_entry(['empty.tif', 200, 1], dst_port = mt_xyz)
    
    run_vdatum('empty.xyz', vc)
    
    if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
        with open('result/empty.xyz') as infile:
            empty_infos = xyz_inf(infile)
        print(empty_infos)

        ll = 'd' if empty_infos[4] < 0 else '0'
        lu = 'd' if empty_infos[5] > 0 else '0'
        wg['datalists'] = ['result/empty.xyz']
        wg['spat'] = False
        wg = waffles_config(**wg)
        out, status = waffles_gmt_surface(wg, tension = 0, upper_limit = lu, lower_limit = ll)
    else:
        out = None
        status = -1
        
    remove_glob('empty.*')
    remove_glob('result/*')
    remove_glob('.mjr.datalist')
    os.removedirs('result')
    return(out, status)

## ==============================================
## Waffles Interpolation Uncertainty module
## make own module!
## ==============================================
_unc_config = {
    'wg': waffles_config_copy(_waffles_grid_info),
    'mod': 'surface',
    'mod_args': (),
    'dem': None,
    'msk': None,
    'prox': None,
    'slp': None,
    'percentile': 95,
    'zones': ['bathy', 'bathy-topo', 'topo'],
    'sims': 10,
    'chnk_lvl': 6,
}

waffles_unc_config = lambda: copy.deepcopy(_unc_config)
waffles_unc_config_copy = lambda uc: copy.deepcopy(uc)

## TODO: update naming for when module is called directly...
def waffles_interpolation_uncertainty(wg = _waffles_grid_info, mod = 'surface', mod_args = (), \
                                      dem = None, msk = None, prox = None, slp = None, \
                                      percentile = 95, zones = ['bathy', 'bathy-topo', 'topo'], \
                                      sims = 10, chnk_lvl = 4):
    '''calculate the interpolation uncertainty
    - as related to distance to nearest measurement.

    returns [[err, dist] ...]'''

    s_dp = None
    s_ds = None

    # ## ==============================================
    # ## set the module and input grids.
    # ## ==============================================
    if mod not in _waffles_modules.keys():
        echo_error_msg('invalid module name `{}`; reverting to `surface`'.format(mod))
        mod = 'surface'
        mod_args = ()
        
    wg['mod'] = mod
    wg['mod_args'] = mod_args
    
    echo_msg('running INTERPOLATION uncertainty module using {}...'.format(wg['mod']))
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)

    # if dem is None or not os.path.exists(dem):
    #     if dem is None: dem = '{}.tif'.format(wg['name'])
    #     tmp_wg = waffles_config_copy(wg)
    #     if dem is None:
    #         dem = '{}.tif'.format(wg['name'])
    #         tmp_wg['name'] = '_{}'.format(wg['name'])
    #     if msk is None or not os.path.exists(msk):
    #         if msk is None: msk = '{}_msk.tif'.format(tmp_wg['name'])
    #         tmp_wg['mask'] = True
    #     else: tmp_wg['mask'] = False
    #     waffles_run(tmp_wg)

    # if msk is None or not os.path.exists(msk):
    #     if msk is None: msk = '_{}_msk.tif'.format(wg['name'])
    #     tmp_wg = waffles_config_copy(wg)
    #     tmp_wg['name'] = '_{}_msk'.format(wg['name'])
    #     tmp_wg['mod'] = 'num'
    #     tmp_wg['mod_args'] = ('mode=k',)
    #     waffles_run(tmp_wg)
        
    # if prox is None or not os.path.exists(prox):
    #     if prox is None: prox = '_{}_prox.tif'.format(wg['name'])
    #     gdal_proximity(msk, prox)
    # if slp is None or not os.path.exists(slp):
    #     if slp is None: slp = '_{}_slp.tif'.format(wg['name'])
    #     gdal_slope(dem, slp)

    ## ==============================================
    ## region analysis
    ## ==============================================
    region_info = {}
        
    ## mask analysis
    num_sum, g_max, num_perc = gdal_mask_analysis(mask = msk)

    ## proximity analysis
    prox_perc_95 = gdal_percentile(prox, 95)
    prox_perc_90 = gdal_percentile(prox, 90)
    prox_percentile = gdal_percentile(prox, percentile)

    region_info[wg['name']] = [wg['region'], g_max, num_sum, num_perc, prox_percentile] 
    for x in region_info.keys():
        echo_msg('region: {}: {}'.format(x, region_info[x]))

    ## ==============================================
    ## chunk region into sub regions
    ## ==============================================
    echo_msg('chunking region into sub-regions using chunk level {}...'.format(chnk_lvl))
    chnk_inc = int(region_info[wg['name']][4] * chnk_lvl)
    sub_regions = region_chunk(wg['region'], wg['inc'], chnk_inc)
    echo_msg('chunked region into {} sub-regions.'.format(len(sub_regions)))

    ## ==============================================
    ## sub-region analysis
    ## ==============================================
    echo_msg('analyzing {} sub-regions...'.format(len(sub_regions)))
    sub_zones = {}    
    for sc, sub_region in enumerate(sub_regions):
        gdal_cut(msk, sub_region, 'tmp_msk.tif')
        gdal_cut(dem, sub_region, 'tmp_dem.tif')
        s_sum, s_g_max, s_perc = gdal_mask_analysis('tmp_msk.tif')
        s_dc = gdal_infos('tmp_dem.tif', True)
        zone = 'Bathy' if s_dc['zr'][1] < 0 else 'Topo' if s_dc['zr'][0] > 0 else 'BathyTopo'
        sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, s_dc['zr'][0], s_dc['zr'][1], zone]
        remove_glob('tmp_*.tif')
        
    s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
    s_5perc = np.percentile(s_dens, 5)
    s_dens = None
    echo_msg('Sampling density for region is: {:.16f}'.format(s_5perc))

    ## ==============================================
    ## zone analysis / generate training regions
    ## ==============================================
    trainers = []
    bathy_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'Bathy']
    bathy_topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'BathyTopo']
    topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'Topo']

    for z, tile_set in enumerate([bathy_tiles, bathy_topo_tiles, topo_tiles]):
        if len(tile_set) > 0:
            t_dens = np.array([x[3] for x in tile_set])
            t_50perc = np.percentile(t_dens, 50)
        else: t_50perc = 0.0
        echo_msg('Minimum sampling for {} tiles: {}'.format(zones[z].upper(), t_50perc))
        t_trainers = [x for x in tile_set if x[3] > t_50perc]
        echo_msg('possible {} training zones: {}'.format(zones[z].upper(), len(t_trainers)))
        trainers.append(t_trainers)
    trains = regions_sort(trainers)
    echo_msg('sorted training tiles.')
    echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))

    ## ==============================================
    ## split-sample simulations and error calculations
    ## ==============================================
    for sim in range(0, sims):
        echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of {} [{:3}%]'.format(sim + 1, sims, 0))
        status = 0
        for z, train in enumerate(trains):
            train_h = train[:25]
            ss_samp = s_5perc

            ## ==============================================
            ## perform split-sample analysis on each training region.
            ## ==============================================
            for n, sub_region in enumerate(train_h):
                perc = int(float(n+(len(train_h) * z))/(len(train_h)*len(trains)) * 100)
                echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of {} [{:3}%]'.format(sim + 1, sims, perc))
                this_region = sub_region[0]
                if sub_region[3] < ss_samp: ss_samp = None

                ## ==============================================
                ## extract the xyz data for the region from the DEM
                ## ==============================================
                o_xyz = '{}_{}.xyz'.format(wg['name'], n)
                ds = gdal.Open(dem)
                with open(o_xyz, 'w') as o_fh:
                    for xyz in gdal_parse(ds, srcwin = gdal_srcwin(ds, region_buffer(this_region, (20 * wg['inc']))), mask = msk):
                        xyz_line(xyz, o_fh)
                ds = None

                if os.stat(o_xyz).st_size == 0:
                    echo_error_msg('no data in sub-region...')
                else:
                    ## ==============================================
                    ## split the xyz data to inner/outer; outer is
                    ## the data buffer, inner will be randomly sampled
                    ## ==============================================
                    s_inner, s_outer = gmt_select_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = False) #verbose = uc['wg']['verbose'])
                    if os.stat(s_inner).st_size != 0:
                        sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                    else: sub_xyz = []
                    ss_len = len(sub_xyz)
                    if ss_samp is not None:
                        sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
                    else: sx_cnt = 1
                    sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                    np.random.shuffle(sub_xyz)
                    np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                    ## ==============================================
                    ## generate the random-sample DEM
                    ## ==============================================
                    wc = waffles_config(name = 'sub_{}'.format(n),
                                        datalists = [s_outer, sub_xyz_head],
                                        region = this_region,
                                        inc = wg['inc'],
                                        mod = wg['mod'],
                                        verbose = False,
                                        mod_args = wg['mod_args'],
                                        mask = True,
                                        )
                    sub_dem = waffles_run(wc)
                    sub_msk = '{}_msk.tif'.format(wc['name'])
                    
                    if os.path.exists(sub_dem) and os.path.exists(sub_msk):
                        ## ==============================================
                        ## generate the random-sample data PROX and SLOPE
                        ## ==============================================        
                        sub_prox = '{}_prox.tif'.format(wc['name'])
                        gdal_proximity(sub_msk, sub_prox)

                        sub_slp = '{}_slp.tif'.format(wc['name'])
                        gdal_slope(sub_dem, sub_slp)

                        ## ==============================================
                        ## Calculate the random-sample errors
                        ## ==============================================
                        sub_xyd = gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                        #sub_dp = gdal_query(sub_xyd, sub_prox, 'zg')
                        sub_dp = gdal_query(sub_xyd, sub_prox, 'xyzg')
                        sub_ds = gdal_query(sub_dp, slp, 'g')
 
                        if len(sub_dp) > 0:
                            if sub_dp.shape[0] == sub_ds.shape[0]:
                                sub_dp = np.append(sub_dp, sub_ds, 1)
                            else:
                                print(n)
                                print(sub_dp.shape)
                                print(sub_dp)
                                print(sub_ds.shape)
                                print(sub_ds)
                                sys.exit()
                                sub_dp = []
                    else: sub_dp = None
                    remove_glob(sub_xyz_head)

                    if s_dp is not None: 
                        if sub_dp is not None and len(sub_dp) > 0:
                            try:
                                s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                            except: s_dp = sub_dp
                    else: s_dp = sub_dp
                remove_glob(o_xyz)
                remove_glob('sub_{}*'.format(n))
    echo_msg('ran INTERPOLATION uncertainty module using {}.'.format(wg['mod']))

    if len(s_dp) > 0:
        ## ==============================================
        ## save err dist data files
        ## ==============================================
        echo_msg('gathered {} error points'.format(len(s_dp)))

        #np.savetxt('{}.err'.format(uc['wg']['name']), s_dp, '%f', ' ')

        d_max = region_info[wg['name']][4]
        s_dp = s_dp[s_dp[:,3] < d_max,:]
        #s_dp = s_dp[s_dp[:,3] > 0,:]

        prox_err = s_dp[:,[2,3]]
        slp_err = s_dp[:,[2,4]]
        
        np.savetxt('{}_prox.err'.format(wg['name']), prox_err, '%f', ' ')
        #np.savetxt('{}_slp.err'.format(wg['name']), slp_err, '%f', ' ')

        ec_d = err2coeff(prox_err[:50000000], dst_name = wg['name'] + '_prox', xa = 'distance')
        #ec_s = err2coeff(slp_err[:50000000], dst_name = wg['name'] + '_slp', xa = 'slope')

        ## ==============================================
        ## apply error coefficient to full proximity grid
        ## ==============================================
        echo_msg('applying coefficient to proximity grid')
        ## USE numpy/gdal instead
        #run_cmd('gdal_calc.py -A {} --outfile {}.tif --calc {}+({}*(A**{}))'.format(prox, ec_d[0], ec_d[1], ec_d[2], wg['name']), verbose = True)
        math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}.tif=gd+n-9999:GTiff\
        '.format(prox, ec_d[2], ec_d[1], 0, wg['name'])
        run_cmd(math_cmd, verbose = wg['verbose'])
        echo_msg('applied coefficient {} to proximity grid'.format(ec_d))

        #remove_glob('_*.tif')
        
        # math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_slp_unc.tif=gd+n-9999:GTiff\
        # # '.format(slp, ec_s[2], ec_s[1], 0, wg['name'])
        # run_cmd(math_cmd, verbose = wg['verbose'])
        # echo_msg('applied coefficient {} to slope grid'.format(ec_s))
        
    return(ec_d, 0)

def waffles_coastline(wg, want_nhd = True, want_gmrt = False):
    '''Generate a coastline polygon from various sources.'''
    
    w_mask = '{}_w.tif'.format(wg['name'])
    
    if wg['datalist'] is not None:    
        ## ==============================================
        ## wet/dry datalist mask
        ## ==============================================
        dly = waffles_yield_datalist(wg)
        if wg['weights']: dly = xyz_block(dly, waffles_dist_region(wg), wg['inc'], weights = True)
        gdal_xyz2gdal(dly, w_mask, waffles_dist_region(wg), wg['inc'], dst_format = wg['fmt'], mode = 'w', verbose = wg['verbose'])
    else:
        ## ==============================================
        ## burn the region
        ## ==============================================
        region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        run_cmd('gdal_rasterize -tr {} {} -te {} -burn -9999 -a_nodata -9999 -ot Int32 -co COMPRESS=DEFLATE region_buff.shp {}'\
                .format(wg['inc'], wg['inc'], region_format(waffles_dist_region(wg), 'te'), w_mask), verbose = False)

    ## ==============================================
    ## load the wet/dry mask array
    ## ==============================================
    ds = gdal.Open(w_mask)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        region = gdal_gt2region(ds_config)
        dst_gt = ds_config['geoT']        
        coast_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
    else: return(-1, -1)
    remove_glob('{}*'.format(w_mask))

    ## ==============================================
    ## Input coastline shapefile `coastpoly`
    ## ==============================================
    
    ## ==============================================
    ## USGS NHD (HIGH-RES U.S. Only)
    ## ==============================================
    if want_nhd:
        u_mask = '{}_u.tif'.format(wg['name'])
        region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        run_cmd('gdal_rasterize -tr {} {} -te {} -burn -9999 -a_nodata -9999 -ot Int32 -co COMPRESS=DEFLATE region_buff.shp {}'\
                .format(wg['inc'], wg['inc'], region_format(waffles_dist_region(wg), 'te'), u_mask), verbose = False)
        remove_glob('region_buff.*')

        fl = fetches.fetch_infos['tnm'][0](waffles_proc_region(wg), [], None)
        r = fl.run(ds = 6, formats = 'FileGDB', extent = 'HU-4 Subregion')
        #fr = fetches.fetch_results(r, waffles_proc_region(wg), fl._outdir, None)
        #fr.start()
        #fr.join()

        if len(r) > 0:
            r_shp = []
            for result in r:
                if fetches.fetch_file(result[0], os.path.join(result[2], result[1]), verbose = True) == 0:
                    gdb_zip = os.path.join(result[2], result[1])
                    gdb_files = unzip(gdb_zip)
                    gdb, gdb_files = procs_unzip(gdb_zip, ['gdb'])
                    gdb_bn = os.path.basename('.'.join(gdb_zip.split('.')[:-1]))
                    gdb = gdb_bn + '.gdb'

                    run_cmd('ogr2ogr {}_NHDArea.shp {} NHDArea -overwrite'.format(gdb_bn, gdb), verbose = True)
                    r_shp.append('{}_NHDArea.shp'.format(gdb_bn))
                    run_cmd('ogr2ogr {}_NHDPlusBurnWaterBody.shp {} NHDPlusBurnWaterBody -overwrite'.format(gdb_bn, gdb), verbose = True)
                    r_shp.append('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn))
                else: echo_error_msg('unable to fetch {}'.format(result))
                    #except: echo_error_msg('unable to process {}'.format(result))

            [run_cmd('ogr2ogr -skipfailures -update -append nhdArea_merge.shp {}'.format(shp), verbose = True) for shp in r_shp]
            run_cmd('gdal_rasterize -burn 1 nhdArea_merge.shp {}'\
                    .format(u_mask), verbose = True)
            remove_glob(gdb_zip)
            remove_glob('{}*'.format(gdb_bn))
            _clean_zips(gdb_files)
            [remove_glob('{}*'.format(shp[:-3])) for shp in r_shp]
            remove_glob('nhdArea_merge.*')

        ## ==============================================
        ## update wet/dry mask with nhd data
        ## ==============================================
        c_ds = gdal.Open(u_mask)
        for this_xyz in gdal_parse(c_ds):
            xpos, ypos = _geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
            try:
                if coast_array[ypos, xpos] == ds_config['ndv']:
                    if this_xyz[2] == 1: coast_array[ypos, xpos] = 0
            except: pass
        c_ds = None            
        remove_glob('{}*'.format(u_mask))
    
    ## ==============================================
    ## GSHHG/GMRT - Global low-res
    ## ==============================================
    g_mask = '{}_g.tif'.format(wg['name'])
    if wg['gc']['GMT'] is not None and not want_gmrt:
        run_cmd('gmt grdlandmask {} -I{} -r -Df -G{}=gd:GTiff -V -N1/0/1/0/1'.format(region_format(waffles_dist_region(wg), 'gmt'), wg['inc'], g_mask), verbose = True)
    else:
        r = fetches.fetch_infos['gmrt'][0](region_buffer(waffles_dist_region(wg), 5, pct = True), [], None).run()
        gmrt_tif = r[0][1]
        if fetches.fetch_file(r[0][0], gmrt_tif, verbose = True) == 0:
            run_cmd('gdalwarp {} {} -tr {} {} -overwrite'.format(gmrt_tif, g_mask, wg['inc'], wg['inc']), verbose = True)
            remove_glob(gmrt_tif)

    ## ==============================================
    ## update wet/dry mask with gsshg/gmrt data
    ## ==============================================
    c_ds = gdal.Open(g_mask)
    for this_xyz in gdal_parse(c_ds):
        xpos, ypos = _geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
        try:
            if coast_array[ypos, xpos] == ds_config['ndv']:
                if this_xyz[2] == 1: coast_array[ypos, xpos] = 0
                elif this_xyz[2] == 0: coast_array[ypos, xpos] = 1
        except: pass
    c_ds = None
    remove_glob('{}*'.format(g_mask))
    gdal_write(coast_array, '{}.tif'.format(wg['name']), ds_config)
    gdal_write(coast_array, 'test.tif', ds_config)

    #run_cmd('gdal_polygonize.py -8 {}.tif {}.shp'.format(wg['name'], wg['name']), verbose = True)
    
    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}.shp'.format(wg['name']))
    if tmp_ds is not None:
        tmp_layer = tmp_ds.CreateLayer('{}'.format(wg['name']), None, ogr.wkbMultiPolygon)
        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
        gdal_polygonize('{}.tif'.format(wg['name']), tmp_layer, verbose = wg['verbose'])
        
    tmp_ds = None

    return(0, 0)

def ogr_clip(src_ogr, dst_ogr, clip_region = None, dn = "ESRI Shapefile"):
    
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    region2ogr(clip_region, 'tmp_clip.shp')
    c_ds = driver.Open('tmp_clip.shp', 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon)

    layer.Clip(c_layer, dst_layer)
    #ogr.Layer.Clip(layer, c_layer, dst_layer)

    ds = c_ds = dst_ds = None

def ogr_empty_p(src_ogr):

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
    else: return(True)


def waffles_gdal_md(wg, cudem = False):
    '''add metadata to the waffles dem'''
    
    ds = gdal.Open('{}.tif'.format(wg['name']), gdal.GA_Update)
    if ds is not None:
        md = ds.GetMetadata()
        if wg['node'] == 'pixel':
            md['AREA_OR_POINT'] = 'Area'
        else: md['AREA_OR_POINT'] = 'Point'
        md['TIFFTAG_DATETIME'] = '{}'.format(this_date())
        if cudem:
            md['TIFFTAG_COPYRIGHT'] = 'DOC/NOAA/NESDIS/NCEI > National Centers for Environmental Information, NESDIS, NOAA, U.S. Department of Commerce'
            md['TIFFTAG_IMAGEDESCRIPTION'] = 'Topography-Bathymetry; NAVD88'
        ds.SetMetadata(md)
        ds = None
    else: echo_error_msg('failed to set metadata')

def waffles_queue(q):
    while True:
        this_wg = q.get()
        dem = waffles_run(this_wg)

        q.task_done()
    
## ==============================================
## Waffles run waffles module via wg_config()
## ==============================================
def waffles_run(wg = _waffles_grid_info):
    '''generate a DEM using wg dict settings
    see waffles_config() to generate a wg config.

    - runs the waffles module to generate the DEM
    - optionally clips the output to shapefile
    - optionally filters the output
    - optionally resamples the output
    - cuts the output to dist-size
    - reformats the output to final format
    - sets metadata in output

    returns dem-fn'''

    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    
    ## ==============================================
    ## validate and/or set the waffles_config
    ## ==============================================
    if wg is None:
        echo_error_msg('invalid configuration, {}'.format(wg))
        sys.exit(-1)

    args_d = {}
    args_d = args2dict(wg['mod_args'], args_d)
    
    if wg['verbose']:
        echo_msg(wg)
        echo_msg('running module {} with {} [{}]...'.format(wg['mod'], wg['mod_args'], args_d))

    dem = '{}.tif'.format(wg['name'])
    vect = True if _waffles_modules[wg['mod']][2] == 'vector' else False
    if wg['mask']: dem_msk = '{}_msk.tif'.format(wg['name'])
    if vect: dem_vect = '{}.shp'.format(wg['name'])
    wg['region'] = waffles_grid_node_region(wg) if wg['node'] == 'grid' else wg['region']
    
    ## ==============================================
    ## optionally generate the DEM in chunks
    ## ==============================================
    if wg['chunk'] is not None:
        xcount, ycount, dst_gt = gdal_region2gt(wg['region'], wg['inc'])
        s_regions = region_chunk(wg['region'], wg['inc'], (xcount/wg['chunk'])+1)
    else: s_regions = [wg['region']]

    chunks = []
    if wg['mask']: chunks_msk = []
    if vect: chunks_vect = []
    for region in s_regions:
        this_wg = waffles_config_copy(wg)
        this_wg['region'] = region
        #this_wg['region'] = waffles_grid_node_region(this_wg) if this_wg['node'] == 'grid' else this_wg['region']
        this_wg['name'] = 'chunk_{}'.format(region_format(region, 'fn'))
        this_dem = this_wg['name'] + '.tif'
        chunks.append(this_dem)
        if this_wg['mask']:
            this_dem_msk = this_wg['name'] + '_msk.tif'
            chunks_msk.append(this_dem_msk)
        if vect: chunks_vect.append('{}.shp'.format(this_wg['name']))
        args_d['wg'] = this_wg

        ## ==============================================
        ## gererate the DEM (run the module)
        ## ==============================================
        #try:
        out, status = _waffles_modules[this_wg['mod']][0](args_d)
        #except KeyboardInterrupt as e:
        #    echo_error_msg('killed by user, {}'.format(e))
        #    sys.exit(-1)
        #except Exception as e:
        #    echo_error_msg('{}'.format(e))
        #    status = -1

        if status != 0: remove_glob(this_dem)
        if not os.path.exists(this_dem): continue
        gdi = gdal_infos(this_dem, scan = True)
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                remove_glob(this_dem)
                if this_wg['mask']: remove_glob(this_dem_msk)
                continue
        else: continue

        gdal_set_epsg(this_dem, this_wg['epsg'])
        waffles_gdal_md(this_wg)
                
        ## ==============================================
        ## optionally filter the DEM 
        ## ==============================================
        if this_wg['fltr'] is not None:
            if this_wg['verbose']: echo_msg('filtering {}...'.format(this_dem))
            fltr_args = {}
            fltr = this_wg['fltr'].split(':')
            fltr_args['fltr'] = gmt_inc2inc(fltr[0])
            fltr_args['use_gmt'] = True
            fltr_args = args2dict(fltr[1:], fltr_args)        
            if fltr_args['use_gmt']: fltr_args['use_gmt'] = True if this_wg['gc']['GMT'] is not None else False
            try:
                gdal_smooth(this_dem, 'tmp_s.tif', **fltr_args)
                os.rename('tmp_s.tif', this_dem)
            except TypeError as e: echo_error_msg('{}'.format(e))

        ## ==============================================
        ## optionally resample the DEM 
        ## ==============================================
        if this_wg['sample'] is not None:
            if this_wg['verbose']: echo_msg('resampling {}...'.format(this_dem))
            if this_wg['gc']['GMT'] is not None:
                gmt_sample_inc(this_dem, inc = this_wg['sample'], verbose = this_wg['verbose'])
                if this_wg['mask']:
                    if this_wg['verbose']: echo_msg('resampling {}...'.format(this_dem_msk))
                    gmt_sample_inc(this_dem_msk, inc = this_wg['sample'], verbose = this_wg['verbose'])
            else:
                out, status = run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} tmp.tif\
                '.format(inc, inc, src_grd, region_format(waffles_proc_region(this_wg)), verbose = verbose))
                if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))

        gdal_set_epsg(this_dem, this_wg['epsg'])
        
        ## ==============================================
        ## optionally clip the DEM to polygon
        ## ==============================================
        if this_wg['clip'] is not None:
            if this_wg['verbose']: echo_msg('clipping {}...'.format(this_dem))
            clip_args = {}
            cp = this_wg['clip'].split(':')
            clip_args['src_ply'] = cp[0]
            clip_args = args2dict(cp[1:], clip_args)
            gdal_clip(this_dem, **clip_args)
            if this_wg['mask']:
                if this_wg['verbose']: echo_msg('clipping {}...'.format(this_dem_msk))
                gdal_clip(this_dem_msk, **clip_args)

        if not os.path.exists(this_dem): continue
        gdi = gdal_infos(this_dem, scan = True)
        if gdi is not None:
            if np.isnan(gdi['zr'][0]):
                remove_glob(this_dem)
                if this_wg['mask']: remove_glob(this_dem_msk)
                continue
        else: continue
                
        ## ==============================================
        ## cut dem to final size -
        ## region buffered by (inc * extend) or
        ## sample * extend) if sample if specified
        ## ==============================================
        try:
            out = gdal_cut(this_dem, waffles_dist_region(this_wg), 'tmp_cut.tif')
            if out is not None: os.rename('tmp_cut.tif', this_dem)
            if this_wg['mask']:
                out = gdal_cut(this_dem_msk, waffles_dist_region(this_wg), 'tmp_cut.tif')
                if out is not None: os.rename('tmp_cut.tif', this_dem_msk)
        except OSError as e:
            remove_glob('tmp_cut.tif')
            echo_error_msg('cut failed, is the dem open somewhere, {}'.format(e)) 
                
    ## ==============================================
    ## merge the chunks and remove
    ## ==============================================
    if len(chunks) > 1:
        out, status = run_cmd('gdal_merge.py -n -9999 -a_nodata -9999 -ps {} -{} -ul_lr {} -o {} {} -co TILED=YES -co COMPRESS=DEFLATE -co PREDICTOR=3\
        '.format(wg['inc'], wg['inc'], waffles_dist_ul_lr(wg), dem, ' '.join(chunks)), verbose = True)
        ## add option to keep chunks.
        [remove_glob(x) for x in chunks]            
    else:
        if os.path.exists(chunks[0]):
            gdal_gdal2gdal(chunks[0], dst_fmt = wg['fmt'], dst_gdal = dem)
            remove_glob(chunks[0])

    if wg['mask']:
        if len(chunks_msk) > 1:
            out, status = run_cmd('gdal_merge.py -n -9999 -a_nodata -9999 -ps {} -{} -ul_lr {} -o {} {}\
            '.format(wg['inc'], wg['inc'], waffles_dist_ul_lr(wg), dem_msk, ' '.join(chunks_msk)), verbose = True)
            ## add option to keep chunks.
            [remove_glob(x) for x in chunks_msk]
        else:
            if os.path.exists(chunks_msk[0]):
                gdal_gdal2gdal(chunks_msk[0], dst_fmt = wg['fmt'], dst_gdal = dem_msk, co = False)
                remove_glob(chunks_msk[0])

    if vect:
        if len(chunks_vect) > 1:
            out, status = run_cmd('ogrmerge.py {} {}'.format(dem_vect, ' '.join(chunks_vect)))
            [remove_glob('{}*'.format(x.split('.')[0])) for x in chunks_vect]
        else:
            remove_glob('{}*'.format(dem_vect.split('.')[0]))
            #ogr_clip(chunks_vect[0], dem_vect, clip_region = wg['region'], dn = "ESRI Shapefile")
            #out, status = run_cmd('ogr2ogr -clipsrc {} {} {}'.format(region_format(wg['region'], 'te'), dem_vect, chunks_vect[0]), verbose = True)
            out, status = run_cmd('ogr2ogr {} {}'.format(dem_vect, chunks_vect[0]), verbose = True)
            remove_glob('{}*'.format(chunks_vect[0].split('.')[0]))

    if os.path.exists(dem):
        ## ==============================================
        ## convert to final format if not geotiff
        ## ==============================================
        if wg['fmt'] != 'GTiff':
            orig_dem = dem
            if wg['gc']['GMT'] is not None:
                dem = gmt_grd2gdal(orig_dem, wg['fmt'])
            else: dem = gdal_gdal2gdal(orig_dem, wg['fmt'])
            remove_glob(orig_dem)

        ## ==============================================
        ## set the projection and other metadata
        ## ==============================================
        gdal_set_epsg(dem, wg['epsg'])
        waffles_gdal_md(wg, True if wg['mod'] == 'cudem' else False)
        if wg['mask']: gdal_set_epsg(dem_msk, wg['epsg'])

    ## ==============================================
    ## optionally generate uncertainty grid
    ## ==============================================
    if wg['unc'] and not vect:
        try:
            if os.path.exists(dem) and os.path.exists(dem_msk):
                echo_msg('generating uncertainty')

                uc = _unc_config
                uc['wg'] = wg
                uc['mod'] = wg['mod']
                uc['mod_args'] = wg['mod_args']
                uc['dem'] = dem
                uc['msk'] = dem_msk

                dem_prox = '{}_prox.tif'.format(wg['name'])
                gdal_proximity(dem_msk, dem_prox)
                uc['prox'] = dem_prox

                dem_slp = '{}_slp.tif'.format(wg['name'])
                gdal_slope(dem, dem_slp)
                uc['slp'] = dem_slp

                echo_msg(uc)
                waffles_interpolation_uncertainty(**uc)
        except Exception as e:
            echo_error_msg('failed to calculate uncertainty, {}'.format(e))

    ## ==============================================
    ## if dem has data, return
    ## ==============================================
    remove_glob('waffles_dem_mjr.datalist')
    remove_glob(wg['datalist'])

    if vect:
        if ogr_empty_p('{}.shp'.format(wg['name'])):
            remove_glob('{}*'.format(wg['name']))
            return(None)
        
    return(dem)

## ==============================================
## Command-line Interfaces (CLIs)
## $ waffles
## $ datalists
## $ regions
## ==============================================

## ==============================================
## waffles cli
## ==============================================
waffles_cli_usage = '''waffles [OPTIONS] <datalist/entry>

Generate DEMs and derivatives and process datalists.

General Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
\t\t\tappend :<inc> to resample the output to the given <inc>: -E.3333333s:.1111111s
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired DEM MODULE and options. (see available Modules below)
\t\t\tsyntax is -M module:mod_opt=mod_val:mod_opt1=mod_val1:...
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
\t\t\tappend :<num> to extend the processing region: -X6:12
  -T, --filter\t\tFILTER the output using a Cosine Arch Filter at -T<dist(km)> search distance.
\t\t\tIf GMT is not available, or if :use_gmt=False, perform a Gaussian filter at -T<factor>. 
\t\t\tAppend :split_value=<num> to only filter values below <num>.
\t\t\te.g. -T10:split_value=0:use_gmt=False to smooth bathymetry using Gaussian filter
  -Z --z-region\t\tRestrict data processing to records that fall within the z-region
\t\t\tUse '-' to indicate no bounding range; e.g. -Z-/0 will restrict processing to data
\t\t\trecords whose z value is below zero.
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk\t\tProcess the region in CHUNKs. -K<chunk-level>
  -W, --w-region\tRestrict data processing to records that fall within the w-region (weight).
\t\t\tUse '-' to indicate no bounding range; e.g. -W1/- will restrict processing to data
\t\t\trecords whose weight value is at least 1.
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tArchive the datalist to the given region.
  -m, --mask\t\tGenerate a data mask raster.
  -u, --uncert\t\tGenerate an associated uncertainty grid.

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, while an entry is a space-delineated line:
  `/path/to/data format weight data,meta,data`

Supported datalist formats: 
  {}

Modules (see waffles --modules for more info):
  {}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_known_datalist_fmts_short_desc(), _waffles_module_short_desc(_waffles_modules))

def waffles_cli(argv = sys.argv):
    '''run waffles from command-line
    e.g. `python waffles.py` 
    generates a waffles_config from the command-line options
    and either outputs the or runs the waffles_config
    on each region supplied (multiple regions can be supplied
    by using a vector file as the -R option.)
    See `waffles_cli_usage` for full cli options.'''
    
    wg = waffles_config_copy(_waffles_grid_info)
    wg_user = None
    dls = []
    region = None
    module = None
    want_prefix = False
    want_verbose = False
    want_config = False
    want_threads = False
    status = 0
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            region = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-R': region = str(arg[2:])
        elif arg == '--module' or arg == '-M':
            module = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-M': module = str(arg[2:])
        elif arg == '--increment' or arg == '-E':
            incs = argv[i + 1].split(':')
            wg['inc'] = gmt_inc2inc(incs[0])
            if len(incs) > 1: wg['sample'] = gmt_inc2inc(incs[1])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            wg['inc'] = gmt_inc2inc(arg[2:].split(':')[0])
            if len(incs) > 1: wg['sample'] = gmt_inc2inc(incs[1])
        elif arg == '--outname' or arg == '-O':
            wg['name'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': wg['name'] = arg[2:]
        elif arg == '--format' or arg == '-F':
            wg['fmt'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': wg['fmt'] = arg[2:]
        elif arg == '--filter' or arg == '-T':
            wg['fltr'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-T': wg['fltr'] = arg[2:]
        elif arg == '--extend' or arg == '-X':
            exts = argv[i + 1].split(':')
            wg['extend'] = int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = int_or(exts[1], 10)
            i += 1
        elif arg[:2] == '-X':
            exts = arg[2:].split(':')
            wg['extend'] = int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = int_or(exts[1], 10)
        elif arg == '--wg-config' or arg == '-G':
            wg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-G': wg_user = arg[2:]
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--chunk' or arg == '-K':
            wg['chunk'] = int_or(argv[i + 1], None)
            i = i + 1
        elif arg[:2] == '-K': wg['chunk'] = int_or(arg[2:], None)
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = int_or(argv[i + 1], 4326)
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = int_or(arg[2:], 4326)
        elif arg == '--z-range' or arg == '-Z':
            zr = argv[i + 1].split('/')
            if len(zr) > 1:
                wg['z_region'] = [None if x == '-' else float(x) for x in zr]
            i = i + 1
        elif arg[:2] == '-Z':
            zr = arg[2:].split('/')
            if len(zr) > 1:
                wg['z_region'] = [None if x == '-' else float(x) for x in zr]
        elif arg == '--w-range' or arg == '-W':
            wr = argv[i + 1].split('/')
            if len(wr) > 1:
                wg['w_region'] = [None if x == '-' else float(x) for x in wr]
            i = i + 1
        elif arg[:2] == '-W':
            wr = arg[2:].split('/')
            if len(wr) > 1:
                wg['w_region'] = [None if x == '-' else float(x) for x in wr]
        elif arg == '-w' or arg == '--weights': wg['weights'] = True
        elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-m' or arg == '--mask': wg['mask'] = True
        elif arg == '-u' or arg == '--uncert':
            wg['mask'] = True
            wg['unc'] = True
        #elif arg == '-s' or arg == 'spat-meta': wg['spat'] = True
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'
        elif arg == '--verbose' or arg == '-V': wg['verbose'] = True
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
            sys.exit(0)
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(_version))
            sys.exit(0)
        else: dls.append(arg)
        i += 1

    ## ==============================================
    ## load the user wg json and run waffles with that.
    ## ==============================================
    if wg_user is not None:
        if os.path.exists(wg_user):
            try:
                with open(wg_user, 'r') as wgj:
                    wg = json.load(wgj)
                    wg = waffles_config(**wg)
                    dem = waffles_run(wg)
                    sys.exit(0)
            except Exception as e:
                wg = waffles_config_copy(wg)
                echo_error_msg(e)
        else:
            echo_error_msg('specified json file does not exist, {}'.format(wg_user))
            sys.exit(0)

    ## ==============================================
    ## Otherwise run from cli options...
    ## set the dem module
    ## ==============================================
    if module is not None:
        mod_opts = {}
        opts = module.split(':')
        if opts[0] in _waffles_modules.keys():
            mod_opts[opts[0]] = list(opts[1:])
        else: echo_error_msg('invalid module name `{}`'.format(opts[0]))
        
        for key in mod_opts.keys():
            mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]
        mod = opts[0]
        mod_args = tuple(mod_opts[mod])
        wg['mod'] = mod
        wg['mod_args'] = mod_args
    else: 
        sys.stderr.write(waffles_cli_usage)
        echo_error_msg('''must specify a waffles -M module.''')
        sys.exit(-1)
        
    #if wg['mod'] != 'vdatum' and wg['mod'] != 'coastline':
    if _waffles_modules[wg['mod']][3]:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
            sys.exit(-1)
            
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['datalists'] = dls
    
    ## ==============================================
    ## reformat and set the region
    ## ==============================================
    if region is not None:
        try:
            these_regions = [[float(x) for x in region.split('/')]]
        except ValueError: these_regions = gdal_ogr_regions(region)
        except Exception as e:
            echo_error_msg('failed to parse region(s), {}'.format(e))
    else: these_regions = [None]
    if len(these_regions) == 0: echo_error_msg('failed to parse region(s), {}'.format(region))
    
    ## ==============================================
    ## run waffles for each input region.
    ## ==============================================
    these_wgs = []
    for this_region in these_regions:
        twg = waffles_config_copy(wg)
        twg['region'] = this_region
        if want_prefix or len(these_regions) > 1:
            twg['name'] = waffles_append_fn(wg['name'], twg['region'], twg['sample'] if twg['sample'] is not None else twg['inc'])
            
        twg = waffles_config(**twg)
        if want_config:
            this_wg = waffles_config_copy(twg)
            if this_wg is not None:
                #echo_msg(json.dumps(this_wg, indent = 4, sort_keys = True))
                echo_msg(this_wg)
                with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                    echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                    echo_msg('generating major datalist: {}_mjr.datalist'.format(this_wg['name']))
                    wg_json.write(json.dumps(this_wg, indent = 4, sort_keys = True))
            else: echo_error_msg('could not parse config.')
        else:
            if want_threads:
                these_wgs.append(twg)
            else: dem = waffles_run(twg)
    if want_threads:
        wq = queue.Queue()
        num_threads = 2 if len(these_wgs) > 1 else len(these_wgs)
        for _ in range(num_threads):
            t = threading.Thread(target = waffles_queue, args = (wq, ))
            t.daemon = True
            t.start()

        [wq.put(x) for x in these_wgs]
        while True:
            #while not wq.empty():
            time.sleep(2)
            echo_msg_inline('generating dem(s) [{}/{}]'.format((len(these_wgs) - wq.qsize()) - num_threads,len(these_wgs)))
            #if wq.qsize == 0:
            if (len(these_wgs) - wq.qsize()) - num_threads == len(these_wgs): break
            #if wq.empty(): break
        echo_msg('queue complete')
        wq.join()

## ==============================================
## dadtalists cli
## ==============================================    
datalists_version = '0.0.1'
datalists_usage = '''{} ({}): Process and generate datalists

usage: {} [ -aghirvFPR [ args ] ] DATALIST ...

Options:
  -R, --region\t\tSpecifies the desired REGION;
  -P, --epsg\t\tSpecify the EPSG of the DATALIST.
  -F, --format\t\tOnly process FORMAT data type.

  -a, --archive\t\tARCHIVE the data from the DATALIST.
  -g, --glob\t\tGLOB FORMAT data into the DATALIST.
  -i, --info-file\tGenerate INF files for the data in the DATALIST
  -l, --list\t\tLIST the datafiles from the DATALIST.
  -m, --mask\t\tMASK the datafiles from the DATALIST.
  -r, --region-info\tReturn the full REGION of the DATALIST.
  -w, --weights\t\toutput weights along with each datalist.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} my_data.datalist -R -90/-89/30/31 -g -i

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            datalists_version, 
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]))

def datalists_cli(argv = sys.argv):

    status = 0
    dls = []
    i_region = None
    i_inc = 0.000277777
    o_bn = None
    epsg = '4326'
    these_regions = []
    want_verbose = False
    want_inf = False
    want_region = False
    want_sm = False
    want_glob = False
    want_list = False
    want_archive = False
    want_mask = False
    want_weights = False
    z_region = None
    w_region = None
    dl_fmt = None
    #dl_fmts = []
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================

    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            i_region = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            i_region = str(arg[2:])

        elif arg == '--output-name' or arg == '-O':
            o_bn = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-O':
            o_bn = str(arg[2:])

        elif arg == '--increment' or arg == '-E':
            try:
                i_inc = float(argv[i + 1])
            except:
                sys.stederr.write('error, -E should be a float value\n')
                sys.exit(1)
            i = i + 1
        elif arg[:2] == '-E':
            try:
                i_inc = float(arg[2:])
            except:
                sys.stederr.write('error, -E should be a float value\n')
                sys.exit(1)

        elif arg == '--epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-P':
            epsg = arg[2:]

        elif arg == '--format' or arg == '-F':
            dl_fmt = int(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-F':
            dl_fmt = int(arg[2:])

        elif arg == '--glob' or arg == '-g':
            want_glob = True

        elif arg == '--spatial-md' or arg == '-s':
            want_sm = True

        elif arg == '--list' or arg == '-l':
            want_list = True

        elif arg == '--mask' or arg == '-m':
            want_mask = True

        elif arg == '--archive' or arg == '-a':
            want_archive = True

        elif arg == '--info-file' or arg == '-i':
            want_inf = True

        elif arg == '--region-info' or arg == '-r':
            want_region = True

        elif arg == '--weights' or arg == '-w':
            want_weights = True

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), datalists_version))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        else:
            dls.append(arg)

        i = i + 1

    if want_glob:
        if dl_fmt is None:
            dl_fmts = list(_known_datalist_fmts.keys())[1:]
        else: dl_fmts = [dl_fmt]
        for key in dl_fmts:
            for f in _known_datalist_fmts[key]:
                globs = glob.glob('*.{}'.format(f))
                [sys.stdout.write('{}\n'.format(' '.join([x, str(key), '1']))) for x in globs]
        sys.exit(0)
        
    if len(dls) == 0:
        print(datalists_usage)
        sys.exit(1)
        
    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if i_region is not None:
        try:
            these_regions = [[float(x) for x in i_region.split('/')]]
        except ValueError: these_regions = gdal_ogr_regions(i_region)
        except Exception as e:
            echo_error_msg('failed to parse region(s), {}'.format(e))
    else: these_regions = [[-180,180,-90,90]]
    if len(these_regions) == 0: echo_error_msg('failed to parse region(s), {}'.format(i_region))

    for rn, this_region in enumerate(these_regions):
        ## ==============================================
        ## Load the input datalist
        ## ==============================================
        dl_m = datalist_major(dls, region = this_region)
        echo_msg('processed datalist')
        dlp_hooks = datalist_default_hooks()
        dlp_hooks.append(lambda e: regions_intersect_ogr_p(this_region, inf_entry(e)))

        if z_region is not None:
            dlp_hooks.append(lambda e: z_region_pass(inf_entry(e), upper_limit = z_region[1], lower_limit = z_region[0]))
        if w_region is not None:
            dlp_hooks.append(lambda e: z_pass(e[2], upper_limit = w_region[1], lower_limit = w_region[0]))
        
        if want_inf:
            datalist_inf_entry([dl_m, -1, 1])
        elif want_region:
            print(datalist_inf(dl_m, False, False))
        elif want_list:
            for this_entry in datalist(dl_m, wt = 1, pass_h = dlp_hooks):
                print(' '.join([','.join(x) if i == 3 else os.path.abspath(str(x)) if i == 0 else str(x) for i,x in enumerate(this_entry[:-1])]))
        elif want_mask: pass
        else:
            datalist_dump_xyz(dl_m, wt = 1 if want_weights else None, pass_h = dlp_hooks, region = this_region)
        remove_glob(dl_m)
        
## ==============================================
## mainline -- run waffles directly...
##
## run waffles:
## % python waffles.py dem <args>
## % python waffles.py <args>
##
## run datalists:
## % python waffles.py datalists <args>
## ==============================================
if __name__ == '__main__': waffles_cli(sys.argv)

### End
