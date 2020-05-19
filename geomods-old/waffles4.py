### waffles.py
##
## Copyright (c) 2013 - 2020 CIRES Coastal DEM Team
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
## Generate Digital Elevation Models and derivatives using a variety of algorithms, etc.
##
## GDAL and gdal-python are required to run waffles.
##
## Recommended external software for full functionality:
## - GMT (FOSS)
## - MBSystem (FOSS)
## - VDatum (US/NOAA)
## - LASTools (Non-Free)
##
### TODO:
##
## switch waffles modules to -A in console and allow multiple input datalists/files
## -- want to be able to run: waffles -A surface *.xyz:168 or some such
## update datalists
## add source uncertainty to uncertainty
## add iso metadata generation to metadata
## add lasf.py for lidar processing
## update confgs for python3, etc.
##
### Code:

import sys
import os

import time
import datetime
import glob
import math
import threading
import subprocess

import numpy as np
import gdal
import ogr
import osr
from gdalconst import *

try:
    import ConfigParser as configparser
except: import configparser

_version = '0.4.0'

_license = '''waffles, version {}
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
## General utility functions - utils.py
##
## =============================================================================

def _con_dec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''

    if x is None:
        return(x)

    fstr = "%." + str(dec) + "f"
    return(fstr % x)

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)'''
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_year():
    '''return the current year'''
    
    return(datetime.datetime.now().strftime('%Y'))

def hav_dst(pnt0, pnt1):
    '''return the distance between pnt0 and pnt1,
    using the haversine formula.'''
    
    x0 = float(pnt0[0])
    y0 = float(pnt0[1])
    x1 = float(pnt1[0])
    y1 = float(pnt1[1])

    ## ==============================================
    ## radians in meters
    ## ==============================================
    rad = 637100
    
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    
    a = math.sin(dx / 2) * math.sin(dx / 2) + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy / 2) * math.sin(dy / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    return(rad * c)

def remove_glob(glob_str):
    '''glob `glob_str` and os.remove results'''

    globs = glob.glob(glob_str)
    if len(globs) > 0:
        for g in globs:
            try:
                os.remove(g)
            except: pass
        return(0)
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

def cmd_check(cmd_str, cmd_vers_str, verbose = False):
    '''check system for availability of 'cmd_str' and return it's version'''
    
    cmd_vers = None
    if verbose: echo_msg('checking for {}...'.format(cmd_str))
    
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str))
        cmd_vers = cmd_vers.rstrip()
    else: cmd_vers = None
    
    if verbose: echo_msg('found {} version {}'.format(cmd_str, cmd_vers))

    return(cmd_vers)

def echo_error_msg(msg):
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('waffles: error, {}\n'.format(msg))

def echo_msg(msg):
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('waffles: {}\n'.format(msg))

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

#CONFIG_FILE = os.path.expanduser('~/waffles.ini')
#_waff_co = configparser.ConfigParser()
#_waff_cf = os.path.expanduser(CONFIG_FILE)

def check_config3(chk_vdatum = False, verbose = False):

    _waff_co = {}
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform

    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers
    
    if host_os == 'win32': ae = '.exe'
    else: ae = ''

    if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal_grid --version')
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version')
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version | grep Version')
    _waff_co['BOUNDS'] = cmd_check('bounds{}'.format(ae), 'bounds --version')
    
    return(_waff_co)
    
## =============================================================================
##
## Region functions and class - regions.py
##
## a region is a geographic bounding-box with 4 corners and the
## region object made from the region_string 'east/west/south/north'
##
## =============================================================================

def regions_intersect_p_depr(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''    
    return(regions_reduce(region_a, region_b)._valid)

def regions_intersect_p(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''
        
    geom_a = _extent2geom(region_a.region)
    geom_b = _extent2geom(region_b.region)

    if geom_a.Intersects(geom_b):
        return(True)
    else: return(False)

def regions_reduce(region_a, region_b):
    '''return the minimum region when combining
    region_a and region_b'''

    region_c = [0, 0, 0, 0]
    if region_a.west > region_b.west: region_c[0] = region_a.west
    else: region_c[0] = region_b.west
    
    if region_a.east < region_b.east:region_c[1] = region_a.east
    else: region_c[1] = region_b.east
    
    if region_a.south > region_b.south: region_c[2] = region_a.south
    else: region_c[2] = region_b.south
    
    if region_a.north < region_b.south: region_c[3] = region_a.north
    else: region_c[3] = region_b.north
    
    return(region(region_c))
    
def regions_merge(region_a, region_b):
    '''merge two regions into a single region'''

    region_c = [0, 0, 0, 0]
    if region_a.west < region_b.west: region_c[0] = region_a.west
    else: region_c[0] = region_b.west
    
    if region_a.east > region_b.east: region_c[1] = region_a.east
    else: region_c[1] = region_b.east
    
    if region_a.south < region_b.south: region_c[2] = region_a.south
    else: region_c[2] = region_b.south
    
    if region_a.north > region_b.north: region_c[3] = region_a.north
    else: region_c[3] = region_b.north
    
    return(region(region_c))
    
class region:
    '''geographic bounding box regions 'w/e/s/n' or [w, e, s, n]'''

    def __init__(self, extent):
        try:
            self.region_string = extent
            self.region = [float(x) for x in extent.split('/')]
        except:
            self.region_string = '/'.join([str(x) for x in extent])
            self.region = extent
            
        self._reset()

    def _reset(self):
        self.west = self.region[0]
        self.east = self.region[1]
        self.south = self.region[2]
        self.north = self.region[3]
        self._format_gmt()
        self._format_bbox()
        self._format_fn()
        self._valid = self._valid_p()        

    def _valid_p(self):
        '''validate region'''

        if self.west < self.east and self.south < self.north: return(True)
        else: return(False)

    def _format_gmt(self):
        '''format region to GMT string'''
        self.gmt = '-R' + '/'.join([str(x) for x in self.region])

    def _format_bbox(self):
        '''format region to bbox string'''

        self.bbox = ','.join([str(self.west), str(self.south), str(self.east), str(self.north)])

    def _format_fn(self):
        '''format region to filename string'''

        if self.north < 0: ns = 's'
        else: ns = 'n'
        if self.west > 0: ew = 'e'
        else: ew = 'w'
        self.fn = ('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(self.north)), abs(int(self.north * 100) % 100), 
                                                            ew, abs(int(self.west)), abs(int(self.west * 100) % 100)))

    def gdal2region(self):
        '''extract the region from a GDAL file.'''

        pass

    def buffer(self, bv, percentage = False):
        '''buffer region'''

        if percentage: bv = self.pct(bv)
        return(region([self.region[0] - bv, self.region[1] + bv, self.region[2] - bv, self.region[3] + bv]))

    def center(self):
        '''return the center point of the region'''
        
        xc = self.west + (self.east - self.west / 2)
        yc = self.south + (self.north - self.south / 2)
        
        return([xc, yc])
    
    def chunk(self, inc, n_chunk = 10):
        '''chunk the region into n_chunk by n_chunk cell regions, given inc.'''

        i_chunk = 0
        x_i_chunk = 0
        x_chunk = n_chunk
        o_chunks = []
        region_x_size = math.floor((self.east - self.west) / inc)
        region_y_size = math.floor((self.north - self.south) / inc)

        while True:
            y_chunk = n_chunk

            while True:
                this_x_chunk = x_chunk
                this_y_chunk = y_chunk
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = this_x_chunk - this_x_origin
                this_y_size = this_y_chunk - this_y_origin
                geo_x_o = self.west + this_x_origin * inc
                geo_x_t = geo_x_o + this_x_size * inc
                geo_y_o = self.south + this_y_origin * inc
                geo_y_t = geo_y_o + this_y_size * inc

                if geo_y_t > self.north: geo_y_t = self.north
                if geo_y_o < self.south: geo_y_o = self.south
                if geo_x_t > self.east: geo_x_t = self.east
                if geo_x_o > self.east: geo_x_o = self.west
                
                this_region = region([geo_x_o, geo_x_t, geo_y_o, geo_y_t])
                o_chunks.append(this_region)

                if y_chunk >= region_y_size:
                    break
                else: 
                    y_chunk += n_chunk
                    i_chunk += 1

            if x_chunk >= region_x_size:
                break
            else:
                x_chunk += n_chunk
                x_i_chunk += 1

        return(o_chunks)

    def pct(self, pctv):
        ewp = (self.east - self.west) * (pctv * .01)
        nsp = (self.north - self.south) * (pctv * .01)

        return((ewp + nsp) / 2)

## =============================================================================
##
## Datalist Class and functions - datalists.py
##
## MBSystem style datalists.
## Recurse through a datalist file and process the results.
##
## a datalist '*.datalits' file should be formatted as in MBSystem:
## ~path ~format ~weight
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
## a format of 200 represents GDAL data
## a format of 300 represents LAS/LAZ data
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
##
## if 'i_region' is specified, will only process data that falls within
## the given region
##
## _dl_fmt: [dump_func, inf_func, [known, extentions]]
##
## todo: make an 'entry' subclass of the datalist subclass
##       to hold a datalist entry object and processing functions.
##
## =============================================================================

_known_delims = [',', ' ', '\t', '/', ':']

_dl_fmts = { 168: [lambda e, d, r, l, v: dump_168(e, d, r, l, v), lambda p: gmt_inf(p), ['xyz', 'dat']],
             200: [lambda e, d, r, l, v: dump_200(e, d, r, l, v), lambda p: gdal_inf(p), ['tif', 'img', 'asc', 'grd']],
             300: [lambda e, d, r, l, v: dump_200(e, d, r, l, v), lambda p: las_inf_lastools(p), ['las', 'laz']],
}

def xyz_parse(src_xyz):
    '''xyz file parsing generator'''
    
    for xyz in src_xyz:
        this_line = xyz.strip()

        for delim in _known_delims:
            this_xyz = this_line.split(delim)
            if len(this_xyz) > 1:
                this_delim = delim
                break

        yield [float(x) for x in this_xyz]

def xyz_line(line, dst_port = sys.stdout, delimiter = ' ', weight = None):
    '''write XYZ `line` to `dst_port` using `delimiter` and `weight`'''
    
    #try:
    if weight is None:
        w_string = ''
    else: w_string = '{}{}'.format(delimiter, weight)

    l = '{}{}\n'.format(delimiter.join([str(x) for x in line]), w_string).encode('utf-8')
    dst_port.write(l)
    #except: sys.exit(1)
        
def xyz_region(src_xyz):
    '''return a region based on a src_xyz file using GMT...'''
    
    out, status = run_cmd('gmt gmtinfo {} -I-'.format(src_xyz), verbose = False)
    o_region = region(out[2:])
    return(o_region)

def xy_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n], else False.'''

    x = src_xy[0]
    y = src_xy[1]

    if x < src_region[0]: return(False)
    elif x > src_region[1]: return(False)
    elif y < src_region[2]: return(False)
    elif y > src_region[3]: return(False)
    else: return(True)

def xyz_inf(src_xyz):
    '''return minmax info from a src_xyz file.'''
    
    minmax = []
    with open(src_xyz, 'r') as in_file:
        for i,l in enumerate(xyz_parse(in_file)):
            if i == 0:
                minmax[0] = l[0], minmax[1] = l[0]
                minmax[2] = l[1], minmax[3] = l[1]
                minmax[4] = l[2], minmax[5] = l[2]
            else:
                if l[0] < minmax[0]: minmax[0] = l[0]
                elif l[0] > minmax[1]: minmax[1] = l[0]
                if l[1] < minmax[2]: minmax[2] = l[1]
                elif l[1] > minmax[3]: minmax[3] = l[1]
                if l[1] < minmax[2]: minmax[4] = l[2]
                elif l[1] > minmax[3]: minmax[5] = l[2]
                
    return(minmax)

def generate_inf(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    
    this_region = parse_or_gen_inf(data_path, data_fmt)
    return(this_region)

def parse_inf(inf_file):
    minmax = [0, 0, 0, 0]
    
    with open(inf_file) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                try: ## GMT/GeoMods inf
                    minmax = [float(x) for x in til]
                except: ## mbsystem inf
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            minmax[0] = til[2]
                            minmax[1] = til[5]
                        elif til[1] == 'Latitude:':
                            minmax[2] = til[2]
                            minmax[3] = til[5]
    return(minmax)

def parse_or_gen_inf(data_path, data_fmt = 168):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file
    or the result of `gmt gmtinfo file.xyz -C`, which is
    a 6 column line with minmax info, etc.
    returns the region of the inf file.'''
    
    minmax = None
    if os.path.exists(data_path):
        path_i = data_path + '.inf'
        if not os.path.exists(path_i):
            _dl_fmts[data_fmt][1](data_path)

        if os.path.exists(path_i): minmax = parse_inf(path_i)[:4]
    try: 
        o_region = region(minmax)
    except: o_region = None

    if o_region is not None:
        if 0. in [float(x) for x in minmax]: o_region = None

    return(o_region)

def datalist_set_weight(data_e, weight = None):
    if weight is not None:
        try:
            dweight = float(weight)
        except: dweight = 1
    else:
        if len(data_e) > 2:
            try:
                dweight = float(data_e[2])
            except: dweight = 1
        else: dweight = 1

    return(dweight)

## TODO: add weighted mean...
def block_entry(entry, region = None, increment = 1, ptArray = [], sumArray = [], dst_gt = None, verbose = False):

    if verbose: echo_msg('using data file {}'.format(entry[0]))        
    if entry[1] == 168:
        with open(entry[0]) as infile:
            for this_xyz in xyz_parse(infile):
                if region is not None:
                    if xy_in_region_p(this_xyz, region.region):
                        xpos, ypos = _geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
                        ptArray[ypos, xpos] += 1
                        sumArray[ypos, xpos] += this_xyz[2]

    elif entry[1] == 200:
        if region is not None:
            srcwin = _srcwin(entry[0], region.region)
        else: srcwin = None
        for this_xyz in gdal_yield(entry[0], dump_nodata = False, srcwin = srcwin):
            if region is not None:
                if xy_in_region_p(this_xyz, region.region):
                    xpos, ypos = _geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
                    ptArray[ypos, xpos] += 1
                    sumArray[ypos, xpos] += this_xyz[2]

def dump_168(entry, dst_port = sys.stdout, region = None, delimiter = ' ', verbose = False):
    with open(entry[0]) as infile:
        for line in xyz_parse(infile):
            if region is not None:
                if xy_in_region_p(line, region.region):
                    xyz_line(line, dst_port, delimiter, entry[2])
            else: xyz_line(line, dst_port, delimiter, entry[2])

def dump_200(entry, dst_port = sys.stdout, region = None, delimiter = ' ', verbose = False):
    if region is not None:
        srcwin = _srcwin(entry[0], region.region)
    else: srcwin = None
    gdal_dump(entry[0], dst_xyz = dst_port, delim = delimiter, weight = entry[2], \
                 dump_nodata = False, srcwin = srcwin, mask = None, warp_to_wgs = False)
                    
def dump_entry(entry, dst_port = sys.stdout, region = None, delimiter = ' ', verbose = False):
    '''dump a datalist entry as xyz.
    a datalist entry is [path, format, weight, ...]'''

    if verbose: echo_msg('using data file {}'.format(entry[0]))
    try:
        _dl_fmts[entry[1]][0](entry, dst_port, region, delimiter, verbose)
    except IOError as e:
        echo_error_msg('data dump broken. {}'.format(e))
        sys.exit(-2)
    
def archive_entry(entry, a_name, dirname = 'archive', region = None, delimiter = ' ', verbose = None):
    '''archive a datalist entry.
    a datalist entry is [path, format, weight, ...]'''
    
    i_dir = os.path.dirname(entry[0])
    i_xyz = os.path.basename(entry[0]).split('.')[0]
    
    a_dir = os.path.join(dirname, a_name, 'data', entry[-1])
    a_xyz_dir = os.path.join(a_dir, 'xyz')
    a_xyz = os.path.join(a_xyz_dir, i_xyz + '.xyz')
    
    if not os.path.exists(a_dir):
        os.makedirs(a_dir)

    if not os.path.exists(a_xyz_dir):
        os.makedirs(a_xyz_dir)

    with open(a_xyz, 'w') as fob:
        dump_entry(entry, dst_port = fob, region = region, delimiter = delimiter, verbose = verbose)

    d = datalist(os.path.join(a_dir, entry[-1] + '.datalist'))
    d._append_entry([os.path.join('xyz', i_xyz + '.xyz') if x[0] == 0 else 168 if x[0] == 1 else x[1] for x in enumerate(entry[:-1])])

    generate_inf(a_xyz, 168)
    
class datalist(object):
    '''Waffles and MBSystem style datalists for elevation data.
    Supports XYZ and GDAL data.'''

    def __init__(self, i_datalist, i_region = None, i_fmt = None, i_inc = None, verbose = False, callback = lambda: False):

        self.verbose = verbose
        self.stop = callback
        
        self.dl_entry = i_datalist.split(':')

        if len(self.dl_entry) > 1:
            if self.dl_entry[1] != -1:
                self._path = self.dl_entry[0].split('.')[0]+'.datalist'
            else: self._path = self.dl_entry[0]
        else: self._path = self.dl_entry[0]

        self._path_dirname = os.path.dirname(self._path)
        self._path_basename = os.path.basename(self._path)
        self._name = os.path.basename(self._path).split('.')[0]
        
        if not os.path.exists(self._path):
            open(self._path, 'a').close()

        if len(self.dl_entry) > 1:
            if self.dl_entry[1] != -1:
                self._append_entry(self.dl_entry)
            
        self.region = i_region
        self.i_fmt = i_fmt
        self.inc = i_inc
        self.want_weights = False
        
        self.delim = ' '
    
        self.datalist = []
        self.datalists = []
        self.datafiles = []
        self.archive_datalist = None

        if self.region is not None:
            self._name_r = self._name + self.region.fn
        else: self._name_r = self._name

        self._parse()
        
    def _reset(self):
        '''reload the datalist'''

        self.datalist = []
        self.datafiles = []
        self._parse()

    def _valid_p(self):
        '''validate the datalist'''

        ## ==============================================
        ## an empty datalist should technically be valid
        ## also check if lines have at least 2 columns
        ## ==============================================
        
        if len(self.datalist) > 0:
            return(True)
        else: return(False)

    def _valid_entry_p(self, data_l):
        '''validate a datafile entry line'''

        if data_l[0] == '#' or data_l[0] == '\n' or data_l[0] == '': return(False)
        dl_cols = [x.strip() for x in data_l.split(' ')]
        if len(dl_cols) < 2: return(False)
        path_d = os.path.join(self._path_dirname, dl_cols[0])
        if not os.path.exists(path_d): return(False)
        try:
            int(dl_cols[1])
        except: return(False)
        return(True)
        
    def _parse(self):
        '''read a datalist and parse the entries.'''
        
        ## ==============================================
        ## datalist file (no file options)
        ## ==============================================
        
        if len(self.dl_entry) == 1:
            with open(self._path, 'r') as fob:
                for dl in fob:
                    if self._valid_entry_p(dl):
                        dl_cols = [x.strip() for x in dl.split(' ')][:3]
                        dl_gmds_cols = [y.strip() for y in [x for x in dl.split(' ')][3:]]
                        dl_gmds_cols = [x.strip() for x in ' '.join(dl_gmds_cols).split(',')]
                        if len(dl_cols) == 3:
                            dl_cols.append(dl_gmds_cols)
                        self.datalist.append([os.path.join(self._path_dirname, x[1]) if x[0] == 0 else int(x[1]) if x[0] < 2 else x[1] for x in enumerate(dl_cols)])

        ## ==============================================
        ## single data file (has file options and there
        ## is no .datalist file to ref)
        ## ==============================================
        
        elif len(self.dl_entry) >= 2:
            if self._valid_entry_p(' '.join(self.dl_entry)):
                dl_cols = [x.strip() for x in self.dl_entry[:3]]
                dl_gmds_cols = [y.strip() for y in [x for x in self.dl_entry][3:]]
                dl_gmds_cols = [x.strip() for x in ' '.join(dl_gmds_cols).split(',')]
                if len(dl_cols) == 3:
                    dl_cols.append(dl_gmds_cols)
                self.datalist.append([self.dl_entry[0] if x[0] == 0 else int(x[1]) if x[0] < 2 else x[1] for x in enumerate(dl_cols)])

    def _proc_data(self, proc = lambda entry: sys.stdout.write('{}\n'.format(entry[0])), weight = None):
        '''Recurse through the datalist and run proc on each data file.'''

        for data_e in self.datalist:
            if self.stop(): break
            
            try:
                dformat = int(data_e[1])
            except: dformat = 168
            
            if weight is None and self.want_weights:
                weight = datalist_set_weight(data_e, weight)
                
            try:
                data_e[2] = weight
            except: data_e.append(weight)
            
            if dformat == -1:
                d = datalist(data_e[0], self.region, self.i_fmt, self.inc, self.verbose, self.stop)
                if self.i_fmt == -1:
                    data_e.append(d._name)
                    data_mb = data_e[:3]
                    data_gm = data_e[3:]
                    proc([x[1] if x[0] == 0 else x[1] for x in enumerate(data_e)])
                d._proc_data(proc, weight)

            else:
                if self.usable_entry_p(data_e):
                    data_e.append(self._name)
                    proc([x[1] if x[0] == 0 else x[1] for x in enumerate(data_e)])

    def _entry_dump(self, entry, dst_port = sys.stdout):
        '''dump a datalist entry as xyz.
        a datalist entry is [path, format, weight, ...]'''

        if self.verbose: echo_msg('using data file {}'.format(entry[0]))
        try:
            _dl_fmts[entry[1]][0](entry, dst_port, self.region, self.delim, self.verbose)
        except IOError as e:
            echo_error_msg('data dump broken. {}'.format(e))
            sys.exit(-2)
        
    def _dump_data(self, dst_port = sys.stdout):
        '''dump the data from the datalist to stdout.'''

        self._proc_data(lambda entry: self._entry_dump(entry, dst_port = dst_port))

    def _load_datalists(self):
        '''load a datalist and gather all datalists found therein.'''

        if self.verbose: echo_msg('scanning datalist `\033[1m{}\033[m` for datalists...'.format(self._name))

        fmt = self.i_fmt
        self.i_fmt = -1
        self._proc_data(lambda entry: self.datalists.append(entry))
        self.i_fmt = fmt

        ## single file or only data in datalist
        if len(self.datalists) == 0:
            self.datalists.append([self._path, -1, 1])

        if self.verbose: echo_msg('scanned datalist `\033[1m{}\033[m` and found {} datalists.'.format(self._name, len(self.datalists)))
        
    def _load_data(self):
        '''load a datalist and gather all datafiles found therein.'''

        if self.verbose: echo_msg('loading data files from datalist `\033[1m{}\033[m`...'.format(self._name))
        if self.i_fmt == -1:
            self._proc_data(lambda entry: self.datalists.append(entry))
        else: self._proc_data(lambda entry: self.datafiles.append(entry))

        if self.verbose: echo_msg('loaded datalist `\033[1m{}\033[m` and found {} data files.'.format(self._name, len(self.datafiles)))

    def _yield_data(self):
        '''Yield the xyz data from the datalist to a generator
        usage: for line in self._yield_data(): proc(line)'''

        self._load_data()
        
        for entry in self.datafiles:
            if entry[1] == 168:
                with open(entry[0], encoding='utf-8') as infile:
                    for line in xyz_parse(infile):
                        yield(line)
            elif entry[1] == 200:
                for this_xyz in gdal_yield(entry[0]):
                    yield(this_xyz)
            
    def _archive(self, dirname = 'archive'):

        self._proc_data(lambda entry: archive_entry(entry, a_name = self._name_r, dirname = dirname, region = self.region, delimiter = self.delim, verbose = self.verbose))
        self.datalists = [[self._path, -1, 1, self._name]]
        self.i_fmt = -1
        
        self._proc_data(lambda entry: self.datalists.append(entry))
        d = datalist(os.path.join(dirname, self._name_r, self._name_r + '.datalist'), verbose = self.verbose)
        self._load_datalists()
        
        for i in self.datalists:
            dl_n = i[-1]
            dl_p = os.path.join('data', dl_n, dl_n + '.datalist')
            if os.path.exists(os.path.join(dirname, self._name_r, dl_p)):
                d._append_entry([dl_p, -1, i[2]])
                
        self.archive_datalist = d
                            
    def _gen_inf(self):
        '''load a dalist and process datafiles'''

        self._proc_data(lambda entry: generate_inf(entry[0], entry[1]))

    def usable_entry_p(self, entry):
        usable = True
        if self.i_fmt is not None and self.i_fmt != entry[1]: usable = False

        if usable:
            if not os.path.exists(entry[0]): usable = False

            if usable:
                if self.region is not None:
                    ## gmt may fail if sep is set to non-space
                    dinf_region = parse_or_gen_inf(entry[0], data_fmt = entry[1])
                    if dinf_region is None: usable = False
                    if usable and not regions_intersect_p(self.region, dinf_region): usable = False
        return(usable)
        
    def _append_entry(self, entry):
        '''append a datalist entry to a datalist, given filename, format and weight'''

        with open(self._path, 'a') as outfile:
            outfile.write('{}\n'.format(' '.join([str(x) for x in entry])))

        if self.verbose:
            echo_msg('adding datalist entry: {}'.format(' '.join([str(x) for x in entry])))
            
        self._reset()

    def to_ogr(self, block = None, o_name = None, fmt = 'ESRI Shapefile'):
        '''dump the data from the DATALIST to OGR-compatible CSV file'''
        
        if o_name is None: o_name = self._name
            
        xyz2ogr(self._block_yield(block), '{}.shp'.format(o_name), dst_fmt = fmt)
                
    def gather_region(self):

        out_regions = []
        self._proc_data(lambda entry: out_regions.append(parse_or_gen_inf(entry[0], entry[1])))
        out_regions = [x for x in out_regions if x is not None]
        out_region = out_regions[0]

        for i in out_regions[1:]:
            out_region = regions_merge(out_region, i)
                
        self.region = out_region

    def _block_yield(self, increment = 1):
        '''block the datalist.'''
        
        if self.region is None: self.gather_region()
            
        xcount, ycount, dst_gt = _extent2gt(self.region.region, cellsize)
        sumArray = np.zeros((ycount, xcount))
        ptArray = np.zeros((ycount, xcount))

        self._proc_data(lambda entry: block_entry(entry, self.region, increment, ptArray, sumArray, dst_gt, self.verbose))

        ptArray[ptArray==0] = np.nan
        outArray = sumArray / ptArray
        for y in range(0, ycount, 1):
            for x_i in range(0, xcount, 1):
                geo_x, geo_y = _pixel2geo(x_i, y, dst_gt)
                z = outArray[y][x_i]
                if not np.isnan(z):
                    yield([geo_x, geo_y, z])
        
    def mask(self, region = None, inc = .000277777, o_name = None, o_fmt = 'GTiff'):
        '''Generate a num-msk with GDAL from a datalist.'''

        if region is None:
            if self.region is None:
                self.gather_region()
            region = self.region.region

        if o_name is None:
            o_name = self._name

        o_mask = '{}_msk.{}'.format(o_name, _fext(o_fmt))
        
        xyz_mask(self._yield_data(), o_mask, region, inc, dst_format = o_fmt, verbose = self.verbose)
        
        return(o_mask)    

# class entry(datalist):

#     _inherited = ['i_datalist', 'i_region', 'i_fmt', 'i_inc', 'verbose', 'callback']
    
#     def __init__(self, entry_line):
        
        
## =============================================================================
##
## GDAL/OGR/OSR - py
## wrapper functions etc. for GDAL
##
## =============================================================================

gdal.PushErrorHandler('CPLQuietErrorHandler')
GDAL_OPTS = ["COMPRESS=LZW", "INTERLEAVE=PIXEL", "TILED=YES",\
             "SPARSE_OK=TRUE", "BIGTIFF=YES" ]

def gdal_inf(src_gdal):
    '''generate an info (.inf) file from a src_gdal file using GDAL.'''
        
    #minmax = (map(str, _extent(src_gdal)))
    minmax = [str(x) for x in _extent(src_gdal)]
    with open('{}.inf'.format(src_gdal), 'w') as inf:
        inf.write('{}\n'.format(' '.join(minmax)))
    return(0, 0)

def _ogr_create_polygon(coords):
    '''convert coords to Wkt'''

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[1], coord[0])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    poly_wkt = poly.ExportToWkt()
    poly = None
    
    return(poly_wkt)

def _extent2geom(extent):
    '''convert an extent [west, east, south, north] to an 
    OGR geometry'''

    eg = [[extent[2], extent[0]], [extent[2], extent[1]],
          [extent[3], extent[1]], [extent[3], extent[0]],
          [extent[2], extent[0]]]

    geom = ogr.CreateGeometryFromWkt(_ogr_create_polygon(eg))
    
    return(geom)

def _ogr_get_fields(src_ogr):
    '''return all fields in src_ogr'''

    schema = []
    source = ogr.Open(src_ogr)
    
    if source is not None:
        layer = source.GetLayer()
        ldefn = layer.GetLayerDefn()

        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
            
    source = None
    return(schema)

def _ogr_get_layer_fields(src_lyr):
    '''return all fields in src_lyr'''

    schema = []
    if src_lyr is not None:
        ldefn = src_lyr.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
    
    return(schema)

def ogr_mask_union(src_layer, src_field, dst_defn = None, callback = lambda: False, verbose = False):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.'''

    if verbose: sys.stderr.write('geomods: unioning polygons...')
    
    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()

    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for f in src_layer:
        if not callback():
            if f.GetField(src_field) == 0:
                src_layer.DeleteFeature(f.GetFID())
            elif f.GetField(src_field) == 1:
                f.geometry().CloseRings()
                wkt = f.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                src_layer.DeleteFeature(f.GetFID())
                        
    union = multi.UnionCascaded()
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometry(union)
    union = multi = None

    if verbose: sys.stderr.write('.ok\n')
    
    return(out_feat)

def _ogr_extents(src_ds):
    '''return the extent(s) of the ogr dataset'''
    
    these_extents = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                these_extents.append(pgeom.GetEnvelope())
                
    return(these_extents)

def _fext(src_drv_name):
    '''return the common file extention given a GDAL driver name'''
    
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER):
            fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)

        if fexts is not None:
            fext = fexts.split()[0]
    except:
        if src_drv_name == 'GTiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        
    return(fext)

def _gather_infos(src_ds):
    '''Gather information from `src_ds` GDAL dataset.'''

    ds_config = {}
    ds_config['nx'] = src_ds.RasterXSize
    ds_config['ny'] = src_ds.RasterYSize
    ds_config['nb'] = src_ds.RasterCount
    ds_config['geoT'] = src_ds.GetGeoTransform()
    ds_config['proj'] = src_ds.GetProjectionRef()
    ds_config['dt'] = src_ds.GetRasterBand(1).DataType
    ds_config['dtn'] = gdal.GetDataTypeName(src_ds.GetRasterBand(1).DataType)
    ds_config['ndv'] = src_ds.GetRasterBand(1).GetNoDataValue()
    ds_config['fmt'] = src_ds.GetDriver().ShortName

    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
    
    return(ds_config)

def _set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''Set a datasource config dictionary'''
    
    ds_config = {}

    ds_config['nx'] = nx
    ds_config['ny'] = ny
    ds_config['nb'] = nb
    ds_config['geoT'] = geoT
    ds_config['proj'] = proj
    ds_config['dt'] = dt
    ds_config['ndv'] = ndv
    ds_config['fmt'] = fmt

    return(ds_config)

def _infos(src_fn, full = False):
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)

        if full:
            t = ds.ReadAsArray()
            max_z = np.max(t)
            min_z = np.min(t)
            ds_config['zmin'] = min_z
            ds_config['zmax'] = max_z
            
        t = ds = None
        
        return(ds_config)
    else: return(None)

def _cpy_infos(src_config):
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''Convert a geographic x,y value to a pixel location of geoTransform'''

    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else:
        pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt( geoTransform))

    return(int(pixel_x), int(pixel_y))

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''Convert a pixel location to geographic coordinates given geoTransform'''

    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)

    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001: return
    invDet = 1.0 / det

    ## ==============================================
    ## compute adjoint and divide by determinate
    ## ==============================================

    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet

    return(outGeoTransform)

def _gt2extent(ds_config, warp_to_wgs = False):
    '''convert a gdal geo-tranform to an extent [w, e, s, n]'''
    
    geoT = ds_config['geoT']
    
    x_origin = geoT[0]
    y_origin = geoT[3]

    if warp_to_wgs:
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt(ds_config['proj'])
        
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(4326)
        
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)

        point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(x_origin, y_origin))
        point.Transform(dst_trans)
        pnt = point.GetPoint()
        x_origin = pnt[0]
        y_origin = pnt[1]
        
    x_inc = geoT[1]
    y_inc = geoT[5]

    return([x_origin, x_origin + x_inc * ds_config['nx'], y_origin + y_inc * ds_config['ny'], y_origin])

## add z
def _extent(src_fn, warp_to_wgs = False):
    '''return the extent of the src_fn gdal file.'''
    
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)
        extent = _gt2extent(ds_config, warp_to_wgs)
        ds = None
        return(extent)
    else: return(None)

def _extent2gt(extent, cellsize):
    '''return a count info and a gdal geotransform based on extent and cellsize
    output is a list (xcount, ycount, geot)'''
    
    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1
    dst_gt = (extent[0], cellsize, 0, extent[3], 0, (cellsize * -1.))

    return(xcount, ycount, dst_gt)
    
def _srcwin(src_fn, extent):
    '''given a gdal file src_fn and an extent [w, e, s, n],
    output the appropriate gdal srcwin.'''
    
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)

        gt = ds_config['geoT']
        this_origin = _geo2pixel(extent[0], extent[3], gt)
        this_end = _geo2pixel(extent[1], extent[2], gt)
        this_size = (int(this_end[0] - this_origin[0]), int(this_end[1] - this_origin[1]))
        this_origin = [0 if x < 0 else x for x in this_origin]
        this_size = [0 if x < 0 else x + 1 for x in this_size]
        if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
        if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
        ds = None
        return(this_origin[0], this_origin[1], this_size[0], this_size[1])
    else: return(None)

def _sr_wkt(epsg, esri = False):
    '''convert an epsg code to wkt'''
    
    wkt = None
    try:
        int(epsg)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        if esri: sr.MorphToESRI()
        wkt = sr.ExportToWkt()
        sr = None
        return(wkt)
    except:
        sys.stderr.write('geomods: error, invalid epsg code\n')
        return(None)

def _prj_file(dst_fn, epsg):
    '''generate a .prj file given an epsg code'''
    
    with open(dst_fn, 'w') as out:
        out.write(_sr_wkt(int(epsg), True))
    return(0)

def _set_gdal_epsg(src_gdal, epsg = 4326):
    '''set the projection of src_gdal to epsg'''
    
    ds = gdal.Open(src_gdal, gdal.GA_Update)
    if ds:
        ds.SetProjection(_sr_wkt(int(epsg)))
        ds = None

        return(0)
    else: return(None)
        
def gdal_write (src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff', verbose = False):
    '''write src_arr Array to gdal file dst_gdal using src_config'''

    if verbose: sys.stderr.write('geomods: writing gdal grid...')

    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        driver.Delete(dst_gdal)
    
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        ds.SetProjection(ds_config['proj'])
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        ds.GetRasterBand(1).WriteArray(src_arr)
        ds = None
        if verbose: sys.stderr.write('ok\n')
        return(0)
    else: return(None)
    
def gdal_chunks(src_fn, n_chunk = 10, verbose = False):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.'''

    band_nums = []
    o_chunks = []
    
    if band_nums == []: band_nums = [1]
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk
    
    src_ds = gdal.Open(src_fn)

    if src_ds is not None:
        ds_config = _gather_infos(src_ds)
        band = src_ds.GetRasterBand(1)
        gt = ds_config['geoT']

        if verbose: sys.stderr.write('geomods: chunking grid...')
        while True:
            if verbose: sys.stderr.write('.')
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
                if this_x_size == 0 or this_y_size == 0:
                    break
                
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)
                dst_gt = [this_geo_x_origin, gt[1], 0.0, this_geo_y_origin, 0.0, gt[5]]
                
                #for band in bands:
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                if not np.all(band_data == band_data[0,:]):
                    o_chunk = '{}_chnk{}x{}.tif'.format(os.path.basename(src_fn).split('.')[0], x_i_chunk, i_chunk)
                    dst_fn = os.path.join(os.path.dirname(src_fn), o_chunk)
                    o_chunks.append(dst_fn)

                    dst_config = _cpy_infos(ds_config)
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
                
        if verbose: sys.stderr.write('ok\n')
        src_ds = None
        
        return(o_chunks)
    else: return(None)

def gdal_crop(src_fn):
    '''Crop `src_fn` GDAL file by it's NoData value.
    Returns cropped array.'''
    
    src_ds = gdal.Open(src_fn)

    if src_ds is not None:
        ds_config = _gather_infos(src_ds)
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray()
        src_ds = None

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
    else: return(None)

def gdal_cut(src_fn, srcwin, dst_fn):
    '''cut src_fn gdal file to srcwin and output dst_fn gdal file'''
    
    src_ds = gdal.Open(src_fn)
    
    if src_ds is not None:
        ds_config = _gather_infos(src_ds)
        gt = ds_config['geoT']
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
        ds_config = _set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, _sr_wkt(4326), ds_config['dt'], ds_config['ndv'], ds_config['fmt'])
        src_ds = None
        return(gdal_write(ds_arr, dst_fn, ds_config))
    
    return(None)

## todo: cleanup this function...
def gdal_dump(src_gdal, dst_xyz = sys.stdout, delim = ' ', weight = None, dump_nodata = False, srcwin = None, mask = None, warp_to_wgs = False):
    '''Dump `src_gdal` GDAL file to ASCII XYZ'''

    status = 0
    band_nums = []
    skip = 1
    
    msk_band = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
    
    if band_nums == []: band_nums = [1]

    src_ds = gdal.Open(src_gdal)

    if weight is None:
        w_string = ''
    else: w_string = '{}{}'.format(delim, weight)

    if src_ds is not None:
        bands = []
        for band_num in band_nums: 
            band = src_ds.GetRasterBand(band_num)
            if band is not None:
                bands.append(band)

        ds_config = _gather_infos(src_ds)
        gt = ds_config['geoT']

        if warp_to_wgs:
            src_srs = osr.SpatialReference()
            src_srs.ImportFromWkt(ds_config['proj'])
            
            dst_srs = osr.SpatialReference()
            dst_srs.ImportFromEPSG(4326)
            
            dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        
        if srcwin is None:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])

        dst_fh = dst_xyz
        zs = len(bands)
        if weight is not None: zs+=1
            
        band_format = (("%g" + delim) * zs).rstrip(delim) + '\n'

        if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
           and abs(ds_config['nx'] * gt[1]) < 180 \
           and abs(ds_config['ny'] * gt[5]) < 180:
            format = '%.10g' + delim + '%.10g' + delim + '%s'
        else: format = '%.3f' + delim + '%.3f' + delim + '%s'

        for y in range(srcwin[1], srcwin[1] + srcwin[3], skip):
            nodata = ['-9999', 'nan']
            data = []
            for band in bands:
                if band.GetNoDataValue() is not None:
                    #nodata.append(('{:.10}'.format(band.GetNoDataValue())))
                    nodata.append('%g' % band.GetNoDataValue())
                band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)

                if msk_band is not None:
                    msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
                    band_data[msk_data==0]=-9999

                band_data = np.reshape(band_data, (srcwin[2], ))
                data.append(band_data)

            for x_i in range(0, srcwin[2], skip):
                x = x_i + srcwin[0]

                #geo_x, geo_y = _pixel2geo(x, y, gt)
                geo_x = gt[0] + (x + 0.5) * gt[1] + (y + 0.5) * gt[2]
                geo_y = gt[3] + (x + 0.5) * gt[4] + (y + 0.5) * gt[5]

                x_i_data = []
                for i in range(len(bands)):
                    x_i_data.append(data[i][x_i])

                z = x_i_data[0]

                if weight is not None:
                    band_str = band_format % (z, weight)
                else: band_str = band_format % (z)
                
                if warp_to_wgs:
                    point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(geo_x, geo_y))
                    point.Transform(dst_trans)
                    pnt = point.GetPoint()
                    line = '{} {} {}{}'.format(pnt[0], pnt[1], band_str, w_string)
                else: line = format % (float(geo_x), float(geo_y), band_str)

                if dst_fh == sys.stdout: line = line.encode(sys.stdout.encoding)

                line = line.encode('utf-8')
                
                if dump_nodata:
                    dst_fh.write(line)
                else:
                    #print band_str.split()[0]
                    if band_str.split()[0] not in nodata:
                        dst_fh.write(line)

        srcds = src_mask = None

def gdal_yield(src_gdal, dump_nodata = False, srcwin = None):
    '''Yield `src_gdal` GDAL file to XYZ.'''

    srcds = gdal.Open(src_gdal)

    if srcds is not None:
        band = srcds.GetRasterBand(1)
        gt = srcds.GetGeoTransform()
        if srcwin is None:
            srcwin = (0,0,srcds.RasterXSize,srcds.RasterYSize)

        for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):

            nodata = ['{:.10f}',format(-9999), 'nan']
            data = []

            if band.GetNoDataValue() is not None: nodata.append('{:.10f}'.format(band.GetNoDataValue()))
            
            band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)    
            band_data = np.reshape(band_data, (srcwin[2], ))
            data.append(band_data)

            for x_i in range(0, srcwin[2], 1):
                x = x_i + srcwin[0]
                geo_x, geo_y = _pixel2geo(x, y, gt)
                x_i_data = []
                for i in range(1):
                    x_i_data.append(data[i][x_i])
                    
                z = x_i_data[0]
                line = [geo_x, geo_y, z]
                
                if dump_nodata:
                    yield(line)
                else:
                    if '{:.10f}'.format(z) not in nodata:
                        yield(line)
        srcds = None
    
def gdal_infos(src_fn, full = False):
    '''return the ds_config and extent of the src_fn gdal file.'''
    
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)
        ds = None
        return(ds_config, _gt2extent(ds_config))
    else: return(None)
        
def gdal_null(dst_fn, extent, cellsize, nodata = -9999, outformat = 'GTiff'):
    '''generate a `null` grid with gdal'''

    xcount, ycount, dst_gt = _extent2gt(extent, cellsize)    
    null_array = np.zeros((ycount, xcount))
    null_array[null_array == 0] = nodata
    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, _sr_wkt(4326), gdal.GDT_Float32, -9999, outformat)
    
    return(gdal_write(null_array, dst_fn, ds_config))
    
def gdal_percentile(src_fn, perc = 95):
    '''calculate the `perc` percentile of src_fn gdal file.'''

    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        x_dim = ds_array.shape[0]
        ds_array_flat = ds_array.flatten()
        p = np.percentile(ds_array_flat, perc)

        if p==p:
            percentile=p
            if percentile < 2:
                percentile = 2
        else: percentile = 1

        ds = ds_array = None

        return(percentile)
    else: return(None)

def gdal_polygonize(src_gdal, dst_layer, verbose = False):
    '''run gdal.Polygonize on src_gdal and add polygon to dst_layer'''
    
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        srcband = src_ds.GetRasterBand(1)

        if verbose: sys.stderr.write('geomods: polygonizing grid...')
        try:
            gdal.Polygonize(srcband, None, dst_layer, 0, [])
        except KeyboardInterrupt as e:
            sys.exit(1)
        if verbose: sys.stderr.write('ok\n')
        src_ds = srcband = None
        return(0)
    else: return(None)

def gdal_proximity(src_fn, dst_fn):
    '''Compute a proximity grid via GDAL'''

    prog_func = None

    src_ds = gdal.Open(src_fn)
    dst_ds = None
    
    if src_ds is not None:
        src_band = src_ds.GetRasterBand(1)
        ds_config = _gather_infos(src_ds)
        
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

def gdal_query(src_xyz, src_grd, out_form):
    '''Query a gdal-compatible grid file with xyz data.'''

    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================

    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = _gather_infos(ds)
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
            except: z = dsnodata

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)

                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd

                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    c = con_dec(math.fabs(float(d / (g + 0.00000001) * 100)), 2)
                    s = con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    d = con_dec(d, 4)

                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    xyzl.append(np.array(outs, dtype = ds_config['dtn']))
        dsband = ds = None
        out_array = np.array(xyzl, dtype = ds_config['dtn'])

    return(out_array)
    
def gdal_split(src_gdal, split_value = 0):
    '''split raster into two based on z value'''

    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_upper.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_lower.tif'.format(os.path.basename(src_gdal)[:-4]))
                                         
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = _gather_infos(src_ds)
        dst_config = _cpy_infos(src_config)
        dst_config['fmt'] = 'GTiff'

        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny']) 

        upper_array = ds_arr
        upper_array[upper_array <= split_value] = src_config['ndv'] 
        gdal_write(upper_array, dst_upper, dst_config)
        upper_array = None

        lower_array = ds_arr
        lower_array[lower_array >= split_value] = src_config['ndv']
        gdal_write(lower_array, dst_lower, dst_config)
        lower_array = src_ds = None
        
        return([dst_upper, dst_lower])
    else: return(None)

def gdal_sum(src_gdal):
    '''sum the z vale of src_gdal'''
    
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_array = ds.GetRasterBand(1).ReadAsArray() 
        sums = np.sum(ds_array)
        ds = ds_array = None
        return(sums)
    else: return(None)

def osr_transformation(src_epsg, dst_epsg):

    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(src_epsg)

    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(dst_epsg)

    dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
    return(dst_trans)
                    
def gdal_transform(src_gdal):

    src_ds = gdal.Open(src_gdal)

    src_config = _gather_infos(src_ds)
    
    ##Getting spatial reference of input raster
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(src_config['proj'])

    print(src_srs.ExportToWkt())
    
    # WGS84 projection reference
    #OSR_WGS84_REF = osr.SpatialReference()
    #OSR_WGS84_REF.ImportFromEPSG(4326)
    
    # OSR transformation
    #wgs84_transformation = osr.CoordinateTransformation(src_srs, OSR_WGS84_REF)

    #return(wgs84_transformation)
    
    #wkt_geom.Transform(wgs84_to_image_trasformation)
        
def xyz2gdal(src_xyz, dst_gdal, extent, cellsize,
             dst_format='GTiff', zvalue='d', xloc=0, yloc=1, zloc=2, 
             delim=' ', verbose=False, overwrite=False):
    '''Create a GDAL supported grid from xyz data
    `zvalue` of `d` generates a num grid'''

    pointnum = 0
    dst_nodata=-9999
    xcount, ycount, dst_gt = _extent2gt(extent, cellsize)    
    if zvalue == 'z': sumArray = np.zeros((ycount, xcount))
    ptArray = np.zeros((ycount, xcount))

    for this_xyz in src_xyz:
        x = this_xyz[xloc]
        y = this_xyz[yloc]
        z = this_xyz[int(zloc)].strip()

        if x > extent[0] and x < extent[1]:
            if y > extent[2] and y < extent[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                ptArray[ypos, xpos] = 1

                if zvalue == 'z': sumArray[ypos, xpos] += float(z)
                if zvalue == 'd' or zvalue == 'z': 
                    ptArray[ypos, xpos] += 1
                else: ptArray[ypos, xpos] = 1
        pointnum += 1

    if zvalue == 'z':
        outarray = sumArray / ptArray
    elif zvalue == 'd': outarray = ptArray
    else: outarray = ptArray

    outarray[np.isnan(outarray)] = dst_nodata
    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, _sr_wkt(4326), gdal.GDT_Int32, dst_nodata, 'GTiff')
    return(gdal_write(outarray, dst_gdal, ds_config))

## add option for grid/pixel node registration...currently pixel node only.
def xyz_mask(src_xyz, dst_gdal, extent, cellsize,
             dst_format='GTiff', xloc=0, yloc=1, zloc=2, 
             delim=' ', verbose=False):
    '''Create a num grid mask of xyz data. The output grid
    will contain 1 where data exists and 0 where no data exists.'''

    xcount, ycount, dst_gt = _extent2gt(extent, cellsize)
    ptArray = np.zeros((ycount, xcount))

    for this_xyz in src_xyz:
        
        x = this_xyz[xloc]
        y = this_xyz[yloc]

        if x > extent[0] and x < extent[1]:
            if y > extent[2] and y < extent[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                ptArray[ypos, xpos] = 1

    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, _sr_wkt(4326), gdal.GDT_Int32, -9999, 'GTiff')
    return(gdal_write(ptArray, dst_gdal, ds_config, verbose = verbose))

## this spits out a lot of errors due to field width with shapefiles when
## sending data from the datalist...update so fields are floats perhaps...
def xyz2ogr(src_xyz, dst_ogr, xloc = 0, yloc = 1, zloc = 2, dst_fmt = 'ESRI Shapefile',\
            overwrite = True, verbose = False):
    '''Make a point vector OGR file from a src_xyz table data'''
    
    driver = ogr.GetDriverByName(dst_fmt)
    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)
    ds = driver.CreateDataSource(dst_ogr)

    layer = ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPoint25D)

    oft_title = "Longitude"
    oft_string = ogr.OFTReal
    fdw = 22
    fdp = 20
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)

    oft_title = "Latitude"
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)

    oft_title = "Elevation"
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)
    
    f = ogr.Feature(feature_def = layer.GetLayerDefn())
    
    #for this_xyz in xyz_parse(src_xyz):
    for this_xyz in src_xyz:
        x = this_xyz[xloc]
        y = this_xyz[yloc]
        z = this_xyz[zloc]

        f.SetField(0, x)
        f.SetField(1, y)
        f.SetField(2, z)

        wkt = 'POINT(%.8f %.8f %.10f)' % (x,y,z)
        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)

def gdal2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''Convert the gdal file to gdal using gdal'''
    
    status = 0
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], _fext(dst_fmt))
        gdal2gdal_cmd = ('gdal_translate {} {} -f {}\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(gdal2gdal_cmd, verbose = verbose)

        if status != 0:
            dst_gdal = None
    else: dst_gdal = None

    return(dst_gdal)        
        
def run_gdal_grid(datalist, region, inc, o_name, module='invdist', node = 'pixel', gmt = None, verbose = False):
    '''run the gdal command gdal_grid and pass the data from the datalist.'''
    
    o_dem = '{}.tif'.format(o_name)
    datalist.to_ogr(block = inc, o_name = o_name)

    out_size_x = int((region.east - region.west) / inc) + 1
    out_size_y = int((region.north - region.south) / inc) + 1

    gg_cmd = 'gdal_grid -zfield "Elevation" -txe {:.10f} {:.10f} -tye {:.10f} {:.10f} -outsize {} {} \
    -a {} -l {} {}.shp {} --config GDAL_NUM_THREADS ALL_CPUS\
    '.format(region.west, region.east, region.north, region.south, out_size_x, out_size_y, module, o_name, o_name, o_dem)
    out, status = run_cmd(gg_cmd, verbose = verbose)

    return(o_dem)
        
## =============================================================================
##
## VDatum - vdatumfun.py
## wrapper functions for NOAA's VDatum
##
## Currently only compatible with VDatum > 4.0
##
## =============================================================================

class vdatum:
    '''vdatum object to communicate with NOAA's VDatum'''

    def __init__(self, vdatum_path = None, verbose = False):
        self.verbose = verbose
        self.status = 0
        
        if vdatum_path is None:
            try:
                self.vdatum_path = _waff_co.get('VDATUM', 'jar')
            except:
                self.vdatum_path = self._find_vdatum()[0]
                if self.status != 0:
                    self.vdatum_path = None
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
            out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, '-'), verbose = self.verbose)
            for i in out.decode('utf-8').split('\n'):
                if '- v' in i.strip():
                    self._version = i.strip().split('v')[-1]
                    break

    def _find_vdatum(self):
        '''Find the VDatum executable on the system and 
        return a list of found vdatum.jar paths'''

        results = []
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.abspath(os.path.join(root, 'vdatum.jar')))
                break
        if len(results) <= 0:
            self.status = -1
        if self.verbose: echo_msg('found vdatum {}'.format(results[0]))

        return(results)

    def run_vdatum(self, src_fn):
        '''Run vdatum on src_fn which is an XYZ file'''
        
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}\
        '.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.fd, src_fn, self.ds_dir, self.region)
        out, status = run_cmd('java -jar {} {}'.format(self.vdatum_path, vdc), verbose = self.verbose)
        #out, status = run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, True)

        return(status)

## =============================================================================
##
## GMT Wrapper Functions - gmtfun.py
## wrapper functions to GMT system commands
##
## GMT must be installed on the system to run these functions
## and commands.
##
## =============================================================================

def gmt_inf(src_xyz):
    '''generate an info (.inf) file from a src_xyz file using GMT.'''
    
    out, status = run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), verbose = False)
    return(out, status)

def gmt_grd_inf(src_gdaxl):
    '''generate an info (.inf) file from a src_gdal file using GMT.'''
        
    out, status = run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_gdal, src_gdal), verbose = False)
    return(out, status)

def gmt_inc2inc(inc_str):
    '''convert an GMT-style inc_str (6s) to native units'''
    
    units = inc_str[-1]

    if units == 'c': inc = float(inc_str[:-1]) / 3600
    elif units == 's': inc = float(inc_str[:-1]) / 3600
    elif units == 'm': inc = float(inc_str[:-1]) / 360
    else: inc = float(inc_str)    
    
    return(inc)

def grd2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''Convert the grd file to tif using GMT'''
    
    status = 0
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], _fext(dst_fmt))
        grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(grd2gdal_cmd, verbose = verbose)

        if status != 0:
            dst_gdal = None
    else: dst_gdal = None

    return(dst_gdal)

def grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''

    status = 0
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    if not os.path.exists(src_grd):
        return([])
    else:
        grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
        out, status = run_cmd(grdinfo_cmd, verbose = verbose)
        try:
            os.remove('gmt.conf')
        except: pass
        
        if status !=0:
            return([])
        else: return(out.split())

def gmtinfo(src_xyz, verbose = False):
    '''Return an info list of `src_xyz`'''

    status = 0
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    if not os.path.exists(src_xyz):
        return([])
    else:
        gmtinfo_cmd = ('gmt gmtinfo {} -C'.format(src_xyz))
        out, status = run_cmd(gmtinfo_cmd, verbose = verbose)
        remove_glob('gmt.conf')
        
        if status !=0:
            return([])
        else: return(out.split())        

def gmt_block(datalist, mode = 'blockmean', inc = '1s', o_name = None, delim = 'SPACE', weights = False, verbose = False):
    '''run block/mean/median on src_xyz'''

    status = 0
    if mode == 'blockmean' or mode == 'blockmean':
        out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = {}'.format(delim.upper()), verbose = verbose)
        if mode == 'blockmean' and weights:
            mode = 'blockmean -Wi'
            datalist.want_weights = True
        if mode == 'blockmedian': mode = 'blockmedian -Q'
        if o_name is None: o_name = datalist._name
        if delim.lower() == 'comma':
            out_ext = 'csv'
            o_vrt = open('{}.vrt'.format(o_name), 'w')
            t = '''<OGRVRTDataSource>
  <OGRVRTLayer name="{}">
    <SrcDataSource>{}.csv</SrcDataSource>
    <GeometryType>wkbPoint</GeometryType>
    <GeometryField encoding="PointFromColumns" x="field_1" y="field_2" z="field_3"/>
  </OGRVRTLayer>
</OGRVRTDataSource>'''.format(o_name, o_name)
            o_vrt.write(t)
            o_vrt.close()

        else: out_ext = 'xyz'
        
        if os.path.exists(datalist._path):
            blk_cmd1 = ('gmt {} -V {} -I{} > {}.{}'.format(mode, datalist.region.gmt, inc, o_name, out_ext))
            out, status = run_cmd(blk_cmd1, verbose = True, data_fun = datalist._dump_data)
        else: status = -1
    else: status = -1
    remove_glob('gmt.conf')
    
    return(status)
        
def gmtselect_split(o_xyz, sub_region, sub_bn, verbose = False):
    '''split an xyz file into an inner and outer region.'''

    status = 0
    out_inner = None
    out_outer = None

    gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
    out, status = run_cmd(gmt_s_inner, verbose = verbose)

    if status == 0: out_inner = '{}_inner.xyz'.format(sub_bn)

    gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
    out, status = run_cmd(gmt_s_outer, verbose = verbose)

    if status == 0:  out_outer = '{}_outer.xyz'.format(sub_bn)

    return([out_inner, out_outer])
        
def grdcut(src_grd, src_region, dst_grd, verbose = False):
    '''Cut `src_grd` to `src_region` '''

    status = 0
    if os.path.exists(src_grd):
        cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
        out, status = run_cmd(cut_cmd1, verbose = True)
    else: status = -1

    return(status)

def grdfilter(src_grd, dst_grd, dist = '3s', verbose = False):
    '''filter `src_grd` '''

    status = 0
    if os.path.exists(src_grd):
        ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1'.format(src_grd, dst_grd, src_grd, dist))
        out, status = run_cmd(ft_cmd1, verbose = verbose)
    else: status = -1

    return(status)

def grd2xyz(src_grd, dst_xyz, region = None, mask = None, verbose = False, want_datalist = False):
    '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region.
    Optionally, generate a datalist and inf file for the resultant xyz data.'''

    status = 0
    if mask:
        grdmask_cmd = ('gmt grdmath -N -V {} {} OR = tmp.grd'.format(src_grd, mask))
        out, status = run_cmd(grdmask_cmd, verbose = verbose)
        if status == 0: 
            src_grd = 'tmp.grd'

    if region and region._valid:
        region_str = region.gmt
    else: region_str = ''

    grd2xyz_cmd = ('gmt grd2xyz -V {} -s {} > {}'.format(src_grd, region_str, dst_xyz))
    out, status = run_cmd(grd2xyz_cmd, verbose = verbose)

    if status == 0:
        if mask:
            if os.path.exists('tmp.grd'):
                os.remove('tmp.grd')

        if want_datalist:
            s_datalist = datalist('{}.datalist'.format(dst_xyz.split('.')[0]))
            s_datalist._append_datafile(['{}'.format(os.path.basename(dst_xyz)), 168, 1])
            s_datalist._reset()

            mb_inf(s_datalist._path, -1)
        
    return(status)

def gmt_nan2zero(src_grd, verbose = False):
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = tmp.grd'.format(src_grd))
    out, status = run_cmd(num_msk_cmd, verbose = verbose)
    if status == 0:
        grd2gdal('tmp.grd')
        remove_glob('tmp.grd')
        os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_grdcut(src_grd, region, verbose = False):
    cut_cmd = ('gmt grdcut -V {} -Gtmp.grd {}\
    '.format(src_grd, region.gmt))
    out, status = run_cmd(cut_cmd, verbose = True)
    if status == 0:
        remove_glob(src_grd)
        os.rename('tmp.grd', '{}'.format(src_grd))
    return(status)

def slope(src_dem, dst_slp, verbose = False):
    '''Generate a Slope grid from a DEM with GMT'''

    status = 0
    o_b_name = '{}'.format(src_dem.split('.')[0])

    slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
    '.format(src_dem, o_name, src_dem))
    out, status = run_cmd(slope_cmd0, verbose = verbose)

    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}\
        '.format(o_b_name, dst_slp))
        out, status = run_cmd(slope_cmd1, verbose = verbose)
        
    if os.path.exists('{}_pslp.grd'.format(o_b_name)):
        os.remove('{}_pslp.grd'.format(o_b_name))

    return(status)

def num_msk(num_grd, dst_msk, verbose = False):
    '''Generate a num-msk from a NUM grid.'''

    status = 0

    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    out, status = run_cmd(num_msk_cmd, verbose = verbose)

    return(status)

def xyz2grd(datalist, region, inc, dst_name, a = 'n', node = 'pixel', verbose = False):
    '''Run the GMT command `xyz2grd` given a datalist, region and increment.'''
   
    status = 0
    if node == 'pixel':
        reg_str = '-r'
    else: reg_str = ''
    
    num_cmd0 = ('gmt xyz2grd -V {} -I{:.10f} -G{} -A{} {}\
    '.format(region.gmt, inc, dst_name, a, reg_str))
    out, status = run_cmd(num_cmd0, verbose = verbose, data_fun = datalist._dump_data)

    return(out, status)

def gmt_sample_gnr(src_grd, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node
    grid registration.'''
    
    out, status = run_cmd('gmt grdsample -T {} -Gtmp.tif=gd+n-9999:GTiff'.format(src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))

    return(status)

## =============================================================================
##
## MB-System Wrapper Functions - mbsfun.py
##
## MS-System must be installed on the system to run
## these functions and commands.
##
## =============================================================================

def mb_inf(src_xyz, src_fmt = 168):
    '''generate an info (.inf) file from a src_xyz file using MBSystem.'''

    out, status = run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz), True)
    return(out, status)

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):
    '''Run the MBSystem command `mbgrid` given a datalist, region and increment.
    The datalist should be an MBSystem-style datalist; if using a waffles datalist, 
    it should be converted first (using datalist._archive().'''
    
    status = 0
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1:
        dist = dist + '/2'
    
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.10f}/{:.10f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} > mb_proc.txt \
    '.format(datalist._path, region.gmt, inc, inc, dst_name, dist, tension, e_switch))
    out, status = run_cmd(mbgrid_cmd, verbose = verbose)

    return(out, status)

## =============================================================================
##
## spatial-metadata module:
##
## Generate spatial metadata from a datalist of xyz elevation data.
## The output is a shapefile with the boundaries of each specific
## datalist as a layer feature.
##
## Specify an input datalist, region and cell-size.
## o_name is the output prefix: `{o_name}_sm.shp`
##
## Uses geomods modules: datalists, gdalfun, utils
##
## Generates the boundaries using BOUNDS > GMT > GDAL
## if BOUNDS is not present on the system, will use
## GDAL to generate the boundaries which is a slower
## on large datasets.
##
## =============================================================================

class spatial_metadata:
    '''Generate spatial metadata from a datalist of xyz elevation data.
    The output is a shapefile with the unioned boundaries of each specific
    datalist as a layer feature.
    Each datalist entry should hold a comma-separated field which holds
    the dataset field information:         
    "Name,Agency,Date,Type,Resolution,HDatum,VDatum,URL"'''
    
    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None, o_extend = 6, callback = lambda: False, verbose = False):

        try:
            import Queue as queue
        except ModuleNotFoundError:
            import queue as queue
        
        self.dl_q = queue.Queue()
        self.datalist = i_datalist
        self.inc = i_inc
        self.region = i_region
        self.extend = int(o_extend)
        self.dist_region = self.region.buffer(self.extend * self.inc)

        self.stop = callback
        self.verbose = verbose
        self.want_queue = True

        self.gc = check_config3(False, self.verbose)
        #self.gc['BOUNDS'] = None
        #self.gc['GMT'] = None
        
        self.v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
        
        if self.gc['BOUNDS'] is not None:
            self.t_fields = ['string', 'string', 'string', 'string', 'string', 'string', 'string', 'string']
        else: self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name
        
    def _gather_from_queue(self):
        '''Gather geometries from a queue of [[datalist, layer], ...].'''

        while True:
            sm_args = self.dl_q.get()
            dl = sm_args[0]
            layer = sm_args[1]            
            if not self.stop():
                self._gather_from_datalist(dl, layer)
            
            self.dl_q.task_done()

    def _gather_from_datalist(self, dl, layer):
        '''gather geometries from datalist `dl` and append
        results to ogr `layer`. Load the datalist, generate
        a NUM-MSK grid, polygonize said NUM-MSK then union
        the polygon and add it to the output layer.'''

        this_datalist = datalist(dl[0], self.region, verbose = self.verbose, callback = self.stop)
        this_o_name = this_datalist._name        
        this_datalist._load_data()

        if len(this_datalist.datafiles) > 0:
            try:
                o_v_fields = entry[3]
                if len(o_v_fields) != 8:
                    o_v_fields = [this_datalist._name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
            except: o_v_fields = [this_datalist._name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
            if self.verbose: echo_msg('gathering geometries from datalist \033[1m{}\033[m...'.format(this_o_name))
            
            if self.gc['BOUNDS'] is not None:
                o_v_fields = ['\\"{}\\"'.format(x) if ' ' in x else x for x in o_v_fields]
                run_cmd('bounds -k {}/{} -n "{}" -gg --verbose >> {}.gmt\
                '.format(self.inc, self.dist_region.region_string, '|'.join(o_v_fields), layer), verbose = self.verbose, data_fun = this_datalist._dump_data)
            else:
                defn = layer.GetLayerDefn()
                if self.gc['GMT'] is None:
                    use_gmt = False
                else: use_gmt = True

                this_dem = waffles(this_datalist, self.region, i_inc = self.inc, o_fmt = 'GTiff', o_name = this_datalist._name, o_extend = self.extend, verbose = self.verbose)
                this_mask = this_dem.run('mask')

                if os.path.exists(this_mask) and not self.stop():
                    ## should use more unique name...crashes when 2 datalists have same name at same time...
                    remove_glob('{}_poly.*'.format(this_o_name))

                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_o_name))
                    tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_o_name), None, ogr.wkbMultiPolygon)
                    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

                    gdal_polygonize(this_mask, tmp_layer, verbose = self.verbose)

                    if len(tmp_layer) > 1:
                        out_feat = ogr_mask_union(tmp_layer, 'DN', defn, self.stop, verbose = self.verbose)
                        for i, f in enumerate(self.v_fields):
                            out_feat.SetField(f, o_v_fields[i])

                        layer.CreateFeature(out_feat)

                    tmp_ds = tmp_layer = out_feat = None
                    remove_glob('{}_poly.*'.format(this_o_name))
                    remove_glob('{}*'.format(this_mask[:-3]))

        echo_msg('gathered geometries from datalist \033[1m{}\033[m.'.format(this_o_name))

    def run(self, epsg = 4269):
        '''Run the spatial-metadata module and Geneate spatial metadata from the datalist
        specify the output project epsg.'''

        echo_msg('generating SPATIAL-METADATA for {}...'.format(self.datalist._name))
        dst_vec = '{}_sm.shp'.format(self.o_name)
        dst_layername = '{}_sm'.format(self.o_name)
        remove_glob('{}.*'.format(dst_layername))
        _prj_file('{}.prj'.format(dst_layername), epsg)

        self.datalist._load_datalists()

        if self.gc['BOUNDS'] is not None:
            with open('{}.gmt'.format(dst_layername), 'w') as gmtf:
                gmtf.write('# @VGMT1.0 @GMULTIPOLYGON\n# @N{}\n# @T{}\n# FEATURE_DATA\n'.format('|'.join(self.v_fields), '|'.join(self.t_fields)))
            for dl in self.datalist.datalists:
                self._gather_from_datalist(dl, dst_layername)
        else:
            ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
            if ds is not None:
                layer = ds.CreateLayer('{}'.format(dst_layername), None, ogr.wkbMultiPolygon)

                for i, f in enumerate(self.v_fields):
                    layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i]))

                for feature in layer:
                    layer.SetFeature(feature)

                if self.want_queue:
                    for _ in range(3):
                           t = threading.Thread(target = self._gather_from_queue, args = ())
                           t.daemon = True
                           t.start()

                    if len(self.datalist.datalists) > 0:
                        for dl in self.datalist.datalists:
                            self.dl_q.put([dl, layer])
                    else:
                        self.dl_q.put([[self.datalist._path, -1, 1], layer])

                    self.dl_q.join()
                else:
                    for dl in self.datalist.datalists:
                        self._gather_from_datalist(dl, layer)
                
        ds = layer = None
        if self.gc['BOUNDS'] is not None:
            run_cmd('ogr2ogr {} {}.gmt'.format(dst_vec, dst_layername, verbose = self.verbose))
            
        if not os.path.exists(dst_vec):
            dst_vec = None

        echo_msg('generated SPATIAL-METADATA for {}.'.format(self.datalist._name))
        return(dst_vec)

## =============================================================================
##
## uncertainty module - uncertainties.py
##
## datasource and interpolation uncertainty.
##
## =============================================================================

def err2coeff(my_data, coeff_guess = [0, 0.1, 0.2], dst_name = None):
    '''data is 2 col file with `err dist`'''

    from scipy import optimize
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText
    except:
        print('you need to install matplotlib to run uncertainty plots...')
    
    #try: 
    #    my_data = np.loadtxt(data, delimiter=' ')
    #except: sys.exit(2)

    error=my_data[:,0]
    distance=my_data[:,1]

    max_int_dist = np.max(distance)
    nbins = 10

    #coeff_guess=[0, 0.1, 0.2]
    n, _ = np.histogram(distance, bins = nbins)

    # want at least 2 values in each bin?
    while 0 or 1 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins = nbins)

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

    # fit
    
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    #plt.show()

    if dst_name is None:
        out_png = 'unc_best_fit.png'
    else: out_png = '{}_bf.png'.format(dst_name)
    plt.savefig(out_png)   # save the figure to file
    plt.close()

    #scatter

    plt.scatter(distance, error)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')

    if dst_name is None:
        out_png = 'unc_scatter.png'
    else: out_png = '{}_scatter.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

    return(out)

class uncertainty:

    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None, o_node = 'pixel', o_extend = 6, callback = lambda: False, verbose = False):

        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.node = o_node

        self.extend = int(o_extend)

        self.proc_region = self.region.buffer((self.extend * 2) * self.inc)
        self.dist_region = self.region.buffer(self.extend * self.inc)
        
        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.gc = check_config3(False, self.verbose)
        
        self.dem_mod = 'mbgrid'
        
        self.dem = { 
            'dem': None,
            'num': None,
            'msk': None,
            'prox': None,
            'int-unc': None,
        }

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name

        self.zones = ['bathy', 'bathy-topo', 'topo']
        self.region_info = {}
        self.sub_zones = {}
        self.trainers = []
        
    def run(self, dem_mod = 'mbgrid', dem = None, msk = None):
        
        i_dp = None
        opts = dem_mod.split(':')
        self.dem_mod = opts[0]
        self.mod_args = list(opts[1:])
        
        if dem is not None: self.dem['dem'] = dem
        if msk is not None: self.dem['msk'] = msk
                
        ## s_dp = self.source()
        ## v_dp = self.vdatum()
        i_dp = self.interpolation()

        ## self.combine(s_dp, v_dp, i_dp)
        
        return(i_dp)
        
    def set_or_make_dem(self):
        '''check if dem dict contains dems, otherwise generate them...'''

        if self.dem['dem'] is None:
            this_dem = waffles(self.datalist, self.region, str(self.inc), o_b_name = self.o_name)
            this_dem.o_fmt = 'GTiff'
            self.dem['dem'] = this_dem.run(self.dem_mod, self.mod_args)
            echo_msg('generated DEM {} using {}.'.format(self.dem['dem'], self.dem_mod))
        else: echo_msg('found DEM {}'.format(self.dem['dem']))

        if self.dem['msk'] is None:
            self.dem['msk'] = self.datalist.mask(region = self.dist_region.region, inc = self.inc, o_name = self.o_name)
            msk_dem = waffles(self.datalist, self.region, str(self.inc))
            msk_dem.o_fmt = 'GTiff'
            self.dem['msk'] = msk_dem.run('mask')
        else: echo_msg('found Data MASK {}'.format(self.dem['msk']))

        if self.dem['prox'] is None:
            self.dem['prox']  = '{}_prox.tif'.format(self.dem['msk'].split('.')[0]) 
            gdal_proximity(self.dem['msk'], self.dem['prox'])
            echo_msg('generated PROXIMITY grid {}.'.format(self.dem['prox']))
        else: echo_msg('found PROXIMITY grid {}.'.format(self.dem['prox']))
            
    def err_plot(self, dp, d_max):
        '''plot a numpy array of 'err dist' values and return the error coefficient.'''
        
        #dp = dp[dp[:,1]<self.region_info[self.o_name][4],:]
        dp = dp[dp[:,1] < d_max,:]
        dp = dp[dp[:,1] > 0,:]
        ec = err2coeff(dp, dst_name = self.o_name)
        echo_msg('error coefficient: {}'.format(ec))

        return(ec)

    def region_analysis(self):
        '''Analyze the input self.region and return infos.'''
        
        region_info = {}
        
        num_sum = gdal_sum(self.dem['msk'])
        gc = _infos(self.dem['msk'])
        g_max = float(gc['nx'] * gc['ny'])
        num_perc = (num_sum / g_max) * 100.
        prox_perc_95 = gdal_percentile(self.dem['prox'], 95)

        region_info[self.o_name] = [self.region.region, g_max, num_sum, num_perc, prox_perc_95]

        return(region_info)
    
    def tile_analysis(self):
        '''Anaylize the chunked regions and return infos about them.'''
        
        sub_count = 0
        sub_zones = {}
        
        for sc, sub_region in enumerate(self.sub_regions):

            gdal_cut(self.dem['msk'], _srcwin(self.dem['msk'], sub_region.region), 'tmp_msk.tif')
            gdal_cut(self.dem['dem'], _srcwin(self.dem['dem'], sub_region.region), 'tmp_dem.tif')

            s_gc = _infos('tmp_msk.tif')
            s_g_max = float(s_gc['nx'] * s_gc['ny'])
            s_sum = gdal_sum('tmp_msk.tif')
            s_perc = (s_sum / s_g_max) * 100
            
            s_dc = _infos('tmp_dem.tif', True)
            
            if s_dc['zmax'] < 0:
                zone = 'Bathy'
            elif s_dc['zmin'] > 0:
                zone = 'Topo'
            else: zone = 'BathyTopo'

            sub_zones[sc+1] = [sub_region.region, s_g_max, s_sum, s_perc, s_dc['zmin'], s_dc['zmax'], zone]

            remove_glob('tmp_msk.tif')
            remove_glob('tmp_dem.tif')
            
        return(sub_zones)

    def zone_analysis(self):
        '''Analyze uncertainty zones and select training tiles.'''
        
        trainers = []
        
        bathy_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Bathy']
        bathy_topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'BathyTopo']
        topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Topo']

        for z, tile_set in enumerate([bathy_tiles, bathy_topo_tiles, topo_tiles]):
            if len(tile_set) > 0:
                t_dens = np.array([x[3] for x in tile_set])
                t_50perc = np.percentile(t_dens, 50)
            else: t_50perc = 0.0
            if self.verbose: echo_msg('Minimum sampling for {} tiles: {}'.format(self.zones[z].upper(), t_50perc))

            t_trainers = [x for x in tile_set if x[3] > t_50perc]
            echo_msg('possible {} training zones: {}'.format(self.zones[z].upper(), len(t_trainers)))
            trainers.append(t_trainers)
                
        return(trainers)

    def zone_sort(self):
        '''sort training tiles by distance'''

        train_sorted = []
        
        for z, train in enumerate(self.trainers):
            train_d = []
            
            np.random.shuffle(train)

            while True:
                if len(train) == 0: break
                
                this_center = region(train[0][0]).center()
                train_d.append(train[0])
                train = train[1:]
                if len(train) == 0: break
                
                dsts = [hav_dst(this_center, region(x[0]).center()) for x in train]
                min_dst = np.percentile(dsts, 50)
                d_t = lambda t: hav_dst(this_center, region(t[0]).center()) > min_dst
                
                np.random.shuffle(train)
                train.sort(reverse=True, key=d_t)
                
            if self.verbose: echo_msg(' '.join([region(x[0]).gmt for x in train_d[:25]]))
            train_sorted.append(train_d)
            
        return(train_sorted)
        
    def split_sample(self, sub_regions, ss_samp, s_dp = None):
        '''perform split-sample analysis on the training tiles and return a list of `error distance` values'''
        
        if self.stop() or self.status !=0: return(s_dp)
        
        for n,sub_region in enumerate(sub_regions):
            #if self.verbose:
            echo_msg('processing sub-region ({}) {}'.format(n, sub_region))
            
            this_region = region(sub_region[0])
            o_xyz = '{}_{}.xyz'.format(self.o_name, n)

            if self.verbose:
                echo_msg('initial sampling density: {}'.format(sub_region[3]))
                echo_msg('desired sampling density: {}'.format(ss_samp))

            if sub_region[3] < ss_samp: ss_samp = None
            
            with open(o_xyz, 'w') as o_fh:
                gdal_dump(self.dem['dem'], o_fh, ' ', None, False, _srcwin(self.dem['dem'], this_region.buffer(20*self.inc).region), self.dem['msk'])

            if os.stat(o_xyz).st_size == 0:
                echo_msg('no data in sub-region...')
                self.status = -1
            else:
                ## use gdal instead
                s_inner, s_outer = gmtselect_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = self.verbose)

                if os.stat(s_inner).st_size != 0:
                    sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                else: sub_xyz = []
                ss_len = len(sub_xyz)

                if ss_samp is not None:
                    sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
                else: sx_cnt = 1
                sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                if self.verbose: echo_msg('withholding {} out of {} points for error sampling'.format(ss_len - sx_cnt, ss_len))

                np.random.shuffle(sub_xyz)
                np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                sub_datalist = datalist('sub_{}.datalist'.format(n), this_region, verbose = self.verbose)
                sub_datalist._append_entry([s_outer, 168, 1])
                sub_datalist._append_entry([sub_xyz_head, 168, 1])

                sub_surf = waffles(sub_datalist, this_region, i_inc = str(self.inc), o_name = 'sub_{}'.format(n), verbose = self.verbose)
                sub_dem = sub_surf.run(self.dem_mod, self.mod_args)

                if sub_dem is not None:
                    sub_msk = sub_datalist.mask(region = this_region.buffer(10*self.inc).region, inc = self.inc)
                    sub_prox = '{}_prox.tif'.format(sub_msk.split('.')[0])
                    gdal_proximity(sub_msk, sub_prox)
                    gdal_infos(sub_prox)
                    sub_xyd = gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                    sub_dp = gdal_query(sub_xyd, sub_prox, 'zg')
                else: sub_dp = None

                remove_glob(sub_xyz_head)
                remove_glob(sub_datalist._path)
                sub_xyz = None

                if s_dp is None:
                    s_dp = sub_dp
                else:
                    try:
                        s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                    except:
                        if self.verbose: echo_error_msg('found no error points...')
                        pass

            remove_glob(o_xyz)
            remove_glob('sub_{}*'.format(n))
            
        return(s_dp)
            
    def interpolation(self):
        '''calculate the interpolation uncertainty.'''

        dp = None
        sims = 10
        sim_loops = 1
        
        self.set_or_make_dem()
        if self.verbose: echo_msg(self.dem)

        self.region_info = self.region_analysis()
        if self.verbose:
            for x in self.region_info.keys():
                echo_msg('region: {}: {}'.format(x, self.region_info[x]))

        ## ==============================================
        ## chunk region into sub regions
        ## ==============================================
        chnk_lvl = 4        
        echo_msg('chunking region into sub-regions using chunk level {}...'.format(chnk_lvl))
        chnk_inc = int(chnk_lvl * self.region_info[self.o_name][4])
        self.sub_regions = self.region.chunk(self.inc, chnk_inc)
        #if self.verbose: print([x.region for x in self.sub_regions])
        echo_msg('chunked region into {} sub-regions.'.format(len(self.sub_regions)))

        echo_msg('analyzing {} sub-regions...'.format(len(self.sub_regions)))
        self.sub_zones = self.tile_analysis()
        if self.verbose:
            for x in self.sub_zones.keys():
                echo_msg('Sub-region {}: {}'.format(x, self.sub_zones[x]))
        
        echo_msg('running \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m...'.format(self.dem_mod))
                    
        s_dens = np.array([self.sub_zones[x][3] for x in self.sub_zones.keys()])
        s_5perc = np.percentile(s_dens, 5)
        s_dens = None
        echo_msg('Sampling density for region is: {:.12f}'.format(s_5perc))
                
        self.trainers = self.zone_analysis()
        trains = self.zone_sort()
        echo_msg('sorted training tiles.')
        ## generate shapefile from tain_h
        echo_msg('analyzed {} sub-regions.'.format(len(self.sub_regions)))
        
        for sim in range(0, sims):
            echo_msg('performing INTERPOLATION UNCERTAINTY simulation {} out of {}...'.format(sim + 1, sims))
            self.status = 0
            
            for s_l in range(0, sim_loops):
                for z, train in enumerate(trains):
                    train_h = train[:25]
                    dp = self.split_sample(train_h, s_5perc, dp)
            if len(dp) == 0:
                self.status = -1
                #break
            echo_msg('performed INTERPOLATION UNCERTAINTY simulation {} out of {}; {} error points accumulated.'.format(sim + 1, sims, len(dp)))

        echo_msg('ran \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m.'.format(self.dem_mod))
        
        if self.status == 0:

            if self.verbose: echo_msg('gathered {} error points'.format(len(dp)))
            np.savetxt('{}.err'.format(self.o_name), dp, '%f', ' ')
            ec = self.err_plot(dp[:50000000], self.region_info[self.o_name][4])
            
            ## ==============================================
            ## apply error coefficient to full proximity grid
            ## ==============================================

            echo_msg('applying coefficient to proximity grid')
            ## USE numpy instead
            
            math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_dst_unc.tif=gd+n-9999:GTiff\
            '.format(self.dem['prox'], ec[2], ec[1], 0, self.o_name)
            utils.run_cmd(math_cmd, self.verbose, self.verbose)
            echo_msg('applyed coefficient to proximity grid')
            
        return(dp)

    def source(self):
        '''source data uncertainty'''
        
        pass
            
## =============================================================================
##
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT and GDAL 
##
## =============================================================================

class waffles(object):
    '''Generate a Digital Elevation Model using one of the dem modules.
    DEM Modules include, `mbgrid`, `surface`, `num`, `mean`, `bathy`'''

    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None,\
                 o_node = 'pixel', o_fmt = 'GTiff', o_extend = 6, fltr = None, clip_ply = None, \
                 callback = lambda: False, verbose = False):

        self.status = 0
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.o_fmt = o_fmt
        self.extend = int(o_extend)
        self.node = o_node

        self.stop = callback
        self.verbose = verbose
        self.o_name = o_name

        self.proc_region = self.region.buffer((self.extend * self.inc) + self.inc * 2)
        self.dist_region = self.region.buffer(self.extend * self.inc)        
        
        if self.node == 'pixel':
            self.reg_str = '-r'
            self.gproc_region = self.proc_region
        else:
            self.reg_str = ''
            self.gproc_region = self.proc_region.buffer(.5 * self.inc)
            
        self.dem = '{}.tif'.format(self.o_name)
        self.gc = check_config3(False, self.verbose)

        self.fltr = fltr
        self.clip_ply = clip_ply
        self.datalist.want_weights = True
        
    def buffer_region(self, buff):

        self.region = self.region.buffer(buff)
        self.proc_region = self.region.buffer((self.extend * self.inc) + self.inc * 2)
        self.dist_region = self.region.buffer(self.extend * self.inc)
        self.datalist.region = self.proc_region
        
    def run(self, dem_mod = 'mbgrid', args = ()):
        '''Run the DEM module `dem_mod` using args `dem_mod_args`.'''

        args_d = {}
        try:
            for arg in args:
                p_arg = arg.split('=')
                args_d[p_arg[0]] = p_arg[1]
        except IndexError as e:
            echo_error_msg('invalid module option: {}, {}'.format(args, e))
            sys.exit(-1)

        ## ==============================================
        ## Generate the DEM using dem_mod and it's args
        ## ==============================================
        
        if self.status == 0:
            try:
                _dem_mods[dem_mod][0](self)(**args_d)
            except TypeError as e:
                echo_error_msg('{}'.format(e))
                self.status = -1
               
        ## ==============================================
        ## optionally filter the DEM
        ## ==============================================

        ## use dem_smooth python/gdal if GMT not avail...
        if self.fltr is not None:
            self.filter_dem(self.fltr)

        ## ==============================================
        ## optionally clip the DEM
        ## ==============================================

        if self.clip_ply is not None:
            clip_args = {}
            cp = self.clip_ply.split(':')
            clip_args['src_ply'] = cp[0]
            cargs = cp[1:]
            for arg in cargs:
                p_arg = arg.split('=')
                clip_args[p_arg[0]] = p_arg[1]

            self.clip_dem(**clip_args)
        
        ## ==============================================
        ## cut the DEM to output dst_region
        ## ==============================================

        self.cut_dem()
            
        ## ==============================================
        ## reformat DEM to final output format
        ## ==============================================
        
        if self.status != 0:
            self.dem = None
        else:
            if self.o_fmt != 'GTiff':
                orig_dem = self.dem
                if self.gc['GMT'] is not None:
                    self.dem = grd2gdal(self.dem, self.o_fmt)
                else: self.dem = gdal2gdal(self.dem, self.o_fmt)
                remove_glob(orig_dem)

        ## ==============================================
        ## set the grid meatadata
        ## grid-node no work.
        ## ==============================================

        ds = gdal.Open(self.dem, gdal.GA_Update)
        if ds is not None:
            md = ds.GetMetadata()
            if self.node == 'pixel':
                md['AREA_OR_POINT'] = 'Area'
            else: md['AREA_OR_POINT'] = 'Point'
            ds.SetMetadata(md)
            ds = None
        else: echo_error_msg('failed to set metadata')

        return(self.dem)

    def filter_dem(self, fltr_dist = 3):
        '''filter self.dem using either GMT grdfilter or numpy'''

        if self.gc['GMT'] is not None:
            status = grdfilter(self.dem, '{}_f.tif=gd+n-9999:GTiff'.format(self.o_name), dist = fltr_dist, verbose = self.verbose)
            if status == 0: os.rename('{}_f.tif'.format(self.o_name), self.dem)

    def cut_dem(self):
        out = gdal_cut(self.dem, _srcwin(self.dem, self.dist_region.region), 'tmp.tif')
        if out is not None: os.rename('tmp.tif', self.dem)
        
    def clip_dem(self, src_ply = None, invert = False):
        '''clip self.dem to either a src_ply polygon or the gsshg via GMT.
        use 'src_ply = None' to use gsshg.'''
        
        if src_ply is not None and os.path.exists(src_ply):
            if invert:
                gr_inv = '-i'
            else: gr_inv = ''
                
            gi = _infos(self.dem)
            gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'\
                     .format(gi['ndv'], gr_inv, os.path.basename(src_ply).split('.')[0], src_ply, self.dem)
            run_cmd(gr_cmd, verbose = self.verbose)
                
        else: ## use gsshg via GMT
            if invert:
                ns = '-N0/1/0/1/0'
            else: ns = '-N1/0/1/0/1'

            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V {} {}\
            '.format(self.inc, self.dist_region.gmt, ns, self.reg_str))
            out, self.status = run_cmd(dem_landmask_cmd, verbose = self.verbose)

            if self.status == 0:
                dem_landmask_cmd1 = ('gmt grdmath -V {} tmp_lm.grd MUL 0 NAN = {}_clp.tif=gd+n-9999:GTiff\
                '.format(self.dem, self.o_name))
                out, self.status = run_cmd(dem_landmask_cmd1, verbose = self.verbose)

                os.rename('{}_clp.tif'.format(self.o_name), self.dem)
                remove_glob('tmp_lm.grd')

    ## ==============================================
    ## run mbgrid on the datalist and generate the DEM
    ## note: mbgrid will cause popen to hang if stdout 
    ## is not cleared...should add output file to send to...
    ## ==============================================
                    
    def mbgrid(self, dist = '10/3', tension = 35):
        '''Generate a DEM with MBSystem's mbgrid program.'''

        import shutil
        
        if self.gc['MBGRID'] is None:
            echo_error_msg('MBSystem must be installed to use the MBGRID module')
            self.status = -1

        if self.gc['GMT'] is None:
            echo_error_msg('GMT must be installed to use the MBGRID module')
            self.status = -1
            
        if self.status == 0:

            self.datalist._archive(dirname = '.mb_tmp_datalist')
            self.datalist = self.datalist.archive_datalist
            
            out, self.status = run_mbgrid(self.datalist, self.gproc_region, self.inc, self.o_name, dist = dist, tension = tension, verbose = self.verbose)
            
            remove_glob('*.cmd')                
            shutil.rmtree('.mb_tmp_datalist')

    ## ==============================================
    ## Run GMT surface on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## ==============================================
    
    def surface(self, tension = .35, relaxation = 1.2, lower_limit = 'd', upper_limit = 'd'):
        '''Generate a DEM with GMT surface'''

        if self.gc['GMT'] is None:
            echo_error_msg('GMT must be installed to use the SURFACE module')
            self.status = -1

        dem_surf_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt surface -V {} -I{:.10f} -G{}=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} {}\
        '.format(self.proc_region.gmt, self.inc, self.reg_str, self.proc_region.gmt, self.inc, self.dem, tension, relaxation, lower_limit, upper_limit, self.reg_str))
        out, self.status = run_cmd(dem_surf_cmd, verbose = self.verbose, data_fun = self.datalist._dump_data)
        
    ## ==============================================
    ## Run GMT triangulate on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## output may have nodata values around the edges.
    ## ==============================================

    def triangulate(self):
        '''Generate a DEM with GMT surface'''

        if self.gc['GMT'] is None:
            echo_error_msg('GMT must be installed to use the TRIANGULATE module')
            self.status = -1
        
        dem_tri_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt triangulate {} -I{:.10f} -V -G{}=gd+n-9999:GTiff {}\
        '.format(self.proc_region.gmt, self.inc, self.reg_str, self.proc_region.gmt, self.inc, self.dem, self.reg_str))
        out, self.status = run_cmd(dem_tri_cmd, verbose = self.verbose, data_fun = self.datalist._dump_data)

    ## ==============================================
    ## Run GMT nearneighbor on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## output may have nodata values around the edges.
    ## ==============================================

    def nearneighbor(self, radius='6s'):
        '''Generate a DEM with GMT surface'''

        ## USE GDAL_GRID NEARNEIGHBOR IF GMT NOT AVAIL...
        if self.gc['GMT'] is None:
            echo_error_msg('GMT must be installed to use the NEARNEIGHBOR module')
            self.status = -1
        
        dem_nn_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}=gd+n-9999:GTiff {}\
        '.format(self.proc_region.gmt, self.inc, self.reg_str, self.proc_region.gmt, self.inc, radius, self.dem, self.reg_str))
        out, self.status = run_cmd(dem_nn_cmd, verbose = self.verbose, data_fun = self.datalist._dump_data)

    ## ==============================================
    ## Generate a NUM/MASK grid
    ##
    ## run GMT xyz2grd on the datalist to generate
    ## an uninterpolated grid
    ## see the -A switch in GMT xyz2grd for mode options
    ## 
    ## `[d|f|l|m|n|r|S|s|u|z]
    ## By  default  we  will calculate mean values if multiple entries
    ## fall on the same node. Use -A to change this behavior, except it
    ## is ignored if -Z is given.
    ## Append f or s to simply keep the first or last data point that
    ## was assigned to each node. Append l or u or d to find the lowest

    ## (minimum) or upper (maximum) value or the difference between the
    ## maximum and miminum value at each node, respectively. Append m
    ## or r or S to compute mean or  RMS  value  or standard deviation at
    ## each node, respectively. Append n to simply count the number of
    ## data points that were assigned to each node (this only requires
    ## two input columns x and y as z is not consulted).
    ## Append z to sum multiple values that belong to the same node.`
    ## ==============================================

    def num(self, mode = 'n'):
        '''Generate a num and num-msk grid with GMT'''
        self.datalist.want_weights = False
        
        ## add xyz2gdal if gmt not avail...
        if self.gc['GMT'] is not None:
            out, self.status = xyz2grd(self.datalist, self.dist_region, self.inc, '{}=gd+n-9999:GTiff'.format(self.dem), mode, self.node, verbose = self.verbose)
            
    def mask(self, use_gmt = True):
        '''Generate a num and num-msk grid'''

        if self.gc['GMT'] is None:
            self.datalist._load_data()
            self.dem = self.datalist.mask(region = self.gproc_region, inc = self.inc, o_name = self.o_name, o_fmt = 'GTiff')
        else:
            self.num('n')
            if self.status == 0:
                self.status = gmt_nan2zero(self.dem)
            else: self.dem = None

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Inverse Distance DEM
    ## ==============================================

    def invdst(self, power = 2.0, smoothing = 0.0, radius1 = 0.1, radius2 = 0.1, angle = 0.0, max_points = 0, min_points = 0, nodata = 0.0):
        '''Generate an inverse distance grid with GDAL'''

        gg_mod = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
                                 .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
        self.dem = run_gdal_grid(self.datalist, self.gproc_region, self.inc, self.o_name, gg_mod, self.node, self.gc['GMT'], self.verbose)
        
    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Moving Average DEM
    ## ==============================================
    
    def m_average(self, radius1 = 0.01, radius2 = 0.01, angle = 0.0, min_points = 0, nodata = 0.0):
        '''Generate a moving average grid with GDAL'''

        gg_mod = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'.format(radius1, radius2, angle, min_points, nodata)
        self.dem = run_gdal_grid(self.datalist, self.gproc_region, self.inc, self.o_name, gg_mod, self.node, self.gc['GMT'], self.verbose)

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## a linear triangulation DEM
    ## ==============================================
    
    def linear(self, radius = -1, nodata = 0.0):
        '''Generate a llinear triangulation grid with GDAL'''

        gg_mod = 'linear:radius={}:nodata={}'.format(radius, nodata)
        self.dem = run_gdal_grid(self.datalist, self.gproc_region, self.inc, self.o_name, gg_mod, self.node, self.gc['GMT'], self.verbose)
        
    ## ==============================================
    ## conversion-grid module:
    ## Generate a VDatum conversion grid
    ##
    ## Requires NOAA's VDatum version >= 4.x
    ## ==============================================

    ## update for gdal only...
    def vdatum(self, ivert = 'navd88', overt = 'mllw', region = '3'):
        '''generate a vertical datum transformation grid with vdatum'''

        if self.gc['VDATUM'] is None:
            echo_error_msg('NOAAs VDATUM must be installed to use the VDATUM module.')
            self.status = -1

        if self.status == 0:

            gdal_null('empty.tif', self.dist_region.region, 0.00083333, nodata = 0)
            with open('empty.xyz', 'w') as mt_xyz: gdal_dump('empty.tif', dst_xyz = mt_xyz, dump_nodata = True)

            this_vd = vdatum(vdatum_path = self.gc['VDATUM'], verbose = self.verbose)
            this_vd.ivert = ivert
            this_vd.overt = overt
            this_vd.run_vdatum('empty.xyz')

            if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
                empty_infos = gmtinfo('result/empty.xyz', self.verbose)

                if empty_infos[4] > 0:
                    ll_switch = '-Lld'
                else: ll_switch = '-Ll0'

                if empty_infos[5] > 0:
                    lu_switch = '-Lud'
                else: lu_switch = '-Lu0'

                gc = 'gmt blockmean result/empty.xyz -V -I{} {} {}\
                | gmt surface -I{} {} -G{} -V -T0 {} {} {}\
                '.format(self.inc, self.dist_region.gmt, reg_str, 
                         self.inc, self.dist_region.gmt, self.dem, ll_switch, lu_switch, reg_str)
                run_cmd(gc, verbose = self.verbose)
            else: self.status = -1

            try:
                remove_glob('empty.*')
                remove_glob('result/*')
                os.removedirs('result')
            except: pass
        else: self.dem = None

    def rem_rep(self):
        '''IBCAO `remove-replace` weighted gridding'''

        ## weight thresh - w <= 1 -> low; w >1 -> high ??
        
        ## Grid low-weights at self.inc / 3 w/ 'surface'
        ## Resample to self.inc to make BASE grid
        ## Grid high-weights at self.inc w/ nearneighbor (possibly invdst) to make HIGH grid
        ## Calculate differences between BASE and HIGH
        ## Grid differences w/ surface/mbgrid/invdst and clip to buffer-zone
        ## Apply differences to BASE grid

    def archive(self):
        self.datalist._archive()
        
    def spatial_metadata(self, inc = None, epsg = 4269):
        '''run the spatial-metadata module class'''

        self.sm = None
        if inc is None: inc = self.inc
        sm().run
        sm = spatial_metadata(self.datalist, self.region, i_inc = inc, o_name = self.o_name, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        sm.want_queue = True
        self.sm = sm.run(epsg)

    def uncertainty(self, dem_mod = 'mbgrid', dem = None, msk = None, prox = None):
        '''run the uncertainty module class'''
        
        unc = uncertainty(self.datalist, self.region, i_inc = self.inc, o_name = self.o_name, o_node = self.node, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        unc.run(dem_mod, dem, msk)
        
## =============================================================================
##
## Console Scripts
##
## Scripts are: 
## waffles
## waffles-datalists (datalists)
##
## =============================================================================

_dem_mods = {
    'mbgrid': [lambda x: x.mbgrid, '''Weighted SPLINE DEM via mbgrid
    \t\t\t< mbgrid:tension=35:dist=10/3:use_datalists=False >
    \t\t\t:tension=[0-100] - Spline tension.
    \t\t\t:dist=[value] - MBgrid -C switch (distance to fill nodata with spline)'''],
    'surface': [lambda x: x.surface, '''SPLINE DEM via GMT surface
    \t\t\t< surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
    \t\t\t:tension=[0-1] - Spline tension.'''],
    'triangulate': [lambda x: x.triangulate, '''TRIANGULATION DEM via GMT triangulate'''],
    'nearneighbor': [lambda x: x.nearneighbor, '''NEARNEIGHBOR DEM via GMT
    \t\t\t< nearneighbor:radius=6s >
    \t\t\t:radius=[value] - Nearest Neighbor search radius'''],
    'num': [lambda x: x.num, '''Uninterpolated DEM via GMT xyz2grd
    \t\t\t< num:mode=n >
    \t\t\t:mode=[key] - xyz2grd -A switch to specify mode of grid population.'''],
    'mask': [lambda x: x.mask, '''Data MASK grid'''],
    'invdst': [lambda x: x.invdst, '''INVERSE DISTANCE DEM via gdal_grid
    \t\t\t< invdst:power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1 >'''],
    'average': [lambda x: x.m_average, '''Moving AVERAGE DEM via gdal_grid
    \t\t\t< average:radius1=0.01:radius2=0.01 >'''],
    'linear': [lambda x: x.linear, '''LINEAR triangulation DEM via gdal_grid
    \t\t\t< linear:radius=-1 >
    \t\t\t:radius=[value] - Linear interpolation search radius.'''],
    'vdatum': [lambda x: x.vdatum, '''VDATUM transformation grid
    \t\t\t< vdatum:ivert=navd88:overt=mhw:region=3 >
    \t\t\t:ivert=[vdatum] - Input VDatum vertical datum.
    \t\t\t:overt=[vdatum] - Output VDatum vertical datum.
    \t\t\t:region=[0-10] - VDatum region (3 is CONUS)'''],
    'uncertainty': [lambda x: x.uncertainty, '''DEM UNCERTAINTY grid <beta>
    \t\t\t< uncertainty:dem_mod=mbgrid:dem=None:msk=None:prox=None >
    \t\t\t:dem_mod=[mbgrid/surface/invdst/triangulate/...] - Waffles DEM interpolation module
    \t\t\t:dem=[dem_path] - Pre-built DEM
    \t\t\t:msk=[msk_path] - Pre-built DEM data Mask grid
    \t\t\t:prox=[prox_path] - Pre-built DEM data proximity grid'''],
    'spatial-metadata': [lambda x: x.spatial_metadata, '''Datalist SPATIAL METADATA <beta>
    \t\t\t< spatial-metadata:inc=None >
    \t\t\t:inc=[increment] - Spatial metadata resolution [default uses -E value]'''],
}

def dem_mod_desc(x):
    fd = ['\033[1m{:18}\033[m{}'.format(key, x[key][-1]) for key in x]
    return('\n  '.join(fd))

_waffles_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

usage: {} [ -hprvCEFIOPRTX [ args ] ] module[:parameter=value]* ...

Modules and their options:
  {}

General Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf - or --, use the data region gathered from the DATALIST.
  -I, --datalist\tThe input DATALIST/DATAFILE:FMT.
\t\t\tThis can either be a waffles/mbsystem style datalist or a
\t\t\tsingle datafile which must include it's format; e.g. data.xyz:168
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired DEM MODULE and options.
\t\t\tsyntax is -M module:mod_opt=mod_val:mod_opt1=mod_val1:...
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION. [6]
  -T, --filter\t\tFILTER the output using a Cosine Arch filter at -T<dist(km)> search distance.
  -C, --clip\t\tCLIP the output to the clip polygon. [clip_ply.shp:invert=False]

  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year/module info to output BASENAME).
  -r, --grid-node\tuse grid-node registration, default is pixel-node

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Examples:
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 -V surface:tension=.7
 % {} -I input.datalist -E .3333333s -X 2 -R input_tiles_ply.shp -V -r -s -u mbgrid
 % {} -R-82.5/-82.25/26.75/27 -E0.0000925 vdatum:ivert=navd88:overt=mhw:region=3 -O ncei -p -r

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format(os.path.basename(sys.argv[0]), 
           _version, 
           os.path.basename(sys.argv[0]), 
           dem_mod_desc(_dem_mods),
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]))

def waffles_cli():
    status = 0
    i_region = None
    i_datalist = None
    i_inc = 0.0000925925
    these_regions = []
    stop_threads = False
    want_verbose = False
    want_archive = False
    want_prefix = False
    want_filter = False
    i_module = None
    mod_opts = {}
    o_bn = None
    o_fmt = 'GTiff'
    node_reg = 'pixel'
    o_extend = 6
    filter_dist = None
    clip_ply = None
    o_epsg = 4326

    argv = sys.argv
        
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

        elif arg == '--datalist' or arg == '-I':
            i_datalist = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-I':
            i_datalist = str(arg[2:])

        elif arg == '--increment' or arg == '-E':
            i_inc = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-E':
            i_inc = arg[2:]

        elif arg == '--module' or arg == '-M':
            i_module = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-M':
            i_module = str(arg[2:])

        elif arg == '--format' or arg == '-F':
            o_fmt = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-F':
            o_fmt = arg[2:]

        elif arg == '--output-name' or arg == '-O':
            o_bn = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-O':
            o_bn = str(arg[2:])

        elif arg == '--epsg' or arg == '-P':
            o_epsg = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-P':
            o_epsg = str(arg[2:])
        
        elif arg == '--extend' or arg == '-X':
            o_extend = int(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-X':
            o_extend = int(arg[2:])

        elif arg == '--filter' or arg == '-T':
            filter_dist = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-T':
            filter_dist = arg[2:]

        elif arg == '--clip' or arg == '-C':
            clip_ply = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-C':
            clip_ply = str(arg[2:])

        elif arg == '--prefix' or arg == '-p':
            want_prefix = True
            
        elif arg == '--grid-node' or arg == '-r':
            node_reg = 'grid'

        elif arg == '--archive' or arg == '-a':
            want_archive = True

        elif arg == '--help' or arg == '-h':
            print(_waffles_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, _license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_waffles_usage)
            sys.exit(0)

        else:
            pass
            # opts = arg.split(':')
            # if opts[0] in _dem_mods.keys():
            #     mod_opts[opts[0]] = list(opts[1:])
            # else: echo_error_msg('invalid module name `{}`'.format(opts[0]))

        i = i + 1

    gc =  check_config3(want_verbose)

    if gc['platform'] == 'win32':
        out, status = run_cmd('py3_env', verbose = True)

    if i_module is not None:
        opts = i_module.split(':')
        if opts[0] in _dem_mods.keys():
            mod_opts[opts[0]] = list(opts[1:])
        else: echo_error_msg('invalid module name `{}`'.format(opts[0]))
        
        for key in mod_opts.keys():
            mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]

    if len(mod_opts) == 0:
        echo_msg('try waffles --help to see full command-line usage')
        while True:
            try:
                opts = raw_input('Waffles Modules:\n  {}\nchoose a waffles module\n<waffles>: '.format(dem_mod_desc(_dem_mods))).split(':')
            except NameError:
                opts = input('Waffles Modules:\n  {}\nchoose a waffles module\n<waffles>: '.format(dem_mod_desc(_dem_mods))).split(':')
            except KeyboardInterrupt as e:
                echo_error_msg('user breakage')
                sys.exit(-1)
                
            if opts[0] in _dem_mods.keys():
                mod_opts[opts[0]] = list(opts[1:])
                break
            
    ## ==============================================
    ## process input region(s) and loop
    ## The input region can either be a GMT region-string;
    ## e.g. w/e/s/n or an OGR compatible vector file, in
    ## which case each feature will be used as an input
    ## region.
    ## ==============================================

    if i_region is None:
        try:
            i_region = raw_input('''Enter region of interest;
This can either be a GMT-style region ( -R xmin/xmax/ymin/ymax ) or an OGR-compatible vector file with regional polygons. 
If a vector file is supplied it will search each region found therein.
If - or -- or ENTER, use the data region gathered from the DATALIST.
<waffles>: ''')
        except NameError:
            i_region = input('''Enter region of interest;
This can either be a GMT-style region ( -R xmin/xmax/ymin/ymax ) or an OGR-compatible vector file with regional polygons. 
If a vector file is supplied it will search each region found therein.
If - or -- or ENTER, use the data region gathered from the DATALIST.
<waffles>: ''')
        except KeyboardInterrupt as e:
            echo_error_msg('user breakage')
            sys.exit(-1)
        if len(i_region) == 0:
            these_regions = [None]
        else: these_regions = i_region

    try: 
        these_regions = [region(i_region)]
    except: these_regions = [region(x) for x in _ogr_extents(i_region)]

    if len(these_regions) == 0:
        these_regions = [None]

    for this_region in these_regions:
        if this_region is not None:
            if not this_region._valid:
                echo_error_msg('failed to load region(s)')
                sys.exit(-1)
    if want_verbose: echo_msg('loaded \033[1m{}\033[m region(s).'.format(len(these_regions)))
            
    ## ==============================================
    ## process the input increment
    ## ==============================================
        
    try:
        i_inc = gmt_inc2inc(i_inc)
    except:
        while True:
            try:
                i_inc = raw_input('Enter gridding inncrement \nGridding CELL-SIZE in native units or GMT-style increments.\n<waffles>: ')
            except NameError:
                i_inc = input('Enter gridding inncrement \nGridding CELL-SIZE in native units or GMT-style increments.\n<waffles>: ')
            except KeyboardInterrupt as e:
                echo_error_msg('user breakage')
                sys.exit(-1)
            try:
                i_inc = gmt_inc2inc(i_inc)
                break
            except: echo_error_msg('the increment value should be a number')

    echo_msg('Increment: {}'.format(i_inc))
        
    for rn, this_region in enumerate(these_regions):

        status = 0
        
        ## ==============================================
        ## Load the input datalist
        ## the input datalist can either be a waffles or
        ## mbsystem style datalist or an individual data
        ## file, which can be either XYZ or a GDAL raster.
        ## if a data file, append the format and weight:
        ## input.dat:168:1
        ## ==============================================

        if stop_threads or status != 0: break

        if i_datalist is not None:
            this_datalist = datalist(i_datalist, this_region, i_inc = i_inc, verbose = want_verbose, callback = lambda: stop_threads)
            if not this_datalist._valid_p():
                echo_error_msg('invalid datalist')
                status = -1
            else: echo_msg('Datalist: {}'.format(this_datalist._path_basename))
        else:
            if 'vdatum' not in mod_opts.keys():
                while True:
                    try:
                        i_datalist = raw_input('''Enter datalist/file path: 
This can either be a waffles/mbsystem style datalist or a
single datafile which must include it's format; e.g. data.xyz:168
<waffles>: ''')
                    except NameError:
                        i_datalist = input('''Enter datalist/file path: 
This can either be a waffles/mbsystem style datalist or a
single datafile which must include it's format; e.g. data.xyz:168
<waffles>: ''')
                    except KeyboardInterrupt as e:
                        echo_error_msg('user breakage')
                        sys.exit(-1)
                    if os.path.exists(i_datalist):
                        this_datalist = datalist(i_datalist, this_region, i_inc = i_inc, verbose = want_verbose, callback = lambda: stop_threads)
                        if this_datalist._valid_p(): break
                    else: echo_error_msg('path doesn\'t exist {}'.format(i_datalist))
                echo_msg('Datalist: {}'.format(this_datalist._path_basename))
        
        ## ==============================================
        ## gather the region from the data/datalist if no
        ## user supplied region is available...
        ## ==============================================
        
        if this_region is None and this_datalist is not None:
            echo_msg('no region supplied, gathering from datalist...')
            this_datalist.gather_region()
            this_region = this_datalist.region

        if this_region is None and this_datalist is None:
            echo_error_msg('no region or datalist; aborting...')
            status = -1
            break

        echo_msg('Region: {}'.format(this_region.region))
        extent = this_region.buffer(o_extend * i_inc).region
        xcount = int((extent[3] - extent[2]) / i_inc) + 1
        ycount = int((extent[1] - extent[0]) / i_inc) + 1
        echo_msg('Output Grid Size: {}/{}'.format(xcount, ycount))
        
        ## ==============================================
        ## buffer the region for the datalist in order to
        ## include data at the edges of grid
        ## ==============================================
        
        this_datalist.region = this_datalist.region.buffer(10 * i_inc)

        ## ==============================================
        ## set the output file name
        ## ==============================================
        
        if o_bn is None:
            if this_datalist is None:
                o_name = 'waffles'
            else: o_name = this_datalist._path_basename.split('.')[0]
        else: o_name = os.path.join(os.path.dirname(o_bn), os.path.basename(o_bn).split('.')[0])
        
        if want_prefix:
            str_inc = inc2str_inc(i_inc)
            o_name = '{}{}_{}_{}'.format(o_name, str_inc, this_region.fn, this_year())

        echo_msg('Output Name: {}'.format(o_name))
            
        ## ==============================================
        ## Run the DEM module(s)
        ## ==============================================
        
        for dem_mod in mod_opts.keys():

            echo_msg('Module: \033[1m{}\033[m'.format(dem_mod))
            args = tuple(mod_opts[dem_mod])
            echo_msg('Module Options: {}'.format(args))

            ## ==============================================
            ## most modules need a datalist and a region
            ## and an increment to be valid, except vdatum doesn't
            ## need a datalist...
            ## ==============================================
            
            if this_datalist is None:
                if dem_mod.lower() != 'vdatum':
                    echo_error_msg('a valid datalist is needed to run module: {}...'.format(dem_mod))
                    break
                else:
                    echo_error_msg('no datalist; aborting...')
                    status = -1
                    break
            
            echo_msg('running geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region))
            dl = waffles(this_datalist, this_region, i_inc = i_inc, o_name = o_name, o_node = node_reg,\
                         o_fmt = o_fmt, o_extend = o_extend, fltr = filter_dist, clip_ply = clip_ply, callback = lambda: stop_threads, verbose = want_verbose)
            t = threading.Thread(target = dl.run, args = (dem_mod, args))
            dl.want_filter = want_filter
            
            try:
                t.start()
                while True:
                    time.sleep(1)
                    if not t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit):
                echo_msg('stopping all threads')
                dl.status = -1
                stop_threads = True

            t.join()

            if dl.dem is not None:
                if os.path.exists(dl.dem):
                    _set_gdal_epsg(dl.dem, o_epsg)

            if dem_mod.lower() == 'spatial-metadata': dl.dem = dl.sm
                
            if dl.status == 0:
                echo_msg('ran geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m and produced {}...\
                '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string, dl.dem))
            else:
                echo_error_msg('failed to run geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
                '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string))

## =============================================================================
##
## waffles-datalists Mainline - run datalists from console.
##
## =============================================================================

datalists_usage = '''{} ({}): Process and generate datalists

usage: {} [ -dghirsvEFOPR [ args ] ] DATALIST ...

Options:
  -R, --region\t\tSpecifies the desired REGION;
  -E, --increment\tThe CELL-SIZE in native units.
  -O, --output-name\tThe OUTPUT file name.
  -P, --epsg\t\tSpecify the EPSG of the DATALIST.
  -F, --format\t\tOnly process FORMAT data type.

  -a, --archive\t\tARCHIVE the data from the DATALIST.
  -g, --glob\t\tGLOB FORMAT data into the DATALIST.
  -i, --info-file\tGenerate INF files for the data in the DATALIST
  -l, --list\t\tLIST the datafiles from the DATALIST.
  -m, --mask\t\tMASK the datafiles from the DATALIST.
  -r, --region-info\tReturn the full REGION of the DATALIST.
  -s, --spatial-md\tGenerate SPATIAL-METADATA from the DATALIST.
  -w, --weights\t\toutput weights along with each datalist.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} my_data.datalist -R -90/-89/30/31 -g -i
 % {} my_data.datalist -R -90/-89/30/31 -E.000925 -O my_data_spatial -s

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def datalists_cli():

    status = 0
    i_datalist = None
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
    dl_fmt = None
    #dl_fmts = []
    
    argv = sys.argv
        
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
            print(datalists_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, _license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(datalists_usage)
            sys.exit(0)

        else: 
            i_datalist = arg

        i = i + 1

    if i_datalist is None:
        print(datalists_usage)
        sys.exit(1)
        
    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if i_region is None:
        these_regions = [None]
    else:
        try: 
            these_regions = [region(i_region)]
        except: these_regions = [regions.region('/'.join([str(r) for r in x])) for x in _ogr_extents(i_region)]

    if len(these_regions) == 0:
        these_regions = [None]

    for this_region in these_regions:
        if this_region is not None:
            if not this_region._valid:
                status = -1
    if want_verbose: echo_msg('loaded \033[1m{}\033[m region(s).'.format(len(these_regions)))
    
    if status == -1:
        echo_error_msg('failed to load region(s)')
        print(datalists_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):
        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        this_datalist = datalist(i_datalist, this_region, dl_fmt, verbose = want_verbose)
        this_datalist.want_weights = want_weights
        
        if want_glob:
            if dl_fmt is None:
                dl_fmts = _dl_fmts.keys()
            else: dl_fmts = [dl_fmt]
            for key in dl_fmts:
                for f in _dl_fmts[key][-1]:
                    globs = glob.glob('*.{}'.format(f))
                    [this_datalist._append_entry([x, key, 1]) for x in globs]

        if this_datalist._valid_p():

            #this_datalist.to_ogr(block=i_inc)            
            if want_inf: this_datalist._gen_inf()
            elif want_region:
                this_datalist.gather_region()
                sys.stdout.write('{}\n'.format(this_datalist.region.gmt))
            elif want_sm:
                import metadata

                if this_datalist.region is None:
                    this_datalist.gather_region()

                sm = spatial_metadata(this_datalist, this_datalist.region, i_inc = i_inc, o_name = o_bn, callback = lambda: False, verbose = want_verbose)
                sm.want_queue = False
                sm.run()
            elif want_mask: this_datalist.mask()
            elif want_list: this_datalist._proc_data()
            elif want_archive: this_datalist._archive()
            else: this_datalist._dump_data()

if __name__ == '__main__':
    waffles_cli()

### End
