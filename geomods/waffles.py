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
## surface (GMT), triangulate (GMT), nearneighbor (GMT)
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
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
## or via `datalists -i`
##
## if 'region' is specified, will only process data that falls within
## the given region
##
### TODO:
## Add uncertainty functions and modules
## Add remove/replace module
## Add source uncertainty to uncertainty module
## Add LAS/LAZ support to datalits
## Merge with fetch and allow for http datatype in datalists?
##
### Code:
import sys
import os
import io
import time
import glob
import math
import subprocess
import copy
## ==============================================
## import gdal, etc.
## ==============================================
import numpy as np
import json
import gdal
import ogr
import osr

## ==============================================
## General utility functions - utils.py
## ==============================================
_version = '0.5.2'

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    returns a str representation of float(inc)'''
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

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
    
    globs = glob.glob(glob_str)
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

def echo_msg2(msg, prefix = 'waffles'):
    '''echo `msg` to stderr using `prefix`
    >> echo_msg2('message', 'test')
    test: message'''
    
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: {}\n'.format(prefix, msg))

## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))

## ==============================================
## regions - regions are a bounding box list:
## [w, e, s, n]
## -- regions.py --
## ==============================================
def region_valid_p(region):
    '''return True if `region` [xmin, xmax, ymin, ymax] appears to be valid'''
    
    if region is not None:
        if region[0] < region[1] and region[2] < region[3]: return(True)
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
    else: return(False)

def region_format(region, t = 'gmt'):
    '''format region to string, defined by `t`
    t = 'str': xmin/xmax/ymin/ymax
    t = 'gmt': -Rxmin/xmax/ymin/ymax
    t = 'bbox': xmin,ymin,xmax,ymax
    t = 'te': xmin ymin xmax ymax
    t = 'fn': ymax_xmin

    returns the formatted region as str'''
    
    if t == 'str': return('/'.join([str(x) for x in region]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'te': return(' '.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'fn':
        if region[3] < 0: ns = 's'
        else: ns = 'n'
        if region[0] > 0: ew = 'e'
        else: ew = 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100) % 100), 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100) % 100)))

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
            if geo_x_o < region[0]: geo_x_0 = region[0]
            o_chunks.append([geo_x_o, geo_x_t, geo_y_o, geo_y_t])
        
            if y_chunk <= ycount:
                y_chunk += n_chunk
                i_chunk += 1
            else: break
        if x_chunk <= xcount:
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
        echo_msg(' '.join([region_format(x[0], 'gmt') for x in train_d[:25]]))
        train_sorted.append(train_d)
    return(train_sorted)

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
    'delim': 'space',
    'result_dir': 'result',
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
    
def run_vdatum(src_fn, vd_config = _vd_config):
    '''run vdatum on src_fn which is an XYZ file
    use vd_config to set vdatum parameters.

    returns [command-output, command-return-code]'''
    
    if vd_config['jar'] is None: vd_config['jar'] = vdatum_locate_jar()[0]
    if vd_config['jar'] is not None:
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}\
        '.format(vd_config['ihorz'], vd_config['ivert'], vd_config['ohorz'], vd_config['overt'], \
                 vd_config['delim'], src_fn, vd_config['result_dir'], vd_config['region'])
        #out, status = run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, True)
        return(run_cmd('java -jar {} {}'.format(vd_config['jar'], vdc), verbose = True))
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
    '.format(src_dem, o_name, src_dem))
    out, status = run_cmd(slope_cmd0, verbose = verbose)
    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}\
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

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):
    '''run the MBSystem command `mbgrid` given a datalist, region and increment.
    the datalist should be an MBSystem-style datalist; if using a waffles datalist, 
    it should be converted first (using datalist_archive()).

    returns [mbgrid-output, mbgrid-return-code]'''
    
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1: dist = dist + '/2'
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.10f}/{:.10f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} > mb_proc.txt \
    '.format(datalist, region_format(region, 'gmt'), inc, inc, dst_name, dist, tension, e_switch))
    return(run_cmd(mbgrid_cmd, verbose = verbose))

## ==============================================
## GDAL Wrappers and Functions - gdalfun.py
## ==============================================
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.UseExceptions()
_gdal_progress = gdal.TermProgress
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
    
def gdal_set_epsg(src_gdal, epsg = 4326):
    '''set the projection of src_gdal to epsg

    returns status-code (0 == success)'''
    
    ds = gdal.Open(src_gdal, gdal.GA_Update)
    if ds is not None:
        ds.SetProjection(gdal_sr_wkt(int(epsg)))
        ds = None
        return(0)
    else: return(None)

def gdal_set_nodata(src_gdal, nodata = -9999):
    '''set the nodata value of src_gdal

    returns 0'''
    
    ds = gdal.Open(src_gdal, gdal.GA_Update)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(float(-9999))
    ds = None
    return(0)

def gdal_infos(src_gdal, scan = False):
    '''scan src_gdal and gather region info.

    returns region dict.'''
    
    if os.path.exists(src_gdal):
        ds = gdal.Open(src_gdal)
        if ds is not None:
            dsc = gdal_gather_infos(ds)
            if scan:
                t = ds.ReadAsArray()
                t = t.astype(float)
                t[t == dsc['ndv']] = np.nan
                t = t[~np.isnan(t)]
                if len(t) > 0:
                    zmin = np.min(t[~np.isnan(t)])
                    zmax = np.max(t[~np.isnan(t)])
                else:
                    zmin = np.nan
                    zmax = np.nan
                dsc['zmin'] = zmin
                dsc['zmax'] = zmax
            ds = None
            return(dsc)
        else: return(None)
    else: return(None)

def gdal_gather_infos(src_ds):
    '''gather information from `src_ds` GDAL dataset

    returns gdal_config dict.'''
    
    ds_config = {
        'nx': src_ds.RasterXSize,
        'ny': src_ds.RasterYSize,
        'nb':src_ds.RasterCount,
        'geoT': src_ds.GetGeoTransform(),
        'proj': src_ds.GetProjectionRef(),
        'dt': src_ds.GetRasterBand(1).DataType,
        'dtn': gdal.GetDataTypeName(src_ds.GetRasterBand(1).DataType),
        'ndv': src_ds.GetRasterBand(1).GetNoDataValue(),
        'fmt': src_ds.GetDriver().ShortName,
    }
    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
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
    if os.path.exists(dst_gdal): driver.Delete(dst_gdal)
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
    
    gi = gdal_infos(src_gdal)
    if gi is not None and src_ply is not None:
        gr_inv = '-i' if invert else ''
        gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'\
                 .format(gi['ndv'], gr_inv, os.path.basename(src_ply).split('.')[0], src_ply, src_gdal)
        return(run_cmd(gr_cmd, verbose = True))
    else: return(None)

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
                    #c = con_dec(math.fabs(d / (g + 0.00000001) * 100), 2)
                    #s = con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    #d = con_dec(d, 4)
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    xyzl.append(np.array(outs, dtype = ds_config['dtn']))
        dsband = ds = None
        out_array = np.array(xyzl, dtype = ds_config['dtn'])
    return(out_array)
    
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
        
def gdal_gdal2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326):
    '''convert the gdal file to gdal using gdal

    return output-gdal-fn'''
    
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
        gdal2gdal_cmd = ('gdal_translate {} {} -f {}\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(gdal2gdal_cmd, verbose = True)
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
    '''convert a gdal geo-tranform to a region [xmin, xmax, ymin, ymax]'''
    
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + geoT[1] * ds_config['nx'], geoT[3] + geoT[5] * ds_config['ny'], geoT[3]])

def gdal_region2gt(region, inc):
    '''return a count info and a gdal geotransform based on extent and cellsize
    output is a list (xcount, ycount, geot)'''
    
    ysize = region[3] - region[2]
    xsize = region[1] - region[0]
    xcount = int(round(xsize / inc)) + 1
    ycount = int(round(ysize / inc)) + 1
    dst_gt = (region[0], inc, 0, region[3], 0, (inc * -1.))
    return(xcount, ycount, dst_gt)

def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.'''
    
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
    out_feat.SetGeometry(multi)
    union = multi = None
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
    '''convert coords to Wkt'''
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[1], coord[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def gdal_region2geom(region):
    '''convert an extent [west, east, south, north] to an OGR geometry'''
    
    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    geom = ogr.CreateGeometryFromWkt(gdal_create_polygon(eg))
    return(geom)

def gdal_region(src_ds):
    '''return the extent of the src_fn gdal file.'''
    
    ds_config = gdal_gather_infos(src_ds)
    return(gdal_gt2region(ds_config))

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''convert a geographic x,y value to a pixel location of geoTransform'''
    
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
    return(int(pixel_x+.5), int(pixel_y+.5))

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''convert a pixel location to geographic coordinates given geoTransform'''
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    '''apply geotransform to in_x,in_y'''
    
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]
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
    output the appropriate gdal srcwin.'''
    
    ds_config = gdal_gather_infos(src_ds)
    this_origin = _geo2pixel(region[0], region[3], ds_config['geoT'])
    this_origin = [0 if x < 0 else x for x in this_origin]
    this_end = _geo2pixel(region[1], region[2], ds_config['geoT'])
    this_size = (int(this_end[0] - this_origin[0]), int(this_end[1] - this_origin[1]))
    this_size = [0 if x < 0 else x for x in this_size]
    if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
    if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
    return(this_origin[0], this_origin[1], this_size[0], this_size[1])

def xyz2gdal_ds(src_xyz, dst_ogr):
    '''Make a point vector OGR DataSet Object from src_xyz'''
    
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
        #print(this_xyz)
        #sys.exit()
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        f.SetField(0, x)
        f.SetField(1, y)
        f.SetField(2, z)
        wkt = 'POINT(%.8f %.8f %.10f)' % (x,y,z)
        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)
    return(ds)

def gdal_xyz2gdal(src_xyz, dst_gdal, region, inc, dst_format='GTiff', mode='n'):
    '''Create a GDAL supported grid from xyz data 
    `mode` of `n` generates a num grid
    `mode` of `m` generates a mean grid
    `mode` of `k` generates a mask grid'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    if mode == 'm':
        sumArray = np.zeros((ycount, xcount))
        gdt = gdal.GDT_Float32
    else: gdt = gdal.GDT_Int32
    ptArray = np.zeros((ycount, xcount))
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(4326), gdt, -9999, dst_format)
    echo_msg('gridding data')
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                if mode == 'm': sumArray[ypos, xpos] += z
                if mode == 'n' or mode == 'm': ptArray[ypos, xpos] += 1
                else: ptArray[ypos, xpos] = 1
    if mode == 'm':
        ptArray[ptArray == 0] = np.nan
        outarray = sumArray / ptArray
    elif mode == 'n': outarray = ptArray
    else: outarray = ptArray
    outarray[np.isnan(outarray)] = -9999
    return(gdal_write(outarray, dst_gdal, ds_config))

def gdal_xyz_mask(src_xyz, dst_gdal, region, inc, dst_format='GTiff', epsg = 4326):
    '''Create a num grid mask of xyz data. The output grid
    will contain 1 where data exists and 0 where no data exists.'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    ptArray = np.zeros((ycount, xcount))
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(epsg), gdal.GDT_Int32, -9999, 'GTiff')
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                except: pass
    return(gdal_write(ptArray, dst_gdal, ds_config))

def np_gaussian_blur(in_array, size):
    '''blur an array using fftconvolve from scipy.signal
    size is the blurring scale-factor.'''
    
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
    only smooth bathymetry (sub-zero)'''
    
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

def gdal_sample_inc(src_grd, inc = 1, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node grid registration.'''
    
    out, status = run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te -R{} -r -Gtmp.tif=gd+n-9999:GTiff'.format(inc, inc, src_grd, src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)
        
def gdal_polygonize(src_gdal, dst_layer):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''
    
    ds = gdal.Open('{}'.format(src_gdal))
    ds_arr = ds.GetRasterBand(1)
    echo_msg('polygonizing {}'.format(src_gdal))
    gdal.Polygonize(ds_arr, None, dst_layer, 0, callback = _gdal_progress)
    ds = ds_arr = None
    return(0, 0)

def gdal_chunks(src_fn, n_chunk = 10):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.'''
    
    band_nums = []
    o_chunks = []
    if band_nums == []: band_nums = [1]
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk
    
    src_ds = gdal.Open(src_fn)
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
                dst_gt = [this_geo_x_origin, gt[1], 0.0, this_geo_y_origin, 0.0, gt[5]]
                
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

## ==============================================
## inf files (data info) inf.py
## mbsystem/gmt/waffles infos
## ==============================================
def inf_generate(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    
    return(inf_entry([data_path, data_fmt]))

def inf_parse(src_inf):
    '''parse an inf file (mbsystem or gmt)
    
    returns region: [xmin, xmax, ymin, ymax]'''
    
    minmax = [0, 0, 0, 0]
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
    return([float(x) for x in minmax])

def inf_entry(src_entry):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file or the 
    result of `gmt gmtinfo file.xyz -C`, which is a 6 column line 
    with minmax info, etc.

    returns the region of the inf file.'''
    
    minmax = None
    if os.path.exists(src_entry[0]):
        path_i = src_entry[0] + '.inf'
        if not os.path.exists(path_i):
            minmax = _dl_inf_h[src_entry[1]](src_entry)
        else: minmax = inf_parse(path_i)
    if not region_valid_p(minmax): minmax = None
    return(minmax)

## ==============================================
## gdal processing (datalist fmt:200)
## ==============================================
def gdal_parse(src_ds, dump_nodata = False, srcwin = None, mask = None):
    '''send the data from gdal file src_gdal to dst_xyz port (first band only)
    optionally mask the output with `mask`.'''
    
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_gather_infos(src_ds)
    gt = ds_config['geoT']
    msk_band = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    nodata = ['{:g}'.format(-9999), 'nan']
    if band.GetNoDataValue() is not None: nodata.append('{:g}'.format(band.GetNoDataValue()))
    if dump_nodata: nodata = []
    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        if msk_band is not None:
            msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            band_data[msk_data==0]=-9999
        band_data = np.reshape(band_data, (srcwin[2], ))
        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            geo_x, geo_y = _pixel2geo(x, y, gt)
            z = band_data[x_i]
            line = [geo_x, geo_y, z]
            if '{:g}'.format(z) not in nodata:
                yield(line)
                
def gdal_inf(src_ds):
    '''generate an info (.inf) file from a src_gdal file using gdal'''
    
    minmax = gdal_region(src_ds)
    with open('{}.inf'.format(src_ds.GetDescription()), 'w') as inf:
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)
                    
def gdal_inf_entry(entry):
    ds = gdal.Open(entry[0])
    minmax = gdal_inf(ds)
    ds = None
    return(minmax)

def gdal_yield_entry(entry, region = None, verbose = False):
    ds = gdal.Open(entry[0])
    if region is not None:
        srcwin = gdal_srcwin(ds, region)
    else: srcwin = None
    for xyz in gdal_parse(ds, dump_nodata = False, srcwin = srcwin):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)
    ds = None
    
def gdal_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False):
    for xyz in gdal_yield_entry(entry, region, verbose):
        xyz_line(xyz, dst_port)
        
## ==============================================
## xyz processing (datalists fmt:168)
## ==============================================
_xyz_config = {'delim':' ', 'xpos': 0, 'ypos': 1, 'zpos': 2}
_known_delims = [',', ' ', '\t', '/', ':']

def xyz_parse(src_xyz):
    '''xyz file parsing generator
    `src_xyz` is a file object or list of xyz data.

    yields each xyz line as a list [x, y, z, ...]'''
    
    this_delim = None
    for xyz in src_xyz:
        this_line = xyz.strip()
        if this_delim is None:
            for delim in _known_delims:
                this_xyz = this_line.split(delim)
                if len(this_xyz) > 1: #break
                    this_delim = delim
                    break
        else: this_xyz = this_line.split(this_delim)
        yield([float(x) for x in this_xyz])

def xyz2py(src_xyz):
    '''return src_xyz as a python list'''
    
    xyzpy = []
    return([xyzpy.append(xyz) for xyz in xyz_parse(src_xyz)])

def xyz_block(src_xyz, region, inc, dst_xyz = sys.stdout, weights = False):
    '''block the src_xyz data to the mean block value

    yields the xyz value for each block with data'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    avg_arr = np.empty((ycount, xcount), object)
    avg_arr.fill([])
    if weights:
        wt_arr = np.empty((ycount, xcount), object)
        wt_arr.fill([])
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if weights:
            w = this_xyz[3]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                avg_arr[ypos, xpos] = avg_arr[ypos, xpos] + [z]
                if weights:
                    wt_arr[ypos, xpos] = wt_arr[ypos, xpos] + [w]
    for y in range(0, ycount):
        for x in range(0, xcount):
            if len(avg_arr[y, x]) > 0:
                geo_x, geo_y = _pixel2geo(x, y, dst_gt)
                if weights:
                    z = np.average(avg_arr[y, x], weights = wt_arr[y, x])
                else: z = np.average(avg_arr[y, x])
                line = [geo_x, geo_y, float(z)]
                yield(line)
    
def xyz_line(line, dst_port = sys.stdout):
    '''write "xyz" `line` to `dst_port`
    `line` should be a list of xyz values [x, y, z, ...].'''
    
    l = '{}\n'.format(_xyz_config['delim'].join([str(x) for x in line]))
    if dst_port != sys.stdout: l = l.encode('utf-8')
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
            if l[0] < minmax[0]: minmax[0] = l[0]
            elif l[0] > minmax[1]: minmax[1] = l[0]
            if l[1] < minmax[2]: minmax[2] = l[1]
            elif l[1] > minmax[3]: minmax[3] = l[1]
            if l[2] < minmax[4]: minmax[4] = l[2]
            elif l[2] > minmax[5]: minmax[5] = l[2]
    with open('{}.inf'.format(src_xyz.name), 'w') as inf:
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)

def xyz_inf_entry(entry):
    with open(entry[0]) as infile:
        try:
            minmax = mb_inf(infile)
        except: minmax = xyz_inf(infile)
    return(minmax)        

def xyz_yield_entry(entry, region = None, verbose = False):
    with open(entry[0]) as infile:
        for line in xyz_parse(infile):
            if region is not None:
                if xyz_in_region_p(line, region):
                    yield(line + [entry[2]] if entry[2] is not None else line)
            else: yield(line + [entry[2]] if entry[2] is not None else line)

def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False):
    for xyz in xyz_yield_entry(entry, region, verbose):
        xyz_line(xyz, dst_port)

## ==============================================
## datalists and entries - datalists.py
##
## datalist processing (datalists fmt:-1)
## entry processing fmt:*
## ==============================================
_known_dl_delims = [' ', ':']
_known_datalist_fmts = {-1: ['datalist', 'mb-1'], 168: ['xyz', 'csv', 'dat'], 200: ['tif', 'img', 'grd', 'nc']}
_dl_inf_h = {
    -1: lambda e: datalist_inf_entry(e),
    168: lambda e: xyz_inf_entry(e),
    200: lambda e: gdal_inf_entry(e)
}
_dl_pass_h = lambda e: os.path.exists(e[0])

def datalist_inf(dl, inf_file = True):
    '''return the region of the datalist and generate
    an associate `.inf` file if `inf_file` is True.'''
    
    out_regions = []
    minmax = None
    for entry in datalist(dl):
        out_regions.append(inf_entry(entry)[:4])
    
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
        echo_msg('generating inf for datalist {}'.format(dl))
        with open('{}.inf'.format(dl), 'w') as inf:
            inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)

def datalist_inf_entry(e):
    '''write an inf file for datalist entry e
    
    returna the region [xmin, xmax, ymin, ymax]'''
    
    return(datalist_inf(e[0]))

def datalist_append_entry(entry, datalist):
    '''append entry to datalist file `datalist`'''
    
    with open(datalist, 'a') as outfile:
        outfile.write('{}\n'.format(' '.join([str(x) for x in entry])))

def archive_entry(entry, dirname = 'archive', region = None, inc = 1, weight = None, verbose = None):
    '''archive a datalist entry.
    a datalist entry is [path, format, weight, ...]'''
    
    if region is None:
        a_name = entry[-1]
    else: a_name = '{}_{}_{}'.format(entry[-1], region_format(region, 'fn'), this_year())
    
    i_dir = os.path.dirname(entry[0])
    i_xyz = os.path.basename(entry[0]).split('.')[0]
    a_dir = os.path.join(dirname, a_name, 'data', entry[-1])
    a_xyz_dir = os.path.join(a_dir, 'xyz')
    a_xyz = os.path.join(a_xyz_dir, i_xyz + '.xyz')
    a_dl = os.path.join(a_xyz_dir, '{}.datalist'.format(entry[-1]))
    
    if not os.path.exists(a_dir): os.makedirs(a_dir)
    if not os.path.exists(a_xyz_dir): os.makedirs(a_xyz_dir)

    dump_xyz = lambda p: xyz_dump_entry(entry, dst_port = p, region = region, verbose = verbose)
    dump_gdal = lambda p: gdal_dump_entry(entry, dst_port = p, region = region, verbose = verbose)
    
    with open(a_xyz, 'w') as fob:
        if entry[1] == 168: dump_xyz(fob)
        elif entry[1] == 200: dump_gdal(fob)

    mb_inf(a_xyz)
    datalist_append_entry([i_xyz + '.xyz', 168, entry[2] if entry[2] is not None else 1], a_dl)

def datalist_archive(wg, arch_dir = 'archive', region = None, verbose = False):
    '''archive the data from wg_config datalist to `arch_dir`'''
    
    if region is not None:
        #dl_p = lambda e: regions_intersect_ogr_p(region, inf_entry(e)) if e[1] != -1 else True
        dl_p = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
    else: dl_p = _dl_pass_h

    for this_entry in datalist(wg['datalist'], wt = 1 if wg['weights'] else None, pass_h = dl_p, verbose = verbose):
        if this_entry[1] == 168 or this_entry[1] == 200:
            archive_entry(this_entry, dirname = arch_dir, region = region, inc = wg['inc'], weight = wg['weights'], verbose = verbose)
    a_dl = os.path.join(arch_dir, '{}.datalist'.format(wg['name']))
            
    for dir_, _, files in os.walk(arch_dir):
        for f in files:
            if '.datalist' in f:
                rel_dir = os.path.relpath(dir_, arch_dir)
                rel_file = os.path.join(rel_dir, f)
                datalist_append_entry([rel_file, -1, 1], a_dl)
    return(a_dl, 0)
        
def datalist_dump(wg, dst_port = sys.stdout, region = None, verbose = False):
    '''dump the xyz data from datalist to dst_port'''
    
    if region is not None:
        #dl_p = lambda e: regions_intersect_ogr_p(region, inf_entry(e)) if e[1] != -1 else True
        dl_p = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
    else: dl_p = _dl_pass_h
    
    for this_entry in datalist(wg['datalist'], wt = 1 if wg['weights'] else None, pass_h = dl_p, verbose = verbose):
        if this_entry[1] == 168:
            xyz_dump_entry(this_entry, region = region, dst_port = dst_port)
        elif this_entry[1] == 200:
            gdal_dump_entry(this_entry, region = region, dst_port = dst_port)

def datalist_echo(entry):
    '''echo datalist entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))
    datalist(entry[0])

def datafile_echo(entry):
    '''echo datafile entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))

def datalist_master(dls, master = '.master.datalist'):
    '''set the master datalist
    `dls` is a list of datalist entries, minimally: ['datafile.xyz']

    returns the master datalist filename'''
    
    with open(master, 'w') as md:        
        for dl in dls:
            if os.path.exists(dl):
                if len(dl.split(':')) == 1:
                    for key in _known_datalist_fmts.keys():
                        if dl.split(':')[0].split('.')[-1] in _known_datalist_fmts[key]:
                            md.write('{} {} 1\n'.format(dl, key))
    if os.stat(master).st_size == 0:
        remove_glob(master)
        echo_error_msg('bad datalist/entry, {}'.format(dls))
        return(None)    
    return(master)

def entry2py(dle):
    '''convert a datalist entry to python

    return the entry as a list [fn, fmt, wt, ...]'''
    
    for delim in _known_dl_delims:
        this_entry = dle.rstrip().split(delim)
        if len(this_entry) > 1: break
    try:
        entry = [x if n == 0 else float(x) if n < 3 else x for n, x in enumerate(this_entry)]
    except ValueError as e: return(None)
    if len(entry) < 2:
        for key in _known_datalist_fmts.keys():
            if entry[0].split('.')[-1] in _known_datalist_fmts[key]:
                entry.append(key)
    if len(entry) < 3: entry.append(1)
    return(entry)

def datalist2py(dl):
    '''convert a datalist to python data
    
    returns a list of datalist entries.'''
    
    these_entries = []
    this_entry = entry2py(dl)
    if this_entry[1] == -1:
        with open(this_entry[0], 'r') as op:
            for this_line in op:
                if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                    these_entries.append(entry2py(this_line.rstrip()))
    else: these_entries.append(this_entry)
    return(these_entries)

def datalist_yield_xyz(dl, fmt = -1, wt = None, pass_h = lambda e: os.path.exits(e), verbose = False):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''
    
    for this_entry in datalist(dl, fmt, wt, pass_h, verbose):
        if this_entry[1] == 168:
            for xyz in xyz_yield_entry(this_entry):
                yield(xyz)
        elif this_entry[1] == 200:
            for xyz in gdal_yield_entry(this_entry):
                yield(xyz)
                
def datalist(dl, fmt = -1, wt = None, pass_h = lambda e: os.path.exists(e[0]), verbose = False):
    '''recurse a datalist/entry
    for entry in datalist(dl): do_something_with entry

    yields entry [path, fmt, wt, ...]'''
    
    this_dir = os.path.dirname(dl)
    these_entries = datalist2py(dl)
    if len(these_entries) == 0: these_entries = [entry2py(dl)]
    for this_entry in these_entries:
        this_entry[0] = os.path.join(this_dir, this_entry[0])
        if wt is not None:
            this_entry[2] = wt * this_entry[2]
        else: this_entry[2] = wt
        this_entry_md = ' '.join(this_entry[3:]).split(',')
        this_entry = this_entry[:3] + [this_entry_md] + [os.path.basename(dl).split('.')[0]]
        if pass_h(this_entry):
            if verbose: echo_msg('{} {}'.format('scanning datalist ({})'.format(this_entry[2]) if this_entry[1] == -1 else 'using datafile', this_entry[0]))
            if this_entry[1] == -1:
                for entry in datalist(this_entry[0], fmt, this_entry[2], pass_h, verbose):
                    yield(entry)
            else: yield(this_entry)
            
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
    'name_prefix': None,
    'node': 'pixel',
    'fmt': 'GTiff',
    'extend': 0,
    'weights': None,
    'fltr': None,
    'sample': None,
    'clip': None,
    'epsg': 4326,
    'mod': 'surface',
    'mod_args': (),
    'verbose': False,
    'gc': config_check()
}

## ==============================================
## the default waffles config dictionary.
## lambda returns dictionary with default waffles
## config data:
## {
##     'datalist': None,
##     'datalists': [],
##     'region': None,
##     'inc': None,
##     'name': 'waffles_dem',
##     'name_prefix': None,
##     'node': 'pixel',
##     'fmt': 'GTiff',
##     'extend': 0,
##     'weights': None,
##     'fltr': None,
##     'sample': None,
##     'clip': None,
##     'epsg': 4326,
##     'mod': 'surface',
##     'mod_args': (),
##     'verbose': False,
##     'gc': config_check()
## }
## ==============================================
waffles_config = lambda: copy.deepcopy(_waffles_grid_info)
    
def waffles_dict2wg(wg = _waffles_grid_info):
    '''copy the `wg` dict and add any missing keys.
    also validate the key values and return the valid waffles_config
    
    returns a complete and validated waffles_config dict.'''
    
    wg = copy.deepcopy(wg)
    keys = wg.keys()

    ## ==============================================
    ## check the waffles config dict and set
    ## missing values and their defaults.
    ## ==============================================
    if 'datalist' not in keys: wg['datalist'] = None
    if 'datalists' not in keys:
        if wg['datalist'] is not None:
            wg['datalists'] = [x[0] for x in datalist2py(wg['datalist'])]
        else: wg['datalists'] = None
    if 'region' not in keys: wg['region'] = None
    if 'inc' not in keys: wg['inc'] = None
    else: wg['inc'] = gmt_inc2inc(str(wg['inc']))
    if 'name' not in keys: wg['name'] = 'waffles_dem'
    if 'name_prefix' not in keys: wg['name_prefix'] = None
    if 'node' not in keys: wg['node'] = 'pixel'
    if 'fmt' not in keys: wg['fmt'] = 'GTiff'
    if 'extend' not in keys: wg['extend'] = 0
    else: wg['extend'] = int_or(wg['extend'], 0)
    if 'weights' not in keys: wg['weights'] = None
    if 'fltr' not in keys: wg['fltr'] = None
    if 'sample' not in keys: wg['sample'] = None
    else: wg['sample'] = gmt_inc2inc(str(wg['sample']))
    if 'clip' not in keys: wg['clip'] = None
    if 'epsg' not in keys: wg['epsg'] = 4326
    else: wg['epsg'] = int_or(wg['epsg'], 4326)
    if 'mod' not in keys: wg['mod'] = 'surface'
    if 'mod_args' not in keys: wg['mod_args'] = ()
    if 'verbose' not in keys: wg['verbose'] = False
    else: wg['verbose'] = False if not wg['verbose'] or str(wg['verbose']).lower() == 'false' or wg['verbose'] is None else True
    wg['gc'] = config_check()
    
    ## ==============================================
    ## set the master datalist to the mentioned
    ## datalists/datasets
    ## note: the vdatum module doesn't need a datalist
    ## ==============================================
    if wg['datalist'] is None and len(wg['datalists']) > 0:
        wg['datalist'] = datalist_master(wg['datalists'], '{}_mstr.datalist'.format(wg['name']))
    if wg['mod'].lower() != 'vdatum':
        if wg['datalist'] is None:
            echo_error_msg('invalid datalist/s entry')
            return(None)
        
        ## ==============================================
        ## set the region to that of the datalist if
        ## the region was not specified.
        ## ==============================================
        if wg['region'] is None or not region_valid_p(wg['region']):
            wg['region'] = datalist_inf(wg['datalist'])

    if wg['region'] is None:
        echo_error_msg('invalid region and/or datalist/s entry')
        return(None)
    if wg['inc'] is None: wg['inc'] = (wg['region'][1] - wg['region'][0]) / 500
    
    ## ==============================================
    ## if `name_prefix` is set; append region/inc/year
    ## to `prefix_name` and set `name` to that.
    ## ==============================================
    if wg['name_prefix'] is not None:
        wg['name'] = waffles_append_fn(wg['name_prefix'], wg['region'], wg['sample'] if wg['sample'] is not None else wg['inc'])        
    return(wg)

_waffles_modules = {
    'surface': [lambda args: waffles_gmt_surface(**args), '''SPLINE DEM via GMT surface
    \t\t\t  < surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
    \t\t\t  :tension=[0-1] - Spline tension.'''],
    'triangulate': [lambda args: waffles_gmt_triangulate(**args), '''TRIANGULATION DEM via GMT triangulate'''],
    'nearest': [lambda args: waffles_nearneighbor(**args), '''NEAREST NEIGHBOR DEM via GMT or GDAL
    \t\t\t  < nearest:radius=6s:use_gdal=False >
    \t\t\t  :radius=[value] - Nearest Neighbor search radius
    \t\t\t  :use_gdal=[True/False] - use gdal grid nearest algorithm'''],
    'num': [lambda args: waffles_num(**args), '''Uninterpolated DEM populated by <mode>.
    \t\t\t  < num:mode=n >
    \t\t\t  :mode=[key] - specify mode of grid population: k (mask), m (mean) or n (num)'''],
    'spat-metadata': [lambda args: waffles_spatial_metadata(**args), '''Datalist SPATIAL METADATA <beta>
    \t\t\t  < spatial-metadata:inc=None >
    \t\t\t  :inc=[increment] - Spatial metadata resolution [default uses -E value]'''],
    'vdatum': [lambda args: waffles_vdatum(**args), '''VDATUM transformation grid
    \t\t\t  < vdatum:ivert=navd88:overt=mhw:region=3:jar=None >
    \t\t\t  :ivert=[vdatum] - Input VDatum vertical datum.
    \t\t\t  :overt=[vdatum] - Output VDatum vertical datum.
    \t\t\t  :region=[0-10] - VDatum region (3 is CONUS).
    \t\t\t  :jar=[/path/to/vdatum.jar] - VDatum jar path - (auto-locates by default)'''],
    'mbgrid': [lambda args: waffles_mbgrid(**args), '''Weighted SPLINE DEM via mbgrid
    \t\t\t  < mbgrid:tension=35:dist=10/3:use_datalists=False >
    \t\t\t  :tension=[0-100] - Spline tension.
    \t\t\t  :dist=[value] - MBgrid -C switch (distance to fill nodata with spline)
    \t\t\t  :use_datalists=[True/False] - use waffles built-in datalists'''],
    'invdst': [lambda args: waffles_invdst(**args), '''INVERSE DISTANCE DEM via GDAL GRID
    \t\t\t  < invdst:power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1 >'''],
    'average': [lambda args: waffles_moving_average(**args), '''Moving AVERAGE DEM via GDAL GRID
    \t\t\t  < average:radius1=0.01:radius2=0.01 >'''],
    'archive': [lambda args: waffles_archive(**args), '''archive the datalist
    \t\t\t  < archive:arch_dir=archive >
    \t\t\t  :arch_dir=[dirname] - archive data to dirname'''],
}

## ==============================================
## module descriptors (used in cli help)
## ==============================================
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['{:22}{}\n'.format(key, x[key][-1]) for key in x]) + '\n'
_waffles_module_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

## ==============================================
## the "proc-region" region_buffer(wg['region'], (wg['inc'] * 20) + (wg['inc'] * wg['extend']))
## ==============================================
waffles_proc_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * 20) + (wg['inc'] * wg['extend']))
waffles_proc_str = lambda wg: region_format(waffles_proc_region(wg), 'gmt')

## ==============================================
## the "dist-region" region_buffer(wg['region'], (wg['inc'] * wg['extend']))
## ==============================================
waffles_dist_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * wg['extend']))

## ==============================================
## the datalist dump function, to use in run_cmd()
## ==============================================
waffles_dl_func = lambda wg: lambda p: datalist_dump(wg, dst_port = p, region = waffles_proc_region(wg), verbose = wg['verbose'])

## ==============================================
## grid registration string for use in GTM programs
## ==============================================
waffles_gmt_reg_str = lambda wg: '-r' if wg['node'] == 'pixel' else ''

## ==============================================
## the 'long-name' used from prefix
## ==============================================
waffles_append_fn = lambda bn, region, inc: '{}{}_{}_{}'.format(bn, inc2str_inc(inc), region_format(region, 'fn'), this_year())
    
def waffles_mbgrid(wg = _waffles_grid_info, dist = '10/3', tension = 35, use_datalists = False):
    '''Generate a DEM with MBSystem's mbgrid program.
    if `use_datalists` is True, will parse the datalist through
    waffles instead of mbsystem.'''
    
    import shutil
    if wg['gc']['MBGRID'] is None:
        echo_error_msg('MBSystem must be installed to use the MBGRID module')
        return(None, -1)
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the MBGRID module')
        return(None, -1)
    if use_datalists:
        datalist_archive(wg, arch_dir = '.mb_tmp_datalist', verbose = True)
        wg['datalist'] = datalist_master(['.mb_tmp_datalist/{}.datalist'.format(wg['name'])])

    out, status = run_mbgrid(wg['datalist'], waffles_dist_region(wg), wg['inc'], wg['name'], dist = dist, tension = tension, verbose = True)
    remove_glob('*.cmd')
    gmt_grd2gdal('{}.grd'.format(wg['name']))
    if wg['node'] == 'pixel': gmt_sample_gnr('{}.tif'.format(wg['name']), verbose = wg['verbose'])
    remove_glob('{}.grd'.format(wg['name']))
    if use_datalists: shutil.rmtree('.mb_tmp_datalist')
    return(0, 0)

def waffles_gmt_surface(wg = _waffles_grid_info, tension = .35, relaxation = 1.2, lower_limit = 'd', upper_limit = 'd'):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the SURFACE module')
        return(None, -1)
    dem_surf_cmd = ('gmt blockmean {} -I{:.10f}{} -V {} | gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} {}\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_gmt_reg_str(wg), waffles_proc_str(wg), \
             wg['inc'], wg['name'], tension, relaxation, lower_limit, upper_limit, waffles_gmt_reg_str(wg)))
    return(run_cmd(dem_surf_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

def waffles_gmt_triangulate(wg = _waffles_grid_info):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the TRIANGULATE module')
        return(None, -1)
    dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V {} | gmt triangulate {} -I{:.10f} -V -G{}.tif=gd+n-9999:GTiff {}\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_gmt_reg_str(wg), waffles_proc_str(wg), \
             wg['inc'], wg['name'], waffles_gmt_reg_str(wg)))
    return(run_cmd(dem_tri_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

def waffles_nearneighbor(wg = _waffles_grid_info, radius = None, use_gdal = False):
    '''genearte a DEM with GMT nearneighbor'''
    
    radius = wg['inc'] * 2 if radius is None else gmt_inc2inc(radius)
    if wg['gc']['GMT'] is not None and not use_gdal:
        dem_nn_cmd = ('gmt blockmean {} -I{:.10f}{} -V {} | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}.tif=gd+n-9999:GTiff {}\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_gmt_reg_str(wg), waffles_proc_str(wg), \
                 wg['inc'], radius, wg['name'], waffles_gmt_reg_str(wg)))
        return(run_cmd(dem_nn_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))
    else: return(waffles_gdal_grid(wg, 'nearest:radius1={}:radius2={}:nodata=-9999'.format(radius, radius)))
    
def waffles_num(wg = _waffles_grid_info, mode = 'n'):
    '''Generate an uninterpolated num grid.
    mode of `k` generates a mask grid
    mode of `m` generates a mean grid
    mode of `n` generates a num grid'''
    
    wg['region'] = region_buffer(wg['region'], wg['inc'] * .5) if wg['node'] == 'grid' else wg['region']
    region = waffles_dist_region(wg)
    #dlh = lambda e: regions_intersect_ogr_p(region, inf_entry(e)) if e[1] != -1 else True
    dlh = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
    if wg['weights']:
        dly = xyz_block(datalist_yield_xyz(wg['datalist'], pass_h = dlh, wt = 1 if wg['weights'] is not None else None, verbose = wg['verbose']), region, wg['inc'], weights = True if wg['weights'] else False)
        #dly = xyz_block(datalist_yield_xyz(wg['datalist'], pass_h = dlh, wt = 1, verbose = wg['verbose']), region, wg['inc'], weights = True)
    else: dly = datalist_yield_xyz(wg['datalist'], pass_h = dlh, verbose = wg['verbose'])
    return(gdal_xyz2gdal(dly, '{}.tif'.format(wg['name']), region, wg['inc'], dst_format = wg['fmt'], mode = mode))
                                
def waffles_spatial_metadata(wg):
    '''generate spatial-metadata for the top-level of the datalist
    top-level entries in the datalist should include a comma-separated list
    of field entries, (see `v_fields`). should be column 3 (from 0) in the 
    datalist entry (e.g."datalist.datalist -1 1 name,agency,date,type,resolution,hdatum,vdatum,url") '''
    
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
        defn = layer.GetLayerDefn()
        entry = datalist2py(wg['datalist'])[0]
        print(entry)
        these_entries = datalist2py(entry[0])
        print(these_entries)
        print(wg['datalist'])
        for this_entry in these_entries:
            this_entry[0] = os.path.join(os.path.dirname(entry[0]), this_entry[0])
            print(this_entry)
            wg['datalist'] = this_entry[0]
            wg['name'] = os.path.basename(this_entry[0]).split('.')[0]
            o_v_fields = [wg['name'], 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
            echo_msg('scanning datalist {}...'.format(wg['datalist']))
            ng, s = waffles_num(wg, mode = 'k')
            
            if gdal_infos(ng, True)['zmax'] == 1:
                tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(wg['name']))
                tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(wg['name']), None, ogr.wkbMultiPolygon)
                tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                gdal_polygonize(ng, tmp_layer)
                
                if len(tmp_layer) > 1:
                    out_feat = gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                    [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(v_fields)]
                    layer.CreateFeature(out_feat)
                    
                tmp_ds = tmp_layer = out_feat = None
                remove_glob('{}_poly.*'.format(wg['name']))
            remove_glob('{}'.format(ng))
        ds = None
    return(dst_vector, 0)

def waffles_archive(wg = _waffles_grid_info, arch_dir = None):
    '''archive the datalist to `arch_dir`
    converts all data to xyz'''
    
    arch_dir = wg['name'] if arch_dir is None else arch_dir
    return(datalist_archive(wg, arch_dir = arch_dir, region = waffles_proc_region(wg), verbose = True))

def waffles_vdatum(wg = _waffles_grid_info, ivert = 'navd88', overt = 'mhw', region = '3', jar = None):
    '''generate a 'conversion-grid' with vdatum.
    output will be the differences (surfaced) between 
    `ivert` and `overt` for the region'''
    
    vc = _vd_config
    vc['jar'] = jar
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
        wg['datalist'] = datalist_master(['result/empty.xyz'])
        out, status = waffles_gmt_surface(wg, tension = 0, upper_limit = lu, lower_limit = ll)
        
    remove_glob('empty.*')
    remove_glob('result/*')
    remove_glob('.master.datalist')
    os.removedirs('result')
    return(out, status)

def waffles_gdal_grid(wg = _waffles_grid_info, alg_str = 'linear:radius=1'):
    '''run gdal grid using alg_str
    parse the data through xyz_block to get weighted mean before
    building the GDAL dataset to pass into gdal_grid'''
    
    wg['region'] = region_buffer(wg['region'], wg['inc'] * .5) if wg['node'] == 'grid' else wg['region']
    region = waffles_proc_region(wg)
    dlh = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
    wt = 1 if wg['weights'] is not None else None
    dly = xyz_block(datalist_yield_xyz(wg['datalist'], pass_h = dlh, wt = wt, verbose = wg['verbose']), region, wg['inc'], weights = False if wg['weights'] is None else True)
    ds = xyz2gdal_ds(dly, '{}'.format(wg['name']))
    xcount, ycount, dst_gt = gdal_region2gt(wg['region'], wg['inc'])
    gd_opts = gdal.GridOptions(outputType = gdal.GDT_Float32, noData = -9999, format = 'GTiff', \
                               width = xcount, height = ycount, algorithm = alg_str, callback = _gdal_progress if wg['verbose'] else None, \
                               outputBounds = [region[0], region[3], region[1], region[2]])
    gdal.Grid('{}.tif'.format(wg['name']), ds, options = gd_opts)
    ds = None
    gdal_set_nodata('{}.tif'.format(wg['name']), -9999)
    return(0, 0)

def waffles_invdst(wg = _waffles_grid_info, power = 2.0, smoothing = 0.0, radius1 = None, radius2 = None, angle = 0.0, \
                   max_points = 0, min_points = 0, nodata = -9999):
    '''Generate an inverse distance grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmt_inc2inc(radius2)
    gg_mod = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
                             .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

def waffles_moving_average(wg = _waffles_grid_info, radius1 = None, radius2 = None, angle = 0.0, min_points = 0, nodata = -9999):
    '''generate a moving average grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmt_inc2inc(radius2)
    gg_mod = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'.format(radius1, radius2, angle, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

def waffles_wg_valid_p(wg = _waffles_grid_info):
    '''return True if wg_config appears valid'''
    
    if wg['datalists'] is None: return(False)
    if wg['mod'] is None: return(False)
    else: return(True)
        
def waffles_gdal_md(wg):
    '''add metadata to the waffles dem'''
    
    ds = gdal.Open('{}.tif'.format(wg['name']), gdal.GA_Update)
    if ds is not None:
        md = ds.GetMetadata()
        if wg['node'] == 'pixel':
            md['AREA_OR_POINT'] = 'Area'
        else: md['AREA_OR_POINT'] = 'Point'
        ds.SetMetadata(md)
        ds = None
    else: echo_error_msg('failed to set metadata')

def waffles_run(wg = _waffles_grid_info):
    '''generate a DEM using wg dict settings
    see waffles_dict2wg() to generate a wg config.

    - runs the waffles module to generate the DEM
    - optionally clips the output to shapefile
    - optionally filters the output
    - optionally resamples the output
    - cuts the output to dist-size
    - reformats the output to final format
    - sets metadata in output

    returns dem-fn'''

    ## ==============================================
    ## validate and/or set the waffles_config
    ## ==============================================
    wg = waffles_dict2wg(wg)
    if wg is None:
        echo_error_msg('invalid configuration, {}'.format(json.dumps(wg, indent=4, sort_keys=True)))
        sys.exit(-1)
    
    args_d = {}
    for arg in wg['mod_args']:
        p_arg = arg.split('=')
        args_d[p_arg[0]] = p_arg[1]
    args_d['wg'] = wg
        
    if wg['verbose']:
        echo_msg(json.dumps(wg, indent = 4, sort_keys = True))
        echo_msg('running module {}...'.format(wg['mod']))

    ## ==============================================
    ## gererate the DEM
    ## ==============================================
    dem = '{}.tif'.format(wg['name'])

    #try:
    out, status = _waffles_modules[wg['mod']][0](args_d)
    #except Exception as e:
    #    echo_error_msg('{}'.format(e))
    #    status = -1
    #except KeyboardInterrupt as e:
    #    echo_error_msg('killed by user, {}'.format(e))
    #    sys.exit(-1)

    if status != 0: remove_glob(dem)
    if not os.path.exists(dem): return(None)
    gdi = gdal_infos(dem, scan=True)
    if gdi is not None:
        if np.isnan(gdi['zmin']):
            remove_glob(dem)
            return(None)
    else: return(None)
    
    ## ==============================================
    ## spat-metadata and archive don't produce grids,
    ## exit here if running one of those mods.
    ## ==============================================
    if wg['mod'] == 'spat-metadata' or wg['mod'] == 'archive': return(0)
    gdal_set_epsg(dem, wg['epsg'])
    waffles_gdal_md(wg)
                
    ## ==============================================
    ## optionally clip the DEM to polygon
    ## ==============================================
    if wg['clip'] is not None:
        if wg['verbose']: echo_msg('clipping {}...'.format(dem))
        clip_args = {}
        cp = wg['clip'].split(':')
        clip_args['src_ply'] = cp[0]
        clip_args = args2dict(cp[1:], clip_args)
        gdal_clip(dem, **clip_args)
        
    ## ==============================================
    ## optionally filter the DEM 
    ## ==============================================
    if wg['fltr'] is not None:
        if wg['verbose']: echo_msg('filtering {}...'.format(dem))
        fltr_args = {}
        fltr = wg['fltr'].split(':')
        fltr_args['fltr'] = gmt_inc2inc(fltr[0])
        fltr_args['use_gmt'] = True
        fltr_args = args2dict(fltr[1:], fltr_args)        
        if fltr_args['use_gmt']: fltr_args['use_gmt'] = True if wg['gc']['GMT'] is not None else False
        try:
            gdal_smooth(dem, 'tmp_s.tif', **fltr_args)
            os.rename('tmp_s.tif', dem)
        except TypeError as e: echo_error_msg('{}'.format(e))

    ## ==============================================
    ## optionally resample the DEM 
    ## ==============================================
    if wg['sample'] is not None:
        if wg['verbose']: echo_msg('resampling {}...'.format(dem))
        if wg['gc']['GMT'] is not None:
            gmt_sample_inc(dem, inc = wg['sample'], verbose = wg['verbose'])
        else:
            out, status = run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} tmp.tif'.format(inc, inc, src_grd, region_format(waffles_proc_region(wg)), verbose = verbose))
            if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    ## ==============================================
    ## cut dem to final size - region buffered by (inc * extend)
    ## ==============================================
    try:
        out = gdal_cut(dem, waffles_dist_region(wg), 'tmp_cut.tif')
        if out is not None: os.rename('tmp_cut.tif', dem)
    except OSError as e:
        remove_glob('tmp_cut.tif')
        echo_error_msg('cut failed, is the dem open somewhere, {}'.format(e))
            
    ## ==============================================
    ## convert to final format
    ## ==============================================
    if wg['fmt'] != 'GTiff':
        orig_dem = dem
        if wg['gc']['GMT'] is not None:
            dem = gmt_grd2gdal(dem, wg['fmt'])
        else: dem = gdal_gdal2gdal(dem, wg['fmt'])
        remove_glob(orig_dem)

    ## ==============================================
    ## set the projection and other metadata
    ## ==============================================
    gdal_set_epsg(dem, wg['epsg'])
    waffles_gdal_md(wg)
    ## ==============================================
    ## if dem has data, return
    ## ==============================================
    return(dem)
        
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
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -W, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -a, --archive\t\tArchive the datalist to the given region.
  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tuse grid-node registration, default is pixel-node

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and master datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Modules (see waffles --modules for more info):
  {}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(_waffles_module_short_desc(_waffles_modules))

def waffles_cli(argv = sys.argv):
    '''run waffles from command-line
    e.g. `python waffles.py` 
    generates a waffles_config from the command-line options
    and either outputs the or runs the waffles_config
    on each region supplied (multiple regions can be supplied
    by using a vector file as the -R option.)
    See `waffles_cli_usage` for full cli options.'''
    
    wg = waffles_config()
    wg_user = None
    dls = []
    region = None
    module = None
    want_prefix = False
    want_verbose = False
    want_config = False
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
            try:
                wg['extend'] = int(argv[i + 1])
            except ValueError as e:
                echo_error_msg('invalid -X option, {}'.format(e))
            i += 1
        elif arg[:2] == '-X':
            try:
                wg['extend'] = int(arg[2:])
            except ValueError as e:
                echo_error_msg('invalid -X option, {}'.format(e))
        elif arg == '--wg-config' or arg == '-W':
            wg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-W': wg_user = arg[2:]
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = arg[2:]
        elif arg == '-w': wg['weights'] = True
        elif arg == '-p': want_prefix = True
        elif arg == '-r': wg['node'] = 'grid'
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
                    dem = waffles_run(wg)
                    sys.exit(0)
            except ValueError as e:
                wg = waffles_config()
                echo_error_msg('could not parse json from {}, {}'.format(wg_user, e))
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

    if wg['mod'] != 'vdatum':
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            echo_error_msg('''must specify a datalist/entry file''')
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
    else: these_regions = [None]
    if len(these_regions) == 0: echo_error_msg('failed to parse region(s)')
    if want_prefix or len(these_regions) > 1: wg['name_prefix'] = wg['name']
    
    ## ==============================================
    ## run waffles for each input region.
    ## ==============================================
    bn = wg['name']
    for this_region in these_regions:
        wg['region'] = this_region
        
        if want_config:
            this_wg = waffles_dict2wg(wg)
            if this_wg is not None:
                echo_msg(json.dumps(this_wg, indent = 4, sort_keys = True))
                with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                    echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                    echo_msg('generating master datalist: {}_mstr.datalist'.format(this_wg['name']))
                    wg_json.write(json.dumps(this_wg, indent = 4, sort_keys = True))
            else: echo_error_msg('could not parse config.')
        else:
            #try:
            dem = waffles_run(wg)
            #except RuntimeError or OSError as e:
            #echo_error_msg('Cannot access {}.tif, may be in use elsewhere, {}'.format(wg['name'], e))

## ==============================================
## datalists cli
## ==============================================
datalists_cli_usage = '''datalists [-ciwR] <datalist/entry>

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''

def datalists_cli(argv = sys.argv):
    dls = []
    region = None
    want_weight = False
    want_infos = False
    rec_infos = False
    want_verbose = False
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            region = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-R': region = str(arg[2:])
        elif arg == '-w': want_weight = True
        elif arg == '-i': want_infos = True
        elif arg == '-c': rec_infos = True
        elif arg == '--verbose' or arg == '-V': want_verbose = True
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(datalists_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(_version))
            sys.exit(0)
        else: dls.append(arg)
        i += 1
    if len(dls) == 0:
        sys.stderr.write(datalists_cli_usage)
        echo_error_msg('''must specify a datalist/entry file''')
        sys.exit(-1)
    dl_p = lambda e: os.path.exists(e[0])
    if region is not None:
        try:
            region = [float(x) for x in region.split('/')]
        except ValueError as e:
            sys.stderr.write(datalists_cli_usage)
            echo_error_msg('bad region, {}'.format(e))
            sys.exit(-1)
        dl_p = lambda e: regions_intersect_ogr_p(region, inf_entry(e)) if e[1] != -1 else True
    wt = 1 if want_weight else None
            
    ## ==============================================
    ## recurse the datalist
    ## ==============================================
    master = datalist_master(dls)
    if master is not None:
        #datalist_archive(master, arch_dir = 'archive', region = region, verbose = want_verbose)
        if want_infos: print(datalist_inf(master, inf_file = False))#sys.stdout.write(' '.format([str(x) for x in datalist_inf(master)]))
        elif rec_infos:
            inf_entry([master, -1, 1])
        else:
            ## ==============================================
            ## dump to stdout
            ## ==============================================
            datalist_dump({'datalist':master, 'weights': want_weight})
        remove_glob(master)

## ==============================================
## mainline -- run waffles directly...
##
## run waffles:
## % python waffles.py dem <args>
## % python waffles.p <args>
##
## run datalists:
## % python waffles.py datalists <args>
## ==============================================
if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'dem':
            waffles_cli(sys.argv[:1] + sys.argv[2:])
        elif sys.argv[1] == 'datalists':
            datalists_cli(sys.argv[:1] + sys.argv[2:])
        else: waffles_cli(sys.argv)
    else: waffles_cli(sys.argv)

### End
