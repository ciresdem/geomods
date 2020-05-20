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
## ~path ~format ~weight
##
## with the additional columns of ~"meta,data" ~etc
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
## or via _datalist_inf_hooks
##
## if 'region' is specified, will only process data that falls within
## the given region
##
## Adjust the _datalist*hooks dictionary to add custom processing/testing
##
### TODO:
## Add vdatum, uncertainty functions and modules
## Add mbgrid module
## Add gdal modules
## Add remove/replace module
## Add source uncertainty to uncertainty module
## Get -w (weight) option working
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

## ==============================================
## import gdal, etc.
## ==============================================
import numpy as np
import gdal
import ogr
import osr

## ==============================================
## General utility functions - utils.py
## ==============================================
_version = '0.5.1'

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
    globs = glob.glob(glob_str)
    for g in globs:
        try:
            os.remove(g)
        except: pass
    return(0)

## ==============================================
## system cmd verification and configs.
## ==============================================
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
    ae = '.exe' if host_os == 'win32' else ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal_grid --version')
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version')
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version | grep Version')
    _waff_co['BOUNDS'] = cmd_check('bounds{}'.format(ae), 'bounds --version')
    return(_waff_co)
    
## ==============================================
## stderr messaging
## ==============================================
def echo_error_msg2(msg, prefix = 'waffles'):
    '''echo error msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: error, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles'):
    '''echo msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: {}\n'.format(prefix, msg))

echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))

## ==============================================
## some known values to make guesses from
## ==============================================
_known_dl_delims = [':', ' ']
_known_delims = [',', ' ', '\t', '/', ':']
_known_datalist_fmts = {-1: ['datalist', 'mb-1'], 168: ['xyz', 'dat'], 200: ['tif', 'img', 'grd', 'nc']}

## ==============================================
## regions - regions are a bounding box list:
## [w, e, s, n]
## -- regions.py --
## ==============================================
def region_valid_p(region):
    '''return True if region appears valid'''
    if region is not None:
        if region[0] < region[1] and region[2] < region[3]: return(True)
        else: return(False)
    else: return(False)

def region_center(region):
    '''return the center point of the region'''
    xc = region[0] + (region[1] - region[0] / 2)
    yc = region[2] + (region[3] - region[2] / 2)
    return([xc, yc])

def region_pct(region, pctv):
    '''return the pctv val of the region'''
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv = 0, pct = False):
    '''return the region buffered by bv'''
    if pct: bv = region_pct(region, bv)
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def regions_reduce(region_a, region_b):
    '''return the minimum region when combining region_a and region_b'''
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] > region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] < region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] > region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] < region_b[3] else region_b[3]
    return(region_c)

def regions_merge(region_a, region_b):
    '''merge two regions into a single region'''
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] < region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] > region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] < region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] > region_b[3] else region_b[3]
    return(region_c)

def regions_intersect_p(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''
    if region_a is not None and region_b is not None:
        return(region_valid_p(regions_reduce(region_a, region_b)))
    else: return(None)
    
def regions_intersect_ogr_p(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''
    if region_a is not None and region_b is not None:
        geom_a = gdal_region2geom(region_a)
        geom_b = gdal_region2geom(region_b)
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(False)

def region_format(region, t = 'gmt'):
    '''format region to string'''
    if t == 'str': return('/'.join([str(x) for x in region]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'fn':
        if region[3] < 0: ns = 's'
        else: ns = 'n'
        if region[0] > 0: ew = 'e'
        else: ew = 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100) % 100), 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100) % 100)))
    
## ==============================================
## GMT Wrapper Functions - gmtfun.py
## wrapper functions to GMT system commands
##
## GMT must be installed on the system to run these
## functions and commands.
## ==============================================
def gmt_inf(src_xyz):
    '''generate an info (.inf) file from a src_xyz file using GMT.'''
    return(run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), verbose = False))

def gmt_grd_inf(src_grd):
    '''generate an info (.inf) file from a src_gdal file using GMT.'''
    return(run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_grd, src_grd), verbose = False))

def gmt_inc2inc(inc_str):
    '''convert an GMT-style inc_str (6s) to native units'''
    if inc_str is None: return(None)
    units = inc_str[-1]
    if units == 'c': inc = float(inc_str[:-1]) / 3600
    elif units == 's': inc = float(inc_str[:-1]) / 3600
    elif units == 'm': inc = float(inc_str[:-1]) / 360
    else: inc = float(inc_str)    
    return(inc)

def gmt_grd2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''Convert the grd file to tif using GMT'''
    dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
    grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
    '.format(src_grd, dst_gdal, dst_fmt))
    out, status = run_cmd(grd2gdal_cmd, verbose = verbose)
    if status == 0:
        return(dst_gdal)
    else: return(None)

def gmt_grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
    out, status = run_cmd(grdinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)

def gmt_gmtinfo(src_xyz, verbose = False):
    '''Return an info list of `src_xyz`'''
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    gmtinfo_cmd = ('gmt gmtinfo {} -C'.format(src_xyz))
    out, status = run_cmd(gmtinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)

def gdal_ogr_csv_vrt(src_fn):
    with open('{}.vrt'.format(src_fn), 'w') as o_vrt:
        o_vrt.write('''<OGRVRTDataSource>
  <OGRVRTLayer name="{}">
    <SrcDataSource>{}.csv</SrcDataSource>
    <GeometryType>wkbPoint</GeometryType>
    <GeometryField encoding="PointFromColumns" x="field_1" y="field_2" z="field_3"/>
  </OGRVRTLayer>
</OGRVRTDataSource>'''.format(src_fn, src_fn))
    
def gmt_block(datalist, mode = 'blockmean', inc = '1s', o_name = None, delim = 'SPACE', weights = False, verbose = False):
    '''run block/mean/median on src_xyz'''
    if mode == 'blockmean' or mode == 'blockmean':
        out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = {}'.format(delim.upper()), verbose = verbose)
        if mode == 'blockmean' and weights:
            mode = 'blockmean -Wi'
            datalist.want_weights = True
        if mode == 'blockmedian': mode = 'blockmedian -Q'
        if o_name is None: o_name = datalist._name
        if delim.lower() == 'comma':
            out_ext = 'csv'
            gdal_ogr_csv_vrt(o_name)
        else: out_ext = 'xyz'
        if os.path.exists(datalist._path):
            blk_cmd1 = ('gmt {} -V {} -I{} > {}.{}'.format(mode, datalist.region.gmt, inc, o_name, out_ext))
            out, status = run_cmd(blk_cmd1, verbose = True, data_fun = datalist._dump_data)
        else: status = -1
    else: status = -1
    remove_glob('gmt.conf')
    return(status)

def gmt_select_split(o_xyz, sub_region, sub_bn, verbose = False):
    '''split an xyz file into an inner and outer region.'''
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
    '''Cut `src_grd` to `src_region` '''
    cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
    return(run_cmd(cut_cmd1, verbose = True))

def gmt_grdfilter(src_grd, dst_grd, dist = '3s', verbose = False):
    '''filter `src_grd` '''
    ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1'.format(src_grd, dst_grd, src_grd, dist))
    return(run_cmd(ft_cmd1, verbose = verbose))

def gmt_grd2xyz(src_grd, dst_xyz, region = None, mask = None, verbose = False, want_datalist = False):
    '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region.
    Optionally, generate a datalist and inf file for the resultant xyz data.'''
    if mask:
        grdmask_cmd = ('gmt grdmath -N -V {} {} OR = tmp.grd'.format(src_grd, mask))
        out, status = run_cmd(grdmask_cmd, verbose = verbose)
        if status == 0: src_grd = 'tmp.grd'

    if region is not None and region_valid_p(region):
        region_str = region_format(region, 'gmt')
    else: region_str = ''

    grd2xyz_cmd = ('gmt grd2xyz -V {} -s {} > {}'.format(src_grd, region_str, dst_xyz))
    out, status = run_cmd(grd2xyz_cmd, verbose = verbose)

    if mask: remove_glob('tmp.grd')
    return(status)

def gmt_nan2zero(src_grd, node = 'pixel', verbose = False):
    '''convert nan values to zero'''
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = tmp.tif=gd+n-9999:GTiff'.format(src_grd))
    out, status = run_cmd(num_msk_cmd, verbose = True)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_grdcut(src_grd, region, verbose = False):
    '''cut a grid to region'''
    cut_cmd = ('gmt grdcut -V {} -Gtmp.grd {}\
    '.format(src_grd, region_format(region, 'gmt')))
    out, status = run_cmd(cut_cmd, verbose = True)
    if status == 0:
        remove_glob(src_grd)
        os.rename('tmp.grd', '{}'.format(src_grd))
    return(status)

def gmt_slope(src_dem, dst_slp, verbose = False):
    '''Generate a Slope grid from a DEM with GMT'''
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
    '''Generate a num-msk from a NUM grid.'''
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    return(run_cmd(num_msk_cmd, verbose = verbose))

def gmt_sample_gnr(src_grd, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node
    grid registration.'''
    out, status = run_cmd('gmt grdsample -T {} -Gtmp.tif=gd+n-9999:GTiff'.format(src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

## ==============================================
## MB-System Wrapper Functions - mbsfun.py
##
## MS-System must be installed on the system to run
## these functions and commands.
## ==============================================
def mb_inf(src_xyz, src_fmt = 168):
    '''generate an info (.inf) file from a src_xyz file using MBSystem.'''
    return(run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz), True))

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):
    '''Run the MBSystem command `mbgrid` given a datalist, region and increment.
    The datalist should be an MBSystem-style datalist; if using a waffles datalist, 
    it should be converted first (using datalist_archive()).'''
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1: dist = dist + '/2'
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.10f}/{:.10f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} > mb_proc.txt \
    '.format(datalist._path, region_format(region, 'gmt'), inc, inc, dst_name, dist, tension, e_switch))
    return(run_cmd(mbgrid_cmd, verbose = verbose))

## ==============================================
## GDAL Wrapper Functions - gdalfun.py
## ==============================================
def gdal_sr_wkt(epsg, esri = False):
    '''convert an epsg code to wkt'''
    try:
        int(epsg)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)

def gdal_fext(src_drv_name):
    '''return the common file extention given a GDAL driver name'''
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

def gdal_write (src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff'):
    '''write src_arr Array to gdal file dst_gdal using src_config'''
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

def gdal_cut(src_gdal, region, dst_fn):
    '''cut src_fn gdal file to srcwin and output dst_fn gdal file'''
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
    '''clip dem to either a src_ply polygon or the gsshg via GMT. use 'src_ply = None' to use gsshg.'''
    gi = gdal_infos(src_gdal)
    if gi is not None:
        gr_inv = '-i' if invert else ''
        gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'\
                 .format(gi['ndv'], gr_inv, os.path.basename(src_ply).split('.')[0], src_ply, src_gdal)
        return(run_cmd(gr_cmd, verbose = True))
    else: return(None)

def np_split(src_arr, sv = 0, nd = -9999):
    '''split numpy `src_arr` by `sv` (turn u/l into `nd`)'''
    u_arr = np.array(src_arr)
    l_arr = np.array(src_arr)
    u_arr[u_arr <= sv] = nd
    l_arr[l_arr >= sv] = nd
    return(u_arr, l_arr)
    
def gdal_split(src_gdal, split_value = 0):
    '''split raster file `src_gdal`into two files based on z value'''
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

def gdal_crop(src_fn):
    '''Crop `src_fn` GDAL file by it's NoData value. Returns cropped array.'''
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
    
def gdal_prj_file(dst_fn, epsg):
    '''generate a .prj file given an epsg code'''
    with open(dst_fn, 'w') as out:
        out.write(gdal_sr_wkt(int(epsg), True))
    return(0)
    
def gdal_set_epsg(src_gdal, epsg = 4326):
    '''set the projection of src_gdal to epsg'''
    ds = gdal.Open(src_gdal, gdal.GA_Update)
    if ds is not None:
        ds.SetProjection(gdal_sr_wkt(int(epsg)))
        ds = None
        return(0)
    else: return(None)
    
def gdal_gdal2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326):
    '''Convert the gdal file to gdal using gdal'''
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
        gdal2gdal_cmd = ('gdal_translate {} {} -f {}\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(gdal2gdal_cmd, verbose = True)
        if status == 0: return(dst_gdal)
        else: return(None)
    else: return(None)

def gdal_infos(src_gdal, scan = False):
    if os.path.exists(src_gdal):
        ds = gdal.Open(src_gdal)
        dsc = None
        if ds is not None:
            dsc = gdal_gather_infos(ds)
            if scan:
                t = ds.ReadAsArray()
                dsc['zmin'] = np.min(t)
                dsc['zmax'] = np.max(t)
        ds = None
        return(dsc)
    else: return(None)

def gdal_gather_infos(src_ds):
    '''gather information from `src_ds` GDAL dataset'''
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
    '''Set a datasource config dictionary'''
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def gdal_cpy_infos(src_config):
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def gdal_gt2region(ds_config):
    '''convert a gdal geo-tranform to an extent [w, e, s, n]'''
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + geoT[1] * ds_config['nx'], geoT[3] + geoT[5] * ds_config['ny'], geoT[3]])

def gdal_region2gt(region, inc):
    '''return a count info and a gdal geotransform based on extent and cellsize
    output is a list (xcount, ycount, geot)'''
    ysize = region[3] - region[2]
    xsize = region[1] - region[0]
    xcount = int(round(xsize / inc))
    ycount = int(round(ysize / inc))
    dst_gt = (region[0], inc, 0, region[3], 0, (inc * -1.))
    return(xcount, ycount, dst_gt)

def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.'''
    if dst_defn is None: dst_defn = src_layer.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for f in src_layer:
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

def gdal_inf(src_ds):
    '''generate an info (.inf) file from a src_gdal file using gdal'''
    minmax = [str(x) for x in gdal_region(src_ds)]
    with open('{}.inf'.format(src_ds.GetDescription()), 'w') as inf:
        inf.write('{}\n'.format(' '.join(minmax)))
    return(minmax)

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''convert a geographic x,y value to a pixel location of geoTransform'''
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt( geoTransform))
    return(int(round(pixel_x)), int(round(pixel_y)))

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
    this_end = _geo2pixel(region[1], region[2], ds_config['geoT'])
    this_size = (int(this_end[0] - this_origin[0]), int(this_end[1] - this_origin[1]))
    this_origin = [0 if x < 0 else x for x in this_origin]
    this_size = [0 if x < 0 else x for x in this_size]
    if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
    if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
    return(this_origin[0], this_origin[1], this_size[0], this_size[1])
                
def gdal_xyz2gdal(src_xyz, dst_gdal, region, inc, dst_format='GTiff', mode='n'):
    '''Create a GDAL supported grid from xyz data 
    `mode` of `n` generates a num grid
    `mode` of `m` generates a mean grid
    `mode` of `k` generates a mask grid'''
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    if mode == 'm': sumArray = np.zeros((ycount, xcount))
    ptArray = np.zeros((ycount, xcount))
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(4326), gdal.GDT_Int32, -9999, dst_format)
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                    if mode == 'm': sumArray[ypos, xpos] += float(z)
                    if mode == 'n' or mode == 'm': 
                        ptArray[ypos, xpos] += 1
                    else: ptArray[ypos, xpos] = 1
                except: pass
    if mode == 'm':
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
    '''blur an array'''
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
    '''gaussian blur on src_gdal'''
    ds = gdal.Open(src_gdal)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        msk_array = np.array(ds_array)
        msk_array[msk_array != ds_config['ndv']] = 1
        msk_array[msk_array == ds_config['ndv']] = np.nan
        ds_array[ds_array == ds_config['ndv']] = 0
        smooth_array = np_gaussian_blur(ds_array, int(sf))
        smooth_array = smooth_array * msk_array
        mask_array = ds_array = None
        smooth_array[np.isnan(smooth_array)] = ds_config['ndv']
        ds = None
        return(gdal_write(smooth_array, dst_gdal, ds_config))
    else: return([], -1)

def gdal_smooth(src_gdal, dst_gdal, fltr = 10, bathy_only = True, use_gmt = False):
    '''smooth `src_gdal` using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero)'''
    if os.path.exists(src_gdal):
        if bathy_only:
            dem_u, dem_l = gdal_split(src_gdal)
        else: dem_l = src_gdal
        if use_gmt:
            out, status = gmt_grdfilter(dem_l, 'tmp_fltr.tif=gd+n-9999:GTiff', dist = fltr, verbose = True)
        out, status = gdal_blur(dem_l, 'tmp_fltr.tif', fltr)
        if bathy_only:
            u_ds = gdal.Open(dem_u)
            if u_ds is not None:
                u_config = gdal_gather_infos(u_ds)
                l_ds = gdal.Open('tmp_fltr.tif')
                if l_ds is not None:
                    l_config = gdal_gather_infos(l_ds)
                    u_arr = u_ds.GetRasterBand(1).ReadAsArray()
                    l_arr = l_ds.GetRasterBand(1).ReadAsArray()
                    u_arr[u_arr == u_config['ndv']] = 0
                    l_arr[l_arr == l_config['ndv']] = 0
                    ds_arr = u_arr + l_arr
                    gdal_write(ds_arr, 'merged.tif', l_config)
                    l_ds = None
                    remove_glob(dem_l)
                u_ds = None
                remove_glob(dem_u)
            remove_glob('tmp_fltr.tif')
            os.rename('merged.tif', dst_gdal)
        else: os.rename('tmp_fltr.tif', dst_gdal)
    
def gdal_polygonize(src_gdal, dst_layer):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''
    ds = gdal.Open('{}'.format(src_gdal))
    ds_arr = ds.GetRasterBand(1)
    gdal.Polygonize(ds_arr, None, dst_layer, 0, [])
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
                if this_x_size == 0 or this_y_size == 0:
                    break
                
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

## ==============================================
## gdal processing (datalist fmt:200)
## ==============================================
def gdal_band2xyz(band_arr, dst_xyz = sys.stdout, gt = None, srcwin = None, nd = -9999, weight = 1, dump_nodata = False):
    if srcwin is None: srcwin = (0, 0, band_arr.size(0), band_arr.size(1))
    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        nodata = ['{:.10f}',format(-9999), 'nan']
        nodata.append('{:.10f}'.format(nd))
        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            geo_x, geo_y = _pixel2geo(x, y, gt)
            x_i_data = [band_arr[0][x_i]]
            z = x_i_data[0]
            line = [geo_x, geo_y, z, weight]

            if dump_nodata:
                xyz_line(line, dst_xyz)
            else:
                if '{:.10f}'.format(z) not in nodata:
                    xyz_line(line, dst_xyz)                    

def gdal_parse(src_ds, dst_xyz = sys.stdout, weight = None, dump_nodata = False, srcwin = None, want_yield = False):
    '''send the data from gdal file src_gdal to dst_xyz port'''
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_gather_infos(src_ds)
    gt = ds_config['geoT']        
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    
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
            x_i_data = [data[0][x_i]]
            z = x_i_data[0]
            line = [geo_x, geo_y, z]
            
            if dump_nodata:
                xyz_line(line, dst_xyz, weight = weight)
            else:
                if '{:.10f}'.format(z) not in nodata:
                    xyz_line(line, dst_xyz, weight = weight)
                    
def gdal_inf_entry(entry):
    ds = gdal.Open(entry[0])
    minmax = gdal_inf(ds)
    ds = None
                    
def gdal_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False):
    ds = gdal.Open(entry[0])
    if region is not None:
        srcwin = gdal_srcwin(ds, region)
    else: srcwin = None
    gdal_parse(ds, dst_xyz = dst_port, weight = entry[2], dump_nodata = False, srcwin = srcwin)
    ds = None
    
## ==============================================
## inf files (data info) inf.py
## mbsystem/gmt/waffles infos
## ==============================================
def inf_generate(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    return(inf_entry([data_path, data_fmt]))

def inf_parse(src_inf):
    '''parse an inf file (mbsystem or gmt) and return minmax'''
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

def inf_entry(src_entry, overwrite = False):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file
    or the result of `gmt gmtinfo file.xyz -C`, which is
    a 6 column line with minmax info, etc.
    returns the region of the inf file.'''
    minmax = None
    if os.path.exists(src_entry[0]):
        path_i = src_entry[0] + '.inf'
        if not os.path.exists(path_i) or overwrite:
            _datalist_inf_hooks[src_entry[1]](src_entry)
        if os.path.exists(path_i): minmax = inf_parse(path_i)[:4]
    if region_valid_p(minmax):
        return(minmax)
    else: return(None)
    
## ==============================================
## xyz processing (datalists fmt:168)
## ==============================================
_xyz_config = {'delim':' ', 'xpos': 0, 'ypos': 1, 'zpos': 2}

def xyz_parse(src_xyz):
    '''xyz file parsing generator'''
    for xyz in src_xyz:
        this_line = xyz.strip()
        for delim in _known_delims:
            this_xyz = this_line.split(delim)
            if len(this_xyz) > 1:
                this_delim = delim
                break
        yield([float(x) for x in this_xyz])

def xyz_block(src_xyz, region, inc):
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    sum_arr = np.zeros((ycount, xcount))
    pt_arr = np.zeros((ycount, xcount))
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    pt_arr[ypos, xpos] = 1
                    sum_arr[ypos, xpos] += float(z)
                    pt_arr[ypos, xpos] += 1
                except: pass
    block_arr = sum_arr / pt_arr
    srcwin = (0, 0, xcount, ycount)
    ## nodata
    gdal_band2xyz(block_arr, sys.stdout, gt, srcwin)
        
def xyz_line(line, dst_port = sys.stdout, weight = None):
    '''write XYZ `line` to `dst_port` using `delimiter` and `weight`'''
    if weight is None:
        w_string = ''
    else: w_string = '{}{}'.format(_xyz_config['delim'], weight)
    l = '{}{}\n'.format(_xyz_config['delim'].join([str(x) for x in line]), w_string)#.encode(dst_port.encoding)
    if dst_port != sys.stdout: l = l.encode('utf-8')
    dst_port.write(l)
    
def xyz_inf(src_xyz):
    '''return minmax info from a src_xyz file.'''
    minmax = []
    for i,l in enumerate(xyz_parse(src_xyz)):
        if i == 0:
            minmax = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            if l[0] < minmax[0]: minmax[0] = l[0]
            elif l[0] > minmax[1]: minmax[1] = l[0]
            if l[1] < minmax[2]: minmax[2] = l[1]
            elif l[1] > minmax[3]: minmax[3] = l[1]
            if l[1] < minmax[2]: minmax[4] = l[2]
            elif l[1] > minmax[3]: minmax[5] = l[2]
    with open('{}.inf'.format(src_xyz.name), 'w') as inf:
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)

def xyz_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n], else False.'''
    if src_xy[0] < src_region[0]: return(False)
    elif src_xy[0] > src_region[1]: return(False)
    elif src_xy[1] < src_region[2]: return(False)
    elif src_xy[1] > src_region[3]: return(False)
    else: return(True)

def xyz_inf_entry(entry):
    with open(entry[0]) as infile:
        minmax = xyz_inf(infile)
        
def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False):
    '''dump ascii xyz data to dst_port'''
    with open(entry[0]) as infile:
        for line in xyz_parse(infile):
            if region is not None:
                if xyz_in_region_p(line, region):
                    xyz_line(line, dst_port, entry[2])
            else: xyz_line(line, dst_port, entry[2])

## ==============================================
## datalist processing (datalists fmt:-1)
## ==============================================
_datalist_pass_hooks = {-1: lambda e: os.path.exists(e[0]), 168: lambda e: os.path.exists(e[0]), 200: lambda e: os.path.exists(e[0])}
_datalist_inf_hooks = {-1: lambda e: inf_datalist(e), 168: lambda e: xyz_inf_entry(e), 200: lambda e: gdal_inf_entry(e)}
_datalist_hooks = {-1: lambda e, v: datalist(e[0], verbose = v), 168: lambda e, v: xyz_dump_entry(e, verbose = v), 200: lambda e, v: gdal_dump_entry(e, verbose = v)}

def datalist_inf(dl, fmt = None):
    '''return the region of the datalist.'''
    out_regions = []
    if fmt is None:
        _datalist_hooks[168] = lambda e, v: out_regions.append(inf_entry(e))
        _datalist_hooks[200] = lambda e, v: out_regions.append(inf_entry(e))
    else:
        try:
            _datalist_hooks[fmt] = lambda e, v: out_regions.append(inf_entry(e))
        except KeyError as e:
            echo_error('bad fmt, {}'.format(e))
    datalist(dl)
    out_regions = [x for x in out_regions if x is not None]
    if len(out_regions) == 0:
        return(None)
    elif len(out_regions) == 1:
        print(out_regions[0])
        return(out_regions[0])
    else:
        out_region = out_regions[0]
        for x in out_regions[1:]:
            out_region = regions_merge(out_region, x)
        return(out_region)

def inf_datalist(e):
    '''write an inf file for datalist entry e'''
    minmax = datalist_inf(e[0])
    if minmax is not None:
        echo_msg('generating inf for datalist {}'.format(e[0]))
        with open('{}.inf'.format(e[0]), 'w') as inf:
            inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))

def datalist_append(entry, datalist):
    '''append entry to datalist'''
    with open(datalist, 'a') as outfile:
        outfile.write('{}\n'.format(' '.join([str(x) for x in entry])))
            
def datalist_dump(dl, dst_port = sys.stdout, region = None, verbose = False):
    '''dump the data from datalist to dst_port'''
    if region is not None:
        #_datalist_pass_hooks[-1] = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
        _datalist_pass_hooks[168] = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
        _datalist_pass_hooks[200] = lambda e: regions_intersect_ogr_p(region, inf_entry(e))
    _datalist_hooks[-1] = lambda e, v: datalist(e[0], verbose = v)
    _datalist_hooks[168] = lambda e, v: xyz_dump_entry(e, region = region, dst_port = dst_port, verbose = v)
    _datalist_hooks[200] = lambda e, v: gdal_dump_entry(e, region = region, dst_port = dst_port, verbose = v)
    datalist(dl, verbose = verbose)

def datalist_echo(entry):
    '''echo datalist entry'''
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))
    datalist(entry[0])

def datafile_echo(entry):
    '''echo datafile entry'''
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))
    
def entry2py(dle):
    '''convert a datalist entry to python'''
    for delim in _known_dl_delims:
        this_entry = dle.rstrip().split(delim)
        if len(this_entry) > 1: break
    try:
        entry = [x if n == 0 else float(x) if n < 3 else x for n, x in enumerate(this_entry)]
    except ValueError as e: return(None)
    return(entry)
    
def datalist2py(dl):
    '''convert a datalist to python data'''
    these_entries = []
    this_entry = entry2py(dl)
    if len(this_entry) < 2: this_entry.append(-1)
    if len(this_entry) < 3: this_entry.append(1)
    if this_entry[1] == -1:
        with open(this_entry[0], 'r') as op:
            while True:
                this_line = op.readline().rstrip()
                if not this_line: break
                if this_line[0] != '#' and this_line[0] != '\n':
                    these_entries.append(entry2py(this_line))
    else: these_entries.append(this_entry)
    return(these_entries)

def datalist_io(dpipe):
    '''`dpipe` is a function to write to `f` bytesio (e.g. `datalist_dump`); we then yield that data through `xyz_parse`
    use as: 
    for xyz in datalist_io(lambda p: datalist_dump('datalist.datalist', dst_port = p)):
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        ...'''
    f = io.BytesIO()
    d = dpipe(f)
    f.seek(0)
    for this_xyz in xyz_parse(f):
        yield(this_xyz)
    f.seek(0)

def datalist(dl, fmt = -1, wt = 1, verbose = False):
    '''recurse a datalist/entry'''
    this_dir = os.path.dirname(dl)
    these_entries = datalist2py(dl)
    if len(these_entries) == 0: these_entries = [entry2py(dl)]
    for this_entry in these_entries:
        this_entry[0] = os.path.join(this_dir, this_entry[0])
        this_entry[2] = wt * this_entry[2]
        this_entry_md = ' '.join(this_entry[3:]).split(',')
        this_entry = this_entry[:3] + [this_entry_md] + [os.path.basename(dl).split('.')[0]]
        if _datalist_pass_hooks[this_entry[1]](this_entry):
            if verbose is True:
                echo_msg('{} {}'.format('scanning datalist' if this_entry[1] == -1 else 'using datafile', this_entry[0]))
            _datalist_hooks[this_entry[1]](this_entry, verbose)

## ==============================================
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT and GDAL
## ==============================================
_waffles_grid_info = {
    'datalist': None,
    'region': None,
    'inc': None,
    'name': 'waffles',
    'node': 'pixel',
    'fmt': 'GTiff',
    'extend': 0,
    'weights': None,
    'fltr': None,
    'clip': None,
    'epsg': 4326,
    'mod': 'surface',
    'mod_args': (),
    'gc': config_check()
}
_waffles_modules = {
    'surface': [lambda args: waffles_gmt_surface(**args), '''SPLINE DEM via GMT surface
    \t\t\t< surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
    \t\t\t:tension=[0-1] - Spline tension.'''],
    'triangulate': [lambda args: waffles_gmt_triangulate(**args), '''TRIANGULATION DEM via GMT triangulate'''],
    'nearneighbor': [lambda args: waffles_gmt_nearneighbor(**args), '''NEAREST NEIGHBOR DEM via GMT nearneighbor
    \t\t\t< nearneighbor:radius=6s >
    \t\t\t:radius=[value] - Nearest Neighbor search radius'''],
    'num': [lambda args: waffles_num(**args), '''Uninterpolated DEM populated by mode.
    \t\t\t< num:mode=n >
    \t\t\t:mode=[key] - specify mode of grid population: k (mask), m (mean) or n (num)'''],
    'spatial-metadata': [lambda args: waffles_spatial_metadata(**args), '''Datalist SPATIAL METADATA <beta>
    \t\t\t< spatial-metadata:inc=None >
    \t\t\t:inc=[increment] - Spatial metadata resolution [default uses -E value]'''],
}
waffles_module_desc = lambda x: '\n  '.join(['{:18}{}'.format(key, x[key][-1]) for key in x])
waffles_proc_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * 10) + (wg['inc'] * wg['extend']))
waffles_dist_region = lambda wg: region_buffer(wg['region'], (wg['inc'] * wg['extend']))
waffles_proc_str = lambda wg: region_format(waffles_proc_region(wg), 'gmt')
waffles_dl_func = lambda wg: lambda p: datalist_dump(wg['datalist'], dst_port = p, region = waffles_proc_region(wg), verbose = True)
waffles_gmt_reg_str = lambda wg: '-r' if wg['node'] == 'pixel' else ''

def waffles_gmt_surface(wg = _waffles_grid_info, tension = .35, relaxation = 1.2, lower_limit = 'd', upper_limit = 'd'):
    '''generate a DEM with GMT surface'''
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the SURFACE module')
        return(None, -1)
    dem_surf_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} {}\
    '.format(waffles_proc_str(wg), wg['inc'], waffles_gmt_reg_str(wg), waffles_proc_str(wg), wg['inc'], wg['name'], tension, relaxation, lower_limit, upper_limit, waffles_gmt_reg_str(wg)))
    return(run_cmd(dem_surf_cmd, verbose = True, data_fun = waffles_dl_func(wg)))

def waffles_gmt_triangulate(wg = _waffles_grid_info):
    '''generate a DEM with GMT surface'''
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the TRIANGULATE module')
        return(None, -1)
    dem_tri_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt triangulate {} -I{:.10f} -V -G{}.tif=gd+n-9999:GTiff {}\
    '.format(waffles_proc_str(wg), wg['inc'], waffles_gmt_reg_str(wg), waffles_proc_str(wg), wg['inc'], wg['name'], waffles_gmt_reg_str(wg)))
    return(run_cmd(dem_tri_cmd, verbose = True, data_fun = waffles_dl_func(wg)))

def waffles_gmt_nearneighbor(wg = _waffles_grid_info, radius = None):
    '''genearte a DEM with GMT nearneighbor'''
    if wg['gc']['GMT'] is None:
        echo_error_msg('GMT must be installed to use the NEARNEIGHBOR module')
        return(None, -1)
    if radius is None: radius = wg['inc'] * 2
    dem_nn_cmd = ('gmt blockmean {} -I{:.10f} -Wi -V {} | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}.tif=gd+n-9999:GTiff {}\
    '.format(waffles_proc_str(wg), wg['inc'], waffles_gmt_reg_str(wg), waffles_proc_str(wg), wg['inc'], radius, wg['name'], waffles_gmt_reg_str(wg)))
    return(run_cmd(dem_nn_cmd, verbose = True, data_fun = waffles_dl_func(wg)))

def waffles_num(wg = _waffles_grid_info, mode = 'n'):
    '''Generate an uninterpolated num grid.
    mode of `k` generates a mask grid
    mode of `m` generates a mean grid
    mode of `n` generates a num grid'''
    wg['region'] = region_buffer(wg['region'], wg['inc'] * .5) if wg['node'] == 'grid' else wg['region']
    return(gdal_xyz2gdal(datalist_io(waffles_dl_func(wg)), '{}.tif'.format(wg['name']), waffles_dist_region(wg), wg['inc'], dst_format = wg['fmt'], mode = mode))
    
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
        these_entries = datalist2py(entry[0])
        
        for this_entry in these_entries:
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

def waffles_wg_valid_p(wg = _waffles_grid_info):
    '''return True if wg appears valid'''
    if wg['datalist'] is None: return(False)
    elif wg['region'] is None: return(False)
    elif wg['inc'] is None: return(False)
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
    '''generate a DEM using wg dict settings'''
    dem = '{}.tif'.format(wg['name'])
    args_d = {}
    for arg in wg['mod_args']:
        p_arg = arg.split('=')
        args_d[p_arg[0]] = p_arg[1]
    args_d['wg'] = wg
    
    ## ==============================================
    ## gererate the DEM
    ## ==============================================
    try:
        out, status = _waffles_modules[wg['mod']][0](args_d)
    except TypeError as e: echo_error_msg('{}'.format(e))
    except KeyboardInterrupt as e:
        echo_error_msg('killed by user, {}'.format(e))
        sys.exit(-1)

    waffles_gdal_md(wg)
        
    if wg['mod'] == 'spatial-metadata': sys.exit(0)
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
    ## optionally filter the DEM batymetry (sub-zero)
    ## ==============================================
    if wg['fltr'] is not None:
        try:
            gdal_smooth(dem, 'tmp_s.tif', fltr = wg['fltr'], bathy_only = True, use_gmt = True if wg['gc']['GMT'] is not None else False)
            os.rename('tmp_s.tif', dem)
        except TypeError as e: echo_error_msg('{}'.format(e))
        
    ## ==============================================
    ## optionally clip the DEM to polygon
    ## ==============================================
    if wg['clip'] is not None:
        clip_args = {}
        cp = wg['clip'].split(':')
        clip_args['src_gdal'] = dem
        clip_args['src_ply'] = cp[0]
        cargs = cp[1:]
        for arg in cargs:
            p_arg = arg.split('=')
            clip_args[p_arg[0]] = p_arg[1]
        gdal_clip(**clip_args)
        
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
## waffles cli
## ==============================================
waffles_cli_usage = '''waffles [OPTIONS] <datalist/entry>

General Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
\t\t\tIf omitted, use the region gathered from the data in DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -M, --module\t\tDesired DEM MODULE and options. (see available Modules below)
\t\t\tsyntax is -M module:mod_opt=mod_val:mod_opt1=mod_val1:...
  -O, --output-name\tBASENAME for all outputs.
  -P, --epsg\t\tHorizontal projection of data as EPSG code [4326]
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION.
\t\t\tappend :<num> to extend the processing region: -X6:12
  -T, --filter\t\tFILTER the output (bathymetry-only) using a Cosine Arch Filter at -T<dist(km)> search distance.
\t\t\tIf GMT is not available, will perform a Gaussian Blur at -T<factor>. 
  -C, --clip\t\tCLIP the output to the clip polygon. [clip_ply.shp:invert=False]

  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tuse grid-node registration, default is pixel-node

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Modules:
  {}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(waffles_module_desc(_waffles_modules))

def waffles_cli(argv = sys.argv):
    wg = _waffles_grid_info
    dls = []
    region = None
    module = None
    want_prefix = False
    want_verbose = False
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
            wg['inc'] = gmt_inc2inc(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-E': wg['inc'] = gmt_inc2inc(arg[2:])
        elif arg == '--outname' or arg == '-O':
            wg['name'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': wg['name'] = arg[2:]
        elif arg == '--format' or arg == '-F':
            wg['fmt'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': wg['fmt'] = arg[2:]
        elif arg == '--filter' or arg == '-T':
            wg['fltr'] = gmt_inc2inc(argv[i + 1])
            i += 1
        elif arg[:2] == '-T': wg['fltr'] = gmt_inc2inc(arg[2:])
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
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = arg[2:]
        elif arg == '-w': wg['weight'] = True
        elif arg == '-p': want_prefix = True
        elif arg == '-r': wg['node'] = 'grid'
        elif arg == '--verbose' or arg == '-V': want_verbose = True
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(waffles_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format(_version))
            sys.exit(0)
        else: dls.append(arg)
        i += 1

    if len(dls) == 0:
        sys.stderr.write(waffles_cli_usage)
        echo_error_msg('''must specify a datalist/entry file''')
        sys.exit(-1)

    ## ==============================================
    ## set the master datalist
    ## ==============================================
    master = '.master.datalist'
    with open(master, 'w') as md:        
        for dl in dls:
            if os.path.exists(dl):
                if len(dl.split(':')) == 1:
                    for key in _known_datalist_fmts.keys():
                        if dl.split(':')[0].split('.')[-1] in _known_datalist_fmts[key]:
                            md.write('{} {} 1\n'.format(dl, key))
    if os.stat(master).st_size == 0:
        remove_glob('.master.datalist')
        sys.stderr.write(waffles_cli_usage)
        echo_error_msg('bad datalist/entry, {}'.format(dls))
        sys.exit(-1)        
    wg['datalist'] = master
    
    ## ==============================================
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
                        
    ## ==============================================
    ## reformat and set the region
    ## ==============================================
    if region is None: these_regions = [datalist_inf(dl)]
    else:
        try:
            these_regions = [[float(x) for x in region.split('/')]]
        except ValueError: these_regions = gdal_ogr_regions(region)

    ## ==============================================
    ## run waffles for each input region.
    ## ==============================================
    bn = wg['name']
    for this_region in these_regions:
        wg['region'] = this_region
        if want_prefix: wg['name'] = '{}{}_{}_{}'.format(bn, inc2str_inc(wg['inc']), region_format(wg['region'], 'fn'), this_year())    
        echo_msg(wg)
        waffles_run(wg)

    ## ==============================================
    ## cleanup waffles detritus
    ## ==============================================
    remove_glob(master)

## ==============================================
## datalists cli
## ==============================================
datalists_cli_usage = '''datalists [-iwR] <datalist/entry>

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''

def datalists_cli(argv = sys.argv):
    dls = []
    region = None
    want_weight = False
    want_infos = False
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
    if region is not None:
        try:
            region = [float(x) for x in region.split('/')]
        except ValueError as e:
            sys.stderr.write(datalists_cli_usage)
            echo_error_msg('bad region, {}'.format(e))
            sys.exit(-1)
        _datalist_pass_hooks[-1] = lambda e: regions_intersect_p(region, inf_entry(e))
        _datalist_pass_hooks[168] = lambda e: regions_intersect_p(region, inf_entry(e))
        _datalist_pass_hooks[200] = lambda e: regions_intersect_p(region, inf_entry(e))
        _datalist_hooks[168] = lambda e, v: xyz_dump_entry(e, region = region, verbose = want_verbose)
        _datalist_hooks[200] = lambda e, v: gdal_dump_entry(e, region = region, verbose = want_verbose)
        
    if want_infos:
        #print(datalist_inf(dl))
        _datalist_hooks[168] = lambda e, v: datafile_echo(e)
        _datalist_hooks[200] = lambda e, v: datafile_echo(e)
        if region is None:
            _datalist_pass_hooks[-1] = lambda e: inf_entry(e, True)
            #_datalist_pass_hooks[168] = lambda e: inf_entry(e)
            #_datalist_pass_hooks[200] = lambda e: inf_entry(e)
        
    ## ==============================================
    ## recurse the datalist
    ## ==============================================
    master = '.geomods-master.datalist'
    with open(master, 'w') as md:
        for dl in dls:
            if os.path.exists(dl):
                if len(dl.split(':')) == 1:
                    for key in _known_datalist_fmts.keys():
                        if dl.split(':')[0].split('.')[-1] in _known_datalist_fmts[key]:
                            md.write('{} {} 1\n'.format(dl, key))
    try:
        #if want_infos:
        #    print(datalist_inf(dl))
        #else:
        datalist(master, verbose = want_verbose)
    except IndexError:
        echo_error_msg('bad datalist file...')
    except KeyboardInterrupt as e:
        echo_error_msg('user killed process')
        remove_glob(master)
        sys.exit(-1)
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
