### gmtfun.py
##
## Copyright (c) 2020 - 2021 CIRES Coastal DEM Team
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

## ==============================================
## GMT Wrapper Functions - gmtfun.py
## wrapper functions to GMT system commands
##
## GMT must be installed on the system to run these
## functions and commands.
## ==============================================

import os
from geomods import utils
from geomods import regions
from geomods import gdalfun

def gmt_inf(src_xyz):
    '''generate an info (.inf) file from a src_xyz file using GMT.

    returns [cmd-output, cmd-return-code]'''
    
    return(utils.run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), verbose = False))

def gmt_inf_parse(src_inf):
    xyzi = {'name': src_inf, 'numpts': 0, 'minmax': [0,0,0,0,0,0], 'wkt': gdalfun.gdal_region2wkt([0,0,0,0,0,0])}
    with open(src_inf) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                xyzi['minmax'] = [float(x) for x in til]
    xyzi['wkt'] = gdalfun.gdal_region2wkt(xyzi['minmax'])
    return(xyzi)

def gmt_grd_inf(src_grd):
    '''generate an info (.inf) file from a src_gdal file using GMT.

    returns [cmd-output, cmd-return-code]'''
    
    return(utils.run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_grd, src_grd), verbose = False))

def gmt_inc2inc(inc_str):
    '''convert a GMT-style `inc_str` (6s) to geographic units
    c/s - arc-seconds
    m - arc-minutes

    return float increment value.'''
    
    if inc_str is None or inc_str.lower() == 'none': return(None)
    units = inc_str[-1]
    if units == 'c': inc = float(inc_str[:-1]) / 3600.
    elif units == 's': inc = float(inc_str[:-1]) / 3600.
    elif units == 'm': inc = float(inc_str[:-1]) / 360.
    else:
        try:
            inc = float(inc_str)
        except ValueError as e:
            utils.echo_error_msg('could not parse increment {}, {}'.format(inc_str, e))
            return(None)
    return(inc)

def gmt_grd2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''convert the grd file to tif using GMT

    returns the gdal file name or None'''
    
    dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdalfun.gdal_fext(dst_fmt))
    grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
    '.format(src_grd, dst_gdal, dst_fmt))
    out, status = utils.run_cmd(grd2gdal_cmd, verbose = verbose)
    if status == 0:
        return(dst_gdal)
    else: return(None)

def gmt_grdinfo(src_grd, verbose = False):
    '''gather infos about src_grd using GMT grdinfo.

    return an info list of `src_grd`'''
    
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
    out, status = utils.run_cmd(grdinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)

def gmt_gmtinfo(src_xyz, verbose = False):
    '''gather infos about src_xyz using GMT gmtinfo

    return an info list of `src_xyz`'''
    
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    gmtinfo_cmd = ('gmt gmtinfo {} -C'.format(src_xyz))
    out, status = utils.run_cmd(gmtinfo_cmd, verbose = verbose)
    remove_glob('gmt.conf')
    if status == 0:
        return(out.split())
    else: return(None)
        
def gmt_select_split(o_xyz, sub_region, sub_bn, verbose = False):
    '''split an xyz file into an inner and outer region.

    returns [inner_region, outer_region]'''
    
    out_inner = None
    out_outer = None
    gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, regions.region_format(sub_region, 'gmt'), sub_bn)
    out, status = utils.run_cmd(gmt_s_inner, verbose = verbose)
    if status == 0: out_inner = '{}_inner.xyz'.format(sub_bn)
    gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, regions.region_format(sub_region, 'gmt'), sub_bn)
    out, status = utils.run_cmd(gmt_s_outer, verbose = verbose)
    if status == 0:  out_outer = '{}_outer.xyz'.format(sub_bn)
    return([out_inner, out_outer])
        
def gmt_grdcut(src_grd, src_region, dst_grd, verbose = False):
    '''cut `src_grd` to `src_region` using GMT grdcut
    
    returns [cmd-output, cmd-return-code]'''
    
    cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
    return(utils.run_cmd(cut_cmd1, verbose = True))

def gmt_grdfilter(src_grd, dst_grd, dist = '3s', node = 'pixel', verbose = False):
    '''filter `src_grd` using GMT grdfilter

    returns [cmd-output, cmd-return-code]'''
    
    #ft_cmd1 = ('gmt grdfilter -V {} -G{} -R{} -Fc{} -D1{}'.format(src_grd, dst_grd, src_grd, dist, ' -r' if node == 'pixel' else ''))
    ft_cmd1 = ('gmt grdfilter -V {} -G{} -Fc{} -D1{}'.format(src_grd, dst_grd, dist, ' -r' if node == 'pixel' else ''))
    return(utils.run_cmd(ft_cmd1, verbose = verbose))

def gmt_nan2zero(src_grd, node = 'pixel', verbose = False):
    '''convert nan values in `src_grd` to zero

    returns status code (0 == success) '''
    
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = _tmp.tif=gd+n-9999:GTiff'.format(src_grd))
    out, status = utils.run_cmd(num_msk_cmd, verbose = True)
    if status == 0: os.rename('_tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_grdcut(src_grd, region, verbose = False):
    '''cut a grid to region using GMT grdcut

    return status code (0 == success)'''
    
    cut_cmd = ('gmt grdcut -V {} -G_tmp.grd {}\
    '.format(src_grd, region_format(region, 'gmt')))
    out, status = utils.run_cmd(cut_cmd, verbose = True)
    if status == 0:
        remove_glob(src_grd)
        os.rename('_tmp.grd', '{}'.format(src_grd))
    return(status)

def gmt_slope(src_dem, dst_slp, verbose = False):
    '''generate a Slope grid from a DEM with GMT

    return status code (0 == success)'''
    
    o_b_name = '{}'.format(src_dem.split('.')[0])
    slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
    '.format(src_dem, o_b_name, src_dem))
    out, status = utils.run_cmd(slope_cmd0, verbose = verbose)
    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}=gd+n-9999:GTiff\
        '.format(o_b_name, dst_slp))
        out, status = utils.run_cmd(slope_cmd1, verbose = verbose)
    remove_glob('{}_pslp.grd'.format(o_b_name))
    return(status)

def gmt_num_msk(num_grd, dst_msk, verbose = False):
    '''generate a num-msk from a NUM grid using GMT grdmath

    returns [cmd-output, cmd-return-code]'''
    
    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    return(utils.run_cmd(num_msk_cmd, verbose = verbose))

def gmt_sample_gnr(src_grd, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node
    grid registration.

    returns status code (0 == success)'''
    
    out, status = utils.run_cmd('gmt grdsample -T {} -G_tmp.tif=gd+n-9999:GTiff'.format(src_grd), verbose = verbose)
    if status == 0: os.rename('_tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_sample_inc(src_grd, inc = 1, verbose = False):
    '''resamele src_grd to increment `inc` using GMT grdsample

    returns status code (0 == success)'''
    
    out, status = utils.run_cmd('gmt grdsample -I{:.10f} {} -R{} -G_tmp.tif=gd+n-9999:GTiff'.format(inc, src_grd, src_grd), verbose = verbose)
    if status == 0: os.rename('_tmp.tif', '{}'.format(src_grd))
    return(status)

def gmt_yield_entry(entry, region = None, verbose = False, z_region = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''
    ln = 0
    delim = ' '
    if z_region is not None: z_region = ['-' if x is None else str(x) for x in z_region]
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    for line in utils.yield_cmd('gmt gmtselect {} {} {}\
    '.format(entry[0], '' if region is None else regions.region_format(region, 'gmt'), '' if z_region is None else '-Z{}'.format('/'.join(z_region))),
                          data_fun = None, verbose = False):
        ln += 1
        yield(line)
    if verbose: utils.echo_msg('read {} data points from {}'.format(ln, entry[0]))

### End
