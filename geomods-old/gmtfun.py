### gmtfun.py
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

from utils import *

def gmt_inc2inc(inc_str):
    units = inc_str[-1]

    if units == 'c': ## arc-seconds (old)
        inc = float(inc_str[:-1]) / 3600
    elif units == 's': ## arc-seconds
        inc = float(inc_str[:-1]) / 3600
    elif units == 'm': ## arc-minutes
        inc = float(inc_str[:-1]) / 360
    else: inc = float(inc_str)    
    
    return(inc)

## =============================================================================
##
## GMT Wrapper Functions - gmtfun.py
## wrapper functions to GMT system commands
##
## =============================================================================

def gmt_inf(src_xyz):
    '''generate an info (.inf) file from a src_xyz file using GMT.'''
    if os.path.exists(src_grd):
        return(run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), verbose = False))
    else: return(None)

def gmt_grd_inf(src_grd):
    '''generate an info (.inf) file from a src_grd file using GMT.'''
    if os.path.exists(src_grd):
        return(run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_grd, src_grd), verbose = False))
    else: return(None)

def gmt_grd2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, verbose = False):
    '''Convert the grd file to tif using GMT'''
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdalfun._fext(dst_fmt))
        grdc_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = run_cmd(grdc_cmd, verbose = verbose)
        if status != 0: dst_gdal = None
    else: dst_gdal = None
    return(dst_gdal)

def grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    if os.path.exists(src_grd):
        grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
        out, status = run_cmd(grdinfo_cmd, verbose = verbose)
        remove_glob('gmt.conf')
        if status == 0:
            return(out.split())
        else: return(None)
    else: return(None)

def gmtinfo(src_xyz, verbose = False):
    '''Return an info list of `src_xyz`'''
    out, status = run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = verbose)
    if os.path.exists(src_xyz):
        gmtinfo_cmd = ('gmt gmtinfo {} -C'.format(src_xyz))
        out, status = run_cmd(gmtinfo_cmd, verbose = verbose)
        remove_glob('gmt.conf')
        if status == 0:
            return(out.split())
        else: return(None)
    else: return(None)

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
        out, status = run_cmd(cut_cmd1, verbose = verbose)
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
### End
