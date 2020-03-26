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
### Code:

import sys
import os

import threading
import time

import regions
import datalists
import gdalfun
import uncertainty
import metadata
import utils
        
_version = '0.2.7'

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)'''
    
    import fractions

    str_inc = str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', '')

    return(str_inc)

def this_year():
    '''return the current year'''

    import datetime

    return(datetime.datetime.now().strftime('%Y'))

## GMT Wrapper Functions

def grd2gdal(src_grd, dst_fmt = 'GTiff', verbose = False):
    '''Convert the grd file to tif using GMT'''
    
    status = 0
    if os.path.exists(src_grd):
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdalfun._fext(dst_fmt))

        grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = utils.run_cmd(grd2gdal_cmd, verbose, verbose)

        if status != 0:
            dst_gdal = None
    else: dst_gdal = None

    return(dst_gdal)

def grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''

    status = 0
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose, True)
    if not os.path.exists(src_grd):
        return([])
    else:
        grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
        out, status = utils.run_cmd(grdinfo_cmd, verbose, False)
        if status !=0:
            return([])
        else: return(out.split())

def grdcut(src_grd, src_region, dst_grd, verbose = False):
    '''Cut `src_grd` to `src_region` '''

    status = 0
    if os.path.exists(src_grd):
        cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
        out, status = utils.run_cmd(cut_cmd1, verbose, verbose)

    return(status)

def grd2xyz(src_grd, dst_xyz, region = None, mask = None, verbose = False, want_datalist = False):
    '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region.
    Optionally, generate a datalist and inf file for the resultant xyz data.'''

    status = 0
    if mask:
        grdmask_cmd = ('gmt grdmath -N -V {} {} OR = tmp.grd'.format(src_grd, mask))
        out, status = utils.run_cmd(grdmask_cmd, verbose, verbose)
        if status == 0: 
            src_grd = 'tmp.grd'

    if region and region._valid:
        region_str = region.gmt
    else: region_str = ''

    grd2xyz_cmd = ('gmt grd2xyz -V {} -s {} > {}'.format(src_grd, region_str, dst_xyz))
    out, status = utils.run_cmd(grd2xyz_cmd, verbose, verbose)

    if status == 0:
        if mask:
            if os.path.exists('tmp.grd'):
                os.remove('tmp.grd')

        if want_datalist:
            xyz_dir = os.path.join(os.getcwd(), 'xyz')
            if not os.path.exists(xyz_dir):
                os.makedirs(xyz_dir)

            s_datalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(dst_xyz.split('.')[0])))
            s_datalist._append_datafile('{}'.format(os.path.basename(dst_xyz)), 168, 1)
            s_datalist._reset()
        
            out, status = utils.run_cmd('mbdatalist -O -I{}'.format(s_datalist._path), False, True)
        
    return(status)

def slope(src_dem, dst_slp, verbose = False):
    '''Generate a Slope grid from a DEM with GMT'''

    status = 0
    o_b_name = '{}'.format(src_dem.split('.')[0])

    slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
    '.format(src_dem, o_name, src_dem))
    out, status = utils.run_cmd(slope_cmd0, verbose, True)

    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}\
        '.format(o_b_name, dst_slp))
        out, status = utils.run_cmd(slope_cmd1, verbose, True)
        
    if os.path.exists('{}_pslp.grd'.format(o_b_name)):
        os.remove('{}_pslp.grd'.format(o_b_name))

    return(status)

def proximity(src_num_msk, dst_prox):
    '''Generate a Proximity grid with GDAL from a NUM-MSK'''

    status = 0
    gdalfun.proximity(src_num_msk, dst_prox)

    if not os.path.exists(dst_prox):
        status = -1

    return(status)

def num_msk(num_grd, dst_msk, verbose = False):
    '''Generate a num-msk from a NUM grid.'''

    status = 0

    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    out, status = utils.run_cmd(num_msk_cmd, verbose, True)

    return(status)

def xyz2grd(datalist, region, inc, dst_name, a = 'n', node = 'pixel', verbose = False):
    '''Run the GMT command `xyz2grd` given a datalist, region and increment.'''
    
    status = 0
    if node == 'pixel':
        reg_str = '-r'
    else: reg_str = ''
    
    num_cmd0 = ('gmt xyz2grd -V {} -I{:.10f} -G{} -A{} {}\
    '.format(region.gmt, inc, dst_name, a, reg_str))
    out, status = utils.run_cmd(num_cmd0, verbose, verbose, datalist._dump_data)

    return(out, status)

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):
    '''Run the MBSystem command `mbgrid` given a datalist, region and increment.'''
    
    status = 0
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1:
        dist = dist + '/2'
    
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.10f}/{:.10f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} > mb_proc.txt \
    '.format(datalist._path, region.gmt, inc, inc, dst_name, dist, tension, e_switch))
    out, status = utils.run_cmd(mbgrid_cmd, verbose, verbose)

    return(out, status)

## =============================================================================
##
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT and GDAL 
##
## =============================================================================

class dem:
    '''Generate a Digital Elevation Model using one of the dem modules.
    DEM Modules include, `mbgrid`, `surface`, `num`, `mean`, `bathy`'''

    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None, o_b_name = None, \
                 o_node = 'pixel', o_fmt = 'GMT', o_extend = 6, callback = lambda: False, verbose = False):
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        
        self.node = o_node
        self.o_fmt = o_fmt
        self.extend = int(o_extend)

        self.proc_region = self.region.buffer((self.extend * 2) * self.inc)
        self.dist_region = self.region.buffer(self.extend * self.inc)
        
        self.status = 0
        self.stop = callback
        self.verbose = verbose

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]

            str_inc = inc2str_inc(self.inc)
            self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        else: self.o_name = o_b_name

        self.dem = '{}.grd'.format(self.o_name)
        
    def run(self, dem_mod = 'mbgrid', args = ()):
        '''Run the DEM module `dem_mod` using args `dem_mod_args`.'''

        #if dem_mod != 'mbgrid' and dem_mod != 'conversion-grid' and dem_mod != 'spatial-metadata' and dem_mod != 'uncertainty':
        #    if len(self.datalist.datafiles) == 0:
        #        self.datalist._load_data()
        #    if len(self.datalist.datafiles) == 0:
        #        self.status = -1

        if self.status == 0:
            _dem_mods[dem_mod][0](self)(*args)

        if self.status == 0:
            if self.o_fmt != 'GMT':
                self.dem = grd2gdal(self.dem)

        if self.status != 0:
            self.dem = None
                
        return(self.dem)

    ## ==============================================
    ## run mbgrid on the datalist and generate the DEM
    ## note: mbgrid will cause popen to hang if stdout 
    ## is not cleared...should add output file to send to...
    ## ==============================================

    def mbgrid(self, dist = '10/3'):
        '''Generate a DEM with MBSystem's mbgrid program.'''

        out, self.status = run_mbgrid(self.datalist, self.dist_region, self.inc, self.o_name, dist = dist, verbose = self.verbose)
        
        if self.status == 0:
            if self.node == 'pixel':
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem), self.verbose, self.verbose)
                os.rename('tmp.grd', self.dem)

            utils.remove_glob('*.cmd')

    ## ==============================================
    ## Run GMT surface on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## ==============================================

    def surface(self):
        '''Generate a DEM with GMT surface'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        self.o_fmt = 'GTiff'

        dem_surf_cmd = ('gmt blockmean {} -I{:.10f} -V {} | gmt surface -V {} -I{:.10f} -G{}_p.grd -T.35 -Z1.2 -Lud -Lld {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd(dem_surf_cmd, self.verbose, True, self.datalist._dump_data)

        if self.status == 0:
            dem_cut_cmd = ('gmt grdcut -V {}_p.grd -G{} {}\
            '.format(self.o_name, self.dem, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, True)

            utils.remove_glob('{}_p.grd'.format(self.o_name))

    ## ==============================================
    ## run GMT xyz2grd on the datalist to generate
    ## an uninterpolated grid
    ## see the -A switch in GMT xyz2grd for mode options
    ## 
    ## `[d|f|l|m|n|r|S|s|u|z]
    ## By  default  we  will calculate mean values if multiple entries fall on the same node. Use -A to change this behavior, except it is ignored if -Z is given. Append f or s to simply keep the first or last data point that
    ## was assigned to each node. Append l or u or d to find the lowest (minimum) or upper (maximum) value or the difference between the maximum and miminum value at each node, respectively. Append m or r or S to compute mean
    ## or  RMS  value  or standard deviation at each node, respectively. Append n to simply count the number of data points that were assigned to each node (this only requires two input columns x and y as z is not consulted).
    ## Append z to sum multiple values that belong to the same node.`
    ## ==============================================

    def num(self, mode = 'n'):
        '''Generate a num and num-msk grid with GMT'''

        self.dem = '{}_{}.grd'.format(self.o_name, mode)
        out, self.status = xyz2grd(self.datalist, self.dist_region, self.inc, self.dem, mode, self.node, verbose = self.verbose)
        
    def mask(self):
        '''Generate a num and num-msk grid'''

        #self.dem = self.datalist.mask(self.inc)
        
        self.num('n')

        if self.status == 0:
            msk = '{}_msk.grd'.format(self.o_name)
            num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
            '.format(self.dem, msk))
            out, status = utils.run_cmd(num_msk_cmd, self.verbose, self.verbose)

            if self.status == 0:
                utils.remove_glob(self.dem)
                self.dem = msk

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Inverse Distance DEM
    ## ==============================================

    def invdst(self):
        '''Generate an inverse distance grid with GDAL'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'
        
        out_size_x = int((self.dist_region.east - self.dist_region.west) / self.inc)
        out_size_y = int((self.dist_region.north - self.dist_region.south) / self.inc)
        
        out, self.status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = COMMA', self.verbose, True)
        bm_cmd = ('gmt blockmean {} -I{:.7f} -V {} > {}.csv\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.o_name))
        out, self.status = utils.run_cmd(bm_cmd, self.verbose, True, self.datalist._dump_data)
        
        o_vrt = open('{}.vrt'.format(self.o_name), 'w')
        t = '''<OGRVRTDataSource>
  <OGRVRTLayer name="{}">
    <SrcDataSource>{}.csv</SrcDataSource>
    <GeometryType>wkbPoint</GeometryType>
    <GeometryField encoding="PointFromColumns" x="field_1" y="field_2" z="field_3"/>
  </OGRVRTLayer>
</OGRVRTDataSource>'''.format(self.o_name, self.o_name)
        o_vrt.write(t)
        o_vrt.close()
        
        gg_cmd = 'gdal_grid -zfield "field_3" -txe {} {} -tye {} {} -outsize {} {} -a invdist -l {} {}.vrt {}.tif --config GDAL_NUM_THREADS ALL_CPUS\
        '.format(self.dist_region.west, self.dist_region.east, self.dist_region.north, self.dist_region.south, out_size_x, out_size_y, self.o_name, self.o_name, self.o_name)
        print gg_cmd
        out, status = utils.run_cmd(gg_cmd, self.verbose, True)

    ## ==============================================
    ## Bathy-Surface module.
    ##
    ## with GMT, gererate a 'bathy-surface'.
    ## specify the coastline mask shapefile with 'mask = '
    ## if mask is not specified, will use gsshg.
    ## returns a netcdf masked grid and an XYZ file
    ## of the bathy surface points.
    ## ==============================================

    def bathy(self, mask = None):
        '''Generate a masked surface with GMT surface and a breakline mask polygon.
        Use a coastline polygon mask to generatre a `bathy-surface` DEM.'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        self.dem = '{}_bs.grd'.format(self.o_name)

        if mask is None:
            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V -N1/0 {}\
            '.format(self.inc, self.proc_region.gmt, reg_str))
            out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, self.verbose)

        dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V {} | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lu0 {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd(dem_surf_cmd, self.verbose, True, self.datalist._dump_data)

        if self.status == 0:
            dem_landmask_cmd1 = ('gmt grdmath -V {}_p.grd tmp_lm.grd MUL 0 NAN = {}_p_bathy.grd\
            '.format(self.o_name, self.o_name))
            out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, True)

            grdcut('{}_p_bathy.grd'.format(self.o_name), self.dist_region, self.dem)
            utils.remove_glob('{}_p*.grd'.format(self.o_name))

            if self.status == 0:
                bathy_xyz = 'xyz/{}_bs.xyz'.format(self.o_name)                
                self.status = grd2xyz(self.dem, bathy_xyz, want_datalist = True)

    ## ==============================================
    ## conversion-grid module:
    ## Generate a VDatum conversion grid
    ##
    ## Requires NOAA's VDatum version >= 4.x
    ## ==============================================

    def vdatum(self, ivert = 'navd88', overt = 'mllw', region = '3'):
        '''generate a vertical datum transformation grid with vdatum'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'
        
        gdalfun.null('empty.tif', self.dist_region.region, 0.00083333, nodata = 0)
        gdalfun.dump('empty.tif', 'empty.xyz', dump_nodata = True)

        this_vd = utils.vdatum(verbose = self.verbose)
        this_vd.ivert = ivert
        this_vd.overt = overt

        this_vd.run_vdatum('empty.xyz')

        if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
            out, status = utils.run_cmd('gmt gmtinfo result/empty.xyz -C')
            empty_infos = out.split()

            if empty_infos[4] > 0:
                ll_switch = '-Lld'
            else: ll_switch = '-Ll0'

            if empty_infos[5] > 0:
                lu_switch = '-Lud'
            else: lu_switch = '-Lu0'

            gc = 'gmt blockmean result/empty.xyz -V -I{} {} {}\
            | gmt surface -I{} {} -G{}_{}to{}.tif=gd:GTiff -V -T0 {} {} {}\
            '.format(self.inc, self.dist_region.gmt, reg_str, 
                     self.inc, self.dist_region.gmt, self.o_name, this_vd.ivert, this_vd.overt, 
                     ll_switch, lu_switch, reg_str)
            utils.run_cmd(gc, self.verbose, True)

            dem_name = '{}to{}'.format(this_vd.ivert, this_vd.overt)
            self.dem[dem_name] = '{}_{}.tif'.format(self.o_name, dem_name)
        else: self.status = -1

        try:
            utils.remove_glob('empty.*')
            utils.remove_glob('result/*')
            os.removedirs('result')
        except: pass

        return(self.dem)

    def spatial_metadata(self, epsg = 4269):
        #self.datalist.i_fmt = -1
        #self.datalist._load_data()
        sm = metadata.spatial_metadata(self.datalist, self.region, i_inc = self.inc, o_name = self.o_name, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        sm.want_queue = True
        sm.run(epsg)

    def uncertainty(self, dem_mod = 'mbgrid', dem = None, msk = None, prox = None):
        unc = uncertainty.uncertainty(self.datalist, self.region, i_inc = self.inc, o_name = self.o_name, o_node = self.node, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        unc.run(dem_mod, dem, msk)
    
## =============================================================================
##
## Console Scripts
##
## Scripts are: 
## waffles
##
## =============================================================================

_dem_mods = {
    'mbgrid': [lambda x: x.mbgrid, 'Weighted SPLINE DEM via mbgrid', ':dist'],
    'surface': [lambda x: x.surface, 'SPLINE DEM via GMT surface', None],
    'invdst': [lambda x: x.invdst, 'Inverse Distance DEM via gdal_grid', None],
    'num': [lambda x: x.num, 'Uninterpolated DEM via GMT xyz2grd', ':mode'],
    'mask': [lambda x: x.mask, 'DEM Data MASK', None],
    'bathy': [lambda x: x.bathy, 'BATHY-DEM', ':coast-poly'],
    'conversion-grid': [lambda x: x.vdatum, 'VDatum CONVERSION grid', ':i_vdatum:o_vdatum:vd_region'],
    'spatial-metadata': [lambda x: x.spatial_metadata, 'DEM SPATIAL METADATA', ':epsg'],
    'uncertainty': [lambda x: x.uncertainty, 'DEM UNCERTAINTY grid <beta>', ':gridding-module'],
}

def dem_mod_desc(x):
    out = 'Modules:\n'
    for key in x:
        out += '  {:16}\t{} [{}]\n'.format(key, x[key][-2], x[key][-1])
    return(out)

_waffles_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

usage: {} [ -hrvAIEFOPR [ args ] ] module[:option1:option2] ...

{}
Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -I, --datalsit\tThe input DATALIST.
  -E, --increment\tThe desired CELL-SIZE in native units.
  -F, --format\t\tThe desired output FORMAT.
  -P, --prefix\t\tThe output naming PREFIX.
  -O, --output-name\tThe output BASENAME; will over-ride any set PREFIX.
  -X, --extend\t\tThe number of cells with which to extend the extent (6)

  -r\t\t\tuse grid-node registration, default is pixel-node

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
--verbose\t\tIncrease the verbosity

 Examples:
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 surface
 % {} --datalist input.datalist --increment 0.0000925925 -X 2 --region input_tiles_ply.shp mbgrid spatial-metadata
 % {} -R-82.5/-82.25/26.75/27 -E0.0000925 conversion-grid:navd88:mhw:3 -P ncei -r

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]), 
            dem_mod_desc(_dem_mods),
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def main():
    status = 0
    i_region = None
    i_datalist = None
    i_inc = 0.0000925925
    these_regions = []
    stop_threads = False
    want_verbose = False
    mod_opts = {}
    o_pre = None
    o_bn = None
    o_fmt = 'GTiff'
    node_reg = 'pixel'
    o_extend = 6

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

        elif arg == '--prefix' or arg == '-P':
            o_pre = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-P':
            o_pre = str(arg[2:])

        elif arg == '--extend' or arg == '-X':
            o_extend = int(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-X':
            o_extend = int(arg[2:])

        elif arg == '-r':
            node_reg = 'grid'

        elif arg == '--help' or arg == '-h':
            print(_waffles_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_waffles_usage)
            sys.exit(0)

        else: 
            opts = arg.split(':')
            if opts[0] in _dem_mods.keys():
                mod_opts[opts[0]] = list(opts[1:])
            else: utils._progress().err_msg('invalid module name `{}`'.format(opts[0]))

        i = i + 1

    if i_region is None:
        utils._error_msg('a region must be specified')
        print(_waffles_usage)
        sys.exit(1)

    if i_datalist is None:
        utils._error_msg('a datalist must be specified')
        print(_waffles_usage)
        sys.exit(1)        

    try:
        #i_inc = float(i_inc)
        i_inc = regions._inc(i_inc)
    except:
        utils._error_msg('the increment value should be a number')
        print(_waffles_usage)
        sys.exit(1)

    for key in mod_opts.keys():
        mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]

    utils.check_config()
        
    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if want_verbose: pb = utils._progress('loading region(s)...')
    try: 
        these_regions = [regions.region(i_region)]
    except: these_regions = [regions.region('/'.join(map(str, x))) for x in gdalfun._ogr_extents(i_region)]

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
          status = -1

    if want_verbose: pb.end(status, 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions)))
            
    if status == -1:
        utils._error_msg('failed to load region(s)')
        print(_waffles_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):

        status = 0
        
        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        if stop_threads: break
        
        this_datalist = datalists.datalist(i_datalist, this_region, verbose = want_verbose)
        if not this_datalist._valid_p():
            utils._error_msg('invalid datalist')
            status = -1

        if status != 0: break

        for dem_mod in mod_opts.keys():

            ## ==============================================
            ## Run the DEM module
            ## ==============================================

            args = tuple(mod_opts[dem_mod])
                        
            pb = utils._progress('running geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string))
            dl = dem(this_datalist, this_region, i_inc = i_inc, o_name = o_pre, o_b_name = o_bn, o_node = node_reg, o_fmt = o_fmt, o_extend = o_extend, callback = lambda: stop_threads, verbose = want_verbose)
            #dl.node = node_reg
            #dl.o_fmt = o_fmt

            t = threading.Thread(target = dl.run, args = (dem_mod, args))
            try:
                t.start()
                while True:
                    time.sleep(1)
                    pb.update()
                    if not t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): 
                dl.status = -1
                stop_threads = True

            t.join()
            pb.end(dl.status, 'ran geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string))
            
if __name__ == '__main__':
    main()

### END
