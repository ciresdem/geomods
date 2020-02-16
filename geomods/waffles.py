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
import ogr

import regions
import datalists
import gdalfun
import uncertainty
import metadata
import utils

import threading
import time
        
_version = '0.2.6'

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

def grd2tif(src_grd, verbose = False):
    '''Convert the grd file to tif using GMT'''
    
    status = 0
    if os.path.exists(src_grd):
        dst_grd = '{}.tif'.format(os.path.basename(src_grd).split('.')[0])

        grd2tif_cmd = ('gmt grdconvert {} {}=gd+n-9999:GTiff -V\
        '.format(src_grd, dst_grd))
        out, status = utils.run_cmd(grd2tif_cmd, verbose, verbose)

        if status != 0:
            dst_grd = None
            #else: dst_grd = None

    return(dst_grd)

def grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''

    status = 0
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
        grdmask_cmd = ('gmt grdmath -V {} {} OR = tmp.grd'.format(src_grd, mask))
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
    out, status = utils.run_cmd(num_msk_cmd, verbose, verbose)

    return(status)

def xyz2grd(datalist, region, inc, dst_name, a = 'n', node = 'pixel', verbose = False):

    status = 0
    if node == 'pixel':
        reg_str = '-r'
    else: reg_str = ''
    
    num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -G{}.grd -A{} {}\
    '.format(region.gmt, inc, dst_name, a, reg_str))
    out, status = utils.run_cmd_with_input(num_cmd0, datalist._cat_port, verbose, True)

    return(out, status)

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):

    status = 0
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1:
        dist = dist + '/2'
    
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.7f}/{:.7f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {}> /dev/null 2> /dev/null \
    '.format(datalist._path, region.gmt, inc, inc, dst_name, dist, tension, e_switch))
    out, status = utils.run_cmd(mbgrid_cmd, verbose, True)

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

    def __init__(self, i_datalist, i_region, i_inc = 0.000277777, o_name = None, o_b_name = None, callback = lambda: False, verbose = False):
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.dem_mods = {
            'mbgrid': self.mbgrid,
            'surface': self.surface,
            'num': self.num,
            'mean': self.mean,
            'bathy': self.bathy,
            'conversion-grid': self.vdatum,
            'spatial-metadata': self.spatial_metadata,
            'uncertainty': self.uncertainty,
        }

        self.dem = { 
            'dem': None,
            'dem-grd': None,
            'num': None,
            'num-grd': None,
            'num-msk': None,
            'num-msk-grd': None,
            'mean': None,
            'mean-grd': None,
            'bathy-grd': None,
            'bathy': None
        }

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]

            str_inc = inc2str_inc(self.inc)
            self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        else: self.o_name = o_b_name

    def run(self, dem_mod = 'mbgrid', args = ()):
        '''Run the DEM module `dem_mod` using args `dem_mod_args`.'''

        dems = self.dem
        if dem_mod != 'mbgrid' and dem_mod != 'vdatum' and dem_mod != 'metadata_spatial':
            self.datalist._load_data()
            if len(self.datalist.datafiles) == 0:
                self.status = -1

        if self.status == 0:
            dems = self.dem_mods[dem_mod](*args)

        return(dems)

    ## ==============================================
    ## run mbgrid on the datalist and generate the DEM
    ## note: mbgrid will cause popen to hang if stdout 
    ## is not cleared...should add output file to send to...
    ## ==============================================

    def mbgrid(self, dist = '10/3'):
        '''Generate a DEM with MBSystem's mbgrid program.'''

        out, self.status = run_mbgrid(self.datalist, self.dist_region, self.inc, self.o_name, dist = dist, verbose = self.verbose)
        
        if self.status == 0:
            self.dem['dem-grd'] = '{}.grd'.format(self.o_name)

            if self.node == 'pixel':
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['dem-grd']), self.verbose, self.verbose)
                os.rename('tmp.grd', self.dem['dem-grd'])

            self.dem['dem'] = grd2tif(self.dem['dem-grd'])
            utils.remove_glob('*.cmd')

        return(self.dem)    

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

        dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V {} | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lud -Lld {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:
            dem_cut_cmd = ('gmt grdcut -V {}_p.grd -G{}.grd {}\
            '.format(self.o_name, self.o_name, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)

            if os.path.exists('{}_p.grd'.format(self.o_name)):
                os.remove('{}_p.grd'.format(self.o_name))            

            if self.status == 0:
                self.dem['dem-grd'] = ('{}.grd'.format(self.o_name))
                self.dem['dem'] = grd2tif(self.dem['dem-grd'])

        return(self.dem)

    ## ==============================================
    ## run GMT xyz2grd on the datalist to generate
    ## a NUM grid
    ## ==============================================

    def num(self):
        '''Generate a num and num-msk grid with GMT'''

        dst_name = '{}_num'.format(self.o_name)
        out, self.status = xyz2grd(self.datalist, self.dist_region, self.inc, dst_name, 'n', self.node, verbose = self.verbose)
        
        if self.status == 0:
            self.dem['num-grd'] = '{}.grd'.format(dst_name)
            self.dem['num'] = grd2tif(self.dem['num-grd'])
            self.dem['num-msk-grd'] = '{}_msk.grd'.format(dst_name)

            self.status = num_msk(self.dem['num-grd'], self.dem['num-msk-grd'])
            if self.status == 0:
                self.dem['num-msk'] = grd2tif(self.dem['num-msk-grd'])

        return(self.dem)

    ## ==============================================
    ## run GMT xyz2grd on the datalist to generate
    ## a MEAN grid (no interpolation)
    ## ==============================================

    def mean(self):
        '''Generate a num and num-msk grid with GMT'''
        dst_name = '{}_mean'.format(self.o_name)
        out, self.status = xyz2grd(self.datalist, self.dist_region, self.inc, dst_name, 'm', self.node, verbose = self.verbose)
        
        if self.status == 0:
            self.dem['mean-grd'] = '{}.grd'.format(dst_name)
            self.dem['mean'] = grd2tif(self.dem['mean-grd'])

        return(self.dem)

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Inverse Distance DEM
    ## ==============================================

    def inv_dist(self):
        '''Generate an inverse distance grid with GDAL'''

        ## datalist ==> point shapefile (block on the way?)
        ## gdal_grid point shpafile to out grid
        pass

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

        bathy_dem = '{}_bs.grd'.format(self.o_name)

        if mask is None:
            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V -N1/0 {}\
            '.format(self.inc, self.proc_region.gmt, reg_str))
            out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, self.verbose)

        dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V {} | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lu0 {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:
            dem_landmask_cmd1 = ('gmt grdmath -V {}_p.grd tmp_lm.grd MUL 0 NAN = {}_p_bathy.grd\
            '.format(self.o_name, self.o_name))
            out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, True)

            grdcut('{}_p_bathy.grd'.format(self.o_name), self.dist_region, bathy_dem)
            utils.remove_glob('{}_p*.grd'.format(self.o_name))

            if self.status == 0:
                self.dem['bathy-grd'] = bathy_dem
                bathy_xyz = 'xyz/{}_bs.xyz'.format(self.o_name)                
                self.status = grd2xyz(self.dem['bathy-grd'], bathy_xyz, want_datalist = True)

                if self.status == 0:
                    self.dem['bathy-xyz'] = bathy_xyz

        return(self.dem)

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
        sm = metadata.spatial_metadata(self.datalist, self.region, self.inc, self.o_name, self.stop, self.verbose)
        sm.spatial_metadata(epsg)

    def uncertainty(self, dem_mod = 'mbgrid', dem = None, num = None, num_msk = None, prox = None):
        unc = uncertainty.uncertainty(self.datalist, self.region, self.inc, self.o_name, self.stop, self.verbose)
        unc.run(dem_mod, dem, num, num_msk, prox)
    
## =============================================================================
##
## Console Scripts
##
## Scripts are: 
## waffles
##
## =============================================================================

_dem_lambda = lambda d, r, i, o, b, c, v: dem(d, r, i, o, b, c, v)

_dem_mods = {
    'conversion-grid': 'generate a VDatum conversion grid [:i_vdatum:o_vdatum:vd_region]',
    'bathy': 'generate a `bathy-surface` masked grid [:coastpoly]',
    'surface': 'generate a DEM via GMT surface',
    'mbgrid': 'generate a DEM via mbgrid [:dist]',
    'spatial-metadata': 'generate DEM spatial metadata [:epsg]',
    'uncertainty': 'calculate DEM uncertainty [:dem_module]',
}

def dem_mod_desc(x):
    dmd = []
    for key in x: 
        dmd.append('{:16}\t{}'.format(key, x[key]))
    return dmd

_waffles_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

usage: {} [ -hvAIER [ args ] ] module[:option1:option2] ...

Modules:
  {}

Options:
  -R, --region\t\tSpecifies the desired region;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -I, --datalsit\tThe input datalist.
  -E, --increment\tThe desired cell-size in native units.
  -P, --prefix\t\tThe output naming prefix.
  -O, --output-name\tThe output basename.

  -r\t\t\tuse grid-node registration, default is pixel-node

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 surface
 % {} --datalist input.datalist --increment 0.000277777 --region input_tiles_ply.shp mbgrid spatial-metadata
 % {} -R-82.5/-82.25/26.75/27 -E0.0000925 conversion-grid:navd88:mhw:3 -P ncei -r

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]), 
            '\n  '.join(dem_mod_desc(_dem_mods)),
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def main():
    status = 0
    i_region = None
    i_datalist = None
    i_inc = 0.0002777
    these_regions = []
    stop_threads = False
    want_verbose = False
    mod_opts = {}
    o_pre = None
    o_bn = None
    node_reg = 'pixel'

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
            try:
                i_inc = float(argv[i + 1])
            except:
                sys.stderr.write('geomods: error, -E should be a float value\n')
                sys.exit(1)
            i = i + 1
        elif arg[:2] == '-E':
            try:
                i_inc = float(arg[2:])
            except:
                sys.stderr.write('geomods: error, -E should be a float value\n')
                sys.exit(1)

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
            mod_opts[opts[0]] = list(opts[1:])

        i = i + 1

    if i_region is None:
        sys.stderr.write('geomods: error, a region must be specified\n')
        print(_waffles_usage)
        sys.exit(1)

    if i_datalist is None:
        sys.stderr.write('geomods: error, a datalist must be specified\n')
        print(_waffles_usage)
        sys.exit(1)        

    for key in mod_opts.keys():
        mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]

    ## ==============================================
    ## check platform and installed software
    ## ==============================================

    if 'mbgrid' in mod_opts.keys():
        mb_vers = utils._cmd_check('mbgrid', 'mbgrid -version | grep Version')

    gmt_vers = utils._cmd_check('gmt', 'gmt --version')
    gdal_vers = utils._cmd_check('gdal-config', 'gdal-config --version')

    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    pb = utils._progress('loading region(s)...')
    try: 
        these_regions = [regions.region(i_region)]
    except:
        if os.path.exists(i_region):
            _poly = ogr.Open(i_region)
            if _poly is not None:
                _player = _poly.GetLayer(0)
                for pf in _player:
                    _pgeom = pf.GetGeometryRef()
                    these_regions.append(regions.region('/'.join(map(str, _pgeom.GetEnvelope()))))

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
          status = -1
    pb.opm = 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions))
    pb.end(status)
            
    if status == -1:
        sys.stderr.write('geomods: error, failed to load region(s).\n')
        print(_waffles_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):

        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        if stop_threads:
            break
        
        if i_datalist is not None:
            this_datalist = datalists.datalist(i_datalist, this_region)
            if not this_datalist._valid:
                sys.stderr.write('geomods: error, invalid datalist\n')
                status = -1

        else: this_datalist = None

        if status != 0: 
            status = 0
            break

        for dem_mod in mod_opts.keys():

            ## ==============================================
            ## Run the DEM module
            ## ==============================================

            args = tuple(mod_opts[dem_mod])
                        
            pb = utils._progress('running geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod, rn + 1, len(these_regions), this_region.region_string))

            dl = dem(this_datalist, this_region, i_inc, o_pre, o_bn, lambda: stop_threads, want_verbose)
            dl.node = node_reg

            #dl.run(dem_mod, args)
            t = threading.Thread(target = dl.run, args = (dem_mod, args))

            try:
                t.start()
                while True:
                    time.sleep(2)
                    pb.update()
                    if not t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): 
                dl.status = -1
                stop_threads = True

            t.join()
            pb.opm = 'ran geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod, rn + 1, len(these_regions), this_region.region_string)
            pb.end(dl.status)
            
if __name__ == '__main__':
    main()

### END
