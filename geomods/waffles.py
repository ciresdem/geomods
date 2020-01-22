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
## Process DEMs
##
### Code:

import sys
import os
import re
import math
import glob

import random
import numpy as np
import fractions
import threading
import datetime
import time

try:
    import osgeo.ogr as ogr
    import osgeo.osr as osr
    import osgeo.gdal as gdal
    has_gdalpy = True
except ImportError:
    try:
        import ogr
        import osr
        import gdal
        has_gdalpy = True
    except ImportError:
        has_gdalpy = False

try:
    import arcpy
    has_arcpy = True
except ImportError:
    has_arcpy = False

import regions
import datalists
import gdalfun
import utils

_version = '0.1.4'

## =============================================================================
##
## DEM Modules
## { mod-name: [module-function, module-description, module-function arguments]}
##
## =============================================================================

_dem_mods = { 'mbgrid': [lambda x: x.mbgrid, 'generate a DEM with mbgrid', 'None'],
              'gmt-surface': [lambda x: x.surface, 'generate a DEM with GMT', 'None'],
              'spatial-metadata': [lambda x: x.spatial_metadata, 'generate spatial-metadata', 'epsg'],
              'conversion-grid': [lambda x: x.conversion_grid, 'generate a conversion grid with vdatum', 'ivert:overt:region'],
              'bathy-surface': [lambda x: x.masked_surface, 'generate a bathymetry surface with GMT', 'coastline']
}

def dem_mod_desc(x):
    dmd = []
    for key in x: 
        dmd.append('{:16}\t{} [{}]'.format(key, x[key][1], x[key][2]))
    return dmd

_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

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
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 gmt-surface
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

## =============================================================================
##
## Uncertainty Analysis
##
## =============================================================================

def err2coeff(data):
    '''data is 2 col file with `err dist`'''
    try: 
        my_data = np.loadtxt(data, delimiter=' ')
    except: sys.exit(2)

    error=my_data[:,0]
    distance=my_data[:,1]
        
    max_int_dist = np.max(distance)
    nbins = 10

    coeff_guess=[0, 0.1, 0.2]
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

    fitfunc = lambda p, x: p[0] + p[1] * (abs(x) ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)
    return(out)

def dem_interpolation_uncertainty(dem_surface):
    status = 0
    sub_count = 0

    #tw = clis.prog_bar('preparing for interpolation uncertainty')._end(status)

    this_region = dem_surface.region
    status = dem_surface.proximity()
    
    for sub_region in this_region.split(4):
        sub_count += 1

        ## ==============================================            
        ## Extract XYZ data for sub-region and 
        ## randomly sample
        ## ==============================================

        status = this_surf.grd2xyz(this_surf.dem['dem-grd'], 
                                   '{}.xyz'.format(this_surf.oname), 
                                   region = sub_region.buffer(10 * float(iinc)), 
                                   mask = this_surf.dem['num-grd'])
        
        sub_xyz = np.loadtxt(this_surf.dem['xyz'], delimiter = ' ')
                
        np.random.shuffle(sub_xyz)
        sx_len = len(sub_xyz)
        sx_len_pct = int(sx_len * .5)

        np.savetxt('sub_{}d_rest.xyz'.format(sub_count), sub_xyz[sx_len_pct:], '%f', ' ')
                
        sub_datalist =  datalists.datalist('sub_{}.datalist'.format(sub_count), sub_region)
        sub_datalist._append_datafile('sub_{}_rest.xyz'.format(sub_count), 168, 1)
        sub_datalist._reset()

        ## ==============================================
        ## Generate the random-sample DEM            
        ## ==============================================

        sub_surf = cudem(sub_datalist, sub_region, iinc)
        sub_surf.mbgrid()
            
        status = sub_surf.proximity()

        ## ==============================================            
        ## Query the random-sample DEM with the 
        ## withheld data
        ## ==============================================

        sub_xyd = gdalfun.query(sub_xyz[:sx_len_pct], sub_surf.dem['tif'], 'xyd')

        ## ==============================================
        ## Query the random-sample proximity grid 
        ## with the withheld data
        ## ==============================================

        sub_dp = gdalfun.query(sub_xyd, sub_surf.dem['prox'], 'zg')
        #print np.max(sub_dp[:,1])
            
        ## ==============================================
        ## Cleanup
        ## ==============================================

        fl = glob.glob('sub_{}*'.format(sub_count))
        for f in fl:
            try: 
                os.remove(f)
            except: pass

## =============================================================================
##
## DEM Class for DEM development et cetra
##
## Generate DEMs with GMT and MBSystem
## idatalist and iregion are datalist and region instances, respectively
## the `dem` dictionary holds the file names of the generated DEM files
## all files are GTiff unless specified as 'grd'
## including 'dem', 'dem-grd', 'num', 'num-grd', 'num-msk', 'prox', 'slope'
##
## =============================================================================

class cudem(threading.Thread):
    def __init__(self, idatalist, iregion, iinc = '0.000277777', oname = None, obname = None, callback = lambda: False, verbose = False):
        '''run a number of modules for DEM generation and analysis'''

        threading.Thread.__init__(self)

        self.status = 0
        self.stop = callback
        self.verbose = verbose
        self.dem = {}
        self.inc = float(iinc)
        self.datalist = idatalist
        self._module = self._valid_p
        self._module_args = ()
        self.node = 'pixel'

        self.region = iregion
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)

        self.max_prox = self.max_num = None

        if obname is None:
            if oname is None: 
                if idatalist is None:
                    oname = 'cgrid'
                else: oname = self.datalist._path_basename.split('.')[0]

            str_inc = str(fractions.Fraction(str(self.inc * 3600)).limit_denominator(10)).replace('/', '')
            self.oname = '{}{}_{}_{}'.format(oname, str_inc, self.region.fn, datetime.datetime.now().strftime('%Y'))
        else: self.oname = obname

    def run(self):
        if len(self._module_args) == 0:
            self._module()
        else: self._module(*self._module_args)

    def _valid_p(self):
        '''make sure all files referenced in dec dict exist'''

        for key in self.dem:
            if not os.path.exists(self.dem[key]):
                return(False)
        return(True)

    def _rectify(self):
        '''reove keys from the dem dictionary if their references
        don't exist'''
 
        for key in self.dem:
            if not os.path.exists(self.dem[key]):
                del self.dem[key]
                
    def _set_dem(self, key, value):
        '''Set a key/value in the dem dictionary'''

        self.dem[key] = value

    def _print_dem(self):
        '''Print out the key/value pairs from the dem dictionary'''

        for key in self.dem:
            print(key, self.dem[key])

    def _dem_remove(self, key):
        '''Remove file associated with key in the dem dictionary'''

        try: 
            os.remove(self.dem[key])
            return(0)
        except: return(-1)

    ## ==============================================
    ## Generic Functions
    ## ==============================================

    def grd2tif(self, src_grd):
        '''Convert the grd file to tif using GMT'''

        if os.path.exists(src_grd):
            grd2tif_cmd = ('gmt grdconvert {} {}.tif=gd+n-9999:GTiff -V\
            '.format(src_grd, os.path.basename(src_grd).split('.')[0]))
            
            out, self.status = utils.run_cmd(grd2tif_cmd, self.verbose)
        else: self.status = -1

        return(self.status)

    def grdinfo(self, src_grd):
        '''Return an info list of `src_grd`'''

        if not os.path.exists(src_grd):
            self.status = -1
            return(self.status)
        else:
            grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
            out, self.status = utils.run_cmd(grdinfo_cmd, self.verbose)
            return(out.split())

    def grdcut(self, src_grd, src_region, dst_grd):
        '''Cut `src_grd` to `src_region` '''

        if os.path.exists(src_grd):
            self.cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))

    def grd2xyz(self, src_grd, dst_xyz, region = None, mask = None):
        '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region'''

        if mask:
            self.grdmask_cmd = ('gmt grdmath -V {} NOT {} = tmp.grd'.format(mask, src_grd))

            out, self.status = utils.run_cmd(self.grdmask_cmd, self.verbose)
            if self.status == 0: src_grd = 'tmp.grd'

        if region and region._valid:
            self.grd2xyz_cmd = ('gmt grd2xyz -V {} {} -s > {}'.format(src_grd, region.gmt, dst_xyz))

        else: self.grd2xyz_cmd = ('gmt grd2xyz {} -V -s > {}'.format(src_grd, dst_xyz))
            
        out, self.status = utils.run_cmd(self.grd2xyz_cmd, self.verbose)

        if self.status == 0:
            if mask:
                if os.path.exists('tmp.grd'):
                    os.remove('tmp.grd')

        return(self.status)

    def slope(self):
        '''Generate a Slope grid from the DEM with GMT'''

        try: self.dem['dem-grd']
        except KeyError: 
            self.status = -1

        if self.status == 0:
            slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
            '.format(self.dem['dem-grd'], self.oname, self.dem['dem-grd']))

            slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}_slp.tif=gd+n-9999:GTiff\
            '.format(self.oname, self.oname))

            out, self.status = utils.run_cmd(slope_cmd0, self.verbose)

            if self.status == 0:
                out, self.status = utils.run_cmd(slope_cmd1, self.verbose)

                if os.path.exists('{}_pslp.grd'.format(self.oname)):
                    os.remove('{}_pslp.grd'.format(self.oname))

                self.dem['slope'] = ('{}_slp.tif'.format(self.oname))

        return(self.status)

    def proximity(self):
        '''Generate a Proximity grid with GDAL'''

        try: self.dem['num-msk']
        except KeyError: 
            self.status = -1

        if self.status == 0:
            self.dem['prox'] = ('{}_prox.tif'.format(self.oname))

            gdalfun.proximity(self.dem['num-msk'], self.dem['prox'])
        
            self.max_prox = self.grdinfo(self.dem['prox'])[6]
            print(self.max_prox)

        return(self.status)
        
    ## ==============================================
    ## Main DEM Modules
    ## ==============================================

    def num(self):
        '''Generate a num and num-msk grid with GMT'''

        #num_cmd0 = ('cat {} | gmt xyz2grd -V {} -I{} -r -G{}_num.grd -An\
        #'.format(self.datalist._echo_datafiles(' '), self.dist_region.gmt, self.inc, self.oname))

        if self.node == 'pixel':
            num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -r -G{}_num.grd -An\
            '.format(self.dist_region.gmt, self.inc, self.oname))
        else:
            num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -G{}_num.grd -An\
            '.format(self.dist_region.gmt, self.inc, self.oname))

        num_cmd1 = ('gmt grdmath -V {}_num.grd 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
        '.format(self.oname, self.oname))

        pb = 'processing {} to NUM grid'.format(self.datalist._path_basename)

        out, self.status = utils.run_cmd_with_input(num_cmd0, self.datalist._cat_port, verbose = self.verbose, prog = pb)

        if self.status == 0:
            if self.verbose:
                pb = 'generating NUM-MSK grid using gmt grdmath 0 MUL 1 ADD 0 AND'
            else: pb = None
            out, self.status = utils.run_cmd(num_cmd1, self.verbose, pb)

            if self.status == 0:
                self.grd2tif('{}_num.grd'.format(self.oname))

                self.dem['num-grd'] = '{}_num.grd'.format(self.oname)
                self.grd2tif(self.dem['num-grd'])
                self.dem['num'] = '{}_num.tif'.format(self.oname)
                self.dem['num-msk'] = '{}_num_msk.tif'.format(self.oname)

                try:
                    self.max_num = int(self.grdinfo(self.dem['num-grd'])[6])
                except: self.max_num = float(self.grdinfo(self.dem['num-grd'])[6])

        return(self.status)

    def num_msk(self):
        gdalfun.xyz_gmask(self.datalist._caty(), 
                                  '{}_num_msk.tif'.format(self.oname), 
                                  self.dist_region.region, 
                                  self.inc, verbose = self.verbose)

        self.dem['num-msk'] = '{}_num_msk.tif'.format(self.oname)

    def masked_surface(self, mask = None):
        '''Generate a masked surface with GMT surface and a breakline mask polygon'''

        ## GMT GRDLANDMASK FOR COASTLINE WHEN COASTLINE IS NONE
        ## OR GRID ALL DATA AND EXTRACT ZERO LINE AS BEST WE CAN
        ## TRY SURFACE -D OPTION FOR BREAKLINE DATA

        if self.node == 'pixel':
            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V -r -N1/0\
            '.format(self.inc, self.proc_region.gmt))

            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -r -V | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -r -Lu0\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))
        else:
            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V -N1/0\
            '.format(self.inc, self.proc_region.gmt))

            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -r -V | gmt surface -V {} -I{} -G{}_p.grd -T.35 -Z1.2 -Lu0\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))

        dem_landmask_cmd1 = ('gmt grdmath -V {}_p.grd tmp_lm.grd MUL 0 NAN = {}_p_bathy.grd\
        '.format(self.oname, self.oname))

        dem_cut_cmd = ('gmt grdcut -V {}_p_bathy.grd -G{}_bs.grd {}\
        '.format(self.oname, self.oname, self.dist_region.gmt))

        if mask is None:
            out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, 'generating landmask from gsshg')
        
        pb = 'generating bathymetry surface using GMT'

        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, pb)

        if self.status == 0:

            #if mask is None:
            out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, 'masking bathy surface')

            pb = 'clipping DEM to final region'
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, pb)
            
            if self.status == 0:

                if os.path.exists('{}_p.grd'.format(self.oname)):
                    os.remove('{}_p.grd'.format(self.oname))

                if os.path.exists('{}_p_bathy.grd'.format(self.oname)):
                    os.remove('{}_p_bathy.grd'.format(self.oname))

                ## DUMP TO XYZ AND ADD TO DATALIST

                xyz_dir = os.path.join(os.getcwd(), 'xyz')
        
                if not os.path.exists(xyz_dir):
                    os.makedirs(xyz_dir)

                self.dem['dem-bathy'] = ('{}_bs.grd'.format(self.oname))
                self.dem['xyz-bathy'] = ('xyz/{}_bs.xyz'.format(self.oname))

                self.grd2xyz(self.dem['dem-bathy'], self.dem['xyz-bathy'])

                ## ==============================================
                ## Add xyz file to datalist
                ## ==============================================

                sdatalist = datalists.datalist(os.path.join(xyz_dir, 'bathy.datalist'))
                sdatalist._append_datafile('{}'.format(os.path.basename(self.dem['xyz-bathy'])), 168, 1)
                sdatalist._reset()

                ## ==============================================
                ## Generate .inf file
                ## ==============================================

                out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, 'bathy.datalist')), False, None)

        return(self.status)

    def surface(self):
        '''Generate a DEM with GMT surface'''

        #dem_cmd = ('cat {} | gmt blockmean {} -I{} -r -V | gmt surface -V {} -I{} -G{}_p.grd -T.35 -Z1.2 -r -Lud -Lld\
        #'.format(self.datalist._echo_datafiles(' '), self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))

        if self.node == 'pixel':
            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -r -V | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -r -Lud -Lld\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))
        else:
            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lud -Lld\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))

        dem_cut_cmd = ('gmt grdcut -V {}_p.grd -G{}.grd {}\
        '.format(self.oname, self.oname, self.dist_region.gmt))

        pb = 'generating DEM using GMT'

        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, pb)

        if self.status == 0:
            pb = 'clipping DEM to final region'
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, pb)
            
            if self.status == 0:

                if os.path.exists('{}_p.grd'.format(self.oname)):
                    os.remove('{}_p.grd'.format(self.oname))

                self.dem['dem-grd'] = ('{}.grd'.format(self.oname))

                self.grd2tif(self.dem['dem-grd'])
                
                self.dem['dem'] = ('{}.tif'.format(self.oname))

        return(self.status)

    def mbgrid(self):
        '''Generate a DEM and num grid with MBSystem'''

        ## mbgrid will cause popen to hang if stdout is not cleared...should add output file to send to...

        pb = 'generating DEM and NUM grid using mbgrid -T35 -X0.1 -C10/3 -M'
        mbgrid_cmd = ('mbgrid -I{} {} -E{:.7f}/{:.7f}/degrees! -O{} -A2 -G100 -F1 -N -C10/3 -S0 -X0.1 -T35 -M > /dev/null 2> /dev/null \
        '.format(self.datalist._path, self.dist_region.gmt, self.inc, self.inc, self.oname))
        out, self.status = utils.run_cmd(mbgrid_cmd, self.verbose, pb)

        if self.status == 0:
        
            self.dem['dem-grd'] = '{}.grd'.format(self.oname)
            self.dem['num-grd'] = '{}_num.grd'.format(self.oname)

            if self.node == 'pixel':
                pb = 'resampling DEM to pixel-node registration'
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['dem-grd']), self.verbose, pb)
                os.rename('tmp.grd', self.dem['dem-grd'])

                pb = 'resampling NUM grid to pixel-node registration'
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['num-grd']), self.verbose, pb)
                os.rename('tmp.grd', self.dem['num-grd'])                

            self.grd2tif(self.dem['dem-grd'])
            self.dem['dem'] = '{}.tif'.format(self.oname)

            self.grd2tif(self.dem['num-grd'])
            self.dem['num'] = '{}_num.tif'.format(self.oname)

            try:
                os.remove('{}.cmd'.format(self.dem['dem-grd']))
                os.remove('{}.cmd'.format(self.dem['num-grd']))
                os.remove('{}.mb-1'.format(self.oname))
                os.remove('{}_sd.grd'.format(self.oname))
                os.remove('{}_sd.grd.cmd'.format(self.oname))
            except: pass

            self.max_num = self.grdinfo(self.dem['num-grd'])[6]

            num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
            '.format(self.dem['num-grd'], self.oname))

            pb = 'generating NUM-MSK grid using gmt grdmath 0 MUL 1 ADD 0 AND'
            out, self.status = utils.run_cmd(num_msk_cmd, self.verbose, pb)

            if self.status == 0:
                self.dem['num-msk'] = ('{}_num_msk.tif'.format(self.oname))

        return(self.status)

    def conversion_grid(self, ivert = 'navd88', overt = 'mllw', region = '3'):
        '''generate a vertical datum transformation grid with vdatum'''
        
        ## ==============================================
        ## Create empty grid and transform to xy0
        ## ==============================================

        gdalfun.null('empty.tif', self.dist_region.region, 0.00083333, nodata = 0)
        gdalfun.dump('empty.tif', 'empty.xyz', dump_nodata = True)

        ## ==============================================
        ## pass empty.xyz through vdatum
        ## ==============================================

        this_vd = utils.vdatum()
        this_vd.ivert = ivert
        this_vd.overt = overt

        this_vd.run_vdatum('empty.xyz')

        ## ==============================================
        ## surface the results
        ## ==============================================

        if os.stat('result/empty.xyz').st_size != 0:
            out, status = utils.run_cmd('gmt gmtinfo result/empty.xyz -C')
            empty_infos = out.split()

            if empty_infos[4] > 0:
                ll_switch = '-Lld'
            else: ll_switch = '-Ll0'

            if empty_infos[5] > 0:
                lu_switch = '-Lud'
            else: lu_switch = '-Lu0'

            if self.node == 'pixel':
                gc = 'gmt blockmean result/empty.xyz -V -I{} {} | gmt surface -I{} {} -G{}_{}to{}.tif=gd:GTiff -V -r -T0 {} {}\
                '.format(self.inc, self.dist_region.gmt, self.inc, self.dist_region.gmt, self.oname, this_vd.ivert, this_vd.overt, ll_switch, lu_switch)
            else:
                gc = 'gmt blockmean result/empty.xyz -V -I{} {} | gmt surface -I{} {} -G{}_{}to{}.tif=gd:GTiff -V -T0 {} {}\
                '.format(self.inc, self.dist_region.gmt, self.inc, self.dist_region.gmt, self.oname, this_vd.ivert, this_vd.overt, ll_switch, lu_switch)
            utils.run_cmd(gc, self.verbose, 'generating conversion grid')

            self.dem['{}to{}'.format(this_vd.ivert, this_vd.overt)] = '{}_{}to{}.tif'.format(self.oname, this_vd.ivert, this_vd.overt)
        else: self.status = -1

        ## ==============================================
        ## Cleanup
        ## ==============================================

        os.remove('empty.tif')
        os.remove('empty.xyz')
        
        vd_results = glob.glob('result/*')
        for vd_r in vd_results:
            os.remove(vd_r)
        os.removedirs('result')

    def spatial_metadata(self, epsg = 4269):
        '''Geneate spatial metadata from a datalist'''

        ## ==============================================
        ## these fields should be found in the datalist 
        ## starting at position 3 (from 0)
        ## ==============================================

        v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
        t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
    
        dst_vec = '{}_sm.shp'.format(self.oname, self.region.fn)
        dst_layername = os.path.basename(dst_vec).split('.')[0]

        shps = glob.glob('{}_sm.*'.format(self.oname))
        if len(shps) > 0:
            for sh in shps:
                os.remove(sh)

        ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
        if ds is None:
            self.status = -1

        if self.status == 0:

            try:
                int(epsg)
            except: epsg = 4326
            gdalfun._prj_file('{}.prj'.format(os.path.basename(dst_vec).split('.')[0]), int(epsg))

            layer = ds.CreateLayer('{}'.format(os.path.basename(dst_vec).split('.')[0]), None, ogr.wkbMultiPolygon)
            defn = layer.GetLayerDefn()

            for i, f in enumerate(v_fields):
                layer.CreateField(ogr.FieldDefn('{}'.format(f), t_fields[i]))
            
            for feature in layer:
                layer.SetFeature(feature)
    
            ## ==============================================
            ## Generate geometry for each datalist 
            ## and add to output layer
            ## ==============================================
            
            for dl in self.datalist.datalist:
                if not self.stop():
                    
                    ## ==============================================
                    ## Load the sub DATALIST
                    ## ==============================================
                    
                    pb = utils._progress('loading datalist...')

                    this_datalist = datalists.datalist(dl[0], self.dist_region)
                    this_dem = cudem(this_datalist, self.region, str(self.inc), verbose = self.verbose)

                    pb.opm = 'loading datalist...{}'.format(this_datalist._path_basename)
                    pb.end(0)

                    try:
                        o_v_fields = [dl[3], dl[4], dl[5], dl[6], dl[7], dl[8], dl[9], dl[10].strip()]
                    except: o_v_fields = [this_datalist._path_dl_name, 'Unknown', 0, 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
                    
                    ## ==============================================
                    ## Gererate the NUM-MSK
                    ## ==============================================
                    
                    #pb1 = utils._progress('generating {} mask'.format(this_datalist._path_basename))
                    #this_dem.num_msk()
                    #pb1.end(self.status)

                    this_dem.num()
                    this_dem._dem_remove('num')
                    this_dem._dem_remove('num-grd')

                    pb = utils._progress('gathering geometries from {}'.format(this_datalist._path_basename))

                    if self.status == 0:
                        src_ds = gdal.Open(this_dem.dem['num-msk'])
                        srcband = src_ds.GetRasterBand(1)

                        shps = glob.glob('{}_poly.*'.format(this_dem.oname))
                        if len(shps) > 0:
                            for sh in shps:
                                os.remove(sh)

                        ## ==============================================
                        ## Process the tmp vector to get the 
                        ## geometries of data cells
                        ## ==============================================
                        
                        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_dem.oname))
                        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_dem.oname), None, ogr.wkbMultiPolygon)
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                        
                        pb1 = utils._progress('polygonizing {} mask'.format(this_datalist._path_basename))
                        result = gdal.Polygonize(srcband, None, tmp_layer, 0, [], callback = None)

                        multi = ogr.Geometry(ogr.wkbMultiPolygon)
                        src_ds = srcband = None

                        if len(tmp_layer) > 1:
                            for i in tmp_layer:
                                if not self.stop():
                                    if i.GetField('DN') == 0:
                                        tmp_layer.DeleteFeature(i.GetFID())
                                    elif i.GetField('DN') == 1:
                                        i.geometry().CloseRings()
                                        wkt = i.geometry().ExportToWkt()
                                        multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                                        tmp_layer.DeleteFeature(i.GetFID())

                            ## ==============================================
                            ## Add geometries to output layer
                            ## ==============================================

                            union = multi.UnionCascaded()
                            out_feat = ogr.Feature(defn)
                            out_feat.SetGeometry(union)

                            for i, f in enumerate(v_fields):
                                out_feat.SetField(f, o_v_fields[i])

                            layer.CreateFeature(out_feat)

                        ## ==============================================
                        ## Cleanup the tmp layer
                        ## ==============================================
                            
                        tmp_ds = tmp_layer = None

                        shps = glob.glob('{}_poly.*'.format(this_dem.oname))
                        if len(shps) > 0:
                            for sh in shps:
                                os.remove(sh)
                        pb1.end(result)

                    this_dem._dem_remove('num-msk')
                    pb.end(self.status)
        ds = layer = None

## =============================================================================
##
## Run cudem.py from command-line
##
## =============================================================================

def main():

    iregion = None
    idatalist = None
    igrid = None
    iinc = '0.000277777'
    status = 0
    these_regions = []
    do_unc = False
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
            iregion = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            iregion = str(arg[2:])

        elif arg == '--datalist' or arg == '-I':
            idatalist = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-I':
            idatalist = str(arg[2:])

        elif arg == '--increment' or arg == '-E':
            iinc = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-E':
            iinc = str(arg[2:])

        elif arg == '--grid' or arg == '-G':
            igrid = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-G':
            igrid = str(arg[2:])

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
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        else: 
            opts = arg.split(':')
            mod_opts[opts[0]] = list(opts[1:])

        i = i + 1

    if iregion is None:
        print(_usage)
        sys.exit(1)

    ## ==============================================
    ## check platform and installed software
    ## ==============================================

    pb = utils._progress('checking platform.')
    platform = sys.platform
    pb.opm = 'checking platform...{}'.format(platform)
    pb.end(0)

    pb = utils._progress('checking for GMT...')
    if utils.cmd_exists('gmt'): 
        gmt_vers, status = utils.run_cmd('gmt --version')
    else: status = -1
    pb.opm = 'checking for GMT...{}'.format(gmt_vers.rstrip())
    pb.end(status)

    pb = utils._progress('checking for MBSystem...')
    if utils.cmd_exists('mbgrid'): 
        mbs_vers, status = utils.run_cmd('mbgrid -version')
        mbs_vers = mbs_vers.split('\n')[3].rstrip().split()[2]
    else: status = -1
    pb.opm = 'checking for MBSystem...{}'.format(mbs_vers)
    pb.end(status)

    pb = utils._progress('checking for GDAL command-line')
    if utils.cmd_exists('gdal-config'): 
        gdal_vers, status = utils.run_cmd('gdal-config --version')
        pb.opm = 'checking for GDAL command-line - {}'.format(gdal_vers.rstrip())
    else: status = -1
    pb.end(status)

    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    pb = utils._progress('processing region(s)...')
    try: 
        these_regions = [regions.region(iregion)]
    except:
        if os.path.exists(iregion):
            _poly = ogr.Open(iregion)
            _player = _poly.GetLayer(0)
            for pf in _player:
                _pgeom = pf.GetGeometryRef()
                these_regions.append(regions.region('/'.join(map(str, _pgeom.GetEnvelope()))))

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
          status = -1
    pb.opm = 'processing region(s)...{}'.format(len(these_regions))
    pb.end(status)
            
    if status == -1:
        print(_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):
        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        if idatalist is not None:
            pb = utils._progress('loading datalist...')
            this_datalist = datalists.datalist(idatalist, this_region)
            if not this_datalist._valid: 
                status = -1

            pb.opm = 'loading datalist...{}'.format(this_datalist._path_basename)
            pb.end(status)
        else: this_datalist = None

        if status != 0: 
            status = 0
            break

        ## ==============================================
        ## Initialize the DEM CLASS
        ## ==============================================

        this_surf = cudem(this_datalist, this_region, iinc, callback = lambda: stop_threads, oname = o_pre, obname = o_bn, verbose = want_verbose)        
        this_surf.node = node_reg

        for dem_mod in mod_opts.keys():

            if this_datalist is None:
                if dem_mod != 'conversion-grid':
                    status = -1
                    break

            ## ==============================================
            ## Run the DEM module
            ## ==============================================

            args = tuple(mod_opts[dem_mod])
                        
            pb = utils._progress('geomods: running {} module'.format(dem_mod))
            pb = utils._progress('running geomods dem module {} on region ({}/{}): {}\
            '.format(dem_mod, rn + 1, len(these_regions), this_region.region_string))

            this_surf._module = _dem_mods[dem_mod][0](this_surf)
            this_surf._module_args = args

            try:
                this_surf.start()
                while True:
                    time.sleep(2)
                    pb.update()
                    if not this_surf.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): 
                this_surf.status = -1
                stop_threads = True

            pb.end(this_surf.status)
            
if __name__ == '__main__':
    main()

### END
