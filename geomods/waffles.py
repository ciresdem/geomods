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
import math
import glob
import threading
import time

import numpy as np
try:
    import osgeo.ogr as ogr
    import osgeo.osr as osr
    import osgeo.gdal as gdal
except ImportError:
    try:
        import ogr
        import osr
        import gdal
    except ImportError:
        sys.exit(-1)

import regions
import datalists
import gdalfun
import utils

_version = '0.2.1'

## =============================================================================
##
## Uncertainty Analysis
##
## =============================================================================

def err2coeff(my_data):
    '''data is 2 col file with `err dist`'''

    from scipy import optimize
    import matplotlib.pyplot as plt
    #try: 
    #    my_data = np.loadtxt(data, delimiter=' ')
    #except: sys.exit(2)

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

    # fit
    
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    #plt.show()

    out_png = 'unc_best_fit.png'
    plt.savefig(out_png)   # save the figure to file
    plt.close()

    #scatter

    plt.scatter(distance, error)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')

    out_png = 'unc_scatter.png'
    plt.savefig(out_png)
    plt.close()

    return(out)

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)'''
    
    import fractions

    str_inc = str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', '')

    return(str_inc)

def this_year():
    '''return the current year'''

    import datetime

    return(datetime.datetime.now().strftime('%Y'))

def remove_glob(glob_str):
    '''glob `glob_str` and os.remove results'''

    globs = glob.glob(glob_str)
    if len(globs) > 0:
        for g in globs:
            try:
                os.remove(g)
            except: pass

## GMT Wrappers

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

def num_msk_gdal(datalist, region, inc, verbose = False):
    '''Generate a num-msk with GDAL from a datalist.'''

    o_mask = '{}_num_msk.tif'.format(datalist._path_dl_name)
    gdalfun.xyz_gmask(datalist._caty(), 
                      o_mask,
                      region.region, 
                      inc, verbose = verbose)
    
    return(o_mask)

def num_msk(num_grd, dst_msk, verbose = False):
    '''Generate a num-msk from a NUM grid.'''

    status = 0

    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    out, status = utils.run_cmd(num_msk_cmd, verbose, verbose)

    return(status)

## =============================================================================
##
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT and GDAL 
##
## =============================================================================

class dem:

    def __init__(self, i_datalist, i_region, i_inc = '0.000277777', o_name = None, o_b_name = None, callback = lambda: False, verbose = False):
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
            'mbgrid': lambda: self.mbgrid(),
            'surface': lambda: self.surface(),
            'num': lambda: self.num(),
            'mean': lambda: self.mean(),
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
        }

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]

            str_inc = inc2str_inc(self.inc)
            self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        else: self.o_name = o_b_name

    def run(self, dem_mod = 'mbgrid'):
        '''Run the DEM module `dem_mod` using args `dem_mod_args`.'''

        dems = self.dem
        if dem_mod != 'mbgrid':
            self.datalist._load_data()
            if len(self.datalist.datafiles) == 0:
                self.status = -1

        if self.status == 0:
            dems = self.dem_mods[dem_mod]()

        return(dems)

    ## ==============================================
    ## run mbgrid on the datalist and generate the
    ## dem and num grids
    ## note: mbgrid will cause popen to hang if stdout 
    ## is not cleared...should add output file to send to...
    ## ==============================================

    def mbgrid(self):
        '''Generate a DEM and num grid with MBSystem's mbgrid program.'''

        mbgrid_cmd = ('mbgrid -I{} {} -E{:.7f}/{:.7f}/degrees! -O{} -A2 -G100 -F1 -N -C10/3 -S0 -X0.1 -T35 > /dev/null 2> /dev/null \
        '.format(self.datalist._path, self.dist_region.gmt, self.inc, self.inc, self.o_name))
        out, self.status = utils.run_cmd(mbgrid_cmd, self.verbose, self.verbose)

        if self.status == 0:
            self.dem['dem-grd'] = '{}.grd'.format(self.o_name)

            if self.node == 'pixel':
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['dem-grd']), self.verbose, self.verbose)
                os.rename('tmp.grd', self.dem['dem-grd'])

            self.dem['dem'] = grd2tif(self.dem['dem-grd'])
            remove_glob('*.cmd')

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

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -G{}_num.grd -An {}\
        '.format(self.dist_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd_with_input(num_cmd0, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:
            self.dem['num-grd'] = '{}_num.grd'.format(self.o_name)
            self.dem['num'] = grd2tif(self.dem['num-grd'])
            self.dem['num-msk-grd'] = '{}_num_msk.grd'.format(self.o_name)

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

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -G{}_mean.grd {}\
        '.format(self.dist_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd_with_input(num_cmd0, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:
            self.dem['mean-grd'] = '{}_mean.grd'.format(self.o_name)
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
## with GMT, gererate a 'bathy-surface'.
## specify the coastline mask shapefile with 'mask = '
## if mask is not specified, will use gsshg.
## returns a netcdf masked grid and an XYZ file
## of the bathy surface points.
## ==============================================

class bathy_surface:

    def __init__(self, i_datalist, i_region, i_inc = '0.000277777', o_name = None, o_b_name = None, callback = lambda: False, verbose = False):
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.dem = { 
            'bathy-grd': None,
            'bathy-xyz': None,
        }

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]

            str_inc = inc2str_inc(self.inc)
            self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        else: self.o_name = o_b_name

    def run(self, mask = None):
        '''Run the bathy-surface module'''

        self.bathy_surface(mask = mask)

        return(self.dem)

    def bathy_surface(self, mask = None):
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
            remove_glob('{}_p*.grd'.format(self.o_name))

            if self.status == 0:
                self.dem['bathy-grd'] = bathy_dem
                bathy_xyz = 'xyz/{}_bs.xyz'.format(self.o_name)                
                self.status = grd2xyz(self.dem['bathy-grd'], bathy_xyz, want_datalist = True)

                if self.status == 0:
                    self.dem['bathy-xyz'] = bathy_xyz

        return(self.dem)

## =============================================================================
##
## conversion-grid module:
## Generate a VDatum conversion grid
##
## Requires NOAA's VDatum version >= 4.x
##
## =============================================================================

class conversion_grid:

    def __init__(self, i_datalist, i_region, i_inc = '0.000277777', o_name = None, o_b_name = None, callback = lambda: False, verbose = False):
        self.region = i_region
        self.inc = float(i_inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.dem = {}

        if o_b_name is None:
            if o_name is None: 
                o_name = 'cgrid'

            self.o_name = '{}{}_{}_{}'.format(o_name, inc2str_inc(self.inc), self.region.fn, this_year())
        else: self.o_name = o_b_name

    def run(self, ivert = 'navd88', overt = 'mllw', region = '3'):
        '''Run the conversion grid module.'''

        self.conversion_grid(ivert, overt, region)

    def conversion_grid(self, ivert = 'navd88', overt = 'mllw', region = '3'):
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
            remove_glob('empty.*')
            remove_glob('result/*')
            os.removedirs('result')
        except: pass

        return(self.dem)

## =============================================================================
##
## spatial-metadata module:
##
## =============================================================================

class spatial_metadata:

    def __init__(self, i_datalist, i_region, i_inc = '0.000277777', o_name = None, o_b_name = None, callback = lambda: False, verbose = False):

        import Queue as queue

        self.dl_q = queue.Queue()
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.v_fields = ['Name', 'Agency', 'Date', 'Type',
                         'Resolution', 'HDatum', 'VDatum', 'URL']
        self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                         ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]

            self.o_name = '{}{}_{}_{}'.format(o_name, inc2str_inc(self.inc), self.region.fn, this_year())
        else: self.o_name = o_b_name

        #self.pb = utils._progress('generating spatial metadata from {}'.format(self.datalist._path_basename))

    def run(self, epsg = 4269):
        '''Run the spatial-metadata module'''

        self.spatial_metadata(epsg)

    def _gather_from_queue(self):
        '''Gather geometries from a queue of [[datalist, layer], ...].'''

        while True:
            sm_args = self.dl_q.get()
            dl = sm_args[0]
            layer = sm_args[1]
            
            self._gather_from_datalist(dl, layer)
            
            self.dl_q.task_done()

    def _gather_from_datalist(self, dl, layer):
        '''gather geometries from datalist `dl` and append
        results to ogr `layer`. Load the datalist, generate
        a NUM-MSK grid, polygonize said NUM-MSK then union
        the polygon and add it to the output layer.'''

        defn = layer.GetLayerDefn()
        this_datalist = datalists.datalist(dl[0], self.dist_region)

        try:
            o_v_fields = [dl[3], dl[4], dl[5], dl[6], dl[7], dl[8], dl[9], dl[10].strip()]
        except: o_v_fields = [this_datalist._path_dl_name, 'Unknown', 0, 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

        this_dem = dem(this_datalist, self.region, str(self.inc), verbose = self.verbose)
        this_dem.run('num')

        if this_dem.status == 0:
            os.remove(this_dem.dem['num'])
            os.remove(this_dem.dem['num-grd'])
            os.remove(this_dem.dem['num-msk-grd'])

            pb = utils._progress('gathering geometries from \033[1m{}\033[m...'.format(this_datalist._path_basename))
            src_ds = gdal.Open(this_dem.dem['num-msk'])
            srcband = src_ds.GetRasterBand(1)

            remove_glob('{}_poly.*'.format(this_dem.o_name))

            tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_dem.o_name))
            tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_dem.o_name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

            pb1 = utils._progress('polygonizing \033[1m{}\033[m mask...'.format(this_datalist._path_basename))
            t = threading.Thread(target = gdal.Polygonize, args = (srcband, None, tmp_layer, 0, [], None))
            t.daemon = True

            t.start()
            while True:
                time.sleep(2)
                pb1.update()
                if not t.is_alive():
                    break
            t.join()

            pb1.opm = 'polygonized \033[1m{}\033[m mask.'.format(this_datalist._path_basename)
            pb1.end(self.status)

            src_ds = srcband = None

            if len(tmp_layer) > 1:
                out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn)
                for i, f in enumerate(self.v_fields):
                    out_feat.SetField(f, o_v_fields[i])

                layer.CreateFeature(out_feat)

            tmp_ds = tmp_layer = None
            remove_glob('{}_poly.*'.format(this_dem.o_name))
            os.remove(this_dem.dem['num-msk'])

            pb.opm = 'gathered geometries from \033[1m{}\033[m.'.format(this_datalist._path_basename)
            pb.end(self.status)

    def spatial_metadata(self, epsg = 4269):
        '''Geneate spatial metadata from the datalist'''

        if self.inc < 0.0000925:
            print 'warning, increments less than 1/3 arc-second may be slow.'
    
        dst_vec = '{}_sm.shp'.format(self.o_name, self.region.fn)
        dst_layername = os.path.basename(dst_vec).split('.')[0]

        remove_glob('{}_sm.*'.format(self.o_name))

        ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
        if ds is not None:
            try:
                int(epsg)
            except: epsg = 4326

            gdalfun._prj_file('{}.prj'.format(os.path.basename(dst_vec).split('.')[0]), int(epsg))
            layer = ds.CreateLayer('{}'.format(os.path.basename(dst_vec).split('.')[0]), None, ogr.wkbMultiPolygon)

            for i, f in enumerate(self.v_fields):
                layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i]))
            
            for feature in layer:
                layer.SetFeature(feature)
    
            for _ in range(3):
                t = threading.Thread(target = self._gather_from_queue, args = ())
                t.daemon = True
                t.start()

            for dl in self.datalist.datalist:
                #self._gather_from_datalist(dl, layer)
                self.dl_q.put([dl, layer])

            self.dl_q.join()

        ds = layer = None
        if not os.path.exists(dst_vec):
            dst_vec = None

        return(dst_vec)

## =============================================================================
##
## uncertainty module:
##
## =============================================================================

class uncertainty:

    def __init__(self, i_datalist, i_region, i_inc = '0.000277777', o_name = None, o_b_name = None, callback = lambda: False, verbose = False):

        import random

        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'
        self.o_b_name = o_b_name

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.dem = { 
            'dem': None,
            'num': None,
            'num-grd': None,
            'num-msk': None,
            'prox': None,
            'int-unc': None,
        }

        if o_b_name is None:
            if o_name is None: 
                o_name = self.datalist._path_basename.split('.')[0]
            str_inc = inc2str_inc(self.inc)
            self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        else: self.o_name = o_b_name

    def run(self, dem_mod = 'mbgrid', dem = None, num = None, num_msk = None, prox = None):

        if dem is not None:
            self.dem['dem'] = dem

        if num is not None:
            self.dem['num'] = num

        if num_msk is not None:
            self.dem['num-msk'] = num_msk

        if prox is not None:
            self.dem['prox'] = prox

        self.datalist._load_data()
        self.interpolation_uncertainty(dem_mod)

    def set_or_make_dem(self):
        '''check if dem dict contains dems, otherwise generate them...'''

        tw = utils._progress('checking for DEMs...')

        if self.dem['dem'] is None:
            tw.err_msg('generating dem...')
            dems = dem(self.datalist, self.region, str(self.inc)).run(dem_mod)
            for key in dems.keys():
                self.dem[key] = dems[key]

        if self.dem['num'] is None:
            tw.err_msg('generating NUM grid...')
            dems = dem(self.datalist, self.region, str(self.inc)).run('num')
            for key in dems.keys():
                self.dem[key] = dems[key]

        if self.dem['num-msk'] is None:
            tw.err_msg('generating NUM mask...')
            self.dem['num-grd-msk'] = '{}_msk.grd'.format(self.dem['num'].split('.')[0]) 
            num_msk(self.dem['num'], self.dem['num-grd-msk'])
            self.dem['num-msk'] = grd2tif(self.dem['num-grd-msk'])

        if self.dem['prox'] is None:
            tw.err_msg('generating proximity grid...')
            self.dem['prox']  = '{}_prox.tif'.format(self.dem['num'].split('.')[0]) 
            proximity(self.dem['num-msk'], self.dem['prox'])

        tw.opm = 'checked for DEMs.'
        tw.end(self.status)        

    def gmtselect_split(self, o_xyz, sub_region, sub_bn):
        '''split an xyz file into an inner and outer region.'''

        out_inner = None
        out_outer = None

        gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
        out, self.status = utils.run_cmd(gmt_s_inner, self.verbose, False)
        
        if self.status == 0:
            out_inner = '{}_inner.xyz'.format(sub_bn)

        gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
        out, self.status = utils.run_cmd(gmt_s_outer, self.verbose, False)

        if self.status == 0:
            out_outer = '{}_outer.xyz'.format(sub_bn)

        return([out_inner, out_outer])

    def interpolation_uncertainty(self, dem_mod = 'mbgrid'):
        '''calculate the interpolation uncertainty.'''

        loop_count = 1
        sub_count = 0        
        this_region = self.region
        dp = None

        self.set_or_make_dem()
        print self.dem

        ## ==============================================
        ## Calculate the percentage of filled cells
        ## and proximity percentiles.
        ## ==============================================
        tw = utils._progress('')

        num_sum = gdalfun.sum(self.dem['num-msk'])
        gi = grdinfo(self.dem['num-msk'])
        g_max = int(gi[9]) * int(gi[10])

        num_perc = (num_sum / g_max) * 100
        prox_perc_95 = gdalfun.percentile(self.dem['prox'], 75)
        tw.err_msg('waffles: proximity 95th perc: {}'.format(prox_perc_95))

        sub_regions = this_region.chunk(self.inc, 500)
        tw.err_msg('waffles: chunking into {} regions.'.format(len(sub_regions)))
        #tw.err_msg('waffles: processing {} regions {} times.'.format(len(sub_regions), loop_count))

        for sub_region in sub_regions:
            if self.stop(): break

            self.status = 0
            sub_count += 1
            o_xyz = '{}.xyz'.format(self.o_name)

            tw = utils._progress('processing sub region \033[1m{}\033[m...'.format(sub_count))

            self.status = grd2xyz(self.dem['dem'], o_xyz, region = sub_region.buffer(10*self.inc), mask = self.dem['num'])
            if os.stat(o_xyz).st_size == 0:
                tw.err_msg('waffles: error, no data in sub-region...')
                self.status = -1
            else:
                s_inner, s_outer = self.gmtselect_split(o_xyz, sub_region, 'sub_{}'.format(sub_count))

                if os.stat(s_inner).st_size != 0:
                    sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                else: self.status = -1

                if self.status == 0 and not self.stop():
                    self.status = grdcut(self.dem['num-msk'], sub_region, 'tmp_sub.grd')

                    gi = grdinfo('tmp_sub.grd')
                    num_max = int(gi[9]) * int(gi[10])

                    tw.err_msg('waffles: total cells in subgrid is {}'.format(num_max))
                    num_sub_sum = gdalfun.sum(grd2tif('tmp_sub.grd'))
                    tw.err_msg('waffles: total filled cells in subgrid is {}'.format(num_sub_sum))

                    try:
                        os.remove('tmp_sub.grd')
                        os.remove('tmp_sub.tif')
                    except: pass

                    num_sub_perc = (num_sub_sum / num_max) * 100
                    tw.err_msg('waffles: {}% of cells have data'.format(num_sub_perc))
                    s_size = 100 - num_sub_perc

                    if s_size >= 100:
                        n_loops = 0
                    else: n_loops = int(((s_size / num_perc) * 2) + 1)

                    tw.err_msg('waffles: loops for this subregion is {}'.format(n_loops))

                    ## ==============================================
                    ## Split Sample n_loops times
                    ## ==============================================

                    for i in range(0, n_loops):
                        if self.stop(): break

                        np.random.shuffle(sub_xyz)
                        sx_len = len(sub_xyz)
                        sx_len_pct = int(sx_len * (num_sub_perc / 100))

                        if sx_len_pct == 0:
                            break

                        tw.err_msg('waffles: extracting {} random data points out of {}'.format(sx_len_pct, sx_len))
                        sub_xyz_head = 'sub_{}_head.xyz'.format(sub_count)

                        np.savetxt(sub_xyz_head, sub_xyz[:sx_len_pct], '%f', ' ')

                        sub_datalist =  datalists.datalist('sub_{}.datalist'.format(sub_count), sub_region)
                        sub_datalist._append_datafile(s_outer, 168, 1)
                        sub_datalist._append_datafile(sub_xyz_head, 168, 1)
                        sub_datalist.reset()

                        pb = utils._progress('generating sub-region {} random-sample grid...{}/{}'.format(sub_count, i, n_loops))

                        sub_surf = dem(sub_datalist, sub_region, str(self.inc))
                        sub_surf.run(dem_mod)
                        sub_surf.run('num')
                        sub_dems = sub_surf.dem

                        pb.opm = 'generated sub-region {} random-sample grid.{}/{}'.format(sub_count, i, n_loops)
                        pb.end(sub_surf.status)

                        if sub_dems['dem'] is not None:

                            sub_prox = '{}_prox.tif'.format(sub_dems['num'].split(',')[0])
                            self.status = proximity(sub_dems['num-msk'], sub_prox)
                            
                            sub_xyd = gdalfun.query(sub_xyz[sx_len_pct:], sub_dems['dem'], 'xyd')
                            sub_dp = gdalfun.query(sub_xyd, sub_prox, 'zg')

                            if len(sub_dp) != 0:
                                if dp is None:
                                    dp = sub_dp
                                else: dp = np.concatenate((dp, sub_dp), axis = 0)
                        os.remove(sub_xyz_head)
                        os.remove(sub_datalist._path)

                remove_glob('sub_{}*'.format(sub_count))

            tw.opm = 'processed sub region \033[1m{}\033[m.'.format(sub_count)
            tw.end(self.status)

        ## ==============================================
        ## calculate the error coefficients and plot results
        ## ==============================================

        if not self.stop():
            dp = dp[dp[:,1]<prox_perc_95,:]
            dp = dp[dp[:,1]>0,:]
            ec = err2coeff(dp)
            print('error coefficient: {}'.format(ec))
            np.savetxt('test.err', dp, '%f', ' ')

        ## ==============================================
        ## apply error coefficient to full proximity grid
        ## ==============================================

## =============================================================================
##
## Console Scripts
##
## Scripts are: 
## waffles
##
## =============================================================================

_dem_mods = { 
    'dem': [lambda d, r, i, o, b, c, v: dem(d, r, i, o, b, c, v), 'generate a DEM [:grid-module]'],
    'bathy-surface': [lambda d, r, i, o, b, c, v: bathy_surface(d, r, i, o, b, c, v), 'generate a `bathy-surface` masked grid [:coastpoly]'],
    'conversion-grid': [lambda d, r, i, o, b, c, v: conversion_grid(d, r, i, o, b, c, v), 'generate a VDatum conversion grid [:i_vdatum:o_vdatum:vd_region]'],
    'spatial-metadata': [lambda d, r, i, o, b, c, v: spatial_metadata(d, r, i, o, b, c, v), 'generate spatial-metadata [:epsg]'],
    'uncertainty': [lambda d, r, i, o, b, c, v: uncertainty(d, r, i, o, b, c, v), 'calculate DEM uncertainty [:grid-module:dem:num]'],
}

def dem_mod_desc(x):
    dmd = []
    for key in x: 
        dmd.append('{:16}\t{}'.format(key, x[key][1]))
    return dmd

_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

usage: {} [ -hvAIER [ args ] ] module[:option1:option2] ...

Modules:
  {}

 note: grid-modules include: `mbgrid`, `surface`, `num`, `mean`

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
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 dem:surface
 % {} --datalist input.datalist --increment 0.000277777 --region input_tiles_ply.shp dem:mbgrid spatial-metadata
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
    iregion = None
    idatalist = None
    igrid = None
    iinc = '0.000277777'
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

    cmd_vers = utils._cmd_check()

    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    pb = utils._progress('loading region(s)...')
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
    pb.opm = 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions))
    pb.end(status)
            
    if status == -1:
        print(_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):

        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        if stop_threads:
            break
        
        if idatalist is not None:
            this_datalist = datalists.datalist(idatalist, this_region)
            if not this_datalist._valid: 
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

            dl = _dem_mods[dem_mod][0](this_datalist, this_region, iinc, o_pre, o_bn, lambda: stop_threads, want_verbose)
            dl.node = node_reg

            t = threading.Thread(target = dl.run, args = args)

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
