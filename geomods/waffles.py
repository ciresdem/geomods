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

_version = '0.2.0'

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
    '''Run a number of modules for DEM generation and analysis'''

    def __init__(self, idatalist, iregion, iinc = '0.000277777', oname = None, obname = None, callback = lambda: False, verbose = False):
        threading.Thread.__init__(self)

        self.status = 0
        self.stop = callback
        self.verbose = verbose
        #self.dem = {}
        self.dem = { 
            'dem': None,
            'dem-grd': None,
            'num': None,
            'num-grd': None,
            'num-msk': None,
            'prox': None,
            'slope': None,
            'xyz': None,
        }

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
        for key in self.dem:
            if self.dem[key] is not None:
                if not os.path.exists(self.dem[key]):
                    return(False)
        return(True)

    def _rectify(self):
        for key in self.dem:
            if self.dem[key] is not None:
                if not os.path.exists(self.dem[key]):
                    self.dem[key] = None
                
    def _set_dem(self, key, value):
        if key in self.dem.keys():
            self.dem[key] = value

    def _print_dem(self):
        for key in self.dem:
            print(key, self.dem[key])

    def _dem_remove(self, key):
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
            out, self.status = utils.run_cmd(grd2tif_cmd, self.verbose, self.verbose)
        else: self.status = -1

        return(self.status)

    def grdinfo(self, src_grd):
        '''Return an info list of `src_grd`'''

        if not os.path.exists(src_grd):
            self.status = -1
            return([])
        else:
            grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
            out, self.status = utils.run_cmd(grdinfo_cmd, self.verbose, False)
            return(out.split())

    def grdcut(self, src_grd, src_region, dst_grd):
        '''Cut `src_grd` to `src_region` '''

        if not os.path.exists(src_grd):
            self.status = -1
        else:
            cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
            out, self.status = utils.run_cmd(cut_cmd1, self.verbose, self.verbose)

        return(self.status)

    def grd2xyz(self, src_grd, dst_xyz, region = None, mask = None):
        '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region'''

        if mask:
            self.grdmask_cmd = ('gmt grdmath -V {} {} OR = tmp.grd'.format(src_grd, mask))

            out, self.status = utils.run_cmd(self.grdmask_cmd, self.verbose, self.verbose)
            if self.status == 0: src_grd = 'tmp.grd'

        if region and region._valid:
            self.grd2xyz_cmd = ('gmt grd2xyz -V {} {} -s > {}'.format(src_grd, region.gmt, dst_xyz))

        else: self.grd2xyz_cmd = ('gmt grd2xyz {} -V -s > {}'.format(src_grd, dst_xyz))
            
        out, self.status = utils.run_cmd(self.grd2xyz_cmd, self.verbose, self.verbose)

        if self.status == 0:
            if mask:
                if os.path.exists('tmp.grd'):
                    os.remove('tmp.grd')
        
        return(self.status)
        
    ## ==============================================
    ## Main DEM Modules
    ## ==============================================

    def slope(self):
        '''Generate a Slope grid from the DEM with GMT'''

        if self.dem['dem-grd'] is None:
            self.surface()

        if self.status == 0:
            slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
            '.format(self.dem['dem-grd'], self.oname, self.dem['dem-grd']))
            out, self.status = utils.run_cmd(slope_cmd0, self.verbose, True)

            if self.status == 0:
                slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}_slp.tif=gd+n-9999:GTiff\
                '.format(self.oname, self.oname))
                out, self.status = utils.run_cmd(slope_cmd1, self.verbose, True)

                if os.path.exists('{}_pslp.grd'.format(self.oname)):
                    os.remove('{}_pslp.grd'.format(self.oname))

                self.dem['slope'] = ('{}_slp.tif'.format(self.oname))

        return(self.status)

    def proximity(self):
        '''Generate a Proximity grid with GDAL from NUM-MSK'''

        if self.dem['num-msk'] is None:
            #self.status = -1
            self.num()

        if self.status == 0:
            self.dem['prox'] = ('{}_prox.tif'.format(self.oname))

            gdalfun.proximity(self.dem['num-msk'], self.dem['prox'])
        
        return(self.status)

    def num(self):
        '''Generate a num and num-msk grid with GMT'''

        if self.node == 'pixel':
            num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -r -G{}_num.grd -An\
            '.format(self.dist_region.gmt, self.inc, self.oname))
        else:
            num_cmd0 = ('gmt xyz2grd -V {} -I{:.7f} -G{}_num.grd -An\
            '.format(self.dist_region.gmt, self.inc, self.oname))

        ## ==============================================
        ## Generate NUM grid with GMT xyz2grd -An
        ## ==============================================

        out, self.status = utils.run_cmd_with_input(num_cmd0, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:

            ## ==============================================
            ## transform num grid to num-msk (0/1)
            ## ==============================================

            num_cmd1 = ('gmt grdmath -V {}_num.grd 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
            '.format(self.oname, self.oname))
            out, self.status = utils.run_cmd(num_cmd1, self.verbose, self.verbose)

            if self.status == 0:

                ## ==============================================
                ## Convert NUM netcdf to GeoTiff
                ## ==============================================

                self.dem['num-grd'] = '{}_num.grd'.format(self.oname)
                self.status = self.grd2tif(self.dem['num-grd'])

                self.dem['num'] = '{}_num.tif'.format(self.oname)
                self.dem['num-msk'] = '{}_num_msk.tif'.format(self.oname)

        return(self.status)

    def num_msk(self):
        '''Generate a num-msk with GDAL'''

        gdalfun.xyz_gmask(self.datalist._caty(), 
                          '{}_num_msk.tif'.format(self.oname), 
                          self.dist_region.region, 
                          self.inc, verbose = self.verbose)

        self.dem['num-msk'] = '{}_num_msk.tif'.format(self.oname)

    def masked_surface(self, mask = None):
        '''Generate a masked surface with GMT surface and a breakline mask polygon.
        Use a coastline polygon mask to generatre a `bathy-surface` DEM.'''

        masked_dem = ('{}_bs.grd'.format(self.oname))

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

        ## ==============================================
        ## Generate the Landmask via gsshg if no mask
        ## was provided...
        ## ==============================================

        if mask is None:
            out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, self.verbose)
        
        ## ==============================================
        ## Generate the bounded surface DEM with GMT
        ## ==============================================

        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:

            ## ==============================================
            ## Mask the bounded surface DEM
            ## ==============================================

            dem_landmask_cmd1 = ('gmt grdmath -V {}_p.grd tmp_lm.grd MUL 0 NAN = {}_p_bathy.grd\
            '.format(self.oname, self.oname))
            out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, True)

            ## ==============================================
            ## cut the DEM to desired output region
            ## ==============================================

            dem_cut_cmd = ('gmt grdcut -V {}_p_bathy.grd -G{} {}\
            '.format(self.oname, masked_dem, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)
            
            if os.path.exists('{}_p.grd'.format(self.oname)):
                os.remove('{}_p.grd'.format(self.oname))

            if os.path.exists('{}_p_bathy.grd'.format(self.oname)):
                os.remove('{}_p_bathy.grd'.format(self.oname))

            if self.status == 0:

                ## ==============================================
                ## Convert masked bounded DEM to XYZ
                ## ==============================================

                xyz_dir = os.path.join(os.getcwd(), 'xyz')
        
                if not os.path.exists(xyz_dir):
                    os.makedirs(xyz_dir)

                #self.dem['dem-bathy'] = ('{}_bs.grd'.format(self.oname))
                #self.dem['xyz-bathy'] = ('xyz/{}_bs.xyz'.format(self.oname))

                masked_xyz = ('xyz/{}_bs.xyz'.format(self.oname))

                #self.status = self.grd2xyz(self.dem['dem-bathy'], self.dem['xyz-bathy'])
                self.status = self.grd2xyz(masked_dem, masked_xyz)

                if self.status == 0:

                    ## ==============================================
                    ## Add xyz file to datalist
                    ## ==============================================

                    sdatalist = datalists.datalist(os.path.join(xyz_dir, 'bathy.datalist'))
                    sdatalist._append_datafile('{}'.format(os.path.basename(masked_xyz)), 168, 1)
                    sdatalist._reset()

                    ## ==============================================
                    ## Generate .inf file
                    ## ==============================================

                    out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, 'bathy.datalist')), False, True)

        return(self.status)

    def surface(self):
        '''Generate a DEM with GMT surface'''

        ## ==============================================
        ## Run GMT surface on the datalist
        ## generate dem at 10x cells to account for
        ## edge effects
        ## ==============================================

        if self.node == 'pixel':
            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -r -V | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -r -Lud -Lld\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))
        else:
            dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lud -Lld\
            '.format(self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))

        out, self.status = utils.run_cmd_with_input(dem_surf_cmd, self.datalist._cat_port, self.verbose, True)

        if self.status == 0:

            ## ==============================================
            ## cut the DEM to desired output region
            ## ==============================================

            dem_cut_cmd = ('gmt grdcut -V {}_p.grd -G{}.grd {}\
            '.format(self.oname, self.oname, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)

            if os.path.exists('{}_p.grd'.format(self.oname)):
                os.remove('{}_p.grd'.format(self.oname))            

            if self.status == 0:
                
                ## ==============================================
                ## convert netcdf DEM to GeoTiff
                ## ==============================================

                self.dem['dem-grd'] = ('{}.grd'.format(self.oname))
                self.grd2tif(self.dem['dem-grd'])
                self.dem['dem'] = ('{}.tif'.format(self.oname))

        return(self.status)

    def mbgrid(self):
        '''Generate a DEM and num grid with MBSystem'''

        ## ==============================================
        ## run mbgrid on the datalist and generate the
        ## dem and num grids
        ## note: mbgrid will cause popen to hang if stdout 
        ## is not cleared...should add output file to send to...
        ## ==============================================

        mbgrid_cmd = ('mbgrid -I{} {} -E{:.7f}/{:.7f}/degrees! -O{} -A2 -G100 -F1 -N -C10/3 -S0 -X0.1 -T35 -M > /dev/null 2> /dev/null \
        '.format(self.datalist._path, self.dist_region.gmt, self.inc, self.inc, self.oname))
        out, self.status = utils.run_cmd(mbgrid_cmd, self.verbose, self.verbose)

        if self.status == 0:
            self.dem['dem-grd'] = '{}.grd'.format(self.oname)
            self.dem['num-grd'] = '{}_num.grd'.format(self.oname)

            ## ==============================================
            ## resample to pixel-node reg. if desired
            ## ==============================================

            if self.node == 'pixel':
                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['dem-grd']), self.verbose, self.verbose)
                os.rename('tmp.grd', self.dem['dem-grd'])

                out, self.status = utils.run_cmd('gmt grdsample -T {} -Gtmp.grd'.format(self.dem['num-grd']), self.verbose, self.verbose)
                os.rename('tmp.grd', self.dem['num-grd'])                

            ## ==============================================
            ## convert netcdf grd files to geotiff
            ## ==============================================

            self.grd2tif(self.dem['dem-grd'])
            self.dem['dem'] = '{}.tif'.format(self.oname)

            self.grd2tif(self.dem['num-grd'])
            self.dem['num'] = '{}_num.tif'.format(self.oname)

            ## ==============================================
            ## Cleanup mbgrid junk
            ## ==============================================

            try:
                os.remove('{}.cmd'.format(self.dem['dem-grd']))
                os.remove('{}.cmd'.format(self.dem['num-grd']))
                os.remove('{}.mb-1'.format(self.oname))
                os.remove('{}_sd.grd'.format(self.oname))
                os.remove('{}_sd.grd.cmd'.format(self.oname))
            except: pass

            ## ==============================================
            ## transform num grid to num-msk (0/1)
            ## ==============================================

            num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
            '.format(self.dem['num-grd'], self.oname))
            out, self.status = utils.run_cmd(num_msk_cmd, self.verbose, self.verbose)

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
                gc = 'gmt blockmean result/empty.xyz -V -I{} {} \
                | gmt surface -I{} {} -G{}_{}to{}.tif=gd:GTiff -V -r -T0 {} {}\
                '.format(self.inc, self.dist_region.gmt, self.inc, self.dist_region.gmt,
                         self.oname, this_vd.ivert,
                         this_vd.overt, ll_switch, lu_switch)
            else:
                gc = 'gmt blockmean result/empty.xyz -V -I{} {} \
                | gmt surface -I{} {} -G{}_{}to{}.tif=gd:GTiff -V -T0 {} {}\
                '.format(self.inc, self.dist_region.gmt, self.inc, self.dist_region.gmt,
                         self.oname, this_vd.ivert, this_vd.overt, ll_switch, lu_switch)
            utils.run_cmd(gc, self.verbose, 'generating conversion grid')

            dem_name = '{}to{}'.format(this_vd.ivert, this_vd.overt)
            self.dem[dem_name] = dem_name
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

        if self.inc < 0.0000925:
            print 'warning, increments less than 1/3 arc-second may be slow.'

        ## ==============================================
        ## these fields should be found in the datalist 
        ## starting at position 3 (from 0)
        ## ==============================================

        v_fields = ['Name', 'Agency', 'Date', 'Type',
                    'Resolution', 'HDatum', 'VDatum', 'URL']
        t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                    ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
    
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
                    this_datalist._load_data()
                    this_dem = cudem(this_datalist, self.region, str(self.inc), verbose = self.verbose)

                    pb.opm = 'loaded datalist <<\033[1m{}\033[m>>'.format(this_datalist._path_basename)
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

                    pb = utils._progress('gathering geometries from \033[1m{}\033[m...'.format(this_datalist._path_basename))

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
                        
                        pb1 = utils._progress('polygonizing \033[1m{}\033[m mask...'.format(this_datalist._path_basename))

                        t = threading.Thread(target = gdal.Polygonize, args = (srcband, None, tmp_layer, 0, [], None))
                        t.daemon = True

                        try:
                            t.start()
                            while True:
                                time.sleep(2)
                                pb1.update()
                                if not t.is_alive():
                                    break
                        except (KeyboardInterrupt, SystemExit): 
                            self.status = -1
                        t.join()

                        pb1.opm = 'polygonized \033[1m{}\033[m mask.'.format(this_datalist._path_basename)
                        pb1.end(self.status)

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

                    this_dem._dem_remove('num-msk')
                    pb.opm = 'gathered geometries from \033[1m{}\033[m.'.format(this_datalist._path_basename)
                    pb.end(self.status)
        ds = layer = None

    def uncertainty(self, grid_mod = None, dem = None):
        sub_count = 0        
        dp = None
        #tw = utils._progress('generating interpolation uncertainty')
        this_region = self.region

        if grid_mod is None:
            grid_mod = 'mbgrid'

        tw = utils._progress('checking for DEMs...')
        self._rectify()

        ## ==============================================
        ## Set or generate the DEM Etc.
        ## ==============================================

        if dem is not None:
            self.dem['dem'] = '{}.tif'.format(dem)
            self.dem['num'] = '{}_num.tif'.format(dem)
            self.dem['num-msk'] = '{}_num_msk.tif'.format(dem)
            self.dem['prox'] = '{}_prox.tif'.format(dem)

        self._rectify()
            
        if self.dem['dem'] is None:
            tw.err_msg('generating grid...')
            self.mbgrid()

        if self.dem['num'] is None or self.dem['num-msk'] is None:
            tw.err_msg('generating NUM grid...')
            self.num()

        if self.dem['prox'] is None:
            tw.err_msg('generating proximity grid...')
            self.proximity()

        tw.opm = 'checked for DEMs.'
        tw.end(self.status)

        ## ==============================================
        ## Calculate the percentage of filled cells
        ## and proximity percentiles.
        ## ==============================================

        num_sum = gdalfun.sum(self.dem['num-msk'])
        gi = self.grdinfo(self.dem['num-msk'])
        g_max = int(gi[9]) * int(gi[10])

        num_perc = (num_sum / g_max) * 100

        prox_perc_95 = gdalfun.percentile(self.dem['prox'], 75)
        tw.err_msg('waffles: proximity 95th perc: {}'.format(prox_perc_95))

        num_perc_5 = gdalfun.percentile(self.dem['num'], 5)
        tw.err_msg('waffles: num 5th perc: {}'.format(num_perc_5))

        ## ==============================================
        ## Split the input region into 200x200 cell chunks
        ## ==============================================

        sub_regions = this_region.chunk(self.inc, 1000)
        tw.err_msg('waffles: chunking into {} regions.'.format(len(sub_regions)))

        for j in range(0, 10):
            if self.stop():
                break

            for sub_region in sub_regions:
                if self.stop():
                    break
                ## ==============================================  
                ## Extract XYZ data for sub-region and 
                ## randomly sample
                ##
                ## Estimate data density and data sample-size
                ## and looping variable
                ##
                ## ==============================================

                self.status = 0
                sub_count += 1
                o_xyz = '{}.xyz'.format(self.oname)

                tw = utils._progress('processing sub region \033[1m{}\033[m...'.format(sub_count))

                ## ==============================================
                ## convert DEM to sub-region XYZ masking out
                ## interpolated cells
                ## ==============================================

                self.status = self.grd2xyz(self.dem['dem'], o_xyz, region = sub_region.buffer(10*self.inc), mask = self.dem['num'])

                if os.stat(o_xyz).st_size == 0:
                    tw.err_msg('waffles: error, no data in sub-region...')
                    self.status = -1
                else:
                    ## ==============================================
                    ## split the XYZ data into inner and outer areas for processing
                    ## and load the inner XYZ for random grid generation
                    ## ==============================================

                    gmt_s_inner = 'gmt gmtselect -V {} {} > sub_{}_inner.xyz'.format(o_xyz, sub_region.gmt, sub_count)
                    out, self.status = utils.run_cmd(gmt_s_inner, self.verbose, False)

                    gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > sub_{}_outer.xyz'.format(o_xyz, sub_region.gmt, sub_count)
                    out, self.status = utils.run_cmd(gmt_s_outer, self.verbose, False)

                    if os.stat('sub_{}_inner.xyz'.format(sub_count)).st_size != 0:
                        sub_xyz = np.loadtxt('sub_{}_inner.xyz'.format(sub_count), ndmin=2, delimiter = ' ')
                    else: self.status = -1
                    
                    if self.status == 0 and not self.stop():
                    
                        ## ==============================================
                        ## Estimate density and sample size for sub-region
                        ## ==============================================
                        self.status = self.grdcut(self.dem['num-msk'], sub_region, 'tmp_sub.grd')

                        gi = self.grdinfo('tmp_sub.grd')
                        #if len(gi) != 0:
                        num_max = int(gi[9]) * int(gi[10])

                        tw.err_msg('waffles: total cells in subgrid is {}'.format(num_max))

                        self.grd2tif('tmp_sub.grd')
                        num_sub_sum = gdalfun.sum('tmp_sub.tif')
                        tw.err_msg('waffles: total filled cells in subgrid is {}'.format(num_sub_sum))

                        os.remove('tmp_sub.grd')
                        os.remove('tmp_sub.tif')

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
                            if self.stop():
                                break

                            ## ==============================================
                            ## randomly shuffle the inner_xyz and extract
                            ## percentage of data points for DEM generation
                            ## ==============================================

                            np.random.shuffle(sub_xyz)
                            sx_len = len(sub_xyz)
                            sx_len_pct = int(sx_len * (num_sub_perc / 100))

                            if sx_len_pct == 0:
                                break

                            tw.err_msg('waffles: extracting {} random data points out of {}'.format(sx_len_pct, sx_len))
                            np.savetxt('sub_{}_rest.xyz'.format(sub_count), sub_xyz[:sx_len_pct], '%f', ' ')

                            sub_datalist =  datalists.datalist('sub_{}.datalist'.format(sub_count), sub_region)
                            sub_datalist._append_datafile('sub_{}_rest.xyz'.format(sub_count), 168, 1)
                            sub_datalist._append_datafile('sub_{}_outer.xyz'.format(sub_count), 168, 1)
                            sub_datalist._reset()
                            sub_datalist._load_data()
                            
                            ## ==============================================
                            ## Generate the random-sample DEM            
                            ## ==============================================

                            sub_surf = cudem(sub_datalist, sub_region, self.inc)
                            sub_surf._module = _dem_mods[grid_mod][0](sub_surf)
                            sub_surf._module_args = ()
                            #self.status = sub_surf.mbgrid()

                            pb = utils._progress('generating sub-region {} random-sample grid...{}/{}'.format(sub_count, i, n_loops))
                            sub_surf.start()
                            while True:
                                time.sleep(2)
                                pb.update()
                                if not sub_surf.is_alive():
                                    break

                            pb.opm = 'generated sub-region {} random-sample grid.{}/{}'.format(sub_count, i, n_loops)
                            pb.end(sub_surf.status)

                            sub_surf._rectify()

                            if sub_surf.dem['dem'] is not None:

                                self.status = sub_surf.proximity()

                                ## ==============================================            
                                ## Query the random-sample DEM with the 
                                ## remaining randomized data
                                ## ==============================================

                                sub_xyd = gdalfun.query(sub_xyz[sx_len_pct:], sub_surf.dem['dem'], 'xyd')

                                ## ==============================================
                                ## Query the random-sample proximity grid 
                                ## with the diff results from above and append
                                ## the results to the error/distance array
                                ## ==============================================

                                sub_dp = gdalfun.query(sub_xyd, sub_surf.dem['prox'], 'zg')

                                #sub_len += len(sub_dp)

                                if len(sub_dp) != 0:
                                    if dp is None:
                                        dp = sub_dp
                                    else: dp = np.concatenate((dp, sub_dp), axis = 0)

                                #if sub_len >= 100000:
                                #    break

                ## ==============================================
                ## Cleanup
                ## ==============================================

                fl = glob.glob('sub_{}*'.format(sub_count))
                for f in fl:
                    try: 
                        os.remove(f)
                    except: pass

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
## Run cudem.py from command-line
##
## DEM Modules are
## { mod-name: [module-function, module-description, module-function arguments]}
##
## =============================================================================

_dem_mods = { 
    'mbgrid': [lambda x: x.mbgrid, 'generate a DEM with mbgrid', 'None'],
    'gmt-surface': [lambda x: x.surface, 'generate a DEM with GMT', 'None'],
    'spatial-metadata': [lambda x: x.spatial_metadata, 'generate spatial-metadata', 'epsg'],
    'conversion-grid': [lambda x: x.conversion_grid, 'generate a conversion grid with vdatum', 'ivert:overt:region'],
    'bathy-surface': [lambda x: x.masked_surface, 'generate a bathymetry surface with GMT', 'coastline'],
    'uncertainty': [lambda x: x.uncertainty, 'generate an uncertainty grid', 'gridding-module:dem']
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
            pb = utils._progress('loading datalist...')
            this_datalist = datalists.datalist(idatalist, this_region)
            if not this_datalist._valid: 
                status = -1
            pb.opm = 'loaded datalist <<\033[1m{}\033[m>>'.format(this_datalist._path_basename)
            pb.end(status)
        else: this_datalist = None

        if status != 0: 
            status = 0
            break

        ## ==============================================
        ## Initialize the DEM CLASS
        ## ==============================================

        this_surf = cudem(this_datalist,
                          this_region,
                          iinc,
                          callback = lambda: stop_threads,
                          oname = o_pre,
                          obname = o_bn,
                          verbose = want_verbose)        
        this_surf.node = node_reg

        for dem_mod in mod_opts.keys():

            if dem_mod == 'gmt-surface' or \
               dem_mod == 'bathy-surface' or \
               dem_mod == 'uncertainty':
                if this_datalist is not None:
                    this_surf.datalist._load_data() # = this_datalist._proc()

            if this_datalist is None:
                if dem_mod != 'conversion-grid':
                    status = -1
                    break

            ## ==============================================
            ## Run the DEM module
            ## ==============================================

            args = tuple(mod_opts[dem_mod])
                        
            pb = utils._progress('running geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
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

            pb.opm = 'ran geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod, rn + 1, len(these_regions), this_region.region_string)
            pb.end(this_surf.status)
            
if __name__ == '__main__':
    main()

### END
