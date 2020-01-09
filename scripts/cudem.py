#!/usr/bin/env python
### cudem.py
##
## Copyright (c) 2013 - 2020 Matthew Love <matthew.love@colorado.edu>
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
import matplotlib.pyplot as plt
import threading
import time

try:
    import osgeo.ogr as ogr
    import osgeo.gdal as gdal
    has_gdalpy = True
except ImportError:
    try:
        import ogr
        import gdal
        has_gdalpy = True
    except ImportError:
        has_gdalpy = False

try:
    import arcpy
    has_arcpy = True
except ImportError:
    has_arcpy = False

import geomods

_version = '0.1.0'

usage = '''cudem.py ({}): Process and generate Digital Elevation Models

usage: cudem.py [ -hvAIER [ args ] ]

Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -I, --datalsit\tThe input datalist.
  -E, --increment\tThe desired cell-size in native units.
      note: use at least 7 digits precision with WGS84 units
  -A, --algoritm\tThe desired gridding algorithm
\t\t\tCurrent options: 'mbgrid', 'gmt-surface'

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

 Examples:
 % cudem.py -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27
 % cudem.py --datalist input.datalist --increment 0.000277777 --region input_tiles_ply.shp

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''.format(_version)

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
    return out

def dem_interpolation_uncertainty(dem_surface):
    status = 0
    sub_count = 0

    tw = geomods.clis.prog_bar('preparing for interpolation uncertainty')._end(status)

    this_region = dem_surface.region
    status = dem_surface.proximity()
    
    for sub_region in this_region.split(4):
        sub_count += 1
            
        ## Extract XYZ data for sub-region and randomly sample
        status = this_surf.grd2xyz(this_surf.dem['dem-grd'], 
                                   '{}.xyz'.format(this_surf.oname), 
                                   region = sub_region.buffer(10 * float(iinc)), 
                                   mask = this_surf.dem['num-grd'])
        
        sub_xyz = np.loadtxt(this_surf.dem['xyz'], delimiter = ' ')
                
        np.random.shuffle(sub_xyz)
        sx_len = len(sub_xyz)
        sx_len_pct = int(sx_len * .5)

        np.savetxt('sub_{}d_rest.xyz'.format(sub_count), sub_xyz[sx_len_pct:], '%f', ' ')
                
        sub_datalist =  geomods.datalists.datalist('sub_{}.datalist'.format(sub_count), sub_region)
        sub_datalist._append_datafile('sub_{}_rest.xyz'.format(sub_count), 168, 1)
        sub_datalist._reset()

        ## Generate the random-sample DEM            
        sub_surf = dem(sub_datalist, sub_region, iinc)
        sub_surf.mbgrid()
            
        status = sub_surf.proximity()
            
        ## Query the random-sample DEM with the withheld data
        sub_xyd = geomods.gdalfun.query(sub_xyz[:sx_len_pct], sub_surf.dem['tif'], 'xyd')

        ## Query the random-sample proximity grid with the withheld data
        sub_dp = geomods.gdalfun.query(sub_xyd, sub_surf.dem['prox'], 'zg')
        #print np.max(sub_dp[:,1])
            
        ## Cleanup
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
class dem:
    def __init__(self, idatalist, iregion, iinc = '0.000277777', oname = None, callback = lambda: False):
        self.status = 0
        self.stop = callback
        self.dem = {}
        self.inc = float(iinc)
        self.datalist = idatalist
        self.region = iregion
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.max_prox = self.max_num = None
        self.oname = oname
        if self.oname is None: 
            self.oname = os.path.basename(self.datalist._path).split('.')[0]

    ## 'Validate' the DEM (check if all files in `dem` exist)
    def _valid_p(self):
        for key in self.dem:
            if not os.path.exists(self.dem[key]):
                return(False)
        return(True)

    ## Remove keys to DEM filenames if they don't exist
    def _rectify(self):
        for key in self.dem:
            if not os.path.exists(self.dem[key]):
                del self.dem[key]
                
    ## Set a key/value in the `dem` dictionary
    def _set_dem(self, key, value):
        self.dem[key] = value

    ## Print out the key/value pairs from the `dem` dictionary
    def _print_dem(self):
        for key in self.dem:
            print(key, self.dem[key])

    ## Remove file associated with key
    def _dem_remove(self, key):
        try: 
            os.remove(self.dem[key])
            return 0
        except: return -1

    ## Generic Functions
    def grd2tif(self, src_grd):
        '''Convert the grd file to tif using GMT'''

        if os.path.exists(src_grd):
            grd2tif_cmd = ('gmt grdconvert {} {}.tif=gd+n-9999:GTiff -V\
            '.format(src_grd, os.path.basename(src_grd).split('.')[0]))
            
            out, self.status = geomods.clis.run_cmd(grd2tif_cmd, False, True)
        else: self.status = -1

        return self.status

    def grdinfo(self, src_grd):
        '''Return an info list of `src_grd`'''

        if not os.path.exists(src_grd):
            self.status = -1
            return self.status
        else:
            grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
            out, self.status = geomods.clis.run_cmd(grdinfo_cmd, False, True)
            return out.split()

    def grdcut(self, src_grd, src_region, dst_grd):
        '''Cut `src_grd` to `src_region` '''

        if os.path.exists(src_grd):
            self.cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))

    def grd2xyz(self, src_grd, dst_xyz, region = None, mask = None):
        '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region'''

        if mask:
            self.grdmask_cmd = ('gmt grdmath -V {} NOT {} = tmp.grd'.format(mask, src_grd))

            out, self.status = geomods.clis.run_cmd(self.grdmask_cmd, False, True)
            if self.status == 0: src_grd = 'tmp.grd'

        if region and region._valid:
            self.grd2xyz_cmd = ('gmt grd2xyz -V {} {} -s > {}'.format(src_grd, region.gmt, dst_xyz))

        else: self.grd2xyz_cmd = ('gmt grd2xyz {} -V -s > {}'.format(igrid, dst_xyz))
            
        out, self.status = geomods.clis.run_cmd(self.grd2xyz_cmd, False, True)

        if self.status == 0:
            if mask:
                if os.path.exists('tmp.grd'):
                    os.remove('tmp.grd')

        return self.status
        
    ## DEM Et Cetra
    def surface(self):
        '''Generate a DEM with GMT surface'''

        dem_cmd = ('cat {} | gmt blockmean {} -I{} -r -V | gmt surface -V {} -I{} -G{}_p.grd -T.35 -Z1.2 -r -Lud -Lld\
        '.format(self.datalist._echo_datafiles(' '), self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname))

        dem_cmd1 = ('gmt grdcut -V {}_p.grd -G{}.grd {}\
        '.format(self.oname, self.oname, self.dist_region.gmt))

        out, self.status = geomods.clis.run_cmd(dem_cmd, False, True)
        if self.status == 0:
            out, self.status = geomods.clis.run_cmd(dem_cmd1, False, True)
            
            if self.status == 0:

                if os.path.exists('{}s_p.grd'.format(self.oname)):
                    os.remove('{}_p.grd'.format(self.oname))

                self.dem['dem-grd'] = ('{}.grd'.format(self.oname))

                self.grd2tif(self.dem['dem-grd'])
                
                self.dem['dem'] = ('{}.tif'.format(self.oname))
                
        return self.status

    def num(self):
        '''Generate a num and num-msk grid with GMT'''

        num_cmd0 = ('cat {} | gmt xyz2grd -V {} -I{} -r -G{}_num.grd -An\
        '.format(self.datalist._echo_datafiles(' '), self.dist_region.gmt, self.inc, self.oname))

        num_cmd1 = ('gmt grdmath -V {}_num.grd 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
        '.format(self.oname, self.oname))

        out, self.status = geomods.clis.run_cmd(num_cmd0, False, True)
        if self.status == 0:
            out, self.status = geomods.clis.run_cmd(num_cmd1, False, True)

            if self.status == 0:
                self.grd2tif('{}_num.grd'.format(self.oname))

                self.dem['num'] = '{}_num.tif'.format(self.oname)
                self.dem['num-grd'] = '{}_num.grd'.format(self.oname)
                self.dem['num-msk'] = '{}_num_msk.tif'.format(self.oname)

                self.max_num = int(self.grdinfo(self.dem['num'])[6])

        return self.status   

    def mbgrid(self, extras = False):
        '''Generate a DEM and num grid with MBSystem'''

        mbgrid_cmd = ('mbgrid -I{} {} -E{}/{}/degrees! -O{} -A2 -G3 -F1 -N -C10/3 -S0 -X0.1 -T35 -M\
        '.format(self.datalist._path, self.dist_region.gmt, self.inc, self.inc, self.oname))

        out, self.status = geomods.clis.run_cmd(mbgrid_cmd, False, True)

        if self.status == 0:
        
            self.dem['dem-grd'] = '{}.grd'.format(self.oname)
            self.dem['num-grd'] = '{}_num.grd'.format(self.oname)

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

            self.max_num = int(self.grdinfo(self.dem['num-grd'])[6])

            num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}_num_msk.tif=gd+n-9999:GTiff\
            '.format(self.dem['num'], self.oname))

            out, self.status = geomods.clis.run_cmd(num_msk_cmd, False, True)

            if self.status == 0:
                self.dem['num-msk'] = ('{}_num_msk.tif'.format(self.oname))

        return self.status

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

            out, self.status = geomods.clis.run_cmd(slope_cmd0, False, True)

            if self.status == 0:
                out, self.status = geomods.clis.run_cmd(slope_cmd1, False, True)

                if os.path.exists('{}_pslp.grd'.format(self.oname)):
                    os.remove('{}_pslp.grd'.format(self.oname))

                self.dem['slope'] = ('{}_slp.tif'.format(self.oname))

        return self.status

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

        return self.status

    def spatial_metadata(self):
        '''Geneate spatial metadata'''

        ## initialize output vector
        ## fields incl: 'name' 'title' 'agency' 'date' 'type' 'res' 'hdatum' 'vdatum' 'url'

        dst_vec = '{}_sm_ply.shp'.format(self.oname)
        dst_layername = os.path.basename(dst_vec).split('.')[0]

        ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
        layer = ds.CreateLayer('{}'.format(os.path.basename(dst_vec).split('.')[0]), None, ogr.wkbMultiPolygon)
        defn = layer.GetLayerDefn()

        layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('Agency', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('Type', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('Resolution', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('HDatum', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('VDatum', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('URL', ogr.OFTString))

        for feature in layer:
            layer.SetFeature(feature)
    
        for dl in self.datalist.datalist:
            print dl
            this_datalist = geomods.datalists.datalist(dl[0], self.dist_region)
            this_dem = dem(this_datalist, self.region, str(self.inc))
            print this_dem.oname

            this_dem.num()

            src_ds = gdal.Open(this_dem.dem['num-msk'])
            srcband = src_ds.GetRasterBand(1)

            ## polygonize num mask.

            tm_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_dem.oname))
            tmp_layer = ds.CreateLayer('{}_poly'.format(this_dem.oname), None, ogr.wkbMultiPolygon)

            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
            
            result = gdal.Polygonize(srcband, None, tmp_layer, 0, [], callback = None)
            src_ds = srcband = None

            multi = ogr.Geometry(ogr.wkbMultiPolygon)

            for i in tmp_layer:
                tmp_layer.SetFeature(i)
                if i.GetField('DN') == 0:
                    tmp_layer.DeleteFeature(i.GetFID())
                elif i.GetField('DN') == 1:
                    i.geometry().CloseRings()
                    wkt = i.geometry().ExportToWkt()
                    multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                    tmp_layer.DeleteFeature(i.GetFID())

            #ds.ExecuteSQL('repack {}'.format(this_dem.oname))

            union = multi.UnionCascaded()

            out_feat = ogr.Feature(defn)
            out_feat.SetGeometry(union)
            out_feat.SetField('Name', this_dem.oname)
            out_feat.SetField('Agency', dl[3])
            out_feat.SetField('Date', int(dl[4]))
            out_feat.SetField('Type', dl[5])
            out_feat.SetField('Resolution', dl[6])
            out_feat.SetField('HDatum', dl[7])
            out_feat.SetField('VDatum', dl[8])
            out_feat.SetField('URL', dl[9])

            layer.CreateFeature(out_feat)
            
            tmp_ds = tmp_layer = None
            
            ogr.GetDriverByName('ESRI Shapefile').DeleteDataSource('{}_poly.shp'.format(this_dem.oname))

            #ds.ExecuteSQL('repack tmp1')
        ds = src_ds = layer = None

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

    dem_alg = 'mbgrid'
    dem_algs = ['mbgrid', 
                'gmt-surface',
                'spatial-metadata']

    argv = sys.argv
        
    ## parse command line arguments.
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

        elif arg == '--alogrithm' or arg == '-A':
            dem_alg = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-A':
            dem_alg = str(arg[2:])

        elif arg == '--uncertainty' or arg == '-u':
            do_unc = True

        elif arg[0] == '-':
            print(usage)
            sys.exit(0)

        elif idatalist is None:
            idatalist = arg

        elif iregion is None:
            iregion = arg

        else:
            print(usage)
            sys.exit(1)

        i = i + 1

    if iregion is None:
        print(usage)
        sys.exit(1)
        #iregion = raw_input('input region: ')

    #if idatalist is None:
        #idatalist = raw_input('input datalist: ')
        #print(usage)
        #sys.exit(1)

    if dem_alg not in dem_algs:
        dem_alg = 'gmt'

    ## check platform
    tw = geomods.clis.prog_bar('checking platform')
    platform = sys.platform
    tw.pm = 'checking platform - {}'.format(platform)
    tw._end(status)

    ## check for installed software
    tw = geomods.clis.prog_bar('checking for GMT')
    if geomods.clis.cmd_exists('gmt'): 
        gmt_vers, status = geomods.clis.run_cmd('gmt --version')
        tw.opm = 'checking for GMT - {}'.format(gmt_vers.rstrip())
    else: status = -1
    tw._end(status)

    tw = geomods.clis.prog_bar('checking for MBSystem')
    if geomods.clis.cmd_exists('mbgrid'): 
        mbs_vers, status = geomods.clis.run_cmd('mbgrid -version')
        tw.opm = 'checking for MBSystem - {}'.format(mbs_vers.split('\n')[3].rstrip().split()[2])
    else: status = -1
    tw._end(status)

    tw = geomods.clis.prog_bar('checking for GDAL command-line')
    if geomods.cmd_exists('gdal-config'): 
        gdal_vers, status = geomods.clis.run_cmd('gdal-config --version')
        tw.opm = 'checking for GDAL command-line - {}'.format(gdal_vers.rstrip())
    else: status = -1
    tw._end(status)

    tw = geomods.clis.prog_bar('checking for GDAL python bindings')
    if has_gdalpy: 
        status = 0
        gdal_vers = gdal.__version__
        tw.opm = 'checking for GDAL python bindings - {}'.format(gdal_vers)
    else: status = -1
    tw._end(status)

    if platform != 'linux2':
        tw = geomods.clis.prog_bar('checking for arcpy python bindings')
        if has_arcpy: 
            status = 0
        else: status = -1
        tw._end(status)

    ## process input region(s)
    tw = geomods.clis.prog_bar('processing region(s)')
    try: 
        these_regions = [geomods.regions.region(iregion)]
    except:
        if os.path.exists(iregion):
            _poly = ogr.Open(iregion)
            _player = _poly.GetLayer(0)
            for pf in _player:
                _pgeom = pf.GetGeometryRef()
                these_regions.append(geomods.regions.region('/'.join(map(str, _pgeom.GetEnvelope()))))

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
            status = -1
    tw.opm = 'processing {} region(s)'.format(len(these_regions))
    tw._end(status)

    if status == -1:
        print(usage)
        sys.exit(1)

    for this_region in these_regions:
        
        ## process input datalist
        this_datalist = geomods.datalists.datalist(idatalist, this_region)
        if not this_datalist._valid: break

        ## ==============================================
        ##
        ## Process Digital Elevation Model
        ## using `dem_alg` gridding function
        ## generate DEM and Num grid using full region
        ##
        ## ==============================================
        tw = geomods.clis.prog_bar("generating Digital Elevation Model using {}".format(dem_alg))

        this_surf = dem(this_datalist, this_region, iinc, callback = lambda: stop_threads)

        ## MBSystem mbgrid
        if dem_alg == 'mbgrid':
            try:
                dem_t = threading.Thread(target = this_surf.mbgrid, args = ())
                dem_t.start()
                while True:
                    time.sleep(5)
                    tw._update()
                    if not dem_t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): this_surf.status = -1

        ## GMT surface
        elif dem_alg == 'gmt-surface':
            try:
                dem_t = threading.Thread(target = this_surf.surface, args = ())
                dem_nt = threading.Thread(target = this_surf.num, args = ())
                dem_nt.start()
                dem_t.start()
                while True:
                    time.sleep(5)
                    tw._update()
                    if not dem_t.is_alive() and not dem_nt.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): this_surf.status = -1

        ## Validate output DEM
        tw = geomods.clis.prog_bar('validating DEM')

        this_surf._rectify()
        if this_surf._valid_p():
            if 'dem' not in this_surf.dem.keys() or 'num' not in this_surf.dem.keys():
                status = -1
        else:
            status = -1
            
        tw._end(status)

        ## Spatial Metadata
        if dem_alg == 'spatial-metadata':
            this_surf.spatial_metadata()

        tw._end(this_surf.status)

        ## Analyze the uncertainty
        if do_unc: dem_interpolation_uncertainty(this_surf)        
            
if __name__ == '__main__':
    main()

### END
