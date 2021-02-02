### waffles.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
## see/set `_waffles_grid_info` dictionary to run a grid, or run waffles_config()
##
## Current DEM modules:
## surface (GMT), triangulate (GMT/GDAL), nearest (GMT/GDAL), mbgrid (MBSYSTEM), num (waffles), average (GDAL), invdst (GDAL), linear (GDAL), spat-meta, vdatum, coastline
##
## optionally, clip, filter, buffer the resulting DEM.
##
## find data to grid with GEOMODS' fetch.py or use fetch modules as datalist entries in waffles.
##
## a datalist '*.datalist' file should be formatted as in MBSystem:
## ~path ~format ~weight ~metadata,list ~etc
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
## a format of 200 represents GDAL data
## a format of 300 represents LAS/LAZ data <not implemented>
## a format of 400 represents a FETCHES module - e.g. `nos:datatype=bag`
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
## or via `datalists -i ~datalist.datalist`
##
## GDAL/LIDAR/FETCHES data don't need inf files; but they may be generated anyway.
##
### TODO:
## Add remove/replace module
## Add source uncertainty to uncertainty module
## -B for 'breakline' (densify line/exract nodes/add to datalist)
## update filters to do in series, add new ones, such as spike filter...
##
### Code:
import sys
import os
import time
import copy
import shutil

try:
    import Queue as queue
except: import queue as queue
import threading

## ==============================================
## import gdal, etc.
## ==============================================
import numpy as np
import math
import json
import gdal
import ogr

## ==============================================
## import geomods
## ==============================================
from geomods import fetches
from geomods import regions
from geomods import utils
from geomods import datalists
from geomods import gmtfun
from geomods import gdalfun
from geomods import xyzfun

_version = '0.6.5'

## ==============================================
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean', etc.
##
## Requires MBSystem, GMT, GDAL and VDatum for full functionality
## ==============================================
_waffles_grid_info = {
    'datalist': None,
    'data': [],
    'region': None,
    'inc': 1,
    'name': 'waffles_dem',
    'node': 'pixel',
    'fmt': 'GTiff',
    'extend': 0,
    'extend_proc': 20,
    'weights': None,
    'z_region': None,
    'w_region': None,
    'fltr': [],
    'sample': None,
    'clip': None,
    'chunk': None,
    'epsg': 4326,
    'mod': 'help',
    'mod_args': (),
    'verbose': False,
    'archive': False,
    'spat': False,
    'mask': False,
    'unc': False,
    'overwrite': True,
    'gc': utils.config_check()
}

## ==============================================
## the default waffles config dictionary.
## lambda returns dictionary with default waffles
## ==============================================
#waffles_config = lambda: copy.deepcopy(_waffles_grid_info)
waffles_config_copy = lambda wg: copy.deepcopy(wg)

def waffles_config(datalist = None, data = [], region = None, inc = None, name = 'waffles_dem',
                   node = 'pixel', fmt = 'GTiff', extend = 0, extend_proc = 20, weights = None,
                   z_region = None, w_region = None, fltr = [], sample = None, clip = None, chunk = None, epsg = 4326,
                   mod = 'surface', mod_args = (), verbose = False, archive = False, spat = False, mask = False,
                   unc = False, gc = None, overwrite = True):
    wg = waffles_config_copy(_waffles_grid_info)
    wg['datalist'] = datalist
    wg['data'] = data
    wg['region'] = region
    wg['inc'] = gmtfun.gmt_inc2inc(str(inc))
    wg['name'] = name
    wg['node'] = node
    wg['fmt'] = 'GTiff'
    wg['extend'] = utils.int_or(extend, 0)
    wg['extend_proc'] = utils.int_or(extend_proc, 20)
    wg['weights'] = weights
    wg['z_region'] = z_region
    wg['w_region'] = w_region
    wg['sample'] = gmtfun.gmt_inc2inc(str(sample))
    wg['clip'] = clip
    wg['chunk'] = chunk
    wg['epsg'] = utils.int_or(epsg, 4326)
    wg['mod'] = mod
    wg['mod_args'] = mod_args
    wg['verbose'] = verbose
    wg['archive'] = archive
    wg['spat'] = spat
    wg['mask'] = mask
    wg['unc'] = unc
    wg['overwrite'] = overwrite
    wg['gc'] = utils.config_check()

    wg['fltr'] = fltr
    
    if wg['data'] is None:
        if wg['datalist'] is not None:
            wg['data'] = [x[0] for x in datalist2py(wg['datalist'])]
        else: wg['data'] = None
    
    if wg['datalist'] is None and len(wg['data']) > 0:
        wg['datalist'] = datalists.datalist_major(wg['data'], region = wg['region'], major = '{}_major.datalist'.format(wg['name']))
        
    if not os.path.exists(wg['datalist']):
        wg['datalist'] = datalists.datalist_major(wg['data'], region = wg['region'], major = '{}_major.datalist'.format(wg['name']))
        
    #if wg['mod'].lower() != 'vdatum' and wg['mod'].lower() != 'coastline':
    if _waffles_modules[wg['mod']][3]:
        if wg['datalist'] is None:
            utils.echo_error_msg('invalid datalist/s entry')
            return(None)

    if wg['region'] is None or not regions.region_valid_p(wg['region']):
        #if wg['mod'].lower() == 'data':
        #    wg['region'] = [-180, 180, -90, 90]
        #else:
        utils.echo_error_msg('invalid region {}'.format(wg['region']))
        return(None)

    if wg['inc'] is None: wg['inc'] = (wg['region'][1] - wg['region'][0]) / 500
    
    return(wg)

## ==============================================
## The waffles modules
## { 'module-name': [module-lambda, module description, raster/vector, requires-datalist], ... }
## the module lambda should point to a function that takes at least
## the waffles config as its first option (e.g. def mod(wg))
## ==============================================
_waffles_modules = {
    'surface': [lambda args: waffles_gmt_surface(**args), '''SPLINE DEM via GMT surface
    Generate a DEM using GMT's surface command

    < surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
     :tension=[0-1] - Spline tension.''', 'raster', True],
    'triangulate': [lambda args: waffles_gmt_triangulate(**args), '''TRIANGULATION DEM via GMT triangulate
    Generate a DEM using GMT's triangulate command

    < triangulate >''', 'raster', True],
    'cudem': [lambda args: waffles_cudem(**args), '''Generate a CUDEM Bathy/Topo DEM <beta>''', 'raster', True],
    'update': [lambda args: waffles_update_dem(**args), '''Update a CUDEM DEM with data from datalist <beta>''', 'raster', True],
    'nearest': [lambda args: waffles_nearneighbor(**args), '''NEAREST NEIGHBOR DEM via GMT or gdal_grid
    Generate a DEM using GDAL's gdal_grid command or GMT's nearest command

    < nearest:radius=6s:use_gdal=False >
     :radius=[value] - Nearest Neighbor search radius
     :use_gdal=[True/False] - use gdal grid nearest algorithm''', 'raster', True],
    'num': [lambda args: waffles_num(**args), '''Uninterpolated DEM populated by <mode>.
    Generate an uninterpolated DEM using <mode> option.
    Using mode of 'A<mode>' uses GMT's xyz2grd command, see gmt xyz2grd --help for more info.

    < num:mode=n >
     :mode=[key] - specify mode of grid population: k (mask), m (mean) or n (num)''', 'raster', True],
    'vdatum': [lambda args: waffles_vdatum(**args), '''VDATUM transformation grid
    Generate a VDatum based vertical datum transformation grid.

    < vdatum:ivert=navd88:overt=mhw:region=3:jar=None >
     :ivert=[vdatum] - Input VDatum vertical datum.
     :overt=[vdatum] - Output VDatum vertical datum.
     :region=[0-10] - VDatum region (3 is CONUS).
     :jar=[/path/to/vdatum.jar] - VDatum jar path - (auto-locates by default)''', 'raster', False],
    'mbgrid': [lambda args: waffles_mbgrid(**args), '''Weighted SPLINE DEM via mbgrid
    Generate a DEM using MBSystem's mbgrid command.

    < mbgrid:tension=35:dist=10/3:use_datalists=False >
     :tension=[0-100] - Spline tension.
     :dist=[value] - MBgrid -C switch (distance to fill nodata with spline)
     :use_datalists=[True/False] - use waffles built-in datalists''', 'raster', True],
    'invdst': [lambda args: waffles_invdst(**args), '''INVERSE DISTANCE DEM via gdal_grid
    Generate a DEM using GDAL's gdal_grid command.

    < invdst:power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1 >''', 'raster', True],
    'average': [lambda args: waffles_moving_average(**args), '''Moving AVERAGE DEM via gdal_grid
    Generate a DEM using GDAL's gdal_grid command.

    < average:radius1=0.01:radius2=0.01 >''', 'raster', True],
    'linear': [lambda args: waffles_linear(**args), '''LINEAR DEM via gdal_grid
    Generate a DEM using GDAL's gdal_grid command.

    < linear:radius=0.01 >''', 'raster', True],
    'spat-meta': [lambda args: waffles_spatial_metadata(**args), '''generate SPATIAL-METADATA
    Generate spatial metadata based on the data in the datalist.

    < spat-meta >''', 'vector', True],
    'coastline': [lambda args: waffles_coastline(**args), '''generate a coastline (landmask)
    Generate a land/sea mask (coastline) based on various datasets.

    < coastline >''', 'vector', False],
    'uncertainty': [lambda args: waffles_interpolation_uncertainty(**args), '''generate DEM UNCERTAINTY
    Calculate the interpolation uncertainty in relation to distance to nearest measurement

    < uncertainty:mod=surface:dem=None:msk=None:prox=None:slp=None:sims=10:chnk_lvl=6 >''', 'raster', False],
}

## ==============================================
## module descriptors (used in cli help)
## ==============================================
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['{:14}{}\n'.format(key, x[key][1]) for key in x]) + '\n'
_waffles_module_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

## ==============================================
## the grid-node region
## ==============================================
waffles_grid_node_region = lambda wg: regions.region_buffer(wg['region'], wg['inc'] * .5)

## ==============================================
## the "proc-region" regions.region_buffer(wg['region'], (wg['inc'] * 20) + (wg['inc'] * wg['extend']))
## ==============================================
waffles_proc_region = lambda wg: regions.region_buffer(wg['region'], (wg['inc'] * wg['extend_proc']) + (wg['inc'] * wg['extend']))
waffles_coast_region = lambda wg: regions.region_buffer(waffles_proc_region(wg), (wg['inc'] * 200))
waffles_proc_str = lambda wg: regions.region_format(waffles_proc_region(wg), 'gmt')
waffles_proc_bbox = lambda wg: regions.region_format(waffles_proc_region(wg), 'bbox')
waffles_proc_ul_lr = lambda wg: regions.region_format(waffles_proc_region(wg), 'ul_lr')

## ==============================================
## the "dist-region" regions.region_buffer(wg['region'], (wg['inc'] * wg['extend']))
## ==============================================
waffles_dist_region = lambda wg: regions.region_buffer(wg['region'], (wg['inc'] * wg['extend']) if wg['sample'] is None else (wg['sample'] * wg['extend']))
waffles_dist_ul_lr = lambda wg: regions.region_format(waffles_dist_region(wg), 'ul_lr')

## ==============================================
## the datalist dump function, to use in utils.run_cmd()
## ==============================================
waffles_dl_func = lambda wg: lambda p: waffles_dump_datalist(wg, dst_port = p)

## ==============================================
## grid registration string for use in GTM programs
## ==============================================
waffles_gmt_reg_str = lambda wg: '-r' if wg['node'] == 'pixel' else ''

## ==============================================
## the 'long-name' used from prefix
## ==============================================
waffles_append_fn = lambda bn, region, inc: '{}{}_{}_{}v1'.format(bn, utils.inc2str_inc(inc), regions.region_format(region, 'fn'), utils.this_year())

def waffles_wg_valid_p(wg = _waffles_grid_info):
    '''return True if wg_config appears valid'''

    try:
        if wg['data'] is None: return(False)
        if wg['mod'] is None: return(False)
        else: return(True)
    except: return(False)

def waffles_dlp_hooks(wg = _waffles_grid_info):
    '''the deafult datalist pass hooks.
    checks for region intersection and possibly z and/or weight range.

    returns a list of hooks'''
    
    region = waffles_proc_region(wg)
    dlp_hooks = datalists.datalist_default_hooks()
    
    if region is not None: dlp_hooks.append(lambda e: regions.regions_intersect_ogr_p(region, datalists.inf_entry(e)))
    if wg['z_region'] is not None:
        dlp_hooks.append(lambda e: regions.z_region_pass(datalists.inf_entry(e), upper_limit = wg['z_region'][1], lower_limit = wg['z_region'][0]))

    if wg['w_region'] is not None:
        dlp_hooks.append(lambda e: regions.z_pass(e[2], upper_limit = wg['w_region'][1], lower_limit = wg['w_region'][0]))

    return(dlp_hooks)

def waffles_yield_datalist(wg = _waffles_grid_info):
    '''recurse the datalist and do things to it.
    
    yield the xyz line data'''
    
    region = waffles_proc_region(wg)
    dlp_hooks = waffles_dlp_hooks(wg)
    
    dly = datalists.datalist_yield_xyz(wg['datalist'], pass_h = dlp_hooks, wt = 1 if wg['weights'] else None, region = region, archive = wg['archive'], verbose = wg['verbose'], z_region = wg['z_region'], epsg = wg['epsg'])
    if wg['mask']: dly = gdalfun.gdal_xyz_mask(dly, '{}_msk.tif'.format(wg['name']), region, wg['inc'], dst_format = wg['fmt'])
    for xyz in dly: yield(xyz)
    
    if wg['archive']:
        a_dl = os.path.join('archive', '{}.datalist'.format(wg['name']))

        for dir_, _, files in os.walk('archive'):
            for f in files:
                if '.datalist' in f:
                    rel_dir = os.path.relpath(dir_, 'archive')
                    rel_file = os.path.join(rel_dir, f)
                    datalists.datalist_append_entry([rel_file, -1, 1], a_dl)

def waffles_dump_datalist(wg = _waffles_grid_info, dst_port = sys.stdout):
    '''dump the xyz data from datalist to `dst_port`.'''

    for xyz in waffles_yield_datalist(wg):
        xyzfun.xyz_line(xyz, dst_port, True)
        
def waffles_datalist_list(wg = _waffles_grid_info):
    '''list the datalist entries in the given region'''
        
    for this_entry in datalist(wg['datalist'], wt = 1, pass_h = waffles_dlp_hooks(wg)):
        print(' '.join([','.join(x) if i == 3 else os.path.abspath(str(x)) if i == 0 else str(x) for i,x in enumerate(this_entry[:-1])]))
        
def waffles_datalists(wg = _waffles_grid_info, dump = False, echo = False, infos = False, recurse = True):
    '''dump the xyz data from datalist and generate a data mask while doing it.'''

    if echo: waffles_datalist_list(wg)
    if infos: print(datalists.datalist_inf(wg['datalist'], inf_file = True))
    if dump:
        recurse = True
        pass_func = lambda xyz: xyzfun.xyz_line(xyz, sys.stdout, True)
    else: pass_func = lambda xyz: None

    if recurse:
        for xyz in waffles_yield_datalist(wg): pass_func(xyz)
    return(0,0)

## ==============================================
## Waffles Filter; filters a DEM
## optionally split by split_value (z) and only filter lt
## 1 = gdalfun.gdalblur (default filter_val = 10)
## 2 = gmt grdfilter (default filter_val = 1s)
## 3 = ?
## ==============================================
def waffles_filter(src_gdal, dst_gdal, fltr = 1, fltr_val = None, split_value = None):
    '''filter `src_gdal` using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero) using a split_value of 0.

    return 0 for success or -1 for failure'''
    
    if os.path.exists(src_gdal):

        ## ==============================================
        ## split the dem by `split_value`
        ## ==============================================
        if split_value is not None:
            dem_u, dem_l = gdal_split(src_gdal, split_value)
        else: dem_l = src_gdal
        
        ## ==============================================
        ## filter the possibly split DEM
        ## ==============================================
        if int(fltr) == 1: out, status = gdalfun.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        elif int(fltr) == 2: out, status = gmtfun.gmt_grdfilter(dem_l, 'tmp_fltr.tif=gd+n-9999:GTiff', dist = fltr_val if fltr_val is not None else '1s', verbose = True)
        else: out, status = gdalfun.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)

        ## ==============================================
        ## merge the split filtered and non filtered DEMs
        ## ==============================================
        if split_value is not None:
            ds = gdal.Open(src_gdal)
            ds_config = gdalfun.gdal_gather_infos(ds)
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
                    gdalfun.gdal_write(ds_arr, 'merged.tif', ds_config)
                    l_ds = None
                    utils.remove_glob(dem_l)
                u_ds = None
                utils.remove_glob(dem_u)
            os.rename('merged.tif', 'tmp_fltr.tif')
        os.rename('tmp_fltr.tif', dst_gdal)
        return(0)
    else: return(-1)

## ==============================================
## Waffles Spatial Metadata
## Polygonize each datalist entry into vector data source
## ==============================================
def waffles_spatial_metadata(wg):
    '''generate spatial-metadata

    returns [output vector fn, status]'''

    dst_vector = '{}_sm.shp'.format(wg['name'])
    dst_layer = '{}_sm'.format(wg['name'])
    v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
    t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
    utils.remove_glob('{}.*'.format(dst_layer))
    gdalfun.gdal_prj_file('{}.prj'.format(dst_layer), wg['epsg'])
    
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer('{}'.format(dst_layer), None, ogr.wkbMultiPolygon)
        [layer.CreateField(ogr.FieldDefn('{}'.format(f), t_fields[i])) for i, f in enumerate(v_fields)]
        [layer.SetFeature(feature) for feature in layer]
    else: layer = None
    defn = layer.GetLayerDefn()

    for this_entry in datalists.datalist(wg['datalist'], wt = 1 if wg['weights'] else None, pass_h = waffles_dlp_hooks(wg), yield_dl_entry = True, verbose = wg['verbose']):
        if this_entry[1] == -1 or this_entry[-1] == wg['datalist'].split('.')[0]:
            utils.echo_msg(this_entry[0])

            defn = None if layer is None else layer.GetLayerDefn()            
            twg = waffles_config_copy(wg)
            twg['datalist'] = this_entry[0]
            twg['name'] = '{}_{}_msk'.format(os.path.basename(this_entry[0]).split('.')[0].split(':')[0], regions.region_format(twg['region'], 'fn'))
            if twg['inc'] < gmtfun.gmt_inc2inc('.3333333s'):
                twg['inc'] = gmtfun.gmt_inc2inc('.3333333s')
                twg['extend'] = twg['extend'] / 3
            twg['verbose'] = False
            twg = waffles_config(**twg)

            if len(this_entry[3]) == 8:
                o_v_fields = this_entry[3]
            else: o_v_fields = [twg['name'], 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

            out, status = waffles_num(twg, mode='k')
            print(out, status)
            if status == 0:
                waffles_gdal_md(twg)
                ng = '{}.tif'.format(twg['name'])
                if gdalfun.gdal_infos(ng, scan = True)['zr'][1] == 1:
                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(twg['name']))
                    if tmp_ds is not None:
                        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(twg['name']), None, ogr.wkbMultiPolygon)
                        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
                        gdalfun.gdal_polygonize(ng, tmp_layer, verbose = twg['verbose'])

                        if len(tmp_layer) > 1:
                            if defn is None: defn = tmp_layer.GetLayerDefn()
                            out_feat = gdalfun.gdal_ogr_mask_union(tmp_layer, 'DN', defn)
                            [out_feat.SetField(f, o_v_fields[i]) for i, f in enumerate(v_fields)]
                            layer.CreateFeature(out_feat)
                    tmp_ds = None
                    utils.remove_glob('{}_poly.*'.format(twg['name']))
                utils.remove_glob(ng)
    ds = None
    return(dst_vector, 0)
    
## ==============================================
## Waffles CUDEM generation module
## ==============================================
def waffles_cudem(wg = _waffles_grid_info, coastline = None):
    '''generate bathy/topo DEM suitable for CUDEM project'''
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the CUDEM module')
        return(None, -1)

    b_wg = waffles_config_copy(wg)
    ul = -0.1

    ## ==============================================
    ## generate coastline
    ## ==============================================
    if coastline is None:
        utils.echo_msg('----------------------------')
        utils.echo_msg('generating coastline')
        utils.echo_msg('----------------------------')
        
        c_wg = waffles_config_copy(wg)
        c_wg['mod'] = 'coastline'
        c_wg['mod_args'] = ()
        c_wg['name'] = 'tmp_coast'
        
        waffles_run(c_wg)

        utils.run_cmd('gdal_polygonize.py {}.tif {}.shp'.format(c_wg['name'], c_wg['name']), verbose = True)
        utils.run_cmd('ogr2ogr -dialect SQLite -sql "select * from {} where DN=0" {}_c.shp {}.shp'.format(c_wg['name'], c_wg['name'], c_wg['name']), verbose = True)

        coastline = '{}_c.shp'.format(c_wg['name'])

    out, status = utils.run_cmd('coastline2xyz.sh -I {} -O {}_coast.xyz -Z {} -W {} -E {} -S {} -N {}'\
                                .format(coastline, b_wg['name'], 0, waffles_coast_region(b_wg)[0], waffles_coast_region(b_wg)[1], waffles_coast_region(b_wg)[2], waffles_coast_region(b_wg)[3]), verbose = True)
    coast_xyz = '{}_coast.xyz'.format(b_wg['name'])
    coast_region = datalists.inf_entry([coast_xyz, 168, 1])

    ## ==============================================
    ## generate the bathy-surface
    ## using 'surface' with upper_limit of -0.1
    ## at 1 arc-second spacing
    ## ==============================================
    utils.echo_msg('----------------------------')
    utils.echo_msg('generating bathy surface')
    utils.echo_msg('----------------------------')
    
    b_wg['z_region'] = [None, ul + .5]
    b_wg['name'] = 'bathy_{}'.format(wg['name'])
    b_wg['spat'] = False
    b_wg['fltr'] = None
    b_wg['datalist'] = None
    b_wg['data'].append('{} 168 .1'.format(coast_xyz))
    b_wg['mod'] = 'surface'
    b_wg['mod_args'] = ('upper_limit={}'.format(ul),)
    #b_wg['mod'] = 'triangulate'
    b_wg['sample'] = wg['inc']
    b_wg['inc'] = wg['inc'] * 3
    b_wg['name'] = 'bathy_{}'.format(wg['name'])
    b_wg['clip'] = '{}:invert=True'.format(coastline)
    b_wg['extend_proc'] = 40
    b_wg['mask'] = False
    
    b_wg = waffles_config(**b_wg)
    
    bathy_surf = waffles_run(b_wg)
    utils.remove_glob(coast_xyz)
    
    ## ==============================================
    ## append the bathy-surface to the datalist and
    ## generate final DEM using 'surface'
    ## ==============================================
    utils.echo_msg('----------------------------')
    utils.echo_msg('generating integrated bathy-topo surface')
    utils.echo_msg('----------------------------')

    wg['datalist'] = None
    wg['data'].append('{} 200 .5'.format(bathy_surf))
    wg['w_region'] = [.4, None]
    wg = waffles_config(**wg)

    return(waffles_gmt_surface(wg))

## ==============================================
## Waffles Update module <beta>
## ==============================================
def waffles_update_dem(wg = _waffles_grid_info, dem = None):
    ## grid datalist with nearneighbor
    nn_wg = waffles_config_copy(wg)
    nn_wg['mod'] = 'surface'
    nn_wg['mod_args'] = ()
    nn_wg['name'] = 'surf_{}'.format(wg['name'])
    nn_wg['spat'] = True
    nn_wg['mask'] = False
    nn_wg['fltr'] = None
    nn_wg['sample'] = wg['inc']
    nn_wg['inc'] = gmtfun.gmt_inc2inc('1s')
    nn_wg['clip'] = '{}_sm.shp:invert=True'.format(wg['name'])

    nn_dem = waffles_run(nn_wg)
    
    ## mask nearneighbor grid
    dst_msk = 'msk_{}.tif'.format(nn_wg['name'])
    gmtfun.gmt_num_msk(nn_dem, '{}=gd:GTiff'.format(dst_msk), verbose = wg['verbose'])
    gdalfun.gdal_set_nodata('msk_{}.tif'.format(nn_wg['name']), -9999)
    msk_inf = gdalfun.gdal_infos(dst_msk, scan = False)
    print(msk_inf)
    
    dem_inf = gdalfun.gdal_infos(dem, scan = False)
    print(dem_inf)
    nd = dem_inf['ndv']
    
    #dst_msk1 = 'msk_{}.tif'.format(nn_wg['name'])
    #out, status = utils.run_cmd('gmt grdmath {} {} MUL = {}=gd:GTiff'.format(dst_msk, nd, dst_msk1), verbose = True)
    #gdal_set_nodata('msk_{}.tif'.format(nn_wg['name']), -9999)
        
    ### clip DEM to nearneighbor mask
    #msk_dem = 'msk_{}'.format(dem)
    ##out, status = utils.run_cmd('gmt grdmath -N {} {} SUM = {}=gd:GTiff'.format(dem, dst_msk1, msk_dem), verbose = True)
    #out, status = utils.run_cmd('gdal_calc.py -A {} -B {} --calc A+B --outfile {}'.format(dem, dst_msk1, msk_dem), verbose = True)
    
    ## add clipped DEM to datalist
    f_wg = waffles_config_copy(wg)
    f_wg['mod'] = 'surface'
    f_wg['mod_args'] = ()
    f_wg['datalist'] = None
    #f_wg['data'] =[]
    f_wg['data'].append('{} 200 10'.format(nn_dem))
    f_wg['data'].append('{} 200 .25'.format(dem)) #msk_dem))
    
    f_wg = waffles_config(**f_wg)
    
    ## grid datalist with surface
    #dem_update = waffles_run(f_wg)
    dem_update = waffles_gmt_surface(f_wg)
    return(0, 0)
    
## ==============================================
## Waffles MBGrid module
## ==============================================
def waffles_mbgrid(wg = _waffles_grid_info, dist = '10/3', tension = 35, use_datalists = False):
    '''Generate a DEM with MBSystem's mbgrid program.
    if `use_datalists` is True, will parse the datalist through
    waffles instead of mbsystem.

    Note: without `use_datalists` as True, only mbsystem supported data formats may be used.'''
    
    if wg['gc']['MBGRID'] is None:
        utils.echo_error_msg('MBSystem must be installed to use the MBGRID module')
        return(None, -1)
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the MBGRID module')
        return(None, -1)

    if use_datalists:
        #datalist_archive(wg, arch_dir = '.mb_tmp_datalist', verbose = True)
        archive = wg['archive']
        wg['archive'] = True
        for xyz in waffles_yield_datalist(wg): pass
        wg['datalist'] = datalists.datalist_major(['archive/{}.datalist'.format(wg['name'])])
        wg['archive'] = archive

    region = waffles_proc_region(wg)
    region = regions.region_buffer(region, wg['inc'] * -.5)
    xsize, ysize, gt = gdalfun.gdal_region2gt(waffles_proc_region(wg), wg['inc'])
    
    if len(dist.split('/')) == 1: dist = dist + '/2'
    mbgrid_cmd = ('mbgrid -I{} {} -D{}/{} -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} \
    '.format(wg['datalist'], regions.region_format(region, 'gmt'), xsize, ysize, wg['name'], dist, tension, '-M' if wg['mask'] else ''))
    for out in utils.yield_cmd(mbgrid_cmd, verbose = wg['verbose']): sys.stderr.write('{}'.format(out))
    #out, status = utils.run_cmd(mbgrid_cmd, verbose = wg['verbose'])

    utils.remove_glob('*.cmd')
    utils.remove_glob('*.mb-1')
    gmtfun.gmt_grd2gdal('{}.grd'.format(wg['name']))
    utils.remove_glob('{}.grd'.format(wg['name']))
    if use_datalists and not wg['archive']: shutil.rmtree('archive')
    if wg['mask']:
        utils.remove_glob('*_sd.grd')
        num_grd = '{}_num.grd'.format(wg['name'])
        dst_msk = '{}_msk.tif=gd+n-9999:GTiff'.format(wg['name'])
        out, status = gmtfun.gmt_num_msk(num_grd, dst_msk, verbose = wg['verbose'])
        utils.remove_glob(num_grd)
    if not use_datalists:
        if wg['spat'] or wg['archive']:
            for xyz in waffles_yield_datalist(wg): pass
    return('{}.tif'.format(wg['name']), 0)

## ==============================================
## Waffles GMT surface module
## ==============================================
def waffles_gmt_surface(wg = _waffles_grid_info, tension = .35, relaxation = 1.2,
                        lower_limit = 'd', upper_limit = 'd'):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the SURFACE module')
        return(None, -1)

    dem_surf_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} -r\
'.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), \
             wg['inc'], wg['name'], tension, relaxation, lower_limit, upper_limit))
    return(utils.run_cmd(dem_surf_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

## ==============================================
## Waffles GMT triangulate module
## ==============================================
def waffles_gmt_triangulate(wg = _waffles_grid_info):
    '''generate a DEM with GMT surface'''
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the TRIANGULATE module')
        return(None, -1)
    dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt triangulate {} -I{:.10f} -V -G{}.tif=gd+n-9999:GTiff -r\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), wg['inc'], wg['name']))
    return(utils.run_cmd(dem_tri_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))

## ==============================================
## Waffles nearest neighbor module
## GMT if available else GDAL
## ==============================================
def waffles_nearneighbor(wg = _waffles_grid_info, radius = None, use_gdal = False):
    '''genearte a DEM with GMT nearneighbor or gdal_grid nearest'''
    
    radius = wg['inc'] * 3 if radius is None else gmtfun.gmt_inc2inc(radius)
    if wg['gc']['GMT'] is not None and not use_gdal:
        dem_nn_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), \
                 wg['inc'], radius, wg['name']))
        return(utils.run_cmd(dem_nn_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))
    else: return(waffles_gdal_grid(wg, 'nearest:radius1={}:radius2={}:nodata=-9999'.format(radius, radius)))

## ==============================================
## Waffles 'NUM grid' module
## Uninterpolated grdi from data;
## num methods include: mask, mean, num, landmask and any gmt grd2xyz -A option.
## ==============================================
def waffles_num(wg = _waffles_grid_info, mode = 'n'):
    '''Generate an uninterpolated num grid.
    mode of `k` generates a mask grid
    mode of `m` generates a mean grid
    mode of `n` generates a num grid
    mode of `w` generates a landmask grid
    mode of `A<mode> ` generates grid using GMT xyz2grd'''

    if mode[0] == 'A':
        if wg['gc']['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the SURFACE module')
            return(None, -1)

        dem_xyz2grd_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt xyz2grd -{} -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', mode, waffles_proc_str(wg), wg['inc'], wg['name']))
        return(utils.run_cmd(dem_xyz2grd_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg)))
    else:
        dly = waffles_yield_datalist(wg)
        if wg['weights']: dly = xyzfun.xyz_block(dly, waffles_proc_region(wg), wg['inc'], weights = True)
        return(gdalfun.gdal_xyz2gdal(dly, '{}.tif'.format(wg['name']), waffles_proc_region(wg), wg['inc'], dst_format = wg['fmt'], mode = mode, verbose = wg['verbose']))

## ==============================================
## Waffles GDAL_GRID module
## ==============================================
def waffles_gdal_grid(wg = _waffles_grid_info, alg_str = 'linear:radius=1'):
    '''run gdal grid using alg_str
    parse the data through xyzfun.xyz_block to get weighted mean before
    building the GDAL dataset to pass into gdal_grid'''

    region = waffles_proc_region(wg)
    dly = xyzfun.xyz_block(waffles_yield_datalist(wg), region, wg['inc'], weights = False if wg['weights'] is None else True)
    ds = gdalfun.xyz2gdal_ds(dly, '{}'.format(wg['name']))
    if ds.GetLayer().GetFeatureCount() == 0: return(-1,-1)
    xcount, ycount, dst_gt = gdalfun.gdal_region2gt(region, wg['inc'])
    gd_opts = gdal.GridOptions(outputType = gdal.GDT_Float32, noData = -9999, format = 'GTiff', \
                               width = xcount, height = ycount, algorithm = alg_str, callback = gdalfun._gdal_progress if wg['verbose'] else None, \
                               outputBounds = [region[0], region[3], region[1], region[2]])
    gdal.Grid('{}.tif'.format(wg['name']), ds, options = gd_opts)
    ds = None
    gdalfun.gdal_set_nodata('{}.tif'.format(wg['name']), -9999)
    return('{}.tif'.format(wg['name']), 0)

## ==============================================
## Waffles GDAL_GRID invdist module
## ==============================================
def waffles_invdst(wg = _waffles_grid_info, power = 2.0, smoothing = 0.0,
                   radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999):
    '''Generate an inverse distance grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmtfun.gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmtfun.gmt_inc2inc(radius2)
    gg_mod = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
                             .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_moving_average(wg = _waffles_grid_info, radius1 = None, radius2 = None,
                           angle = 0.0, min_points = 0, nodata = -9999):
    '''generate a moving average grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmtfun.gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmtfun.gmt_inc2inc(radius2)
    gg_mod = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'.format(radius1, radius2, angle, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_linear(wg = _waffles_grid_info, radius = None, nodata = -9999):
    '''generate a moving average grid with GDAL'''
    
    radius1 = wg['inc'] * 2 if radius is None else gmtfun.gmt_inc2inc(radius)
    gg_mod = 'linear:radius={}:nodata={}'.format(radius, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles VDATUM 'conversion grid' module
## ==============================================
def waffles_vdatum(wg = _waffles_grid_info, ivert = 'navd88', overt = 'mhw',
                   region = '3', jar = None):
    '''generate a 'conversion-grid' with vdatum.
    output will be the differences (surfaced) between 
    `ivert` and `overt` for the region'''
    
    vc = _vd_config
    if jar is None:
        vc['jar'] = vdatum_locate_jar()[0]
    else: vc['jar'] = jar
    vc['ivert'] = ivert
    vc['overt'] = overt
    vc['region'] = region

    gdalfun.gdal_null('empty.tif', waffles_proc_region(wg), 0.00083333, nodata = 0)
    with open('empty.xyz', 'w') as mt_xyz:
        gdalfun.gdal_dump_entry(['empty.tif', 200, 1], dst_port = mt_xyz)
    
    run_vdatum('empty.xyz', vc)
    
    if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
        with open('result/empty.xyz') as infile:
            empty_infos = xyz_inf(infile)
        print(empty_infos)

        ll = 'd' if empty_infos[4] < 0 else '0'
        lu = 'd' if empty_infos[5] > 0 else '0'
        wg['data'] = ['result/empty.xyz']
        wg['spat'] = False
        wg = waffles_config(**wg)
        out, status = waffles_gmt_surface(wg, tension = 0, upper_limit = lu, lower_limit = ll)
    else:
        out = None
        status = -1
        
    utils.remove_glob('empty.*')
    utils.remove_glob('result/*')
    utils.remove_glob('.mjr.datalist')
    os.removedirs('result')
    return(out, status)

## ==============================================
## Waffles Interpolation Uncertainty module
## make own module!
## ==============================================
_unc_config = {
    'wg': waffles_config_copy(_waffles_grid_info),
    'mod': 'surface',
    'mod_args': (),
    'dem': None,
    'msk': None,
    'prox': None,
    'slp': None,
    'percentile': 95,
    'zones': ['low-dens', 'mid-dens', 'high-dens'],
    'sims': 10,
    'chnk_lvl': 6,
}

waffles_unc_config = lambda: copy.deepcopy(_unc_config)
waffles_unc_config_copy = lambda uc: copy.deepcopy(uc)

def waffles_interpolation_uncertainty(wg = _waffles_grid_info, mod = 'surface', mod_args = (), \
                                      dem = None, msk = None, prox = None, slp = None, \
                                      percentile = 95, zones = ['low-dens', 'mid-dens', 'high-dens', 'low-slp', 'mid-slp', 'high-slp'], \
                                      sims = 10, chnk_lvl = 6):
    '''calculate the interpolation uncertainty
    - as related to distance to nearest measurement.

    returns [[err, dist] ...]'''

    s_dp = None
    s_ds = None

    #zones = ['low-dens-low-slp', 'low-dens-mid-slp', 'low-dens-high-slp', 'mid-dens-low-slp', 'mid-dens-mid-slp', 'mid-dens-high-slp', 'high-dens-low-slp', 'high-dens-mid-slp', 'high-dens-high-slp']
    zones = ['bathy', 'bathy-topo', 'topo']
    
    ## ==============================================
    ## set the module and input grids.
    ## generate any necessary grids that don't exists/
    ## aren't specified...
    ## ==============================================
    if mod not in _waffles_modules.keys():
        utils.echo_error_msg('invalid module name `{}`; reverting to `surface`'.format(mod))
        mod = 'surface'
        mod_args = ()

    wg['mod'] = mod
    wg['mod_args'] = mod_args
    
    utils.echo_msg('running INTERPOLATION uncertainty module using {}...'.format(wg['mod']))
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)

    if dem is None:
        dem = 'dem_{}.tif'.format(wg['mod'])
        tmp_wg = waffles_config_copy(wg)
        tmp_wg['datalist'] = None
        tmp_wg['name'] = 'dem_{}'.format(wg['mod'])
        
        if msk is None:
            msk = 'dem_{}_msk.tif'.format(tmp_wg['mod'])
            tmp_wg['mask'] = True
        else: tmp_wg['mask'] = False
        tmp_wg = waffles_config(**tmp_wg)
        waffles_run(tmp_wg)

    if msk is None:
        if msk is None: msk = '{}_msk.tif'.format(wg['name'])
        tmp_wg = waffles_config_copy(wg)
        tmp_wg['name'] = '{}_msk'.format(wg['name'])
        tmp_wg['mod'] = 'num'
        tmp_wg['mod_args'] = ('mode=k',)
        tmp_wg['datalist'] = None
        tmp_wg = waffles_config(**tmp_wg)
        waffles_run(tmp_wg)
        
    if prox is None:
        if prox is None: prox = '{}_prox.tif'.format(wg['name'])
        gdalfun.gdal_proximity(msk, prox)
        if wg['epsg'] is not None: gdalfun.gdal_set_epsg(prox, wg['epsg'])
    # if slp is None:
    #      if slp is None: slp = '{}_slp.tif'.format(wg['name'])
    #      gdalfun.gdal_slope(dem, slp, 1)
    #      if wg['epsg'] is not None: gdalfun.gdal_set_epsg(slp, wg['epsg'])

    ## ==============================================
    ## region analysis
    ## ==============================================
    region_info = {}

    ## ==============================================
    ## mask analysis
    ## ==============================================
    num_sum, g_max, num_perc = gdalfun.gdal_mask_analysis(mask = msk)

    ## ==============================================
    ## proximity analysis
    ## ==============================================
    prox_percentile = gdalfun.gdal_percentile(prox, percentile)
    prox_perc_33 = gdalfun.gdal_percentile(prox, 33)
    prox_perc_66 = gdalfun.gdal_percentile(prox, 66)
    prox_perc_100 = gdalfun.gdal_percentile(prox, 100)

    # ## ==============================================
    # ## slope analysis
    # ## ==============================================
    # slp_percentile = gdalfun.gdal_percentile(slp, percentile)
    # slp_perc_33 = gdalfun.gdal_percentile(slp, 25)
    # slp_perc_66 = gdalfun.gdal_percentile(slp, 75)
    # slp_perc_100 = gdalfun.gdal_percentile(slp, 100)
    slp_percentile = 0
    
    region_info[wg['name']] = [wg['region'], g_max, num_sum, num_perc, prox_percentile, slp_percentile] 
    for x in region_info.keys():
        utils.echo_msg('region: {}: {}'.format(x, region_info[x]))

    ## ==============================================
    ## chunk region into sub regions
    ## ==============================================
    #chnk_lvl = 2
    if int(region_info[wg['name']][3] > region_info[wg['name']][4]):
        chnk_lvl = int(region_info[wg['name']][3] / region_info[wg['name']][4])
    else: chnk_lvl = int(region_info[wg['name']][4] / region_info[wg['name']][3])
    chnk_lvl = 4
    #chnk_inc = int(region_info[wg['name']][4] * int(chnk_lvl))
    utils.echo_msg('chunking region into sub-regions using chunk level {}...'.format(chnk_lvl))
    #print(chnk_inc)
    #if chnk_inc
    chnk_inc = int(region_info[wg['name']][4] * int(chnk_lvl))
    #chnk_inc = int((int(region_info[wg['name']][1]) / math.sqrt(g_max)) / int(region_info[wg['name']][4] * int(chnk_lvl)))
    sub_regions = regions.region_chunk(wg['region'], wg['inc'], chnk_inc)
    utils.echo_msg('chunked region into {} sub-regions @ {}x{} cells.'.format(len(sub_regions), chnk_inc, chnk_inc))

    ## ==============================================
    ## sub-region analysis
    ## ==============================================
    utils.echo_msg('analyzing {} sub-regions...'.format(len(sub_regions)))
    sub_zones = {}
    dem_ds = gdal.Open(dem)
    msk_ds = gdal.Open(msk)
    prox_ds = gdal.Open(prox)
    #slp_ds = gdal.Open(slp)
    for sc, sub_region in enumerate(sub_regions):
        utils.echo_msg_inline('analyzing sub-regions [{}]'.format(sc))
        ## s_sum, s_g_max, s_perc = gdalfun.gdal_mask_analysis('tmp_msk.tif')
        s_sum, s_g_max, s_perc = gdalfun.gdal_mask_analysis2(msk_ds, region = sub_region)
        p_perc = gdalfun.gdal_prox_analysis2(prox_ds, region = sub_region)
        #slp_perc = gdalfun.gdal_prox_analysis2(slp_ds, region = sub_region)
        slp_perc = 0
        s_dc = gdalfun.gdal_gather_infos(dem_ds, region = sub_region, scan = True)
        if p_perc < prox_perc_33 or abs(p_perc - prox_perc_33) < 0.01: zone = zones[2]
        elif p_perc < prox_perc_66 or abs(p_perc - prox_perc_66) < 0.01: zone = zones[1]
        else: zone = zones[0]
        #if slp_perc < slp_perc_33 or abs(slp_perc - slp_perc_33) < 0.01: zone = zones[6]
        #elif slp_perc < slp_perc_66 or abs(slp_perc - slp_perc_66) < 0.01: zone = zones[7]
        #else: zone = zones[8]
        #if slp_perc < slp_perc_33 or abs(slp_perc - slp_perc_33) < 0.01: zone = zones[3]
        #elif slp_perc < slp_perc_66 or abs(slp_perc - slp_perc_66) < 0.01: zone = zones[4]
        #else: zone = zones[5]
        #if slp_perc < slp_perc_33 or abs(slp_perc - slp_perc_33) < 0.01: zone = zones[0]
        #elif slp_perc < slp_perc_66 or abs(slp_perc - slp_perc_66) < 0.01: zone = zones[1]
        #else: zone = zones[2]

        sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, slp_perc, s_dc['zr'][0], s_dc['zr'][1], zone]
    #dem_ds = msk_ds = prox_ds = slp_ds = None
    dem_ds = msk_ds = prox_ds = slp_ds = None
    utils.echo_msg_inline('analyzing sub-regions [OK]\n')
    #print(sub_zones)
    s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
    s_5perc = np.percentile(s_dens, 5)
    s_dens = None
    utils.echo_msg('Sampling density for region is: {:.16f}'.format(s_5perc))

    ## ==============================================
    ## zone analysis / generate training regions
    ## ==============================================
    trainers = []
    # bathy_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[0]]
    # bathy_topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[1]]
    # topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[2]]
    # flat_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][9] == zones[3]]
    # mid_slp_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][9] == zones[4]]
    # steep_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][9] == zones[5]]
    t_perc = 50
    s_perc = 50

    # #for z, tile_set in enumerate([bathy_tiles, bathy_topo_tiles, topo_tiles]):
    # for z, this_zone in enumerate(zones):
    #     tile_set = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[z]]
    #     if len(tile_set) > 0:
    #         t_dens = np.array([x[4] for x in tile_set])
    #         t_dens = t_dens[t_dens != 0]
    #         if len(t_dens) == 0: continue
    #         t_50perc = np.percentile(t_dens, t_perc)
    #         while t_50perc == 100:
    #             t_perc = t_perc - 5
    #             t_50perc = np.percentile(t_dens, t_perc)
    #     else: continue
    #     utils.echo_msg('Maximum gap for {} tiles: {} cells'.format(zones[z].upper(), int(t_50perc)))
    #     t_trainers = [x for x in tile_set if x[4] < t_50perc or abs(x[4] - t_50perc) < 0.01]
    #     utils.echo_msg('possible {} training zones: {}'.format(zones[z].upper(), len(t_trainers)))
    #     trainers.append(t_trainers)
        
    #for z, tile_set in enumerate([bathy, bathy-topo, topo]):
    samp_percs = {}
    for z, this_zone in enumerate(zones):
        tile_set = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[z]]
        if len(tile_set) > 0:
            dens = np.array([x[3] for x in tile_set])
            d_50perc = np.percentile(dens, 50)
            d_5perc = np.percentile(dens, 5)
            t_dens = np.array([x[4] for x in tile_set])
            t_dens = t_dens[t_dens != 0]
            if len(t_dens) == 0: continue
            t_50perc = np.percentile(t_dens, t_perc)
            while t_50perc == 100:
                t_perc = t_perc - 5
                t_50perc = np.percentile(t_dens, t_perc)
        else: continue
        utils.echo_msg('Maximum proximity, sampling for {} tiles: {}, {}'.format(zones[z].upper(), t_50perc, d_50perc))
        samp_percs[this_zone] = d_5perc
        t_trainers = [x for x in tile_set if x[4] > t_50perc or abs(x[4] - t_50perc) < 0.01]
        #t_trainers = [x for x in tile_set if x[3] < d_50perc or abs(x[3] - d_50perc) < 0.01]
        utils.echo_msg('possible {} training zones: {}'.format(zones[z].upper(), len(t_trainers)))
        trainers.append(t_trainers)
        
    utils.echo_msg('sorting training tiles by distance...')
    trains = regions.regions_sort(trainers, verbose = False)
    utils.echo_msg('sorted sub-regions into {} training tiles.'.format(len([x for s in trains for x in s])))
    utils.echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))

    ## ==============================================
    ## split-sample simulations and error calculations
    ## sims = max-simulations
    ## ==============================================
    utils.echo_msg('performing SPLIT-SAMPLE simulations...')
    #utils.echo_msg('simulation\terrors\tproximity\tp_diff\tslope\ts_diff')
    utils.echo_msg('simulation\terrors\tproximity-coeff\tp_diff')
    last_ec_d = None
    #last_ec_s = None
    sim = 0
    while True:
        #trains = regions.regions_sort(trainers, verbose = False)
        utils.echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of {} [{:3}%]'.format(sim, sims, 0))
        status = 0
        sim += 1
        for z, train in enumerate(trains):
            train_h = train[:25]
            ss_samp = s_5perc
            
            ## ==============================================
            ## perform split-sample analysis on each training region.
            ## ==============================================
            for n, sub_region in enumerate(train_h):
                ss_samp = s_5perc
                #ss_samp = samp_percs[sub_region[8]]
                perc = int(float(n+(len(train_h) * z))/(len(train_h)*len(trains)) * 100)
                utils.echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of {} [{:3}%]'.format(sim, sims, perc))
                #utils.echo_msg_inline('performing SPLIT-SAMPLE ({} - {}) simulation {} out of {} [{:3}%]'.format(ss_samp, sub_region[8], sim, sims, perc))
                #print(sub_region)
                this_region = sub_region[0]
                if sub_region[3] < ss_samp: ss_samp = None

                ## ==============================================
                ## extract the xyz data for the region from the DEM
                ## ==============================================
                o_xyz = '{}_{}.xyz'.format(wg['name'], n)
                ds = gdal.Open(dem)
                with open(o_xyz, 'w') as o_fh:
                    for xyz in gdalfun.gdal_parse(ds, srcwin = gdalfun.gdal_srcwin(ds, regions.region_buffer(this_region, (20 * wg['inc']))), mask = msk):
                        xyzfun.xyz_line(xyz, o_fh)
                ds = None

                if os.stat(o_xyz).st_size != 0:
                    ## ==============================================
                    ## split the xyz data to inner/outer; outer is
                    ## the data buffer, inner will be randomly sampled
                    ## ==============================================
                    s_inner, s_outer = gmtfun.gmt_select_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = False) #verbose = uc['wg']['verbose'])
                    if os.stat(s_inner).st_size != 0:
                        sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                    else: sub_xyz = []
                    ss_len = len(sub_xyz)
                    if ss_samp is not None:
                        sx_cnt = int(sub_region[1] * (ss_samp / 100.)) + 1
                    else: sx_cnt = 1

                    sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                    np.random.shuffle(sub_xyz)
                    np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                    ## ==============================================
                    ## generate the random-sample DEM
                    ## ==============================================
                    wc = waffles_config(name = 'sub_{}'.format(n),
                                        data = [s_outer, sub_xyz_head],
                                        region = this_region,
                                        inc = wg['inc'],
                                        mod = wg['mod'],
                                        verbose = False,
                                        mod_args = wg['mod_args'],
                                        mask = True,
                                        overwrite = True,
                                        epsg = wg['epsg'],
                                        )
                    sub_dem = waffles_run(wc)
                    sub_msk = '{}_msk.tif'.format(wc['name'])

                    if os.path.exists(sub_dem) and os.path.exists(sub_msk):
                        ## ==============================================
                        ## generate the random-sample data PROX and SLOPE
                        ## ==============================================        
                        sub_prox = '{}_prox.tif'.format(wc['name'])
                        gdalfun.gdal_proximity(sub_msk, sub_prox)

                        # sub_slp = '{}_slp.tif'.format(wc['name'])
                        # gdalfun.gdal_slope(sub_dem, sub_slp)
                        
                        ## ==============================================
                        ## Calculate the random-sample errors
                        ## ==============================================
                        sub_xyd = gdalfun.gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                        #sub_dp = gdalfun.gdal_query(sub_xyd, sub_prox, 'zg')
                        sub_dp = gdalfun.gdal_query(sub_xyd, sub_prox, 'xyzg')
                        # sub_ds = gdalfun.gdal_query(sub_dp, slp, 'g')
                        
                        # if len(sub_dp) > 0:
                        #     if sub_dp.shape[0] == sub_ds.shape[0]:
                        #         sub_dp = np.append(sub_dp, sub_ds, 1)
                        #     else: sub_dp = []
                    else: sub_dp = None
                    utils.remove_glob(sub_xyz_head)

                    if s_dp is not None: 
                        if sub_dp is not None and len(sub_dp) > 0:
                            try:
                                s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                            except: s_dp = sub_dp
                    else: s_dp = sub_dp
                utils.remove_glob(o_xyz)
                utils.remove_glob('sub_{}*'.format(n))

        if len(s_dp) > 0:
            d_max = region_info[wg['name']][4]
            #s_max = region_info[wg['name']][5]
            s_dp = s_dp[s_dp[:,3] < d_max,:]
            # s_dp = s_dp[s_dp[:,4] < s_max,:]
            # s_dp = s_dp[s_dp[:,3] > 0,:]
            
            prox_err = s_dp[:,[2,3]]
            #slp_err = s_dp[:,[2,4]]
            if last_ec_d is None: last_ec_d = [0, 0.1, 0.2]
            ec_d = utils.err2coeff(prox_err[:50000000], dst_name = wg['name'] + '_prox', xa = 'distance')
            #ec_s = utils.err2coeff(slp_err[:50000000], dst_name = wg['name'] + '_slp', xa = 'slope')
            #utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_d[2] - ec_d[1], ec_s, ec_s[2] - ec_s[1]))
            utils.echo_msg('{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_d[2] - ec_d[1]))
            
            # if last_ec_d is None: last_ec_d = ec_d
            # if last_ec_s is None: last_ec_s = ec_s
            if ec_d[2] < 0.001: continue
            if sim >= int(sims): break
            if len(s_dp) >= int(region_info[wg['name']][1] / 10): break
            # if abs(ec_d[2] - ec_d[1]) > abs(last_ec_d[2] - last_ec_d[1]):
            #     ec_d = last_ec_d
            #     ec_s = last_ec_s
            #     break
            last_ec_d = ec_d
            #last_ec_s = ec_s
        else: utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None, None, None))

    ## ==============================================
    ## Save/Output results
    ## apply error coefficient to full proximity grid
    ## ==============================================

    utils.echo_msg('applying coefficient to proximity grid')
    ## USE numpy/gdal instead
    #utils.run_cmd('gdal_calc.py -A {} --outfile {}_prox_unc.tif --calc "{}+({}*(A**{}))"'.format(prox, wg['name'], 0, ec_d[1], ec_d[2]), verbose = True)
    math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_prox_unc.tif=gd+n-9999:GTiff\
    '.format(prox, ec_d[2], ec_d[1], 0, wg['name'])
    utils.run_cmd(math_cmd, verbose = wg['verbose'])
    if wg['epsg'] is not None: status = gdalfun.gdal_set_epsg('{}_prox_unc.tif'.format(wg['name']), epsg = wg['epsg'])
    utils.echo_msg('applied coefficient {} to proximity grid'.format(ec_d))

    want_full_out = False
    
    # if want_full_out:
    #     #np.savetxt('{}_prox.err'.format(wg['name']), prox_err, '%f', ' ')
    #     #np.savetxt('{}_slp.err'.format(wg['name']), slp_err, '%f', ' ')
        
    #     utils.run_cmd('gdal_calc.py -A {} --outfile _tmp_msk.tif --calc "1/A"'.format(msk), verbose = True)
    #     gdalfun.gdal_mask('_tmp_msk.tif', '_tmp_msk_inverted.tif', invert = True)
    #     utils.run_cmd('gdal_calc.py -A {} -B {} --outfile {}_slp_unc.tif --calc "B*({}+({}*(A**{})))"'.format(slp, '_tmp_msk_inverted.tif', wg['name'], 0, ec_s[1], ec_s[2]), verbose = True)
    #     utils.remove_glob('_tmp_msk*')
    #     #math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_slp_unc.tif=gd+n-9999:GTiff\
    #         #'.format(slp, ec_s[2], ec_s[1], 0, wg['name'])
    #     #utils.run_cmd(math_cmd, verbose = wg['verbose'])
    #     if wg['epsg'] is not None: gdalfun.gdal_set_epsg('{}_slp_unc.tif'.format(wg['name']), epsg = wg['epsg'])
    #     utils.echo_msg('applied coefficient {} to slope grid'.format(ec_s))

    #     utils.run_cmd('gdal_calc.py -A {}_prox_unc.tif -B {}_slp_unc.tif --outfile={}_unc.tif --calc="((A*A)+(B*B))**(1/2.0)" --overwrite'.format(wg['name'], wg['name'], wg['name']))
    
    #     unc_out = {
    #     'unc': '{}_unc.tif'.format(wg['name']),
    #     'prox_unc': '{}_prox_unc.tif'.format(wg['name']),
    #     'prox_bf': '{}_prox_bf.png'.format(wg['name']),
    #     'prox_scatter': '{}_prox_scatter.png'.format(wg['name']),
    #     'prox': '{}_prox.tif'.format(wg['name']),
    #     'slp_unc': '{}_slp_unc.tif'.format(wg['name']),
    #     'slp_bf': '{}_slp_bf.png'.format(wg['name']),
    #     'slp_scatter': '{}_slp_scatter.png'.format(wg['name']),
    #     'slp': '{}_slp.tif'.format(wg['name']),
    #     'prox_coeff': ec_d,
    #     'slp_coeff': ec_s,
    #     }
    # else:
    unc_out = [{
        'prox_unc': '{}_prox_unc.tif'.format(wg['name']),
        'prox_bf': '{}_prox_bf.png'.format(wg['name']),
        'prox_scatter': '{}_prox_scatter.png'.format(wg['name']),
    }, ec_d]

    return(unc_out, 0)

## ==============================================
## Waffles Coastline module
## generate a coastline (wet/dry mask)
## ==============================================
def waffles_coastline(wg, want_nhd = True, want_gmrt = False):
    '''Generate a coastline polygon from various sources.'''
    
    w_mask = '{}_w.tif'.format(wg['name'])
    
    if wg['datalist'] is not None:    
        ## ==============================================
        ## wet/dry datalist mask
        ## ==============================================
        dly = waffles_yield_datalist(wg)
        if wg['weights']: dly = xyzfun.xyz_block(dly, waffles_dist_region(wg), wg['inc'], weights = True)
        gdalfun.gdal_xyz2gdal(dly, w_mask, waffles_dist_region(wg), wg['inc'], dst_format = wg['fmt'], mode = 'w', verbose = wg['verbose'])
    else:
        ## ==============================================
        ## burn the region
        ## ==============================================
        gdalfun.gdal_region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        utils.run_cmd('gdal_rasterize -tr {} {} -te {} -burn -9999 -a_nodata -9999 -ot Int32 -co COMPRESS=DEFLATE region_buff.shp {}'\
                .format(wg['inc'], wg['inc'], regions.region_format(waffles_dist_region(wg), 'te'), w_mask), verbose = False)

    ## ==============================================
    ## load the wet/dry mask array
    ## ==============================================
    ds = gdal.Open(w_mask)
    if ds is not None:
        ds_config = gdalfun.gdal_gather_infos(ds)
        region = gdalfun.gdal_gt2region(ds_config)
        dst_gt = ds_config['geoT']        
        coast_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
    else: return(-1, -1)
    utils.remove_glob('{}*'.format(w_mask))

    ## ==============================================
    ## Input coastline shapefile `coastpoly`
    ## ==============================================
    
    ## ==============================================
    ## USGS NHD (HIGH-RES U.S. Only)
    ## ==============================================
    if want_nhd:
        u_mask = '{}_u.tif'.format(wg['name'])
        gdalfun.gdal_region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        utils.run_cmd('gdal_rasterize -tr {} {} -te {} -burn -9999 -a_nodata -9999 -ot Int32 -co COMPRESS=DEFLATE region_buff.shp {}'\
                .format(wg['inc'], wg['inc'], regions.region_format(waffles_dist_region(wg), 'te'), u_mask), verbose = False)
        utils.remove_glob('region_buff.*')

        fl = fetches.fetch_infos['tnm'][0](waffles_proc_region(wg), [], None)
        r = fl.run(ds = 6, formats = 'FileGDB', extent = 'HU-4 Subregion')

        if len(r) > 0:
            r_shp = []
            for result in r:
                if fetches.fetch_file(result[0], os.path.join(result[2], result[1]), verbose = True) == 0:
                    gdb_zip = os.path.join(result[2], result[1])
                    gdb_files = utils.unzip(gdb_zip)
                    gdb, gdb_files = utils.procs_unzip(gdb_zip, ['gdb'])
                    gdb_bn = os.path.basename('.'.join(gdb_zip.split('.')[:-1]))
                    gdb = gdb_bn + '.gdb'

                    utils.run_cmd('ogr2ogr {}_NHDArea.shp {} NHDArea -overwrite'.format(gdb_bn, gdb), verbose = False)
                    r_shp.append('{}_NHDArea.shp'.format(gdb_bn))
                    utils.run_cmd('ogr2ogr {}_NHDPlusBurnWaterBody.shp {} NHDPlusBurnWaterBody -overwrite'.format(gdb_bn, gdb), verbose = False)
                    r_shp.append('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn))
                else: utils.echo_error_msg('unable to fetch {}'.format(result))

            [utils.run_cmd('ogr2ogr -skipfailures -update -append nhdArea_merge.shp {}'.format(shp), verbose = False) for shp in r_shp]
            utils.run_cmd('gdal_rasterize -burn 1 nhdArea_merge.shp {}'\
                    .format(u_mask), verbose = True)
            utils.remove_glob(gdb_zip)
            utils.remove_glob('{}*'.format(gdb_bn))
            utils._clean_zips(gdb_files)
            [utils.remove_glob('{}*'.format(shp[:-3])) for shp in r_shp]
            utils.remove_glob('nhdArea_merge.*')

        ## ==============================================
        ## update wet/dry mask with nhd data
        ## ==============================================
        c_ds = gdal.Open(u_mask)
        for this_xyz in gdalfun.gdal_parse(c_ds):
            xpos, ypos = gdalfun._geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
            try:
                if coast_array[ypos, xpos] == ds_config['ndv']:
                    if this_xyz[2] == 1: coast_array[ypos, xpos] = 0
            except: pass
        c_ds = None            
        utils.remove_glob('{}*'.format(u_mask))
    
    ## ==============================================
    ## GSHHG/GMRT - Global low-res
    ## ==============================================
    g_mask = '{}_g.tif'.format(wg['name'])
    if wg['gc']['GMT'] is not None and not want_gmrt:
        utils.run_cmd('gmt grdlandmask {} -I{} -r -Df -G{}=gd:GTiff -V -N1/0/1/0/1'.format(regions.region_format(waffles_dist_region(wg), 'gmt'), wg['inc'], g_mask), verbose = True)
    else:
        r = fetches.fetch_infos['gmrt'][0](regions.region_buffer(waffles_dist_region(wg), 5, pct = True), [], None).run()
        gmrt_tif = r[0][1]
        if fetches.fetch_file(r[0][0], gmrt_tif, verbose = True) == 0:
            utils.run_cmd('gdalwarp {} {} -tr {} {} -overwrite'.format(gmrt_tif, g_mask, wg['inc'], wg['inc']), verbose = True)
            utils.remove_glob(gmrt_tif)

    ## ==============================================
    ## update wet/dry mask with gsshg/gmrt data
    ## ==============================================
    c_ds = gdal.Open(g_mask)
    for this_xyz in gdalfun.gdal_parse(c_ds):
        xpos, ypos = gdalfun._geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
        try:
            if coast_array[ypos, xpos] == ds_config['ndv']:
                if this_xyz[2] == 1: coast_array[ypos, xpos] = 0
                elif this_xyz[2] == 0: coast_array[ypos, xpos] = 1
        except: pass
    c_ds = None
    utils.remove_glob('{}*'.format(g_mask))

    ## ==============================================
    ## write coast_array to file
    ## 1 = land
    ## 0 = water
    ## ==============================================
    #coast_array[coast_array == 0] = ds_config['ndv']
    gdalfun.gdal_write(coast_array, '{}.tif'.format(wg['name']), ds_config)

    #utils.run_cmd('gdal_polygonize.py -8 {}.tif {}.shp'.format(wg['name'], wg['name']), verbose = True)

    ## ==============================================
    ## convert to vector
    ## ==============================================
    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp_c_{}.shp'.format(wg['name']))
    if tmp_ds is not None:
        tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(wg['name']), None, ogr.wkbMultiPolygon)
        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
        gdalfun.gdal_polygonize('{}.tif'.format(wg['name']), tmp_layer, verbose = wg['verbose'])        
        tmp_ds = None

    utils.run_cmd('ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 order by ST_AREA(geometry) desc limit 8" {}.shp tmp_c_{}.shp'\
                  .format(wg['name'], wg['name'], wg['name']), verbose = False)
        
    return(0, 0)

def waffles_gdal_md(wg, cudem = False):
    '''add metadata to the waffles dem'''
    
    ds = gdal.Open('{}.tif'.format(wg['name']), gdal.GA_Update)
    if ds is not None:
        md = ds.GetMetadata()
        if wg['node'] == 'pixel':
            md['AREA_OR_POINT'] = 'Area'
        else: md['AREA_OR_POINT'] = 'Point'
        md['TIFFTAG_DATETIME'] = '{}'.format(utils.this_date())
        if cudem:
            md['TIFFTAG_COPYRIGHT'] = 'DOC/NOAA/NESDIS/NCEI > National Centers for Environmental Information, NESDIS, NOAA, U.S. Department of Commerce'
            md['TIFFTAG_IMAGEDESCRIPTION'] = 'Topography-Bathymetry; NAVD88'
        ds.SetMetadata(md)
        ds = None
    else: utils.echo_error_msg('failed to set metadata')

def waffles_queue(q):
    while True:
        this_wg = q.get()
        dem = waffles_run(this_wg)

        q.task_done()
    
## ==============================================
## Waffles run waffles module via wg_config()
## ==============================================
def waffles_run(wg = _waffles_grid_info):
    '''generate a DEM using wg dict settings
    see waffles_config() to generate a wg config.

    - runs the waffles module to generate the DEM
    - optionally clips the output to shapefile
    - optionally filters the output
    - optionally resamples the output
    - cuts the output to dist-size
    - reformats the output to final format
    - sets metadata in output

    returns dem-fn'''

    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    no_clobber = False
    
    ## ==============================================
    ## validate and/or set the waffles_config
    ## ==============================================
    if wg is None:
        utils.echo_error_msg('invalid configuration, {}'.format(wg))
        sys.exit(-1)

    args_d = {}
    args_d = utils.args2dict(wg['mod_args'], args_d)
    
    if wg['verbose']:
        #utils.echo_msg(json.dumps(wg, indent=4))
        utils.echo_msg(wg)
        utils.echo_msg('running module {} with {} [{}]...'.format(wg['mod'], wg['mod_args'], args_d))

    dem = '{}.tif'.format(wg['name'])
    vect = True if _waffles_modules[wg['mod']][2] == 'vector' else False
    rast = True if _waffles_modules[wg['mod']][2] == 'raster' else False
    if wg['mask']:
        dem_msk = '{}_msk.tif'.format(wg['name'])
        if os.path.exists(dem_msk) and not wg['overwrite']: no_clobber = True
    if vect:
        if wg['mod'] == 'spat-meta':
            dem_vect = '{}_sm.shp'.format(wg['name'])
        else: dem_vect = '{}.shp'.format(wg['name'])
        if os.path.exists(dem_vect) and not wg['overwrite']: no_clobber = True
    #if rast:
    if os.path.exists(dem) and not wg['overwrite']: no_clobber = True
        
    wg['region'] = waffles_grid_node_region(wg) if wg['node'] == 'grid' else wg['region']

    ## ==============================================
    ## optionally generate the DEM in chunks
    ## skip if wg['overwrite'] is False and file exists
    ## ==============================================
    if not no_clobber:
        if wg['chunk'] is not None:
            xcount, ycount, dst_gt = gdalfun.gdal_region2gt(wg['region'], wg['inc'])
            s_regions = regions.region_chunk(wg['region'], wg['inc'], (xcount/wg['chunk'])+1)

        else: s_regions = [wg['region']]

        chunks = []
        if wg['mask']: chunks_msk = []
        if vect: chunks_vect = []
        for region in s_regions:
            this_wg = waffles_config_copy(wg)
            this_wg['region'] = region
            #this_wg['extend'] = wg['extend'] + 5
            #this_wg['region'] = waffles_grid_node_region(this_wg) if this_wg['node'] == 'grid' else this_wg['region']
            this_wg['name'] = 'chunk_{}'.format(regions.region_format(region, 'fn'))
            this_dem = this_wg['name'] + '.tif'
            chunks.append(this_dem)
            if this_wg['mask']:
                this_dem_msk = this_wg['name'] + '_msk.tif'
                chunks_msk.append(this_dem_msk)
            if vect:
                if this_wg['mod'] == 'spat-meta':
                    chunks_vect.append('{}_sm.shp'.format(this_wg['name']))
                else: chunks_vect.append('{}.shp'.format(this_wg['name']))
            args_d['wg'] = this_wg

            ## ==============================================
            ## gererate the DEM (run the module)
            ## ==============================================
            #try:
            waffles_out, status = _waffles_modules[this_wg['mod']][0](args_d)
            #except KeyboardInterrupt as e:
            #    utils.echo_error_msg('killed by user, {}'.format(e))
            #    sys.exit(-1)
            #except Exception as e:
            #    utils.echo_error_msg('{}'.format(e))
            #    status = -1

            if status != 0: utils.remove_glob(this_dem)
            if not os.path.exists(this_dem): continue
            gdi = gdalfun.gdal_infos(this_dem, scan = True)
            if gdi is not None:
                if np.isnan(gdi['zr'][0]):
                    utils.remove_glob(this_dem)
                    if this_wg['mask']: utils.remove_glob(this_dem_msk)
                    continue
            else: continue

            gdalfun.gdal_set_epsg(this_dem, this_wg['epsg'])
            waffles_gdal_md(this_wg)

            ## ==============================================
            ## optionally filter the DEM
            ## this_wg['fltr'] should hold a list of filters
            ## and optional arguments; e.g. = ['1:10', '2:1s']
            ## ==============================================
            if this_wg['fltr'] is not None:
                for fltr in this_wg['fltr']:
                    fltr_args = {}
                    fltrs = fltr.split(':')
                    fltr_args['fltr'] = fltrs[0]
                    fltr_args['fltr_val'] = gmtfun.gmt_inc2inc(fltrs[1])
                    fltr_args = utils.args2dict(fltrs[2:], fltr_args)        

                    if this_wg['verbose']: utils.echo_msg('filtering {} using filter {}@{}...'.format(this_dem, fltr_args['fltr'], fltr_args['fltr_val']))
                    if fltr_args['fltr'] == 2:
                        if this_wg['gc']['GMT'] is None:
                            continue
                    try:
                        waffles_filter(this_dem, 'tmp_s.tif', **fltr_args)
                        os.rename('tmp_s.tif', this_dem)
                    except TypeError as e: utils.echo_error_msg('{}'.format(e))

            ## ==============================================
            ## optionally resample the DEM 
            ## ==============================================
            if this_wg['sample'] is not None:
                if this_wg['verbose']: utils.echo_msg('resampling {}...'.format(this_dem))
                if this_wg['gc']['GMT'] is not None:
                    gmtfun.gmt_sample_inc(this_dem, inc = this_wg['sample'], verbose = this_wg['verbose'])
                    if this_wg['mask']:
                        if this_wg['verbose']: utils.echo_msg('resampling {}...'.format(this_dem_msk))
                        gmtfun.gmt_sample_inc(this_dem_msk, inc = this_wg['sample'], verbose = this_wg['verbose'])
                else:
                    out, status = utils.run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} tmp.tif\
                    '.format(inc, inc, src_grd, regions.region_format(waffles_proc_region(this_wg)), verbose = verbose))
                    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))

            gdalfun.gdal_set_epsg(this_dem, this_wg['epsg'])

            ## ==============================================
            ## optionally clip the DEM to polygon
            ## ==============================================
            if this_wg['clip'] is not None:
                if this_wg['verbose']: utils.echo_msg('clipping {}...'.format(this_dem))
                clip_args = {}
                cp = this_wg['clip'].split(':')
                clip_args['src_ply'] = cp[0]
                clip_args = utils.args2dict(cp[1:], clip_args)
                print(clip_args)
                gdalfun.gdal_clip(this_dem, **clip_args)
                if this_wg['mask']:
                    if this_wg['verbose']: utils.echo_msg('clipping {}...'.format(this_dem_msk))
                    gdalfun.gdal_clip(this_dem_msk, **clip_args)

            if not os.path.exists(this_dem): continue
            gdi = gdalfun.gdal_infos(this_dem, scan = True)
            if gdi is not None:
                if np.isnan(gdi['zr'][0]):
                    utils.remove_glob(this_dem)
                    if this_wg['mask']: utils.remove_glob(this_dem_msk)
                    continue
            else: continue

            ## ==============================================
            ## cut dem to final size -
            ## region buffered by (inc * extend) or
            ## sample * extend) if sample if specified
            ## ==============================================
            try:
                out = gdalfun.gdal_cut(this_dem, waffles_dist_region(this_wg), 'tmp_cut.tif')
                if out is not None: os.rename('tmp_cut.tif', this_dem)
                if this_wg['mask']:
                    out = gdalfun.gdal_cut(this_dem_msk, waffles_dist_region(this_wg), 'tmp_cut.tif')
                    if out is not None: os.rename('tmp_cut.tif', this_dem_msk)
            except OSError as e:
                utils.remove_glob('tmp_cut.tif')
                utils.echo_error_msg('cut failed, is the dem open somewhere, {}'.format(e)) 

        ## ==============================================
        ## merge the chunks and remove
        ## ==============================================

        if len(chunks) > 1:
            out, status = utils.run_cmd('gdal_merge.py -n -9999 -a_nodata -9999 -ps {} -{} -ul_lr {} -o {} {} -co TILED=YES -co COMPRESS=DEFLATE -co PREDICTOR=3\
            '.format(wg['inc'], wg['inc'], waffles_dist_ul_lr(wg), dem, ' '.join(chunks)), verbose = True)
            [utils.remove_glob(x) for x in chunks]            
        else:
            if os.path.exists(chunks[0]):
                gdalfun.gdal_gdal2gdal(chunks[0], dst_fmt = wg['fmt'], dst_gdal = dem)
                utils.remove_glob(chunks[0])

        if wg['mask']:
            if len(chunks_msk) > 1:
                out, status = utils.run_cmd('gdal_merge.py -n -9999 -a_nodata -9999 -ps {} -{} -ul_lr {} -o {} {}\
                '.format(wg['inc'], wg['inc'], waffles_dist_ul_lr(wg), dem_msk, ' '.join(chunks_msk)), verbose = True)
                [utils.remove_glob(x) for x in chunks_msk]
            else:
                if os.path.exists(chunks_msk[0]):
                    gdalfun.gdal_gdal2gdal(chunks_msk[0], dst_fmt = wg['fmt'], dst_gdal = dem_msk, co = False)
                    utils.remove_glob(chunks_msk[0])

        if vect:
            if len(chunks_vect) > 1:
                out, status = utils.run_cmd('ogrmerge.py {} {}'.format(dem_vect, ' '.join(chunks_vect)))
                [utils.remove_glob('{}*'.format(x.split('.')[0])) for x in chunks_vect]
            else:
                utils.remove_glob('{}*'.format(dem_vect.split('.')[0]))
                out, status = utils.run_cmd('ogr2ogr {} {}'.format(dem_vect, chunks_vect[0]), verbose = True)
                #gdalfun.ogr_clip(chunks_vect[0], dem_vect, clip_region = wg['region'], dn = "ESRI Shapefile")
                #out, status = utils.run_cmd('ogr2ogr -clipsrc {} {} {}'.format(regions.region_format(wg['region'], 'te'), dem_vect, chunks_vect[0]), verbose = True)

                utils.remove_glob('{}*'.format(chunks_vect[0].split('.')[0]))

    if os.path.exists(dem):
        ## ==============================================
        ## convert to final format if not geotiff
        ## ==============================================
        if wg['fmt'] != 'GTiff':
            orig_dem = dem
            if wg['gc']['GMT'] is not None:
                dem = gmtfun.gmt_grd2gdal(orig_dem, wg['fmt'])
            else: dem = gdalfun.gdal_gdal2gdal(orig_dem, wg['fmt'])
            utils.remove_glob(orig_dem)

        ## ==============================================
        ## set the projection and other metadata
        ## ==============================================
        try:
            gdalfun.gdal_set_epsg(dem, wg['epsg'])
            waffles_gdal_md(wg, True if wg['mod'] == 'cudem' else False)
            if wg['mask']: gdalfun.gdal_set_epsg(dem_msk, wg['epsg'])
        except: pass

    ## ==============================================
    ## optionally generate uncertainty grid
    ## ==============================================
    if wg['unc'] and not vect and wg['mod'] != 'uncertainty':
        try:
            if os.path.exists(dem) and os.path.exists(dem_msk):
                utils.echo_msg('generating uncertainty')

                uc = _unc_config
                uc['wg'] = wg
                uc['mod'] = wg['mod']
                uc['mod_args'] = wg['mod_args']
                uc['dem'] = dem
                uc['msk'] = dem_msk
                uc['chnk_lvl'] = 6

                dem_prox = '{}_prox.tif'.format(wg['name'])
                gdalfun.gdal_proximity(dem_msk, dem_prox)
                uc['prox'] = dem_prox

                dem_slp = '{}_slp.tif'.format(wg['name'])
                gdalfun.gdal_slope(dem, dem_slp)
                uc['slp'] = dem_slp

                utils.echo_msg(uc)
                waffles_interpolation_uncertainty(**uc)
        except Exception as e:
            utils.echo_error_msg('failed to calculate uncertainty, {}'.format(e))

    if wg['mod'] == 'uncertainty':
        for out_key in waffles_out[0].keys():
            os.rename(waffles_out[0][out_key], '{}_{}.{}'.format(wg['name'], out_key, waffles_out[0][out_key].split('.')[-1]))
        # os.rename(waffles_out['unc'], '{}_unc.tif'.format(wg['name']))
        # os.rename(waffles_out['prox_unc'], '{}_prox_unc.tif'.format(wg['name']))
        # os.rename(waffles_out['prox_bf'], '{}_prox_bf.png'.format(wg['name']))
        # os.rename(waffles_out['prox_scatter'], '{}_prox_scatter.png'.format(wg['name']))
        # os.rename(waffles_out['prox'], '{}_prox.tif'.format(wg['name']))
        # os.rename(waffles_out['slp_unc'], '{}_slp_unc.tif'.format(wg['name']))
        # os.rename(waffles_out['slp'], '{}_slp.tif'.format(wg['name']))
        # os.rename(waffles_out['slp_bf'], '{}_slp_bf.png'.format(wg['name']))
        # os.rename(waffles_out['slp_scatter'], '{}_slp_scatter.png'.format(wg['name']))

    ## ==============================================
    ## if dem has data, return
    ## ==============================================
    utils.remove_glob('waffles_dem_mjr.datalist')
    utils.remove_glob(wg['datalist'])

    if vect:
        if gdalfun.ogr_empty_p(dem_vect):
            utils.remove_glob('{}.*'.format(dem_vect.split('.')[0]))
            return(None)
        
    return(dem)

## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
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
  -T, --filter\t\tFILTER the output DEM using one or multiple filters. <fltr:fltr_val:split_value=z>
\t\t\tAvailable filters:
\t\t\t1: perform a Gaussian filter at -T1:<factor>.
\t\t\t2: use a Cosine Arch Filter at -T2:<dist(km)> search distance.
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAppend :split_value=<num> to only filter values below z-value <num>.
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry using Gaussian filter
  -Z --z-region\t\tRestrict data processing to records that fall within the z-region
\t\t\tUse '-' to indicate no bounding range; e.g. -Z-/0 will restrict processing to data
\t\t\trecords whose z value is below zero.
  -C, --clip\t\tCLIP the output to the clip polygon -C<clip_ply.shp:invert=False>
  -K, --chunk\t\tProcess the region in CHUNKs. -K<chunk-level>
  -W, --w-region\tRestrict data processing to records that fall within the w-region (weight).
\t\t\tUse '-' to indicate no bounding range; e.g. -W1/- will restrict processing to data
\t\t\trecords whose weight value is at least 1.
  -G, --wg-config\tA waffles config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a waffles_config JSON file using the --config flag.

  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year info to output BASENAME).
  -r, --grid-node\tUse grid-node registration, default is pixel-node
  -w, --weights\t\tUse weights provided in the datalist to weight overlapping data.

  -a, --archive\t\tArchive the datalist to the given region.
  -m, --mask\t\tGenerate a data mask raster.
  -u, --uncert\t\tGenerate an associated uncertainty grid.
  -c, --continue\tDon't clobber existing files.

  --help\t\tPrint the usage text
  --config\t\tSave the waffles config JSON and major datalist
  --modules\t\tDisply the module descriptions and usage
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Datalists and data formats:
  A datalist is a file that contains a number of datalist entries, while an entry is a space-delineated line:
  `/path/to/data format weight data,meta,data`

Supported datalist formats: 
  {}

Modules (see waffles --modules <module-name> for more info):
  {}

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(datalists._known_datalist_fmts_short_desc(), _waffles_module_short_desc(_waffles_modules))

def waffles_cli(argv = sys.argv):
    '''run waffles from command-line
    e.g. `python waffles.py` 
    generates a waffles_config from the command-line options
    and either outputs the or runs the waffles_config
    on each region supplied (multiple regions can be supplied
    by using a vector file as the -R option.)
    See `waffles_cli_usage` for full cli options.'''
    
    wg = waffles_config_copy(_waffles_grid_info)
    wg_user = None
    dls = []
    region = None
    module = None
    want_prefix = False
    want_verbose = False
    want_config = False
    want_threads = False
    status = 0
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
            wg['inc'] = gmtfun.gmt_inc2inc(incs[0])
            if len(incs) > 1: wg['sample'] = gmtfun.gmt_inc2inc(incs[1])
            i = i + 1
        elif arg[:2] == '-E':
            incs = arg[2:].split(':')
            wg['inc'] = gmtfun.gmt_inc2inc(arg[2:].split(':')[0])
            if len(incs) > 1: wg['sample'] = gmtfun.gmt_inc2inc(incs[1])
        elif arg == '--outname' or arg == '-O':
            wg['name'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-O': wg['name'] = arg[2:]
        elif arg == '--format' or arg == '-F':
            wg['fmt'] = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': wg['fmt'] = arg[2:]
        elif arg == '--filter' or arg == '-T':
            wg['fltr'].append(argv[i + 1])
            i += 1
        elif arg[:2] == '-T': wg['fltr'].append(arg[2:])
        elif arg == '--extend' or arg == '-X':
            exts = argv[i + 1].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
            i += 1
        elif arg[:2] == '-X':
            exts = arg[2:].split(':')
            wg['extend'] = utils.int_or(exts[0], 0)
            if len(exts) > 1: wg['extend_proc'] = utils.int_or(exts[1], 10)
        elif arg == '--wg-config' or arg == '-G':
            wg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-G': wg_user = arg[2:]
        elif arg == '--clip' or arg == '-C':
            wg['clip'] = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-C': wg['clip'] = arg[2:]
        elif arg == '--chunk' or arg == '-K':
            wg['chunk'] = utils.int_or(argv[i + 1], None)
            i = i + 1
        elif arg[:2] == '-K': wg['chunk'] = utils.int_or(arg[2:], None)
        elif arg == '--epsg' or arg == '-P':
            wg['epsg'] = utils.int_or(argv[i + 1], 4326)
            i = i + 1
        elif arg[:2] == '-P': wg['epsg'] = utils.int_or(arg[2:], 4326)
        elif arg == '--z-range' or arg == '-Z':
            zr = argv[i + 1].split('/')
            if len(zr) > 1:
                wg['z_region'] = [None if x == '-' else float(x) for x in zr]
            i = i + 1
        elif arg[:2] == '-Z':
            zr = arg[2:].split('/')
            if len(zr) > 1:
                wg['z_region'] = [None if x == '-' else float(x) for x in zr]
        elif arg == '--w-range' or arg == '-W':
            wr = argv[i + 1].split('/')
            if len(wr) > 1:
                wg['w_region'] = [None if x == '-' else float(x) for x in wr]
            i = i + 1
        elif arg[:2] == '-W':
            wr = arg[2:].split('/')
            if len(wr) > 1:
                wg['w_region'] = [None if x == '-' else float(x) for x in wr]
        elif arg == '-w' or arg == '--weights': wg['weights'] = True
        elif arg == '-t' or arg == '--threads': want_threads = True
        elif arg == '-p' or arg == '--prefix': want_prefix = True
        elif arg == '-a' or arg == '--archive': wg['archive'] = True
        elif arg == '-m' or arg == '--mask': wg['mask'] = True
        elif arg == '-u' or arg == '--uncert':
            wg['mask'] = True
            wg['unc'] = True
        elif arg == '-c' or arg == '--continue': wg['overwrite'] = False
        #elif arg == '-s' or arg == 'spat-meta': wg['spat'] = True
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'
        elif arg == '--verbose' or arg == '-V': wg['verbose'] = True
        elif arg == '--config': want_config = True
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in _waffles_modules.keys():
                    sys.stderr.write(_waffles_module_long_desc({k: _waffles_modules[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
            except: sys.stderr.write(_waffles_module_long_desc(_waffles_modules))
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
                    wg = waffles_config(**wg)
                    dem = waffles_run(wg)
                    sys.exit(0)
            except Exception as e:
                wg = waffles_config_copy(wg)
                utils.echo_error_msg(e)
        else:
            utils.echo_error_msg('specified json file does not exist, {}'.format(wg_user))
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
        else: utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
        
        for key in mod_opts.keys():
            mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]
        mod = opts[0]
        mod_args = tuple(mod_opts[mod])
        wg['mod'] = mod
        wg['mod_args'] = mod_args
    else: 
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a waffles -M module.''')
        sys.exit(-1)
        
    if _waffles_modules[wg['mod']][3]:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
            sys.exit(-1)
            
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls
    
    ## ==============================================
    ## reformat and set the region
    ## ==============================================
    if region is not None:
        try:
            these_regions = [[float(x) for x in region.split('/')]]
        except ValueError: these_regions = gdalfun.gdal_ogr_regions(region)
        except Exception as e:
            utils.echo_error_msg('failed to parse region(s), {}'.format(e))
    else: these_regions = [None]
    if len(these_regions) == 0: utils.echo_error_msg('failed to parse region(s), {}'.format(region))
    
    ## ==============================================
    ## run waffles for each input region.
    ## ==============================================
    these_wgs = []
    for this_region in these_regions:
        twg = waffles_config_copy(wg)
        twg['region'] = this_region
        if want_prefix or len(these_regions) > 1:
            twg['name'] = waffles_append_fn(wg['name'], twg['region'], twg['sample'] if twg['sample'] is not None else twg['inc'])
            
        twg = waffles_config(**twg)
        if want_config:
            this_wg = waffles_config_copy(twg)
            if this_wg is not None:
                utils.echo_msg(json.dumps(this_wg, indent = 4, sort_keys = True))
                #utils.echo_msg(this_wg)
                with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                    utils.echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                    utils.echo_msg('generating major datalist: {}_mjr.datalist'.format(this_wg['name']))
                    wg_json.write(json.dumps(this_wg, indent = 4, sort_keys = True))
            else: utils.echo_error_msg('could not parse config.')
        else:
            if want_threads:
                these_wgs.append(twg)
            else: dem = waffles_run(twg)
    if want_threads:
        wq = queue.Queue()
        num_threads = 2 if len(these_wgs) > 1 else len(these_wgs)
        for _ in range(num_threads):
            t = threading.Thread(target = waffles_queue, args = (wq, ))
            t.daemon = True
            t.start()

        [wq.put(x) for x in these_wgs]
        while True:
            #while not wq.empty():
            time.sleep(2)
            utils.echo_msg_inline('generating dem(s) [{}/{}]'.format((len(these_wgs) - wq.qsize()) - num_threads,len(these_wgs)))
            #if wq.qsize == 0:
            if (len(these_wgs) - wq.qsize()) - num_threads == len(these_wgs): break
            #if wq.empty(): break
        utils.echo_msg('queue complete')
        wq.join()

### End
