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
## run waffles_config() to generate a waffles config dictionary to run a grid with waffle()
##
## Current DEM modules:
## surface (GMT), triangulate (GMT/GDAL), nearest (GMT/GDAL), mbgrid (MBSYSTEM), num (waffles),
## average (GDAL), invdst (GDAL), linear (GDAL), spat-meta, vdatum, coastline, uncertainty
##
## optionally, clip, filter, buffer the resulting DEM.
##
## find data to grid with GEOMODS' fetch.py or use fetch modules as datalist entries in waffles.
##
## a datalist '*.datalist' file should be formatted as in MBSystem, with additional metadata comma separated list:
## ~path ~format ~weight ~metadata,list ~etc
##
## a format of -1 represents a datalist
## a format of -2 represents an archived dataset
## a format of -4 represents a FETCHES module - e.g. `nos:datatype=bag`
## a format of 168 represents XYZ data
## a format of 200 represents GDAL data
## a format of 300 represents LAS/LAZ data <not implemented>
##
## each xyz file in a datalist should have an associated '*.inf' file for faster processing
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
##
### Code:
import sys
import os
import time
import copy

## ==============================================
## queues and threading
## ==============================================
import threading
try:
    import Queue as queue
except: import queue as queue

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
from geomods import vdatumfun

_version = '0.7.0'

## ==============================================
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean', etc.
##
## Requires MBSystem, GMT, GDAL and VDatum for full functionality
## the default waffles config dictionary.
## lambda returns dictionary with default waffles
## ==============================================
waffles_config_copy = lambda wg: copy.deepcopy(wg)

def waffles_config(datalist = None, data = [], region = None, inc = None, name = 'waffles_dem',
                   node = 'pixel', fmt = 'GTiff', extend = 0, extend_proc = 20, weights = None,
                   z_region = None, w_region = None, fltr = [], sample = None, clip = None, chunk = None, epsg = 4326,
                   mod = 'surface', mod_args = (), verbose = False, archive = False, spat = False, mask = False,
                   unc = False, gc = None, clobber = True, init = False):
    if not init:
        if data is None:
            if datalist is not None:
                data = [x[0] for x in datalist2py(datalist)]
            else: data = None
        if datalist is None and len(data) > 0:
            datalist = datalists.datalist_major(data, region = region, major = '{}_major.datalist'.format(name))
        if _waffles_modules[mod]['datalist-p']:
            if datalist is None:
                utils.echo_error_msg('invalid datalist/s entry')
                return(None)
        if region is None or not regions.region_valid_p(region):
            utils.echo_error_msg('invalid region {}'.format(region))
            return(None)
        if inc is None: inc = (region[1] - region[0]) / 500
    return({'datalist': datalist, 'data': data, 'region': region, 'inc': gmtfun.gmt_inc2inc(str(inc)), 'name': name,
            'node': node, 'fmt': 'GTiff', 'extend': utils.int_or(extend, 0), 'extend_proc': utils.int_or(extend_proc, 20),
            'weights': weights, 'z_region': z_region, 'w_region': w_region, 'sample': gmtfun.gmt_inc2inc(str(sample)),
            'clip': clip, 'chunk': chunk, 'epsg': utils.int_or(epsg, 4326), 'mod': mod, 'mod_args': mod_args, 'verbose': verbose,
            'archive': archive, 'spat': spat, 'mask': mask, 'unc': unc, 'clobber': clobber, 'gc': utils.config_check(), 'fltr': fltr})

## ==============================================
## The waffles modules
## the module lambda should point to a function that takes at least
## the waffles config as its first option (e.g. def mod(wg))
## and return a dict with 'dem' = output raster
## ==============================================
_waffles_modules = {
    'surface': {
        'run': lambda args: waffles_gmt_surface(**args),
        'description': '''SPLINE DEM via GMT surface\n
Generate a DEM using GMT's surface command
        
< surface:tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d >
 :tension=[0-1] - Spline tension.''',
        'dem-p': True,
        'datalist-p': True,
    },
    'triangulate': {
        'run': lambda args: waffles_gmt_triangulate(**args),
        'description': '''TRIANGULATION DEM via GMT triangulate\n
Generate a DEM using GMT's triangulate command.
        
< triangulate >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'cudem': {
        'run': lambda args: waffles_cudem(**args),
        'description': '''Generate a CUDEM Bathy/Topo DEM <beta>\n
Generate a DEM using CUDEM methods.
Generates a pre-surface bathymetric grid at a lower resolution, clipped to 
a coastline, for use in final DEM generation.

< cudem:coastline=None:spat=False:upper_limit=-0.1 >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'nearest': {
        'run': lambda args: waffles_nearneighbor(**args),
        'description': '''NEAREST NEIGHBOR DEM via GMT or gdal_grid\n
Generate a DEM using GDAL's gdal_grid command or GMT's nearest command
        
< nearest:radius=6s:use_gdal=False >
 :radius=[value] - Nearest Neighbor search radius
 :use_gdal=[True/False] - use gdal grid nearest algorithm''',
        'dem-p': True,
        'datalist-p': True,
    },
    'num': {
        'run': lambda args: waffles_num(**args),
        'description': '''Uninterpolated DEM populated by <mode>.\n
Generate an uninterpolated DEM using <mode> option.
Using mode of 'A<mode>' uses GMT's xyz2grd command, 
see gmt xyz2grd --help for more info.

< num:mode=n >
 :mode=[key] - specify mode of grid population: 
k (mask), m (mean), n (num), w (wet)''',
        'dem-p': False,
        'datalist-p': True,
    },
    'vdatum': {
        'run': lambda args: waffles_vdatum(**args),
        'description': '''VDATUM transformation grid\n
Generate a VDatum based vertical datum transformation grid.
VDatum coverage areas only!
        
< vdatum:ivert=navd88:overt=mhw:region=3:jar=None >
 :ivert=[vdatum] - Input VDatum vertical datum.
 :overt=[vdatum] - Output VDatum vertical datum.
 :region=[0-10] - VDatum region (3 is CONUS).
 :jar=[/path/to/vdatum.jar] - VDatum jar path - (auto-locates by default)''',
        'dem-p':False,
        'datalist-p': False,
    },
    'mbgrid': {
        'run': lambda args: waffles_mbgrid(**args),
        'description': '''Weighted SPLINE DEM via mbgrid\n
Generate a DEM using MBSystem's mbgrid command.
        
< mbgrid:tension=35:dist=10/3:use_datalists=False >
 :tension=[0-100] - Spline tension.
 :dist=[value] - MBgrid -C switch (distance to fill nodata with spline)
 :use_datalists=[True/False] - use waffles built-in datalists''',
        'dem-p': True,
        'datalist-p': True,
    },
    'IDW': {
        'run': lambda args: waffles_idw(**args),
        'description': '''INVERSE DISTANCE WEIGHTED DEM \n

< IDW:power=2.0:radius=1s >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'invdst': {
        'run': lambda args: waffles_invdst(**args),
        'description': '''INVERSE DISTANCE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.

< invdst:power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1 >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'average': {
        'run': lambda args: waffles_moving_average(**args),
        'description': '''Moving AVERAGE DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.
        
< average:radius1=0.01:radius2=0.01 >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'linear': {
        'run': lambda args: waffles_linear(**args),
        'description': '''LINEAR DEM via gdal_grid\n
Generate a DEM using GDAL's gdal_grid command.

< linear:radius=0.01 >''',
        'dem-p': True,
        'datalist-p': True,
    },
    'spat-meta': {
        'run': lambda args: waffles_spatial_metadata(**args),
        'description': '''generate SPATIAL-METADATA\n
Generate spatial metadata based on the data in the datalist.
Append metadata to the end of each datalist entry to apply in the output 
vector; e.g. /path/to/data/datalist.datalist -1 10 
'Name','Agency','Date','Type','Resolution','HDatum','VDatum','URL'
        
< spat-meta >''',
        'dem-p': False,
        'datalist-p': True,
    },
    'coastline': {
        'run': lambda args: waffles_coastline(**args),
        'description': '''generate a coastline (landmask)\n
Generate a land/sea mask (coastline) based on various datasets.
Include a datalist to use in forming initial coastline.

< coastline:want_nhd=True:want_gmrt=False >
 :want_nhd=[True/False] - Use the USGS NHD database (US Only)
 :want_gmrt=[True/False] - USE the GMRT to fill unknown coastal areas.''',
        'dem-p': False,
        'datalist-p': False,
    },
    'uncertainty': {
        'run': lambda args: waffles_interpolation_uncertainty(**args),
        'description': '''generate DEM UNCERTAINTY\n
Calculate the interpolation uncertainty in relation to distance to nearest measurement
A DEM, Mask and Proximity gird will be generated if not supplied.

< uncertainty:mod=surface:dem=None:msk=None:prox=None:slp=None:sims=10 >
 :mod=[value] - The waffles module to use for interpolation.
 :dem=[path] - Interpolated DEM
 :msk=[path] - Data mask in same size as DEM (Mask is 1 for data, 0 for nodata)
 :prox=[path] - The proximity grid of the data mask
 :sims=[value] - The maximum number of split-sample simulations to perform.''',
        'dem-p': False,
        'datalist-p': False,
    },
}

## ==============================================
## module descriptors (used in cli help)
## ==============================================
_waffles_module_long_desc = lambda x: 'waffles modules:\n% waffles ... -M <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['\033[1m{:14}\033[0m{}\n'.format(key, x[key]['description']) for key in x]) + '\n'
_waffles_module_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

## ==============================================
## the grid-node region
## ==============================================
waffles_grid_node_region = lambda wg: regions.region_buffer(wg['region'], wg['inc'] * .5)

## ==============================================
## the "proc-region" regions.region_buffer(wg['region'],
## (wg['inc'] * 20) + (wg['inc'] * wg['extend']))
## ==============================================
waffles_proc_region = lambda wg: regions.region_buffer(wg['region'], (wg['inc'] * wg['extend_proc']) + (wg['inc'] * wg['extend']))
waffles_coast_region = lambda wg: regions.region_buffer(waffles_proc_region(wg), (wg['inc'] * 200))
waffles_proc_str = lambda wg: regions.region_format(waffles_proc_region(wg), 'gmt')
waffles_proc_bbox = lambda wg: regions.region_format(waffles_proc_region(wg), 'bbox')
waffles_proc_ul_lr = lambda wg: regions.region_format(waffles_proc_region(wg), 'ul_lr')

## ==============================================
## the "dist-region" regions.region_buffer(wg['region'],
## (wg['inc'] * wg['extend']))
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

def waffles_wg_valid_p(wg):
    """check if the waffles config dictionary can be used.

    Checks if wg_config appears valid, in that it has data and module 
    defined.

    Args:
      wg (dict): a waffles config dictionary

    Returns:
      boolean: True if wg is valid else False
    """
    
    try:
        if wg['data'] is None: return(False)
        if wg['mod'] is None: return(False)
        else: return(True)
    except: return(False)
    
def waffles_dlp_hooks(wg):
    """the deafult datalist pass hooks.

    checks for region intersection and possibly z and/or weight range.

    Args:
      wg (dict): a waffles config dictionary

    Returns:
      list: a list of hooks
    """
    
    region = waffles_proc_region(wg)
    #dlp_hooks = datalists.datalist_default_hooks()
    dlp_hooks = []    
    if region is not None:
        dlp_hooks.append(lambda e: datalists.intersect_p(region, e))
    if wg['z_region'] is not None:
        dlp_hooks.append(lambda e: regions.z_region_pass(datalists.inf_entry(e)['minmax'],
                                                         upper_limit = wg['z_region'][1],
                                                         lower_limit = wg['z_region'][0]))
    if wg['w_region'] is not None:
        dlp_hooks.append(lambda e: regions.z_pass(e[2],
                                                  upper_limit = wg['w_region'][1],
                                                  lower_limit = wg['w_region'][0]))
    return(dlp_hooks)

def waffles_yield_datalist(wg):
    """recurse the datalist and do things to it.

    Yield the XYZ data from the datalist, based on info in the
    waffles config dictionary (wg). Uses the default waffles_dlp_hooks.
    

    Args:
      wg (dict): a waffles config dictionary

    Yields:
      list: the xyz line data as list
    """
    
    region = waffles_proc_region(wg)
    dlp_hooks = waffles_dlp_hooks(wg)
    
    dly = datalists.datalist_yield_xyz(wg['datalist'], pass_h = dlp_hooks, wt = 1 if wg['weights'] else None,
                                       region = region, archive = wg['archive'], verbose = wg['verbose'],
                                       z_region = wg['z_region'], epsg = wg['epsg'])
    if wg['mask']: dly = gdalfun.gdal_xyz_mask(dly, '{}_msk.tif'.format(wg['name']), region, wg['inc'],
                                               dst_format = wg['fmt'])
    for xyz in dly: yield(xyz)
    
    if wg['archive']:
        a_dl = os.path.join('archive', '{}.datalist'.format(wg['name']))

        for dir_, _, files in os.walk('archive'):
            for f in files:
                if '.datalist' in f:
                    rel_dir = os.path.relpath(dir_, 'archive')
                    rel_file = os.path.join(rel_dir, f)
                    datalists.datalist_append_entry([rel_file, -1, 1], a_dl)

def waffles_dump_datalist(wg, dst_port = sys.stdout):
    """dump the xyz data from datalist to `dst_port`.

    Args:
      wg (dict): a waffles config dictionary
      dst_port (port): an open output port to write to
    """
    
    for xyz in waffles_yield_datalist(wg):
        xyzfun.xyz_line(xyz, dst_port, True)
        
def waffles_datalist_list(wg):
    """list the datalist entries in the given region.

    Args:
      wg (dict): a waffles config dictionary
    """
    
    for this_entry in datalist(wg['datalist'], wt = 1, pass_h = waffles_dlp_hooks(wg)):
        print(' '.join([','.join(x) if i == 3 else os.path.abspath(str(x)) \
                        if i == 0 else str(x) for i,x in enumerate(this_entry[:-1])]))
        
def waffles_datalists(wg, dump = False, echo = False, infos = False, recurse = True):
    """dump the xyz data from datalist and generate a data mask while doing it.

    Args:
      wg (dict): a waffles config dictionary
      dump (bool): dump the data
      echo (bool): echo the data-entry
      infos (bool): generate info file for entry
      recurse (bool): recurse the datalist

    Returns:
      list: 0,0
    """
    
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
## 3 = spike (outlier) filter (default filter_val = 10) (st.d.)
## ==============================================
def waffles_filter(src_gdal, dst_gdal, fltr = 1, fltr_val = None, split_value = None, mask = None, node = 'pixel'):
    """filter `src_gdal` using smoothing factor `fltr`; optionally
    only smooth bathymetry (sub-zero) using a split_value of 0.

    Args: 
      src_gdal (path): path to source gdal file
      dst_gdal (path): path to output gdal file
      fltr (int): the filter to use, 1, 2 or 3
      flt_val (varies): the filter value, varies by filter.
      split_value (float): an elevation value (only filter below this value)
      node (str): pixel or grid

    Returns:
      0 for success or -1 for failure
    """
    
    if os.path.exists(src_gdal):
        if split_value is not None:
            dem_u, dem_l = gdalfun.gdal_split(src_gdal, split_value)
        else: dem_l = src_gdal

        if int(fltr) == 1: out, status = gdalfun.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        elif int(fltr) == 2: out, status = gmtfun.gmt_grdfilter(dem_l, 'tmp_fltr.tif=gd+n-9999:GTiff',
                                                                dist = fltr_val if fltr_val is not None else '1s',
                                                                node = node, verbose = False)
        elif int(fltr) == 3: out, status = gdalfun.gdal_filter_outliers(dem_l, 'tmp_fltr.tif',
                                                                        fltr_val if fltr_val is not None else 10)
        else: out, status = gdalfun.gdal_blur(dem_l, 'tmp_fltr.tif', fltr_val if fltr_val is not None else 10)
        if status != 0: return(status)
        
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
def waffles_spatial_metadata(wg, geojson = False):
    """generate spatial-metadata

    Args:
      wg (dict): a waffles config dictionary
      geojson(bool): generate a geojson output

    Returns:
      list: [output-vector-fn, status]
    """
    
    dst_layer = '{}_sm'.format(wg['name'])
    dst_vector = dst_layer + '.shp'
    v_fields = ['Name', 'Agency', 'Date', 'Type', 'Resolution', 'HDatum', 'VDatum', 'URL']
    t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]
    utils.remove_glob('{}.*'.format(dst_layer))
    gdalfun.gdal_prj_file('{}.prj'.format(dst_layer), wg['epsg'])
    
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vector)
    if ds is not None: 
        layer = ds.CreateLayer('{}'.format(dst_layer), None, ogr.wkbMultiPolygon)
        [layer.CreateField(ogr.FieldDefn('{}'.format(f), t_fields[i])) for i, f in enumerate(v_fields)]
        [layer.SetFeature(feature) for feature in layer]
    else: layer = None

    for this_entry in datalists.datalist(wg['datalist'], wt = 1 if wg['weights'] else None,
                                         pass_h = waffles_dlp_hooks(wg), yield_dl_entry_p = True,
                                         verbose = wg['verbose']):
        if this_entry[1] == -1 or this_entry[-1] == wg['datalist'].split('.')[0]:
            defn = None if layer is None else layer.GetLayerDefn()            
            twg = waffles_config_copy(wg)
            twg['datalist'] = this_entry[0]
            twg['name'] = '{}_{}_msk'.format(os.path.basename(this_entry[0]).split('.')[0].split(':')[0], regions.region_format(twg['region'], 'fn'))
            if twg['inc'] < gmtfun.gmt_inc2inc('.3333333s'):
                twg['inc'] = gmtfun.gmt_inc2inc('.3333333s')
                twg['extend'] = twg['extend'] / 3
            twg['verbose'] = True
            twg = waffles_config(**twg)

            if len(this_entry[3]) == 8:
                o_v_fields = this_entry[3]
            else: o_v_fields = [twg['name'], 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

            num_out, status = waffles_num(twg, mode='k')
            if status != 0: continue
            ng = num_out['dem'][0]
            waffles_gdal_md(ng, twg)
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
                utils.remove_glob('{}_poly.*'.format(twg['name']), ng)
                #utils.remove_glob(ng)
    ds = None

    gdal_v = utils.int_or(wg['gc']['GDAL'].split('.')[0])
    if gdal_v is not None and gdal_v >= 3:
        utils.run_cmd('ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}\
        '.format(dst_layer, dst_vector))
    sm_out = {'sm': [dst_vector, 'vector']}
    if geojson:
        dst_gj = '{}_sm.inf'.format(wg['name'])
        utils.run_cmd('ogr2ogr {} {} -f GeoJSON'.format(dst_gj, dst_vector), verbose = True)
        sm_out['geojson'] = [dst_gj, 'infos']
    return(sm_out, 0)
    
## ==============================================
##
## Waffles CUDEM generation module
##
## generate CUDEM DEM products, including the integrated
## bathy/topo DEM, Data Mask, Spatial Metadata and Uncertainty Grid.
##
## TODO: Report and Metadta
##
## ==============================================
def waffles_cudem(wg, coastline = None, spat = False, upper_limit = -0.1, fltr = ['1:10'], inc = None):
    """generate bathy/topo DEM suitable for CUDEM project

    generate spatial-metadata as well as a DEM by first
    generating a 'bathy-surface' and possibly a coastline.

    Args:
      wg (dict): a waffles config dictionary
      coastline (path): path to coastline polygon; None to auto-gen
      spat (bool): generate spatial-metadata
      upper_limit (float): bathy/topo zero value
      fltr (list): list of filters to perform
      inc (float): output increment

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the CUDEM module')
        return(None, -1)
    upper_limit = -0.1

    ## ==============================================
    ## Generate the spatial-metadata
    ## ==============================================
    if spat:
        spat_out = waffle(waffles_config(
            region = wg['region'],
            inc = wg['inc'],
            name = wg['name'],
            mod = 'spat-meta',
            verbose = True,
            extend = wg['extend'],
            mask = True,
            epsg = wg['epsg'],
            datalist = None,
            data = wg['data']))
    
    ## ==============================================
    ## generate/process coastline
    ## ==============================================
    _prog = utils._progress('generating coastline XYZ..')
    
    if coastline is None:
        c_wg = waffles_config_copy(wg)
        c_wg['datalist'] = None
        c_wg['data'] = None
        c_wg['mod'] = 'coastline'
        c_wg['mod_args'] = ()
        c_wg['name'] = 'tmp_coast'
        
        coast_out = waffle(c_wg)
        coast_ply = coast_out['cst_ply'][0]
        coast_msk = coast_out['cst'][0]
    else: coast_ply = coastline

    coast_xyz = '{}_coast.xyz'.format(wg['name'])
    c_cmd = 'coastline2xyz.sh -I {} -O {} -Z 0 -W {} -E {} -S {} -N {}'\
        .format(coast_ply, coast_xyz,
                waffles_coast_region(wg)[0], waffles_coast_region(wg)[1],
                waffles_coast_region(wg)[2], waffles_coast_region(wg)[3])
    out, status = utils.run_cmd(c_cmd, verbose = True)
    
    coast_region = datalists.inf_entry([coast_xyz, 168, 1])['minmax']
    wg['data'].append('{} 168 .1'.format(coast_xyz)),

    _prog.end(0, 'generated coastline XYZ.')
    
    ## ==============================================
    ## generate the bathy-surface
    ## using 'surface' with upper_limit of -0.1
    ## at inc*3 spacing
    ## ==============================================
    bathy_out = waffle(waffles_config(
        region = wg['region'],
        z_region = [None, upper_limit + .5],
        name = 'bathy_{}'.format(wg['name']),
        spat = False,
        datalist = None,
        data = wg['data'],
        mod = 'surface',
        mod_args = ('upper_limit={}'.format(upper_limit),),
        sample = wg['inc'],
        fltr = fltr,
        weights = wg['weights'],
        inc = wg['inc'] * 3 if inc is None else inc,
        epsg = wg['epsg'],
        clip = '{}:invert=True'.format(coast_ply),
        extend = wg['extend'],
        extend_proc = 40,
        mask = False,
        verbose = True))

    wg['data'].append('{} 200 .5'.format(bathy_out['dem'][0]))
    
    ## ==============================================
    ## append the bathy-surface to the datalist and
    ## generate final DEM using 'surface'
    ## utils.remove_glob(coast_xyz)
    ## ==============================================
    wg['datalist'] = None
    wg['w_region'] = [.4, None]
    wg = waffles_config(**wg)
    _prog = utils._progress('generating final DEM...')
    surf_dem = waffles_gmt_surface(wg)
    _prog.end(surf_dem[1], 'generated final DEM')
    return(surf_dem)
    
## ==============================================
## Waffles MBGrid module
## ==============================================
def waffles_mbgrid(wg, dist = '10/3', tension = 35, use_datalists = False):
    """Generate a DEM with MBSystem's mbgrid program.
    
    Args:
      wg (dict): a waffles config dictionary
      dist (str): mbgrid -C switch
      tension (float): spline tension while gridding
      use_datalists (bool): use geomods datalists to parse data
    
    if `use_datalists` is True, will parse the datalist through
    waffles instead of mbsystem.
    Note: without `use_datalists` as True, only mbsystem supported data formats may be used.

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]    
    """
    
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
    xsize, ysize, gt = regions.region2gt(waffles_proc_region(wg), wg['inc'])
    
    if len(dist.split('/')) == 1: dist = dist + '/2'
    mbgrid_cmd = ('mbgrid -I{} {} -D{}/{} -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} \
    '.format(wg['datalist'], regions.region_format(region, 'gmt'), xsize, ysize, wg['name'], dist, tension, '-M' if wg['mask'] else ''))
    for out in utils.yield_cmd(mbgrid_cmd, verbose = wg['verbose']): sys.stderr.write('{}'.format(out))
    #out, status = utils.run_cmd(mbgrid_cmd, verbose = wg['verbose'])

    gmtfun.gmt_grd2gdal('{}.grd'.format(wg['name']))
    utils.remove_glob('*.cmd', '*.mb-1', '{}.grd'.format(wg['name']))
    if use_datalists and not wg['archive']: utils.remove_glob('archive')

    if wg['mask']:
        num_grd = '{}_num.grd'.format(wg['name'])
        dst_msk = '{}_msk.tif=gd+n-9999:GTiff'.format(wg['name'])
        out, status = gmtfun.gmt_num_msk(num_grd, dst_msk, verbose = wg['verbose'])
        utils.remove_glob(num_grd, '*_sd.grd')
    if not use_datalists:
        if wg['spat'] or wg['archive']:
            for xyz in waffles_yield_datalist(wg): pass
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, 0)

## ==============================================
## Waffles GMT surface module
## ==============================================
def waffles_gmt_surface(wg, tension = .35, relaxation = 1.2,
                        lower_limit = 'd', upper_limit = 'd'):
    """generate a DEM with GMT surface
    
    Args: 
      wg (dict): a waffles config dictionary
      tension (float): spline tension while gridding
      relaxation (float): spline relaxation while gridding
      lower_limit (float): constrain interpolation to lower limit
      upper_limit (float): constrain interpolation to upper limit

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the SURFACE module')
        return(None, -1)

    dem_surf_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | \
    gmt surface -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -T{} -Z{} -Ll{} -Lu{} -r\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg),
             wg['inc'], wg['name'], tension, relaxation, lower_limit, upper_limit))
    out, status = utils.run_cmd(dem_surf_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg))
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, status)

## ==============================================
## Waffles GMT triangulate module
## ==============================================
def waffles_gmt_triangulate(wg):
    """generate a DEM with GMT triangulate

    Args: 
      wg (dict): a waffles config dictionary    

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    if wg['gc']['GMT'] is None:
        utils.echo_error_msg('GMT must be installed to use the TRIANGULATE module')
        return(None, -1)
    dem_tri_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt triangulate {} -I{:.10f} -V -G{}.tif=gd:GTiff -r\
    '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), wg['inc'], wg['name']))
    out, status = utils.run_cmd(dem_tri_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg))
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, status)

## ==============================================
## Waffles nearest neighbor module
## GMT if available else GDAL
## ==============================================
def waffles_nearneighbor(wg, radius = None, use_gdal = False):
    """genearte a DEM with GMT nearneighbor or gdal_grid nearest

    Args: 
      wg (dict): a waffles config dictionary    
      radius (float): radius to search for neighbors
      use_gdal (bool): use gdal_grid instead of gmt nearneighbor

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    radius = wg['inc'] * 3 if radius is None else gmtfun.gmt_inc2inc(radius)
    if wg['gc']['GMT'] is not None and not use_gdal:
        dem_nn_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '', waffles_proc_str(wg), \
                 wg['inc'], radius, wg['name']))
        out, status = utils.run_cmd(dem_nn_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg))
    else: out, status = waffles_gdal_grid(wg, 'nearest:radius1={}:radius2={}:nodata=-9999'.format(radius, radius))
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, status)

## ==============================================
## Waffles 'NUM grid' module
## Uninterpolated grdi from data;
## num methods include: mask, mean, num, landmask and any gmt grd2xyz -A option.
## ==============================================
def waffles_num(wg, mode = 'n'):
    """Generate an uninterpolated num grid.
    Args: 
      wg (dict): a waffles config dictionary    
      mode (str): num-gridding mode.

    mode of `k` generates a mask grid
    mode of `m` generates a mean grid
    mode of `n` generates a num grid
    mode of `w` generates a landmask grid
    mode of `A<mode> ` generates grid using GMT xyz2grd

    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    if mode[0] == 'A':
        if wg['gc']['GMT'] is None:
            utils.echo_error_msg('GMT must be installed to use the SURFACE module')
            return(None, -1)

        dem_xyz2grd_cmd = ('gmt blockmean {} -I{:.10f}{} -V -r | \
        gmt xyz2grd -{} -V {} -I{:.10f} -G{}.tif=gd+n-9999:GTiff -r\
        '.format(waffles_proc_str(wg), wg['inc'], ' -Wi' if wg['weights'] else '',
                 mode, waffles_proc_str(wg), wg['inc'], wg['name']))
        out, status = utils.run_cmd(dem_xyz2grd_cmd, verbose = wg['verbose'], data_fun = waffles_dl_func(wg))
    else:
        dly = waffles_yield_datalist(wg)
        if wg['weights']: dly = xyzfun.xyz_block(dly, waffles_proc_region(wg), wg['inc'], weights = True)
        out, status = gdalfun.gdal_xyz2gdal(dly, '{}.tif'.format(wg['name']), waffles_proc_region(wg),
                                            wg['inc'], dst_format = wg['fmt'], mode = mode, verbose = wg['verbose'])
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, status)

## ==============================================
## Waffles IDW
## with -w (weights) UIDW
## ==============================================
def waffles_idw(wg, radius='1s', power=2):
    
    def distance(self, pnt0, pnt1):
        return(math.sqrt(sum([(a-b) ** 2 for a, b in zip(pnt0, pnt1)])))

    radius = wg['inc'] * 2 if radius is None else gmtfun.gmt_inc2inc(radius)
    xcount, ycount, dst_gt = regions.region2gt(waffles_proc_region(wg), wg['inc'])
    ds_config = gdalfun.gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt,
                                       gdalfun.gdal_sr_wkt(wg['epsg']), gdal.GDT_Float32,
                                       -9999, wg['fmt'])

    outArray = np.empty((ycount, xcount))
    outArray[:] = np.nan

    if wg['verbose']:
        progress = utils._progress('generating IDW grid @ {} and {}/{}'.format(radius, ycount, xcount))
        i=0

    hash_t, xyz_t = xyzfun.xyz_block_t(waffles_yield_datalist(wg), waffles_proc_region(wg), wg['inc'])

    for y in range(0, ycount):
        if wg['verbose']:
            i+=1
            progress.update_perc((i, ycount))

        for x in range(0, xcount):
            z_list = []
            dw_list = []
            if wg['weights']:
                ww_list = []

            xyz_bucket = []

            xg, yg = utils._pixel2geo(x, y, dst_gt)
            
            block_region = [xg-radius, xg+radius, yg-radius, yg+radius]
            srcwin = regions.region2srcwin(block_region, dst_gt, xcount, ycount)
            
            for y_i in range(srcwin[1], srcwin[1] + srcwin[3], 1):
                for x_i in range(srcwin[0], srcwin[0] + srcwin[2], 1):
                    xyz_s = hash_t[y_i, x_i]
                    [xyz_bucket.append(xyz_t[b]) for b in xyz_s]

            for this_xyz in xyz_bucket:
                d = utils.euc_dst([this_xyz[0], this_xyz[1]], [xg, yg])
                z_list.append(this_xyz[2])
                dw_list.append(1./(d**power))

                if wg['weights']:
                    w = this_xyz[3]
                    ww_list.append(1./(w**power))

            if len(dw_list) > 0:
                dwt = np.transpose(dw_list)
                if wg['weights']:
                    wwt = np.transpose(ww_list)
                    outArray[y,x] = np.dot(z_list, (np.array(dwt)*np.array(wwt)))/sum(np.array(dw_list)*np.array(ww_list))
                else: outArray[y,x] = np.dot(z_list, dwt)/sum(dw_list)
    ds = None

    if wg['verbose']:
        progress.end(0, 'generated IDW grid {}/{}'.format(ycount, xcount))

    outArray[np.isnan(outArray)] = -9999
    out, status = gdalfun.gdal_write(outArray, '{}.tif'.format(wg['name']), ds_config)

    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, status)

## ==============================================
## Waffles GDAL_GRID module
## ==============================================
def waffles_gdal_grid(wg, alg_str = 'linear:radius=1'):
    """run gdal grid using alg_str

    parse the data through xyzfun.xyz_block to get weighted mean before
    building the GDAL dataset to pass into gdal_grid

    Args: 
      wg (dict): a waffles config dictionary
      alg_str (str): the gdal_grid algorithm string
    
    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    _prog = utils._progress('running GDAL GRID {} algorithm @ {}...\
    '.format(alg_str.split(':')[0], regions.region_format(wg['region'], 'fn')))
    _prog_update = lambda x, y, z: _prog.update()
    region = waffles_proc_region(wg)
    dly = xyzfun.xyz_block(waffles_yield_datalist(wg), region, wg['inc'], weights = False if wg['weights'] is None else True)
    ds = gdalfun.xyz2gdal_ds(dly, '{}'.format(wg['name']))
    if ds.GetLayer().GetFeatureCount() == 0: return({},-1)
    xcount, ycount, dst_gt = regions.region2gt(region, wg['inc'])
    gd_opts = gdal.GridOptions(outputType = gdal.GDT_Float32, noData = -9999, format = 'GTiff', width = xcount,
                               height = ycount, algorithm = alg_str, callback = _prog_update if wg['verbose'] else None,
                               outputBounds = [region[0], region[3], region[1], region[2]])
    gdal.Grid('{}.tif'.format(wg['name']), ds, options = gd_opts)
    ds = None
    gdalfun.gdal_set_nodata('{}.tif'.format(wg['name']), -9999)
    _prog.end(0, 'ran GDAL GRID {} algorithm @ {}.'.format(alg_str.split(':')[0], regions.region_format(wg['region'], 'fn')))
    return({'dem': ['{}.tif'.format(wg['name']), 'raster']}, 0)

## ==============================================
## Waffles GDAL_GRID invdist module
## ==============================================
def waffles_invdst(wg, power = 2.0, smoothing = 0.0,
                   radius1 = None, radius2 = None, angle = 0.0,
                   max_points = 0, min_points = 0, nodata = -9999):
    """Generate an inverse distance grid with GDAL

    Args: 
      wg (dict): a waffles config dictionary
    
    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmtfun.gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmtfun.gmt_inc2inc(radius2)
    gg_mod = 'invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={}'\
                             .format(power, smoothing, radius1, radius2, angle, max_points, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_moving_average(wg, radius1 = None, radius2 = None,
                           angle = 0.0, min_points = 0, nodata = -9999):
    """generate a moving average grid with GDAL

    Args: 
      wg (dict): a waffles config dictionary
    
    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    radius1 = wg['inc'] * 2 if radius1 is None else gmtfun.gmt_inc2inc(radius1)
    radius2 = wg['inc'] * 2 if radius2 is None else gmtfun.gmt_inc2inc(radius2)
    gg_mod = 'average:radius1={}:radius2={}:angle={}:min_points={}:nodata={}'.format(radius1, radius2, angle, min_points, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles GDAL_GRID average module
## ==============================================
def waffles_linear(wg, radius = None, nodata = -9999):
    """generate a moving average grid with GDAL

    Args: 
      wg (dict): a waffles config dictionary
      radius (float): search radius
    
    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    radius = wg['inc'] * 4 if radius is None else gmtfun.gmt_inc2inc(radius)
    gg_mod = 'linear:radius={}:nodata={}'.format(radius, nodata)
    return(waffles_gdal_grid(wg, gg_mod))

## ==============================================
## Waffles VDATUM 'conversion grid' module
## U.S. Only
## ==============================================
def waffles_vdatum(wg, ivert = 'navd88', overt = 'mhw',
                   region = '3', jar = None):
    """generate a 'conversion-grid' with vdatum.
    
    output will be the differences (surfaced) between 
    `ivert` and `overt` for the region

    Args: 
      wg (dict): a waffles config dictionary
      ivert (str): input vertical datum string
      overt (str): output vertical datum string
      region (str): vdatum grid region
      jar (path): path to vdatum .jar file
    
    Returns:
      list: [{'dem': ['dem-fn', 'raster']}, status]
    """
    
    vc = vdatumfun._vd_config
    if jar is None:
        vc['jar'] = vdatumfun.vdatum_locate_jar()[0]
    else: vc['jar'] = jar
    vc['ivert'] = ivert
    vc['overt'] = overt
    vc['region'] = region

    gdalfun.gdal_null('empty.tif', waffles_proc_region(wg), 0.00083333, nodata = 0)
    with open('empty.xyz', 'w') as mt_xyz:
        for xyz in gdalfun.gdal_yield_entry(['empty.tif', 200, 1]):
            xyzfun.xyz_line(xyz, mt_xyz, False)
    
    vdatumfun.run_vdatum('empty.xyz', vc)
    
    if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
        with open('result/empty.xyz') as infile:
            empty_infos = xyzfun.xyz_inf(infile)
        print(empty_infos)

        ll = 'd' if empty_infos['minmax'][4] < 0 else '0'
        lu = 'd' if empty_infos['minmax'][5] > 0 else '0'
        wg['data'] = ['result/empty.xyz']
        wg['spat'] = False
        wg['unc'] = False
        wg = waffles_config(**wg)
        vd_out, status = waffles_gmt_surface(wg, tension = 0, upper_limit = lu, lower_limit = ll)
    else:
        utils.echo_error_msg('failed to generate VDatum grid, check settings')
        vd_out = {}
        status = -1
        
    utils.remove_glob('empty.*', 'result/*', '.mjr.datalist', 'result')
    #os.removedirs('result')
    return(vd_out, status)

## ==============================================
## Waffles Interpolation Uncertainty module
## ==============================================
def waffles_interpolation_uncertainty(wg, mod = 'surface', mod_args = (), \
                                      dem = None, msk = None, prox = None, slp = None, \
                                      percentile = 95, zones = ['low-dens', 'mid-dens', 'high-dens'], \
                                      sims = None, chnk_lvl = 6):
    """calculate the interpolation uncertainty
    - as related to distance to nearest measurement.

    Args: 
      wg (dict): a waffles config dictionary
      mod (str): the waffles gridding module
      mod_args (tuple): the waffles module arguments
      dem (path): path to a dem (auto-gen on None)
      msk (path): path to a data mask raster (auto-gen on None)
      prox (path): path to a data proximity raster (auto-gen on None)
      slp (path): path to DEM slope raster (auto-gen on None)
      percentile (int): percentile to calculate
      zones (list): list of density zones
      sims (int): simulations to perform (auto-gen on None)
      chnk_lvl (int): chunking level (auto-gen)

    Returns:
      list: [{output-dictionary}, status]
    """
    
    s_dp = s_ds = None
    unc_out = {}
    zones = ['low-dens','mid-dens','high-dens','low-slp','mid-slp','high-slp']
    
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
    #print(wg)
    if dem is None:
        dem = '{}.tif'.format(wg['mod'])
        unc_out['dem'] = [dem, 'raster']
        tmp_wg = waffles_config_copy(wg)
        tmp_wg['name'] = wg['mod']
        tmp_wg['datalist'] = None
        
        if msk is None:
            msk = '{}_msk.tif'.format(tmp_wg['mod'])
            unc_out['msk'] = [msk, 'raster']
            tmp_wg['mask'] = True
        else: tmp_wg['mask'] = False
        tmp_wg = waffles_config(**tmp_wg)
        waffle(tmp_wg)

    if msk is None:
        
        if msk is None: msk = '{}_msk.tif'.format(wg['mod'])
        tmp_wg = waffles_config_copy(wg)
        tmp_wg['name'] = '{}_msk'.format(wg['name'])
        tmp_wg['mod'] = 'num'
        tmp_wg['mod_args'] = ('mode=k',)
        tmp_wg['datalist'] = None
        tmp_wg = waffles_config(**tmp_wg)
        waffle(tmp_wg)
        
    if prox is None:
        prox = '{}_prox.tif'.format(wg['mod'])
        utils.echo_msg('generating proximity grid {}...'.format(prox))
        gdalfun.gdal_proximity(msk, prox)
        if wg['epsg'] is not None: gdalfun.gdal_set_epsg(prox, wg['epsg'])
    # if slp is None:
    #      slp = '{}_slp.tif'.format(wg['name'])
    #      slp_msk = '{}_slp_msk.tif'.format(wg['name'])
    #      utils.echo_msg('generating slope grid {}...'.format(slp))
    #      gdalfun.gdal_slope(dem, slp, 1)
    #      #utils.run_cmd('gdal_calc.py -A {} -B {} --outfile {} --calc "(1/(A*B))*(A*B)"'.format(slp, msk, slp_msk))
    #      # gdalfun.gdal_slope(dem, '_tmp_slp.tif', 1)
    #      #utils.run_cmd('gdal_calc.py -A {} --outfile _tmp_msk.tif --calc "1/A"'.format(msk), verbose = True)
    #      #gdalfun.gdal_mask('_tmp_msk.tif', '_tmp_msk_inverted.tif', invert = True)
    #      #utils.run_cmd('gdal_calc.py -A {} -B {} --outfile {} --calc "A*B"'.format(slp, '_tmp_msk_inverted.tif', slp_msk), verbose = True)
    #      #utils.remove_glob('_tmp*')
    #      if wg['epsg'] is not None: gdalfun.gdal_set_epsg(slp, wg['epsg'])
    #      #if wg['epsg'] is not None: gdalfun.gdal_set_epsg(slp_msk, wg['epsg'])

    ## ==============================================
    ## region and der. analysis
    ## ==============================================
    region_info = {}
    msk_ds = gdal.Open(msk)
    num_sum, g_max, num_perc = gdalfun.gdal_mask_analysis(msk_ds)
    msk_ds = None

    prox_percentile = gdalfun.gdal_percentile(prox, percentile)
    prox_perc_33 = gdalfun.gdal_percentile(prox, 25)
    prox_perc_66 = gdalfun.gdal_percentile(prox, 75)
    prox_perc_100 = gdalfun.gdal_percentile(prox, 100)

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
    chnk_inc = int((region_info[wg['name']][1] / math.sqrt(g_max)) / region_info[wg['name']][3])
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
        s_sum, s_g_max, s_perc = gdalfun.gdal_mask_analysis(msk_ds, region = sub_region)
        p_perc = gdalfun.gdal_prox_analysis(prox_ds, region = sub_region)
        #slp_perc = gdalfun.gdal_prox_analysis(slp_ds, region = sub_region)
        slp_perc = 0
        s_dc = gdalfun.gdal_gather_infos(dem_ds, region = sub_region, scan = True)
        if p_perc < prox_perc_33 or abs(p_perc - prox_perc_33) < 0.01: zone = zones[2]
        elif p_perc < prox_perc_66 or abs(p_perc - prox_perc_66) < 0.01: zone = zones[1]
        else: zone = zones[0]
        #if slp_perc < slp_perc_33 or abs(slp_perc - slp_perc_33) < 0.01: zone = zones[3]
        #elif slp_perc < slp_perc_66 or abs(slp_perc - slp_perc_66) < 0.01: zone = zones[4]
        #else: zone = zones[5]

        sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, p_perc, slp_perc, s_dc['zr'][0], s_dc['zr'][1], zone]
    dem_ds = msk_ds = prox_ds = slp_ds = None
    utils.echo_msg_inline('analyzing sub-regions [OK]\n')
    
    ## ==============================================
    ## sub-region density and percentiles
    ## ==============================================
    s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
    s_5perc = np.percentile(s_dens, 5)
    s_dens = None
    utils.echo_msg('Sampling density for region is: {:.16f}'.format(s_5perc))

    ## ==============================================
    ## zone analysis / generate training regions
    ## ==============================================
    trainers = []
    t_perc = 95
    s_perc = 50
    
    for z, this_zone in enumerate(zones):
        tile_set = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][8] == zones[z]]
        if len(tile_set) > 0:
            d_50perc = np.percentile(np.array([x[3] for x in tile_set]), 50)
        else: continue
        t_trainers = [x for x in tile_set if x[3] < d_50perc or abs(x[3] - d_50perc) < 0.01]
        utils.echo_msg('possible {} training zones: {} @ MAX {}'.format(zones[z].upper(), len(t_trainers), d_50perc))
        trainers.append(t_trainers)
        
    utils.echo_msg('sorting training tiles by distance...')
    trains = regions.regions_sort(trainers, verbose = False)
    tot_trains = len([x for s in trains for x in s])
    utils.echo_msg('sorted sub-regions into {} training tiles.'.format(tot_trains))
    utils.echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))
    
    ## ==============================================
    ## split-sample simulations and error calculations
    ## sims = max-simulations
    ## ==============================================
    if sims is None: sims = int(len(sub_regions)/tot_trains)
    utils.echo_msg('performing MAX {} SPLIT-SAMPLE simulations...'.format(sims))
    utils.echo_msg('simulation\terrors\tproximity-coeff\tpc_delta')
    #utils.echo_msg('simulation\terrors\tproximity-coeff\tp_diff\tslp-coeff\tslp_diff')
    sim = 0
    last_ec_d = None
    while True:
        utils.echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of MAX {} [{:3}%]'.format(sim, sims, 0))
        status = 0
        sim += 1
        trains = regions.regions_sort(trainers, verbose = False)
        for z, train in enumerate(trains):
            train_h = train[:25]
            ss_samp = s_5perc
            
            ## ==============================================
            ## perform split-sample analysis on each training region.
            ## ==============================================
            for n, sub_region in enumerate(train_h):
                ss_samp = s_5perc
                perc = int(float(n+(len(train_h) * z))/(len(train_h)*len(trains)) * 100)
                utils.echo_msg_inline('performing SPLIT-SAMPLE simulation {} out of MAX {} [{:3}%]'.format(sim, sims, perc))
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
                                        datalist = None,
                                        data = [s_outer, sub_xyz_head],
                                        region = this_region,
                                        inc = wg['inc'],
                                        mod = wg['mod'],
                                        mod_args = wg['mod_args'],
                                        epsg = wg['epsg'],
                                        verbose = False,
                                        mask = True,
                                        chunk = None,
                                        clobber = True)
                    sub_dems = waffle(wc)

                    #if 'dem' in sub_dems.keys() and 'msk' in sub_dem.keys():
                    if os.path.exists(sub_dems['dem'][0]) and os.path.exists(sub_dems['msk'][0]):
                        ## ==============================================
                        ## generate the random-sample data PROX and SLOPE
                        ## ==============================================        
                        sub_prox = '{}_prox.tif'.format(wc['name'])
                        gdalfun.gdal_proximity(sub_dems['msk'][0], sub_prox)
                        
                        # sub_slp = '{}_slp.tif'.format(wc['name'])
                        # gdalfun.gdal_slope(sub_dem, sub_slp)
                        
                        ## ==============================================
                        ## Calculate the random-sample errors
                        ## ==============================================
                        sub_xyd = gdalfun.gdal_query(sub_xyz[sx_cnt:], sub_dems['dem'][0], 'xyd')
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
                utils.remove_glob(o_xyz, 'sub_{}*'.format(n))

        if len(s_dp) > 0:
            d_max = region_info[wg['name']][4]
            s_max = region_info[wg['name']][5]
            s_dp = s_dp[s_dp[:,3] < d_max,:]
            s_dp = s_dp[s_dp[:,3] > 0,:]
            prox_err = s_dp[:,[2,3]]

            if last_ec_d is None:
                last_ec_d = [0, 0.1, 0.2]
                last_ec_diff = 10
            else: last_ec_diff = abs(last_ec_d[2] - last_ec_d[1])

            ec_d = utils.err2coeff(prox_err[:50000000], coeff_guess = last_ec_d, dst_name = wg['name'] + '_prox', xa = 'distance')
            
            ec_diff = abs(ec_d[2] - ec_d[1])
            ec_l_diff = abs(last_ec_diff - ec_diff)
            
            #s_dp = s_dp[s_dp[:,4] < s_max,:]
            #slp_err = s_dp[:,[2,4]]
            #ec_s = utils.err2coeff(slp_err[:50000000], dst_name = wg['name'] + '_slp', xa = 'slope')
            #utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_d[2] - ec_d[1], ec_s, ec_s[2] - ec_s[1]))
            utils.echo_msg('{}\t{}\t{}\t{}'.format(sim, len(s_dp), ec_d, ec_l_diff))
            
            if ec_d[2] < 0.0001: continue
            if abs(ec_d[2] - ec_d[1]) > 2: continue
            if sim >= int(sims): break
            if abs(last_ec_diff - ec_diff) < 0.001: break
            if len(s_dp) >= int(region_info[wg['name']][1] / 10): break
            last_ec_d = ec_d
            #else: utils.echo_msg('{}\t{}\t{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None, None, None))
        else: utils.echo_msg('{}\t{}\t{}\t{}'.format(sim, len(s_dp), None, None))

    ## ==============================================
    ## Save/Output results
    ## apply error coefficient to full proximity grid
    ## ==============================================

    utils.echo_msg('applying coefficient to proximity grid')
    ## USE numpy/gdal instead
    if wg['gc']['GMT'] is None:
        utils.run_cmd('gdal_calc.py -A {} --outfile {}_prox_unc.tif --calc "{}+({}*(A**{}))"'.format(prox, wg['name'], 0, ec_d[1], ec_d[2]), verbose = True)
    else:
        math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_prox_unc.tif=gd+n-9999:GTiff\
        '.format(prox, ec_d[2], ec_d[1], 0, wg['name'])
        utils.run_cmd(math_cmd, verbose = wg['verbose'])
    if wg['epsg'] is not None: status = gdalfun.gdal_set_epsg('{}_prox_unc.tif'.format(wg['name']), epsg = wg['epsg'])
    utils.echo_msg('applied coefficient {} to proximity grid'.format(ec_d))

    # if wg['gc']['GMT'] is None:
    #     utils.run_cmd('gdal_calc.py -A {} --outfile {}_slp_unc.tif --calc "{}+({}*(A**{}))"'.format(slp, wg['name'], 0, ec_s[1], ec_s[2]), verbose = True)
    # else:
    #     math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_slp_unc.tif=gd+n-9999:GTiff\
    #     '.format(slp, ec_s[2], ec_s[1], 0, wg['name'])
    #     utils.run_cmd(math_cmd, verbose = wg['verbose'])
    # if wg['epsg'] is not None: gdalfun.gdal_set_epsg('{}_slp_unc.tif'.format(wg['name']), epsg = wg['epsg'])
    # utils.echo_msg('applied coefficient {} to slope grid'.format(ec_s))
    # utils.run_cmd('gdal_calc.py -A {}_prox_unc.tif -B {}_slp_unc.tif --outfile={}_unc.tif --calc="((A*A)+(B*B))**(1/2.0)" --overwrite'.format(wg['name'], wg['name'], wg['name']))

    utils.remove_glob(prox)
    #utils.remove_glob(slp)

    unc_out['prox_unc'] = ['{}_prox_unc.tif'.format(wg['name']), 'raster']
    unc_out['prox_bf'] = ['{}_prox_bf.png'.format(wg['name']), 'image']
    unc_out['prox_scatter'] = ['{}_prox_scatter.png'.format(wg['name']), 'image']

    return(unc_out, 0)

## ==============================================
## Waffles Coastline module
## generate a coastline (wet/dry mask)
##
## Using various sources (local datalists/USGS NHD/GMRT/ETc.), generate
## a detailed coastline at the given resolution.
## using a local datalist can take a long time if there is a lot
## of lidar to process. Speed things up by limiting the data
## by weight or z (w-range, z-range), and only including the best
## coastal data.
##
## GMT or GMRT will fill in the gaps where local-data and NHD can't fill;
## hopefully this is just to fill areas off-shore and far-inland.
## ==============================================
def waffles_coastline(wg, want_nhd = True, want_gmrt = False):
    """Generate a coastline polygon from various sources.
    
    Args: 
      wg (dict): a waffles config dictionary
      want_nhd (bool): use USGS NHD database
      want_gmrt (bool): use GMRT to infill else use GMT
    
    Returns:
      list: [{output-dictionary}, status]    
    """
    
    w_mask = '{}_w.tif'.format(wg['name'])

    ## ==============================================
    ## wet/dry datalist mask or burn region
    ## ==============================================
    if wg['datalist'] is not None:
        cwg = waffles_config_copy(wg)
        cwg['z-region'] = [-1, 1]
        dly = waffles_yield_datalist(cwg)
        if wg['weights']: dly = xyzfun.xyz_block(dly, waffles_dist_region(cwg), cwg['inc'], weights = True)
        gdalfun.gdal_xyz2gdal(dly, w_mask, waffles_dist_region(cwg), cwg['inc'],
                              dst_format = cwg['fmt'], mode = 'w', verbose = cwg['verbose'])
    else:
        regions.region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        xsize, ysize, gt = regions.region2gt(waffles_dist_region(wg), wg['inc'])
        utils.run_cmd('gdal_rasterize -ts {} {} -te {} -burn -9999 -a_nodata -9999 \
        -ot Int32 -co COMPRESS=DEFLATE -a_srs EPSG:{} region_buff.shp {}\
        '.format(xsize, ysize, regions.region_format(waffles_dist_region(wg), 'te'), wg['epsg'], w_mask), verbose = False)

    ## ==============================================
    ## load the wet/dry mask array
    ## ==============================================
    ds = gdal.Open(w_mask)
    if ds is not None:
        ds_config = gdalfun.gdal_gather_infos(ds)
        region = regions.gt2region(ds_config)
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
    ## Fetch NHD (NHD High/Plus) data from TNM
    ## to fill in near-shore areas. High resoultion data
    ## varies by location...
    ## ==============================================
    if want_nhd:
        u_mask = '{}_u.tif'.format(wg['name'])
        regions.region2ogr(waffles_dist_region(wg), 'region_buff.shp')
        xsize, ysize, gt = regions.region2gt(waffles_proc_region(wg), wg['inc'])
        utils.run_cmd('gdal_rasterize -ts {} {} -te {} -burn -9999 -a_nodata -9999 \
        -ot Int32 -co COMPRESS=DEFLATE -a_srs EPSG:{} region_buff.shp {}\
        '.format(xsize, ysize, regions.region_format(waffles_dist_region(wg), 'te'), wg['epsg'], u_mask), verbose = False)
        utils.remove_glob('region_buff.*')

        fl = fetches._fetch_modules['tnm'](waffles_proc_region(wg), ["Name LIKE '%Hydro%'"], None, True)
        r_shp = []
        for result in fl._parse_results(e = 'HU-2 Region,HU-4 Subregion,HU-8 Subbasin'):
            if fetches.fetch_file(result[0], os.path.join(result[2], result[1]), verbose = True) == 0:
                gdb_zip = os.path.join(result[2], result[1])
                gdb_files = utils.unzip(gdb_zip)
                gdb_bn = os.path.basename('.'.join(gdb_zip.split('.')[:-1]))
                gdb = gdb_bn + '.gdb'

                utils.run_cmd('ogr2ogr {}_NHDArea.shp {} NHDArea -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, regions.region_format(wg['region'], 'ul_lr')), verbose = False)
                if os.path.exists('{}_NHDArea.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDArea.shp'.format(gdb_bn))
                utils.run_cmd('ogr2ogr {}_NHDPlusBurnWaterBody.shp {} NHDPlusBurnWaterBody -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, regions.region_format(wg['region'], 'ul_lr')), verbose = False)
                if os.path.exists('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn))
                utils.run_cmd('ogr2ogr {}_NHDWaterBody.shp {} NHDWaterBody -where "FType = 390" -clipdst {} -overwrite 2>&1\
                '.format(gdb_bn, gdb, regions.region_format(wg['region'], 'ul_lr')), verbose = False)
                if os.path.exists('{}_NHDWaterBody.shp'.format(gdb_bn)):
                    r_shp.append('{}_NHDWaterBody.shp'.format(gdb_bn))
                utils.remove_glob(gbd)
            else: utils.echo_error_msg('unable to fetch {}'.format(result))

            [utils.run_cmd('ogr2ogr -skipfailures -update -append nhdArea_merge.shp {} 2>&1\
            '.format(shp), verbose = False) for shp in r_shp]
            utils.run_cmd('gdal_rasterize -burn 1 nhdArea_merge.shp {}'.format(u_mask), verbose = True)
            utils.remove_glob('nhdArea_merge.*', 'NHD_*', *r_shp)

        ## ==============================================
        ## update wet/dry mask with nhd data
        ## ==============================================
        utils.echo_msg('filling the coast mask with NHD data...')
        c_ds = gdal.Open(u_mask)
        c_ds_arr = c_ds.GetRasterBand(1).ReadAsArray()
        c_ds = gdal.Open(u_mask)
        for this_xyz in gdalfun.gdal_parse(c_ds):
            xpos, ypos = utils._geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
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
        utils.run_cmd('gmt grdlandmask {} -I{} -r -Df -G{}=gd:GTiff -V -N1/0/1/0/1\
        '.format(regions.region_format(waffles_dist_region(wg), 'gmt'), wg['inc'], g_mask), verbose = wg['verbose'])
    else:
        fl = fetches._fetch_modules['gmrt'](regions.region_buffer(waffles_dist_region(wg), 5, pct = True), [], None, True)
        r = fl._parse_results(fl._filter_results(), layer = 'topo-mask')
        gmrt_tif = r[0][1]
        if fetches.fetch_file(r[0][0], gmrt_tif, verbose = True) == 0:
            utils.run_cmd('gdalwarp {} {} -tr {} {} -overwrite'.format(gmrt_tif, g_mask, wg['inc'], wg['inc']), verbose = True)
            utils.remove_glob(gmrt_tif)

    ## ==============================================
    ## update wet/dry mask with gsshg/gmrt data
    ## speed up!
    ## ==============================================
    utils.echo_msg('filling the coast mask with gsshg/gmrt data...')
    c_ds = gdal.Open(g_mask)
    c_ds_arr = c_ds.GetRasterBand(1).ReadAsArray()
    c_ds = gdal.Open(g_mask)
    for this_xyz in gdalfun.gdal_parse(c_ds):
        xpos, ypos = utils._geo2pixel(this_xyz[0], this_xyz[1], dst_gt)
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
    gdalfun.gdal_write(coast_array, '{}.tif'.format(wg['name']), ds_config)

    ## ==============================================
    ## convert to vector
    ## ==============================================
    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('tmp_c_{}.shp'.format(wg['name']))
    if tmp_ds is not None:
        tmp_layer = tmp_ds.CreateLayer('tmp_c_{}'.format(wg['name']), None, ogr.wkbMultiPolygon)
        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
        gdalfun.gdal_polygonize('{}.tif'.format(wg['name']), tmp_layer, verbose = wg['verbose'])        
        tmp_ds = None
    utils.run_cmd('ogr2ogr -dialect SQLITE -sql "SELECT * FROM tmp_c_{} WHERE DN=0 order by ST_AREA(geometry) desc limit 8"\
    {}.shp tmp_c_{}.shp'.format(wg['name'], wg['name'], wg['name']), verbose = True)
    utils.remove_glob('tmp_c_{}.*'.format(wg['name']))
    utils.run_cmd('ogrinfo -dialect SQLITE -sql "UPDATE {} SET geometry = ST_MakeValid(geometry)" {}.shp\
    '.format(wg['name'], wg['name']))
                      
    coast_out = {
        'cst': ['{}.tif'.format(wg['name']), 'raster'],
        'cst_ply': ['{}.shp'.format(wg['name']), 'vector'],
    }
    
    return(coast_out, 0)

def waffles_gdal_md(in_gdal, wg, cudem = False):
    """add metadata to the waffles dem
    
    Args: 
      in_gdal (path): path to input gdal file
      wg (dict): a waffles config dictionary
      cudem (bool): add CUDEM metadata    
    """
    
    try:
        ds = gdal.Open(in_gdal, gdal.GA_Update)
    except: ds = None
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
    """a queue to run waffle"""
    
    while True:
        this_wg = q.get()
        dem = waffle(this_wg)

        q.task_done()
        
## ==============================================
## Waffles run waffles module via wg_config()
## ==============================================
def waffle(wg):
    """generate a DEM using wg dict settings

    note: see waffles_config() to generate a wg config.

    - runs the waffles module to generate the DEM
    - optionally clips the output to shapefile
    - optionally filters the output
    - optionally resamples the output
    - cuts the output to dist-size
    - reformats the output to final format
    - sets metadata in output

    Args:
      wg (dict): a waffles config dictionary

    Returns:
      dict: {output-dictionary}    
    """
    
    if wg is None:
        utils.echo_error_msg('invalid configuration, {}'.format(wg))
        sys.exit(-1)
        
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    dems = {}
    args_d = {}
    args_d = utils.args2dict(wg['mod_args'], args_d)
    wg['region'] = waffles_grid_node_region(wg) if wg['node'] == 'grid' else wg['region']
    if wg['verbose']: _prog = utils._progress('running WAFFLES module {} with {} [{}]...'.format(wg['mod'], wg['mod_args'], args_d))

    ## ==============================================
    ## optionally generate the DEM in chunks
    ## skip if wg['overwrite'] is False and file exists
    ## ==============================================
    if wg['chunk'] is not None:
        xcount, ycount, dst_gt = regions.region2gt(wg['region'], wg['inc'])
        s_regions = regions.region_chunk(wg['region'], wg['inc'], n_chunk=(xcount/wg['chunk']) + 1, buff=100)
    else: s_regions = [wg['region']]

    chunks = []
    for region in s_regions:
        this_wg = waffles_config_copy(wg)
        this_wg['region'] = region
        this_wg['name'] = 'chunk_{}'.format(regions.region_format(region, 'fn'))
        args_d['wg'] = this_wg
        waffles_out = {}
        
        ## ==============================================
        ## gererate the DEM (run the module)
        ## ==============================================
        #try:
        waffles_out, status = _waffles_modules[this_wg['mod']]['run'](args_d)
        if wg['mask']: waffles_out['msk'] = ['{}_msk.tif'.format(this_wg['name']), 'raster']
        chunks.append(waffles_out)
        #except KeyboardInterrupt as e:
        #    utils.echo_error_msg('killed by user, {}'.format(e))
        #    sys.exit(-1)
        #except Exception as e:
        #    utils.echo_error_msg('{}'.format(e))
        #    [utils.remove_glob('{}.*'.format(waffles_out[x][0].split('.')[0])) for x in waffles_out.keys()]
        ##    status = -1
        #    continue
        
        for out_key in waffles_out.keys():
            #print(waffles_out[out_key])
            if waffles_out[out_key][1] == 'raster':

                ## ==============================================
                ## check DEM for content and set projection, etc.
                ## ==============================================
                this_dem = waffles_out[out_key][0]
                #gdalfun.gdal_set_nodata(this_dem)
                if not os.path.exists(this_dem): continue
                gdi = gdalfun.gdal_infos(this_dem, scan = True)
                #print(gdi)
                if gdi is not None:
                    if np.isnan(gdi['zr'][0]):
                        #print(gdi['zr'][0])
                        utils.echo_warning_msg('no data found11 , skipping')
                        #sys.exit()
                        utils.remove_glob(this_dem)
                        continue
                else: continue

                gdalfun.gdal_set_epsg(this_dem, this_wg['epsg'])
                waffles_gdal_md(this_dem, this_wg)

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
                        fltr_args['node'] = this_wg['node']
                        fltr_args = utils.args2dict(fltrs[2:], fltr_args)        

                        if this_wg['verbose']: utils.echo_msg('filtering {} using filter {}@{}...\
                        '.format(this_dem, fltr_args['fltr'], fltr_args['fltr_val']))
                        if fltr_args['fltr'] == 2:
                            if this_wg['gc']['GMT'] is None:
                                continue
                        try:
                            status = waffles_filter(this_dem, 'tmp_s.tif', **fltr_args)
                            if status == 0: os.rename('tmp_s.tif', this_dem)
                        except TypeError as e: utils.echo_error_msg('{}'.format(e))
                        except Exception as e: utils.echo_error_msg(e)

                ## ==============================================
                ## optionally resample the DEM 
                ## ==============================================
                if this_wg['sample'] is not None:
                    if this_wg['verbose']: utils.echo_msg('resampling {}...'.format(this_dem))
                    # if this_wg['gc']['GMT'] is not None:
                    #     gmtfun.gmt_sample_inc(this_dem, inc = this_wg['sample'], verbose = this_wg['verbose'])
                    #     if this_wg['mask']:
                    #         if this_wg['verbose']: utils.echo_msg('resampling {}...'.format(this_dem_msk))
                    #         gmtfun.gmt_sample_inc(this_dem_msk, inc = this_wg['sample'], verbose = this_wg['verbose'])
                    # else:
                    out, status = utils.run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} _tmp.tif\
                    '.format(this_wg['sample'], this_wg['sample'], this_dem,
                             regions.region_format(waffles_proc_region(this_wg), 'te'), verbose = this_wg['verbose']))
                    if status == 0: os.rename('_tmp.tif', '{}'.format(this_dem))
                    if this_wg['mask']:
                        out_status = utils.run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te {} _tmp_msk.tif\
                        '.format(this_wg['sample'], this_wg['sample'], this_dem_msk,
                                 regions.region_format(waffles_proc_region(this_wg), 'te'), verbose = this_wg['verbose']))
                        if status == 0: os.rename('_tmp_msk.tif', '{}'.format(this_dem_msk))
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
                        utils.echo_warning_msg('no data found, skipping')
                        utils.remove_glob(this_dem, this_dem_msk)
                        #if this_wg['mask']: utils.remove_glob(this_dem_msk)
                        continue
                else: continue

                ## ==============================================
                ## cut dem to final size -
                ## region buffered by (inc * extend) or
                ## sample * extend) if sample if specified
                ## ==============================================
                #utils.echo_msg('cutting {}: {}'.format(out_key, this_dem))
                try:
                    out = gdalfun.gdal_cut(this_dem, waffles_dist_region(this_wg), 'tmp_cut.tif')
                    if out is not None: os.rename('tmp_cut.tif', this_dem)
                except OSError as e:
                    utils.remove_glob('tmp_cut.tif')
                    utils.echo_error_msg('cut failed, is the dem open somewhere, {}'.format(e)) 

    ## ==============================================
    ## merge the chunks and remove any remnants
    ## ==============================================
    if len(chunks) > 0 and len(chunks[0].keys()) > 0:
        for out_key in chunks[0].keys():
            if out_key == 'dem': out_dem = '{}.{}'.format(wg['name'], chunks[0][out_key][0].split('.')[-1])
            else: out_dem = '{}_{}.{}'.format(wg['name'], out_key, chunks[0][out_key][0].split('.')[-1])

            if len(chunks) > 1:
                if chunks[0][out_key][1] == 'raster':
                    #print(chunks)
                    
                    out, status = utils.run_cmd('gdal_merge.py -n -9999 -a_nodata -9999 \
                    -ps {inc} -{inc} -ul_lr {ullr} -o {o_dem} {chunks} -co TILED=YES -co COMPRESS=DEFLATE -co PREDICTOR=3\
                    '.format(inc=wg['inc'], ullr=waffles_dist_ul_lr(wg), o_dem=out_dem, chunks=' '.join([x[out_key][0] for x in chunks])),
                                                verbose = True)
                elif chunks[0][out_key][1] == 'vector':
                    out, status = utils.run_cmd('ogrmerge.py {} {}'.format(dem_vect, ' '.join(chunks_vect)))
                [utils.remove_glob('{}.*'.format(x[out_key][0].split('.')[0])) for x in chunks]            
            else:
                if os.path.exists(chunks[0][out_key][0]):
                    if chunks[0][out_key][1] == 'raster':
                        gdalfun.gdal_gdal2gdal(chunks[0][out_key][0], dst_fmt = wg['fmt'], dst_gdal = out_dem)
                        utils.remove_glob(chunks[0][out_key][0])
                    elif chunks[0][out_key][1] == 'vector':
                        out, status = utils.run_cmd('ogr2ogr {} {} -overwrite\
                        '.format(out_dem, chunks[0][out_key][0]), verbose = True)
                        gdalfun.ogr_remove_ds(chunks[0][out_key][0])
                        #utils.remove_glob('{}.*'.format(chunks[0][out_key][0].split('.')[0]))
                    else: os.rename(chunks[0][out_key][0], '{}_{}.{}\
                    '.format(wg['name'], out_key, chunks[0][out_key][0].split('.')[-1]))
                else:
                    utils.echo_error_msg('{} not found.'.format(chunks[0][out_key][0]))
                    #continue

            ## ==============================================
            ## convert to final format if not geotiff and
            ## set the projection and other metadata
            ## ==============================================
            if os.path.exists(out_dem):
                if chunks[0][out_key][1] == 'raster':
                    if wg['fmt'] != 'GTiff':
                        orig_dem = out_dem
                        if wg['gc']['GMT'] is not None:
                            out_dem = gmtfun.gmt_grd2gdal(orig_dem, wg['fmt'])
                            utils.remove_glob(orig_dem)
                        else: out_dem = gdalfun.gdal_gdal2gdal(orig_dem, wg['fmt'])
                    try:
                        gdalfun.gdal_set_epsg(out_dem, wg['epsg'])
                        waffles_gdal_md(out_dem, wg, True if wg['mod'] == 'cudem' else False)
                    except: pass
                    dems[out_key] = [out_dem, 'raster']
                elif chunks[0][out_key][1] == 'vector':
                    if gdalfun.ogr_empty_p(out_dem):
                        utils.echo_error_msg('{} is empty, something may have gone wrong...'.format(out_dem))
                        utils.remove_glob('{}.*'.format(out_dem.split('.')[0]))
                    else: dems[out_key] = [out_dem, 'vector']
                else: dems[out_key] = [out_dem, chunks[0][out_key][1]]

    ## ==============================================
    ## optionally generate associated uncertainty grid
    ## ==============================================
    if wg['unc'] and _waffles_modules[wg['mod']]['dem-p']:
        try:
            if os.path.exists(dems['dem']) and os.path.exists(dems['msk']):
                utils.echo_msg('generating uncertainty')
                out_unc = waffles_interpolation_uncertainty(wg = wg, mod = wg['mod'],
                                                            mod_args = wg['mod_args'],
                                                            dem = dems['dem'], msk = dems['msk'])
        except Exception as e: utils.echo_error_msg('failed to calculate uncertainty, {}'.format(e))
                
    utils.remove_glob('waffles_dem_mjr.datalist', wg['datalist'])
    if wg['verbose']: _prog.end(status, 'ran WAFFLES module {} with {} [{}].'.format(wg['mod'], wg['mod_args'], args_d))
    return(dems)

## ==============================================
## Command-line Interface (CLI)
## $ waffles
##
## waffles cli
## ==============================================
waffles_cli_usage = """waffles [OPTIONS] <datalist/entry>

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
\t\t\t3: Spike Filter at -T3:<stand-dev. threshhold>.
\t\t\tThe -T switch may be set multiple times to perform multiple filters.
\t\t\tAppend :split_value=<num> to only filter values below z-value <num>.
\t\t\te.g. -T1:10:split_value=0 to smooth bathymetry using Gaussian filter
  -Z, --z-region\t\tRestrict data processing to records that fall within the z-region
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
  -q, --quiet\t\tLower verbosity to a quiet. (overrides --verbose)

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
""".format(datalists._datalist_fmts_short_desc(), _waffles_module_short_desc(_waffles_modules))

def waffles_cli(argv = sys.argv):
    """run waffles from command-line

    generates a waffles_config from the command-line options
    and either outputs the or runs the waffles_config
    on each region supplied (multiple regions can be supplied
    by using a vector file as the -R option.)
    See `waffles_cli_usage` for full cli options.
    """
    
    wg = waffles_config(init = True)
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
        elif arg == '-c' or arg == '--continue': wg['clobber'] = False
        #elif arg == '-s' or arg == 'spat-meta': wg['spat'] = True
        elif arg == '-r' or arg == '--grid-node': wg['node'] = 'grid'
        elif arg == '--verbose' or arg == '-V': wg['verbose'] = True
        elif arg == '--quiet' or arg == '-q': wg['verbose'] = False
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
        elif arg[0] == '-':
            print(waffles_cli_usage)
            utils.echo_error_msg('{} is not a valid waffles cli switch'.format(arg))
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
                    dem = waffle(wg)
                    sys.exit(0)
            except Exception as e:
                wg = waffles_config_copy(wg)
                utils.echo_error_msg(e)
                sys.exit(-1)
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
        else:
            utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
            sys.exit(-1)

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
        
    if _waffles_modules[wg['mod']]['datalist-p']:
        if len(dls) == 0:
            sys.stderr.write(waffles_cli_usage)
            utils.echo_error_msg('''must specify a datalist/entry, try gmrt or srtm for global data.''')
            sys.exit(-1)

    ## ==============================================
    ## check the increment
    ## ==============================================
    if wg['inc'] is None:
        sys.stderr.write(waffles_cli_usage)
        utils.echo_error_msg('''must specify a gridding increment.''')
        sys.exit(-1)
    
    ## ==============================================
    ## set the datalists and names
    ## ==============================================
    wg['data'] = dls
    
    ## ==============================================
    ## reformat and set the region
    ## ==============================================
    # if region is not None:
    #     try:
    #         these_regions = [[float(x) for x in region.split('/')]]
    #     except ValueError: these_regions = regions.gdal_ogr_regions(region)
    #     except Exception as e:
    #         utils.echo_error_msg('failed to parse region(s), {}'.format(e))
    # else: these_regions = [None]
    # if len(these_regions) == 0: utils.echo_error_msg('failed to parse region(s), {}'.format(region))

    if region is not None:
        this_region = regions.str2region(region)
        these_regions = [this_region]
    else: these_regions = []
    if len(these_regions) == 0: utils.echo_error_msg('failed to parse region(s), {}'.format(region))
        
    ## ==============================================
    ## run waffles for each input region.
    ## ==============================================
    these_wgs = []
    if want_threads:
        _prog = utils._progress('Amassing WAFFLES config info for {} regions...'.format(len(these_regions)))
    for i, this_region in enumerate(these_regions):
        this_region = regions.region_region(this_region)
        twg = waffles_config_copy(wg)
        twg['region'] = this_region
        if want_prefix or len(these_regions) > 1:
            twg['name'] = waffles_append_fn(wg['name'], twg['region'], twg['sample'] if twg['sample'] is not None else twg['inc'])
            
        twg = waffles_config(**twg)

        if want_config:
            this_wg = waffles_config_copy(twg)
            if this_wg is not None:
                utils.echo_msg(json.dumps(this_wg, indent = 4, sort_keys = True))
                with open('{}.json'.format(this_wg['name']), 'w') as wg_json:
                    utils.echo_msg('generating waffles config file: {}.json'.format(this_wg['name']))
                    utils.echo_msg('generating major datalist: {}_mjr.datalist'.format(this_wg['name']))
                    wg_json.write(json.dumps(this_wg, indent = 4, sort_keys = True))
            else: utils.echo_error_msg('could not parse config.')
        else:
            if want_threads:
                _prog.update_perc((i, len(these_regions)))
                these_wgs.append(twg)
            else: dems = waffle(twg)
    if want_threads:
        _prog.end(0, 'Amassed WAFFLES configs for {} regions.'.format(len(these_wgs)))
        wq = queue.Queue()
        num_threads = 2 if len(these_wgs) > 1 else len(these_wgs)
        for _ in range(num_threads):
            t = threading.Thread(target = waffles_queue, args = (wq, ))
            t.daemon = True
            t.start()
        [wq.put(x) for x in these_wgs]
        wq.join()
### End
