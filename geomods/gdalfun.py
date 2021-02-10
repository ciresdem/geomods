### gdalfun.py
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
## General GDAL functions and processing
##
### Code:

import os

import gdal
import ogr
import osr
import sys
import numpy as np
import affine
import json

from geomods import utils

## ==============================================
## GDAL Wrappers and Functions - gdalfun.py
## ==============================================
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.UseExceptions()
_gdal_progress = gdal.TermProgress #crashes on osgeo4w
_gdal_progress_nocb = gdal.TermProgress_nocb

def gdal_fext(src_drv_name):
    '''find the common file extention given a GDAL driver name
    older versions of gdal can't do this, so fallback to known standards.

    returns list of known file extentions or None'''
    
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER): fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None: return(fexts.split()[0])
        else: return(None)
    except:
        if src_drv_name == 'GTiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        return(fext)

def _geo2pixel_old(geo_x, geo_y, geoTransform):
   '''convert a geographic x,y value to a pixel location of geoTransform'''
   if geoTransform[2] + geoTransform[4] == 0:
       pixel_x = ((geo_x - geoTransform[0]) / geoTransform[1])# + .5
       pixel_y = ((geo_y - geoTransform[3]) / geoTransform[5])# + .5
   else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
   return(round(pixel_x), round(pixel_y))

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''convert a geographic x,y value to a pixel location of geoTransform'''
    forward_transform = affine.Affine.from_gdal(*geoTransform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x + 0.5), int(pixel_y + 0.5)
    return(pixel_x, pixel_y)

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''convert a pixel location to geographic coordinates given geoTransform'''
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    '''apply geotransform to in_x,in_y'''
    
    out_x = geoTransform[0] + (in_x + 0.5) * geoTransform[1] + (in_y + 0.5) * geoTransform[2]
    out_y = geoTransform[3] + (in_x + 0.5) * geoTransform[4] + (in_y + 0.5) * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    '''invert the geotransform'''
    
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs(det) < 0.000000000000001: return
    invDet = 1.0 / det
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet
    return(outGeoTransform)
    
## ==============================================
## GDAL projections metadata and regions
## ==============================================
def gdal_sr_wkt(epsg, esri = False):
    '''convert an epsg code to wkt

    returns the sr Wkt or None'''
    
    try:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(epsg))
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)

def gdal_prj_file(dst_fn, epsg):
    '''generate a .prj file given an epsg code

    returns 0'''
    
    with open(dst_fn, 'w') as out:
        out.write(gdal_sr_wkt(int(epsg), True))
    return(0)
    
def gdal_set_epsg(src_fn, epsg = 4326):
    '''set the projection of gdal file src_fn to epsg

    returns status-code (0 == success)'''

    try:
        ds = gdal.Open(src_fn, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        ds.SetProjection(gdal_sr_wkt(int(epsg)))
        ds = None
        return(0)
    else: return(None)

def gdal_getEPSG(src_ds):
    '''returns the EPSG of the given gdal data-source'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_region = gdal_gt2region(ds_config)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)

    return(srs_auth)
    
def gdal_set_nodata(src_fn, nodata = -9999):
    '''set the nodata value of gdal file src_fn

    returns 0'''

    try:
        ds = gdal.Open(src_fn, gdal.GA_Update)
    except: ds = None
    if ds is not None:
        band = ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        ds = None
        return(0)
    else: return(None)

def gdal_gt2region(ds_config):
    '''convert a gdal geo-tranform to a region [xmin, xmax, ymin, ymax] via a data-source config dict.

    returns region of gdal data-source'''
    
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + geoT[1] * ds_config['nx'], geoT[3] + geoT[5] * ds_config['ny'], geoT[3]])

def gdal_region2gt(region, inc):
    '''return a count info and a gdal geotransform based on extent and cellsize

    returns a list [xcount, ycount, geot]'''

    dst_gt = (region[0], inc, 0, region[3], 0, (inc * -1.))
    
    this_origin = _geo2pixel(region[0], region[3], dst_gt)
    this_end = _geo2pixel(region[1], region[2], dst_gt)
    #this_size = ((this_end[0] - this_origin[0]) + 1, (this_end[1] - this_origin[1]) + 1)
    this_size = (this_end[0] - this_origin[0], this_end[1] - this_origin[1])
    
    return(this_size[0], this_size[1], dst_gt)

def gdal_region_warp(region, s_warp = 4326, t_warp = 4326):
    '''warp region from source `s_warp` to target `t_warp`, using EPSG keys

    returns the warped region'''

    if s_warp is None or s_warp == t_warp: return(region)
    
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(int(s_warp))

    if t_warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(t_warp))
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)        
        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[0], region[2]))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[1], region[3]))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(region)

def gdal_ogr_regions(src_ds):
    '''return the region(s) of the ogr dataset'''
    
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                penv = ogr.CreateGeometryFromWkt(pwkt).GetEnvelope()
                these_regions.append(penv)
        poly = None
    return(these_regions)
    
def gdal_create_polygon(coords, xpos = 1, ypos = 0):
    '''convert coords to Wkt

    returns polygon as wkt'''
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[xpos], coord[ypos])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def gdal_wkt2geom(wkt):
    return(ogr.CreateGeometryFromWkt(wkt))

def gdal_region2wkt(region):

    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    return(gdal_create_polygon(eg))
    
def gdal_region2geom(region):
    '''convert an extent [west, east, south, north] to an OGR geometry

    returns ogr geometry'''
    
    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    geom = ogr.CreateGeometryFromWkt(gdal_create_polygon(eg))
    return(geom)

def gdal_region2ogr(region, dst_ogr, append = False):
    '''convert a region string to an OGR vector'''

    dst_wkt = gdal_region2wkt(region)
    
    driver = ogr.GetDriverByName('ESRI Shapefile')

    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)
        
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPolygon)

    dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
    dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(dst_wkt))
    dst_feat.SetField('id', 1)
    dst_lyr.CreateFeature(dst_feat)
    dst_feat = None

    dst_ds = None

def gdal_region(src_ds, warp = None):
    '''return the extent of the src_fn gdal file.
    warp should be an epsg to warp the region to.

    returns the region of the gdal data-source'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_region = gdal_gt2region(ds_config)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)
    
    if srs_auth is None or srs_auth == warp: warp = None

    if warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(warp))
        #dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)

        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(ds_region[0], ds_region[2]))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(ds_region[1], ds_region[3]))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        ds_region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(ds_region)

def gdal_srcwin(src_ds, region):
    '''given a gdal file src_fn and a region [w, e, s, n],
    output the appropriate gdal srcwin.

    returns the gdal srcwin'''

    ds_config = gdal_gather_infos(src_ds)
    geoT = ds_config['geoT']

    this_origin = [0 if x < 0 else x for x in _geo2pixel(region[0], region[3], geoT)]
    this_end = [0 if x < 0 else x for x in _geo2pixel(region[1], region[2], geoT)]
    this_size = [0 if x < 0 else x for x in ((this_end[0] - this_origin[0]), (this_end[1] - this_origin[1]))]
    #this_size = [0 if x < 0 else x for x in ((this_end[0] - this_origin[0]), (this_end[1] - this_origin[1]))]
    if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
    if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
    #this_size = [0 if x < 0 else x for x in this_size]
    return(this_origin[0], this_origin[1], this_size[0], this_size[1])

def gdal_sample_inc(src_grd, inc = 1, verbose = False):
    '''resamele src_grd to toggle between grid-node and pixel-node grid registration.'''
    
    out, status = utils.run_cmd('gdalwarp -tr {:.10f} {:.10f} {} -r bilinear -te -R{} -r -Gtmp.tif=gd+n-9999:GTiff'.format(inc, inc, src_grd, src_grd), verbose = verbose)
    if status == 0: os.rename('tmp.tif', '{}'.format(src_grd))
    return(status)

## ==============================================
## GDAL gather/set infos on a gdal file/datasource
## ==============================================
def gdal_infos(src_fn, region = None, scan = False):
    '''scan gdal file src_fn and gather region info.

    returns region dict.'''
    
    if os.path.exists(src_fn):
        try:
            ds = gdal.Open(src_fn)
        except: ds = None
        if ds is not None:
            dsc = gdal_gather_infos(ds, region = region, scan = scan)
            ds = None
            return(dsc)
        else: return(None)
    else: return(None)

def gdal_gather_infos(src_ds, region = None, scan = False):
    '''gather information from `src_ds` GDAL dataset

    returns gdal_config dict.'''

    gt = src_ds.GetGeoTransform()
    if region is not None:
        srcwin = gdal_srcwin(src_ds, region)
    else: srcwin = (0, 0, src_ds.RasterXSize, src_ds.RasterYSize)
    src_band = src_ds.GetRasterBand(1)
    dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])

    ds_config = {
        'nx': srcwin[2],
        'ny': srcwin[3],
        'nb': srcwin[2] * srcwin[3],
        'geoT': dst_gt,
        'proj': src_ds.GetProjectionRef(),
        'dt': src_band.DataType,
        'dtn': gdal.GetDataTypeName(src_band.DataType),
        'ndv': src_band.GetNoDataValue(),
        'fmt': src_ds.GetDriver().ShortName,
        'zr': None,
    }
    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
    if scan:
        src_arr = src_band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        ds_config['zr'] = (np.amin(src_arr), np.amax(src_arr))
        #ds_config['zr'] = src_band.ComputeRasterMinMax()
        src_arr = None
    return(ds_config)

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''set a datasource config dictionary

    returns gdal_config dict.'''
    
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def gdal_cpy_infos(src_config):
    '''copy src_config

    returns copied src_config dict.'''
    
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

## ==============================================
## Write an array to a gdal file
## ==============================================
def gdal_write (src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff'):
    '''write src_arr to gdal file dst_gdal using src_config

    returns [output-gdal, status-code]'''
    
    driver = gdal.GetDriverByName(dst_fmt)
    if os.path.exists(dst_gdal):
        try:
            driver.Delete(dst_gdal)
        except Exception as e:
            utils.echo_error_msg(e)
            utils.remove_glob(dst_gdal)
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    if ds is not None:
        ds.SetGeoTransform(ds_config['geoT'])
        ds.SetProjection(ds_config['proj'])
        ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
        ds.GetRasterBand(1).WriteArray(src_arr)
        ds = None
        return(dst_gdal, 0)
    else: return(None, -1)

def gdal_null(dst_fn, region, inc, nodata = -9999, outformat = 'GTiff'):
    '''generate a `null` grid with gdal

    returns [output-gdal, status-code]'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    null_array = np.zeros((ycount, xcount))
    null_array[null_array == 0] = nodata
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(4326), gdal.GDT_Float32, -9999, outformat)
    return(gdal_write(null_array, dst_fn, ds_config))

## ==============================================
## Manipulate GDAL file
## ==============================================
def gdal_cut(src_gdal, region, dst_fn):
    '''cut src_fn gdal file to srcwin and output dst_fn gdal file

    returns [output-gdal, status-code]'''

    try:
        src_ds = gdal.Open(src_gdal)
    except: src_ds = None
    if src_ds is not None:
        ds_config = gdal_gather_infos(src_ds)
        srcwin = gdal_srcwin(src_ds, region)
        gt = ds_config['geoT']
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        dst_gt = (gt[0] + (srcwin[0] * gt[1]), gt[1], 0., gt[3] + (srcwin[1] * gt[5]), 0., gt[5])
        out_ds_config = gdal_set_infos(srcwin[2], srcwin[3], srcwin[2] * srcwin[3], dst_gt, ds_config['proj'], ds_config['dt'], ds_config['ndv'], ds_config['fmt'])
        src_ds = None
        return(gdal_write(ds_arr, dst_fn, out_ds_config))
    else: return(None, -1)

def gdal_clip(src_gdal, src_ply = None, invert = False):
    '''clip dem to polygon `src_ply`, optionally invert the clip.

    returns [gdal_raserize-output, gdal_rasterize-return-code]'''

    # clip src_ply to src_gdal extent
    gi = gdal_infos(src_gdal)
    g_region = gdal_gt2region(gi)
    tmp_ply = 'tmp_clp_ply.shp'
    utils.run_cmd('ogr2ogr {} {} -clipsrc {} {} {} {}'.format(tmp_ply, src_ply, g_region[0], g_region[3], g_region[1], g_region[2]), verbose = True)
    if gi is not None and src_ply is not None:
        gr_inv = '-i' if invert else ''
        gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'\
                 .format(gi['ndv'], gr_inv, os.path.basename(tmp_ply).split('.')[0], tmp_ply, src_gdal)
        out, status = utils.run_cmd(gr_cmd, verbose = True)
    else: return(None)
    utils.remove_glob('tmp_clp_ply.*')
    return(out, status)

def np_split(src_arr, sv = 0, nd = -9999):
    '''split numpy `src_arr` by `sv` (turn u/l into `nd`)

    returns [upper_array, lower_array]'''
    
    try:
        sv = int(sv)
    except: sv = 0
    u_arr = np.array(src_arr)
    l_arr = np.array(src_arr)
    u_arr[u_arr <= sv] = nd
    l_arr[l_arr >= sv] = nd
    return(u_arr, l_arr)

def gdal_split(src_gdal, split_value = 0):
    '''split raster file `src_gdal`into two files based on z value

    returns [upper_grid-fn, lower_grid-fn]'''
    
    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_u.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_l.tif'.format(os.path.basename(src_gdal)[:-4]))
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = gdal_gather_infos(src_ds)
        dst_config = gdal_cpy_infos(src_config)
        dst_config['fmt'] = 'GTiff'
        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny'])
        ua, la = np_split(ds_arr, split_value, src_config['ndv'])
        gdal_write(ua, dst_upper, dst_config)
        gdal_write(la, dst_lower, dst_config)
        ua = la = ds_arr = src_ds = None
        return([dst_upper, dst_lower])
    else: return(None)

def gdal_crop(src_ds):
    '''crop `src_ds` GDAL datasource by it's NoData value. 

    returns [cropped array, cropped_gdal_config].'''
    
    ds_config = gdal_gather_infos(src_ds)
    ds_arr = src_ds.GetRasterBand(1).ReadAsArray()

    src_arr[elev_array == ds_config['ndv']] = np.nan
    nans = np.isnan(src_arr)
    nancols = np.all(nans, axis=0)
    nanrows = np.all(nans, axis=1)

    firstcol = nancols.argmin()
    firstrow = nanrows.argmin()        
    lastcol = len(nancols) - nancols[::-1].argmin()
    lastrow = len(nanrows) - nanrows[::-1].argmin()

    dst_arr = src_arr[firstrow:lastrow,firstcol:lastcol]
    src_arr = None

    dst_arr[np.isnan(dst_arr)] = ds_config['nv']
    GeoT = ds_config['geoT']
    dst_x_origin = GeoT[0] + (GeoT[1] * firstcol)
    dst_y_origin = GeoT[3] + (GeoT[5] * firstrow)
    dst_geoT = [dst_x_origin, GeoT[1], 0.0, dst_y_origin, 0.0, GeoT[5]]
    ds_config['geoT'] = dst_geoT
    return(dst_arr, ds_config)

def gdal_mask(src_gdal, dst_gdal, invert = False):
    '''transform src_gdal to a raster mask (1 = data; 0 = nodata)
    if invert is True, 1 = nodata, 0 = data.'''
    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    if ds is not None:
        ds_band = ds.GetRasterBand(1)
        ds_array = ds_band.ReadAsArray()
        ds_config = gdal_gather_infos(ds)
        ds_config['dtn'] = 'Int32'
        ndv = ds_band.GetNoDataValue()
        #print(ndv)
        if not invert:
            ds_array[ds_array != ndv] = 1
            ds_array[ds_array == ndv] = 0
        else:
            ds_array[ds_array != ndv] = 0
            ds_array[ds_array == ndv] = 1
            
        return(gdal_write(ds_array, dst_gdal, ds_config))
    else: return(-1, -1)
    
def gdal_mask_analysis(mask = None):
    '''mask is a GDAL mask grid of 0/1

    returns [sum, max, percentile]'''
    
    msk_sum = gdal_sum(mask)
    msk_gc = gdal_infos(mask)
    msk_max = float(msk_gc['nx'] * msk_gc['ny'])
    msk_perc = float((msk_sum / msk_max) * 100.)
    return(msk_sum, msk_max, msk_perc)

def gdal_mask_analysis2(src_gdal, region = None):

    ds_config = gdal_gather_infos(src_gdal)
    if region is not None:
        srcwin = gdal_srcwin(src_gdal, region)
    else: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
        
    msk_sum = np.sum(ds_arr)
    msk_max = float(srcwin[2] * srcwin[3])
    msk_perc = float((msk_sum / msk_max) * 100.)
    dst_arr = None
    
    return(msk_sum, msk_max, msk_perc)

def gdal_prox_analysis2(src_gdal, region = None):

    ds_config = gdal_gather_infos(src_gdal)
    if region is not None:
        srcwin = gdal_srcwin(src_gdal, region)
    else: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    ds_arr = src_gdal.GetRasterBand(1).ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])

    prox_perc = np.percentile(ds_arr, 95)
    dst_arr = None
    
    return(prox_perc)
        
def gdal_proximity(src_fn, dst_fn):
    '''compute a proximity grid via GDAL

    return 0 if success else None'''
    
    prog_func = None
    src_ds = gdal.Open(src_fn)
    dst_ds = None
    if src_ds is not None:
        src_band = src_ds.GetRasterBand(1)
        ds_config = gdal_gather_infos(src_ds)

        if dst_ds is None:
            drv = gdal.GetDriverByName('GTiff')
            dst_ds = drv.Create(dst_fn, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'], [])
        dst_ds.SetGeoTransform(ds_config['geoT'])
        dst_ds.SetProjection(ds_config['proj'])
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(ds_config['ndv'])
        gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)
        dst_band = src_band = dst_ds = src_ds = None
        return(0)
    else: return(None)

## ==============================================
## query a GDAL file
## ==============================================
def gdal_query(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    returns array of values'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    try:
        ds = gdal.Open(src_grd)
    except: ds = None
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()

        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                #print(g)
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    xyzl.append(np.array(outs))
        dsband = ds = None
        out_array = np.array(xyzl)
    return(out_array)

def gdal_yield_query(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()
        dsband = ds = None
        
        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    yield(outs)

def gdal_query2(src_xyz, src_grd, out_form):
    '''query a gdal-compatible grid file with xyz data.
    out_form dictates return values

    yields out_form results'''
    
    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================
    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_band = ds.GetRasterBand(1)
        ds_gt = ds_config['geoT']
        ds_nd = ds_config['ndv']
        tgrid = ds_band.ReadAsArray()

        ## ==============================================   
        ## Process the src xyz data
        ## ==============================================
        for xyz in src_xyz:
            x = xyz[0]
            y = xyz[1]
            try: 
                z = xyz[2]
            except: z = ds_nd

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)
                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd
                #print(g)
                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    #c = con_dec(math.fabs(d / (g + 0.00000001) * 100), 2)
                    #s = con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    #d = con_dec(d, 4)
                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    yield(outs)
        dsband = ds = None
        
## ==============================================
## GDAL command-line wrappers
## ==============================================
def gdal_gdal2gdal(src_grd, dst_fmt = 'GTiff', epsg = 4326, dst_gdal = None, co = True):
    '''convert the gdal file to gdal using gdal

    return output-gdal-fn'''
    
    if os.path.exists(src_grd):
        if dst_gdal is None:
            dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdal_fext(dst_fmt))
        if not co:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {}\
            '.format(src_grd, dst_gdal, dst_fmt))
        else:
            gdal2gdal_cmd = ('gdal_translate {} {} -f {} -co TILED=YES -co COMPRESS=DEFLATE\
            '.format(src_grd, dst_gdal, dst_fmt))
        out, status = utils.run_cmd(gdal2gdal_cmd, verbose = False)
        if status == 0: return(dst_gdal)
        else: return(None)
    else: return(None)

def gdal_slope(src_gdal, dst_gdal, s = 111120):
    '''generate a slope grid with GDAL

    return cmd output and status'''
    
    gds_cmd = 'gdaldem slope {} {} {} -compute_edges'.format(src_gdal, dst_gdal, '' if s is None else '-s {}'.format(s))
    return(utils.run_cmd(gds_cmd))

def gdal_polygonize(src_gdal, dst_layer, verbose = False):
    '''run gdal.Polygonize on src_ds and add polygon to dst_layer'''

    try:
        ds = gdal.Open('{}'.format(src_gdal))
    except: ds = None
    if ds is not None:
        ds_arr = ds.GetRasterBand(1)
        if verbose: utils.echo_msg('polygonizing {}'.format(src_gdal))
        status = gdal.Polygonize(ds_arr, None, dst_layer, 0, callback = _gdal_progress if verbose else None)
        ds = ds_arr = None
        return(0, 0)
    else: return(-1, -1)

## ==============================================
## numpy and arrays
## ==============================================
def gdal_sum(src_gdal):
    '''sum the z vale of src_gdal

    return the sum'''

    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    if ds is not None:
        ds_array = ds.GetRasterBand(1).ReadAsArray() 
        sums = np.sum(ds_array)
        ds = ds_array = None
        return(sums)
    else: return(None)

def gdal_percentile(src_gdal, perc = 95):
    '''calculate the `perc` percentile of src_fn gdal file.

    return the calculated percentile'''
    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    if ds is not None:
        ds_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        x_dim = ds_array.shape[0]
        ds_array_flat = ds_array.flatten()
        ds_array = ds_array_flat[ds_array_flat != 0]
        if len(ds_array) > 0:
            p = np.percentile(ds_array, perc)
            #percentile = 2 if p < 2 else p
        else: p = 2
        ds = ds_array = ds_array_flat = None
        return(p)
    else: return(None)

def np_gaussian_blur(in_array, size):
    '''blur an array using fftconvolve from scipy.signal
    size is the blurring scale-factor.

    returns the blurred array'''
    
    from scipy.signal import fftconvolve
    from scipy.signal import convolve
    import scipy.fftpack._fftpack as sff
    padded_array = np.pad(in_array, size, 'symmetric')
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(size)))
    g = (g / g.sum()).astype(in_array.dtype)
    in_array = None
    #try:
    out_array = fftconvolve(padded_array, g, mode = 'valid')
    #except:
    #print('switching to convolve')
    #out_array = convolve(padded_array, g, mode = 'valid')
    return(out_array)

def gdal_blur(src_gdal, dst_gdal, sf = 1):
    '''gaussian blur on src_gdal using a smooth-factor of `sf`
    runs np_gaussian_blur(ds.Array, sf)'''

    try:
        ds = gdal.Open(src_gdal)
    except: ds = None
    
    if ds is not None:
        ds_config = gdal_gather_infos(ds)
        ds_array = ds.GetRasterBand(1).ReadAsArray(0, 0, ds_config['nx'], ds_config['ny'])
        ds = None
        msk_array = np.array(ds_array)
        msk_array[msk_array != ds_config['ndv']] = 1
        msk_array[msk_array == ds_config['ndv']] = np.nan
        ds_array[ds_array == ds_config['ndv']] = 0
        smooth_array = np_gaussian_blur(ds_array, int(sf))
        smooth_array = smooth_array * msk_array
        mask_array = ds_array = None
        smooth_array[np.isnan(smooth_array)] = ds_config['ndv']
        return(gdal_write(smooth_array, dst_gdal, ds_config))
    else: return([], -1)

## ==============================================
## OGR functions
## ==============================================
def gdal_ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.

    returns the output feature class'''
    
    if dst_defn is None: dst_defn = src_layer.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    feats = len(src_layer)
    utils.echo_msg('unioning {} features'.format(feats))
    for n, f in enumerate(src_layer):
        _gdal_progress_nocb((n+1 / feats) * 100)
        if f.GetField(src_field) == 0:
            src_layer.DeleteFeature(f.GetFID())
        elif f.GetField(src_field) == 1:
            f.geometry().CloseRings()
            wkt = f.geometry().ExportToWkt()
            multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
            src_layer.DeleteFeature(f.GetFID())
    #union = multi.UnionCascaded() ## slow on large multi...
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometryDirectly(multi)
    #union = multi = None
    return(out_feat)

def ogr_clip(src_ogr, dst_ogr, clip_region = None, dn = "ESRI Shapefile"):
    
    driver = ogr.GetDriverByName(dn)
    ds = driver.Open(src_ogr, 0)
    layer = ds.GetLayer()

    gdal_region2ogr(clip_region, 'tmp_clip.shp')
    c_ds = driver.Open('tmp_clip.shp', 0)
    c_layer = c_ds.GetLayer()
    
    dst_ds = driver.CreateDataSource(dst_ogr)
    dst_layer = dst_ds.CreateLayer(dst_ogr.split('.')[0], geom_type=ogr.wkbMultiPolygon)

    layer.Clip(c_layer, dst_layer)
    #ogr.Layer.Clip(layer, c_layer, dst_layer)

    ds = c_ds = dst_ds = None

def ogr_empty_p(src_ogr):

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(src_ogr, 0)

    if ds is not None:
        layer = ds.GetLayer()
        fc = layer.GetFeatureCount()
        if fc == 0:
            return(True)
        else: return(False)
    else: return(True)


def ogr_remove_ds(src_ds, src_fmt = 'ESRI Shapefile'):
    drv = ogr.GetDriverByName(src_fmt)
    drv.DeleteDataSource(src_ds)
    
## ==============================================
## GDAL XYZ functions
##
## gdal processing (datalist fmt:200)
## ==============================================
def xyz2gdal_ds(src_xyz, dst_ogr):
    '''Make a point vector OGR DataSet Object from src_xyz

    returns the in-memory GDAL data-source'''
    
    ds = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
    layer = ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPoint25D)
    fd = ogr.FieldDefn('long', ogr.OFTReal)
    fd.SetWidth(10)
    fd.SetPrecision(8)
    layer.CreateField(fd)
    fd = ogr.FieldDefn('lat', ogr.OFTReal)
    fd.SetWidth(10)
    fd.SetPrecision(8)
    layer.CreateField(fd)
    fd = ogr.FieldDefn('elev', ogr.OFTReal)
    fd.SetWidth(12)
    fd.SetPrecision(12)
    layer.CreateField(fd)
    f = ogr.Feature(feature_def = layer.GetLayerDefn())
    for this_xyz in src_xyz:
        f.SetField(0, this_xyz[0])
        f.SetField(1, this_xyz[1])
        f.SetField(2, this_xyz[2])
        wkt = 'POINT({:.8f} {:.8f} {:.10f})'.format(this_xyz[0], this_xyz[1], this_xyz[2])
        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)
    return(ds)

def gdal_xyz2gdal(src_xyz, dst_gdal, region, inc, dst_format = 'GTiff', mode = 'n', epsg = 4326, verbose = False):
    '''Create a GDAL supported grid from xyz data 
    `mode` of `n` generates a num grid
    `mode` of `m` generates a mean grid
    `mode` of `k` generates a mask grid
    `mode` of 'w' generates a wet/dry mask grid

    returns output, status'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    if verbose:
        utils.echo_msg('gridding data with mode: {} to {}'.format(mode, dst_gdal))
        utils.echo_msg('grid size: {}/{}'.format(ycount, xcount))
    if mode == 'm' or mode == 'w':
        sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    #else: gdt = gdal.GDT_Int32
    ptArray = np.zeros((ycount, xcount))
    #if mode == 'w': ptArray[ptArray == 0] = np.nan
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(epsg), gdt, -9999, dst_format)
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    if mode == 'm' or mode == 'w':
                        sumArray[ypos, xpos] += z
                    if mode == 'n' or mode == 'm':
                        ptArray[ypos, xpos] += 1
                    else: ptArray[ypos, xpos] = 1
                except Exception as e:
                    if verbose: utils.echo_error_msg(e)
                    pass
    if mode == 'm' or mode == 'w':
        ptArray[ptArray == 0] = np.nan
        outarray = sumArray / ptArray
        if mode == 'w':
            outarray[outarray >= 0] = 1
            outarray[outarray < 0] = 0
    elif mode == 'n': outarray = ptArray
    else: outarray = ptArray
    outarray[np.isnan(outarray)] = -9999
    return(gdal_write(outarray, dst_gdal, ds_config))

def gdal_xyz_mask(src_xyz, dst_gdal, region, inc, dst_format='GTiff', epsg = 4326):
    '''Create a num grid mask of xyz data. The output grid
    will contain 1 where data exists and 0 where no data exists.

    yields the xyz data'''
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    ptArray = np.zeros((ycount, xcount))
    ds_config = gdal_set_infos(xcount, ycount, xcount * ycount, dst_gt, gdal_sr_wkt(epsg), gdal.GDT_Float32, -9999, 'GTiff')
    for this_xyz in src_xyz:
        yield(this_xyz)
        x = this_xyz[0]
        y = this_xyz[1]
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                try:
                    ptArray[ypos, xpos] = 1
                except: pass
    out, status = gdal_write(ptArray, dst_gdal, ds_config)
    
def gdal_chunks(src_fn, n_chunk = 10):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.

    returns a list of chunked filenames or None'''
    
    band_nums = []
    o_chunks = []
    if band_nums == []: band_nums = [1]
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk

    try:
        src_ds = gdal.Open(src_fn)
    except: src_ds = None
    if src_ds is not None:
        ds_config = gdal_gather_infos(src_ds)
        band = src_ds.GetRasterBand(1)
        gt = ds_config['geoT']

        while True:
            y_chunk = n_chunk
            while True:
                if x_chunk > ds_config['nx']:
                    this_x_chunk = ds_config['nx']
                else: this_x_chunk = x_chunk

                if y_chunk > ds_config['ny']:
                    this_y_chunk = ds_config['ny']
                else: this_y_chunk = y_chunk

                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = this_x_chunk - this_x_origin
                this_y_size = this_y_chunk - this_y_origin

                ## chunk size aligns with grid cellsize
                if this_x_size == 0 or this_y_size == 0: break
                
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)
                dst_gt = [this_geo_x_origin, float(gt[1]), 0.0, this_geo_y_origin, 0.0, float(gt[5])]
                
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                if not np.all(band_data == band_data[0,:]):
                    o_chunk = '{}_chnk{}x{}.tif'.format(os.path.basename(src_fn).split('.')[0], x_i_chunk, i_chunk)
                    dst_fn = os.path.join(os.path.dirname(src_fn), o_chunk)
                    o_chunks.append(dst_fn)

                    dst_config = gdal_cpy_infos(ds_config)
                    dst_config['nx'] = this_x_size
                    dst_config['ny'] = this_y_size
                    dst_config['geoT'] = dst_gt                    
                    gdal_write(band_data, dst_fn, dst_config)

                band_data = None

                if y_chunk > ds_config['ny']:
                    break
                else: 
                    y_chunk += n_chunk
                    i_chunk += 1
            if x_chunk > ds_config['nx']:
                break
            else:
                x_chunk += n_chunk
                x_i_chunk += 1
        src_ds = None
        return(o_chunks)
    else: return(None)

def gdal_parse(src_ds, dump_nodata = False, srcwin = None, mask = None, warp = None, verbose = False, z_region = None, step = 1):
    '''parse the data from gdal dataset src_ds (first band only)
    optionally mask the output with `mask` or transform the coordinates to `warp` (epsg-code)

    yields the parsed xyz data'''

    #if verbose: sys.stderr.write('waffles: parsing gdal file {}...'.format(src_ds.GetDescription()))
    ln = 0
    band = src_ds.GetRasterBand(1)
    ds_config = gdal_gather_infos(src_ds)
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds_config['proj'])
    src_srs.AutoIdentifyEPSG()
    srs_auth = src_srs.GetAuthorityCode(None)
    
    if srs_auth is None or srs_auth == warp: warp = None

    if warp is not None:
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(warp))
        ## GDAL 3+
        #dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        
    gt = ds_config['geoT']
    msk_band = None
    if mask is not None:
        src_mask = gdal.Open(mask)
        msk_band = src_mask.GetRasterBand(1)
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])
    nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
    if band.GetNoDataValue() is not None: nodata.append('{:g}'.format(band.GetNoDataValue()))
    if dump_nodata: nodata = []
    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        if z_region is not None:
            if z_region[0] is not None:
                band_data[band_data < z_region[0]] = -9999
            if z_region[1] is not None:
                band_data[band_data > z_region[1]] = -9999
        if msk_band is not None:
            msk_data = msk_band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
            band_data[msk_data==0]=-9999
        band_data = np.reshape(band_data, (srcwin[2], ))
        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            z = band_data[x_i]
            if '{:g}'.format(z) not in nodata:
                ln += 1
                geo_x,geo_y = _pixel2geo(x, y, gt)
                if warp is not None:
                    point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(geo_x, geo_y))
                    point.Transform(dst_trans)
                    pnt = point.GetPoint()
                    line = [pnt[0], pnt[1], z]
                else: line = [geo_x, geo_y, z]
                yield(line)
    band = None
    src_mask = None
    msk_band = None
    if verbose: utils.echo_msg('parsed {} data records from {}'.format(ln, src_ds.GetDescription()))

def xyz_transform(src_fn, xyz_c = None, inc = None, vdatum = None, region = None, datalist = None, verbose = False):
    '''transform xyz data file src_fn to xyz data'''
    
    src_xyz, src_zips = procs_unzip(src_fn, _known_datalist_fmts[168])
    if xyz_c is None: xyz_c = copy.deepcopy(_xyz_config)
    if not os.path.exists(os.path.join(os.path.dirname(src_fn), 'xyz')): os.mkdir(os.path.join(os.path.dirname(src_fn), 'xyz'))

    if vdatum is not None:
        vds = vdatum.split(',')
        iv = vds[0]
        ov = vds[1]
        vdc = _vd_config
        if vdc['jar'] is None: vdc['jar'] = vdatum_locate_jar()[0]
        vdc['ivert'] = iv
        vdc['overt'] = ov
        vdc['delim'] = 'comma' if xyz_c['delim'] == ',' else 'space' if xyz_c['delim'] == ' ' else xyz_c['delim']
        vdc['skip'] = xyz_c['skip']
        vdc['xyzl'] = ','.join([str(x) for x in [xyz_c['xpos'], xyz_c['ypos'], xyz_c['zpos']]])
        out, status = run_vdatum(src_xyz, vdc)
        src_result = os.path.join('result', os.path.basename(src_xyz))
    else: src_result = src_xyz
    xyz_final = os.path.join('xyz', os.path.basename(src_xyz))

    if datalist is not None:
        datalist = os.path.join(os.path.dirname(src_fn), 'xyz', datalist)
    
    if os.path.exists(src_result):
        with open(src_result, 'r') as in_n, open(xyz_final, 'w') as out_n:
            for xyz in xyz_parse(in_n, xyz_c = xyz_c, region = region, verbose = verbose):
                xyz_line(xyz, out_n)

    if os.path.exists(xyz_final):
        if datalist is not None:
            mbsfun.mb_inf(xyz_final)
            datalist_append_entry([os.path.basename(xyz_final), 168, 1], datalist)
            if verbose: utils.echo_msg('appended xyz file {} to datalist {}'.format(xyz_final, datalist))
                
    if src_xyz != src_fn: utils.remove_glob(src_xyz)
    
def gdal2xyz_chunks(src_fn, chunk_value = 1000, inc  = None, epsg = None, vdatum = None, datalist = None, verbose = False):
    '''chunk and transform gdal file src_fn to xyz file'''
    
    if verbose:
        utils.echo_msg('------------------------------------')
        utils.echo_msg('input gdal:\t\t{}'.format(src_fn))
        utils.echo_msg('chunk size:\t\t{}'.format(chunk_value))
        utils.echo_msg('output epsg:\t\t{}'.format(epsg))
        utils.echo_msg('output increment:\t{}'.format(inc))
        utils.echo_msg('vdatum string:\t{}'.format(vdatum))
        utils.echo_msg('output datalist:\t{}'.format(datalist))
        utils.echo_msg('------------------------------------')

    if not os.path.exists('xyz'): os.mkdir('xyz')

    src_gdal, src_zips = procs_unzip(src_fn, _known_datalist_fmts[200])
    src_c = gdal_infos(src_gdal)
    if verbose:
        utils.echo_msg('{}'.format(src_c))
        utils.echo_msg('chunking grid file {}...'.format(src_gdal))
    chunks = gdal_chunks(src_gdal, chunk_value)
    if verbose: utils.echo_msg('generated {} chunks from {}'.format(len(chunks), src_gdal))
    if src_gdal != src_fn: utils.remove_glob(src_gdal)

    if vdatum is not None:
        vds = vdatum.split(',')
        if len(vds) < 2:
            utils.echo_error_msg('bad vdatum string {}'.format(vdatum))
            vdatum = None
        else:
            iv = vds[0]
            ov = vds[1]
            vdc = _vd_config
            if vdc['jar'] is None: vdc['jar'] = vdatum_locate_jar()[0]
            vdc['ivert'] = iv
            vdc['overt'] = ov
            vdc['delim'] = 'space'
            vdc['skip'] = '0'
            vdc['xyzl'] = '0,1,2'

    if datalist is not None:
        datalist = os.path.join('xyz', datalist)
        
    for i,chunk in enumerate(chunks):
        utils.echo_msg('* processing chunk {} [{}/{}]...'.format(chunk, i+1, len(chunks)))
        xyz_chunk = '{}.xyz'.format(chunk.split('.')[0])
        xyz_chunk_final = os.path.join('xyz', os.path.basename(xyz_chunk))

        if epsg is not None or inc is not None:
            tif_chunk = '{}_warp.tif'.format(chunk.split('.')[0])
            gdw = 'gdalwarp {} -dstnodata -9999 -overwrite {}'.format(chunk, tif_chunk)
            if epsg is not None: gdw += ' -t_srs EPSG:{}'.format(epsg)
            if inc is not None: gdw += ' -tr {} {}'.format(inc, inc)
            out, status = run_cmd(gdw, verbose = verbose)
            utils.remove_glob(chunk)
        else: tif_chunk = chunk

        with open(xyz_chunk, 'w') as xyz_c:
            gdal_dump_entry([tif_chunk, 200, None], dst_port = xyz_c, verbose = verbose)
        utils.remove_glob(tif_chunk)

        if vdatum is not None:
            out, status = run_vdatum(xyz_chunk, vdc)
            utils.remove_glob(xyz_chunk)
            xyz_chunk = os.path.join('result', os.path.basename(xyz_chunk)) 
            os.rename(xyz_chunk, xyz_chunk_final)
            vdatum_clean_result()

            if verbose: utils.echo_msg('transformed {} chunk to {}'.format(iv, ov))
        else: os.rename(xyz_chunk, xyz_chunk_final)

        if datalist is not None:
            mbsfun.mb_inf(xyz_chunk_final)
            datalist_append_entry([os.path.basename(xyz_chunk_final), 168, 1], datalist)
            if verbose: utils.echo_msg('appended xyz chunk {} to datalist {}'.format(xyz_chunk_final, datalist))
    
def gdal_inf(src_ds, warp = None, overwrite = False):
    '''generate an info (.inf) file from a src_gdal file using gdal

    returns the region [xmin, xmax, ymin, ymax] of src_ds'''
    
    utils.echo_msg('generating inf file for {}'.format(src_ds.GetDescription()))
    gdali = {}
    gdali['name'] = src_ds.GetDescription()
    
    minmax = gdal_region(src_ds, warp)
    try: zr = src_ds.GetRasterBand(1).ComputeRasterMinMax()
    except: zr = [None, None]
    minmax = minmax + list(zr)
    gdali['minmax'] = minmax

    ds_band = src_ds.GetRasterBand(1)
    ds_array = ds_band.ReadAsArray()
    ds_config = gdal_gather_infos(src_ds)
    ds_config['dtn'] = 'Int32'
    ndv = ds_band.GetNoDataValue()
    ds_array[ds_array != ndv] = 1
    ds_array[ds_array == ndv] = 0

    gdali['numpts'] = ds_config['nb']
    
    driver = gdal.GetDriverByName('MEM')
    ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    ds.SetGeoTransform(ds_config['geoT'])
    ds.SetProjection(ds_config['proj'])
    ds_band = ds.GetRasterBand(1)
    ds_band.SetNoDataValue(ds_config['ndv'])
    ds_band.WriteArray(ds_array)
    
    tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
    tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    
    gdal.Polygonize(ds_band, None, tmp_layer, 0, callback = _gdal_progress)

    #print(len(tmp_layer))
    feat = tmp_layer.GetFeature(0)
    geom = feat.GetGeometryRef()
    wkt = geom.ExportToWkt()
    gdali['wkt'] = wkt
    tmp_ds = ds = None
    
    with open('{}.inf'.format(src_ds.GetDescription()), 'w') as inf:
        inf.write(json.dumps(gdali))
        
    return(gdali)
                    
def gdal_inf_entry(entry, warp = None):
    ''' scan a gdal entry and find it's region

    returns the region [xmin, xmax, ymin, ymax] of the gdal entry'''

    try:
        ds = gdal.Open(entry[0])
    except: ds = None
    if ds is not None:
        minmax = gdal_inf(ds, warp)
        ds = None
    else: minmax = None
    return(minmax)

def gdal_yield_entry(entry, region = None, verbose = False, epsg = None, z_region = None):
    '''yield the xyz data from the datalist entry.

    yields [x, y, z, <w, ...>]'''

    try:
        ds = gdal.Open(entry[0])
    except: ds = None
    if ds is not None:
        if region is not None:
            srcwin = gdal_srcwin(ds, gdal_region_warp(region, s_warp = epsg, t_warp = gdal_getEPSG(ds)))
        else: srcwin = None
        for xyz in gdal_parse(ds, dump_nodata = False, srcwin = srcwin, warp = epsg, verbose = verbose, z_region = z_region):
            yield(xyz + [entry[2]] if entry[2] is not None else xyz)
        ds = None
        
### End
