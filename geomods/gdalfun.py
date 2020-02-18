### gdalfun.py
##
## Copyright (c) 2012 - 2020 CIRES Coastal DEM Team
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
import numpy as np
import ogr
import gdal
import osr
from gdalconst import *

GDAL_OPTS = ["COMPRESS=LZW", "INTERLEAVE=PIXEL", "TILED=YES",\
        "SPARSE_OK=TRUE", "BIGTIFF=YES" ]

_known_delims = [',', ' ', '\t', '/', ':']

def xyz_parse(src_xyz):
    '''xyz file parsing generator'''
    
    for xyz in src_xyz:
        this_line = xyz.strip()

        for delim in _known_delims:
            this_xyz = this_line.split(delim)
            if len(this_xyz) > 1:
                this_delim = delim
                break

        yield this_xyz

def _ogr_create_polygon(coords):
    '''convert coords to Wkt'''

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[1], coord[0])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    poly_wkt = poly.ExportToWkt()
    poly = None
    
    return(poly_wkt)

def bounds2geom(bounds):
    '''convert a bounds [west, east, south, north] to an 
    OGR geometry'''

    b1 = [[bounds[2], bounds[0]],
          [bounds[2], bounds[1]],
          [bounds[3], bounds[1]],
          [bounds[3], bounds[0]],
          [bounds[2], bounds[0]]]

    geom = ogr.CreateGeometryFromWkt(_ogr_create_polygon(b1))
    
    return(geom)

def _ogr_get_fields(src_ogr):
    '''return all fields in src_ogr'''

    schema = []
    source = ogr.Open(src_ogr)
    layer = source.GetLayer()
    ldefn = layer.GetLayerDefn()

    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        schema.append(fdefn.name)
    
    return(schema)

def _ogr_get_layer_fields(src_lyr):
    '''return all fields in src_lyr'''

    schema = []
    if src_lyr is not None:
        ldefn = src_lyr.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
    
    return(schema)

def ogr_mask_union(src_layer, src_field, dst_defn = None, callback = lambda: False):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.'''

    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()

    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for f in src_layer:
        if not callback():
            if f.GetField(src_field) == 0:
                src_layer.DeleteFeature(f.GetFID())
            elif f.GetField(src_field) == 1:
                f.geometry().CloseRings()
                wkt = f.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
                src_layer.DeleteFeature(f.GetFID())
                        
    union = multi.UnionCascaded()
    out_feat = ogr.Feature(dst_defn)
    out_feat.SetGeometry(union)
    union = multi = None

    return(out_feat)

def _ogr_extents(src_ds):
    '''return the extent(s) of the ogr dataset'''
    these_extents = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                these_extents.append(pgeom.GetEnvelope())
                
    return(these_extents)


def _gather_infos(src_ds):
    '''Gather information from `src_ds` GDAL dataset.'''

    ds_config = {}

    ds_config['nx'] = src_ds.RasterXSize
    ds_config['ny'] = src_ds.RasterYSize
    ds_config['nb'] = src_ds.RasterCount
    ds_config['geoT'] = src_ds.GetGeoTransform()
    ds_config['proj'] = src_ds.GetProjectionRef()
    ds_config['dt'] = src_ds.GetRasterBand(1).DataType
    ds_config['dtn'] = gdal.GetDataTypeName(src_ds.GetRasterBand(1).DataType)
    ds_config['ndv'] = src_ds.GetRasterBand(1).GetNoDataValue()
    ds_config['fmt'] = src_ds.GetDriver().ShortName

    if ds_config['ndv'] is None:
        ds_config['ndv'] = -9999
    
    return(ds_config)

def _set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''Set a datasource config dictionary'''
    
    ds_config = {}

    ds_config['nx'] = nx
    ds_config['ny'] = ny
    ds_config['nb'] = nb
    ds_config['geoT'] = geoT
    ds_config['proj'] = proj
    ds_config['dt'] = dt
    ds_config['ndv'] = ndv
    ds_config['fmt'] = fmt

    return(ds_config)

def _cpy_infos(src_config):
    dst_config = {}
    for dsc in src_config.keys():
        dst_config[dsc] = src_config[dsc]
    return(dst_config)

def _con_dec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''

    if x is None:
        return

    fstr = "%." + str( dec ) + "f"
    return(fstr % x)

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''Convert a geographic x,y value to a pixel location of geoTransform'''

    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else:
        pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt( geoTransform))

    return(int(pixel_x), int(pixel_y))

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''Convert a pixel location to geographic coordinates given geoTransform'''

    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)

    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs( det ) < 0.000000000000001:
        return

    invDet = 1.0 / det

    ## ==============================================
    ## compute adjoint and divide by determinate
    ## ==============================================

    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet

    return(outGeoTransform)

def _gt2extent(ds_config):
    geoT = ds_config['geoT']
    
    x_origin = geoT[0]
    y_origin = geoT[3]
    x_inc = geoT[1]
    y_inc = geoT[5]

    return([x_origin, x_origin + x_inc * ds_config['nx'], y_origin + y_inc * ds_config['ny'], y_origin])

def _extent(src_fn):
    extent = None
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)
        extent = _gt2extent(ds_config)
        
    return(extent)

def sr_wkt(epsg, esri = False):
    '''convert an epsg code to wkt'''
    
    wkt = None
    try:
        int(epsg)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        if esri: sr.MorphToESRI()
        
        wkt = sr.ExportToWkt()
    except: sys.stderr.write('geomods: error, invalid epsg code\n')

    sr = None
    return(wkt)

def _prj_file(dst_fn, epsg):
    '''generate a .prj file given an epsg code'''
    
    with open(dst_fn, 'w') as out:
        out.write(sr_wkt(epsg, True))

def _write_gdal(src_arr, dst_gdal, ds_config, dst_fmt = 'GTiff', verbose = False):
    '''write src_arr Array to gdal file dst_gdal using src_config'''

    if verbose: sys.stderr.write('geomods: writing gdal grid...')
    
    driver = gdal.GetDriverByName(ds_config['fmt'])
    if os.path.exists(dst_gdal):
        driver.Delete(dst_gdal)
    
    ds = driver.Create(dst_gdal, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    ds.SetGeoTransform(ds_config['geoT'])
    ds.SetProjection(ds_config['proj'])
    ds.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
    ds.GetRasterBand(1).WriteArray(src_arr)
    ds = None

    if verbose: sys.stderr.write('ok\n')
    
## =============================================================================
##
## GDAL Functions for external use
##
## TODO: Allow for specified output format
##       Allow for export to XYZ
##       Allow for return of chunks list
##
## =============================================================================

def chunks(src_fn, n_chunk = 10, verbose = False):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.'''

    band_nums = []
    o_chunks = []
    
    if band_nums == []: band_nums = [1]

    src_ds = gdal.Open(src_fn)

    if src_ds is not None:
        ds_config = _gather_infos(src_ds)

        #bands = []
        #for band_num in band_nums:
        band = src_ds.GetRasterBand(1)
        #if band is not None:
        #bands.append(band)
                
        gt = ds_config['geoT']

        i_chunk = 0
        x_i_chunk = 0
        x_chunk = n_chunk

        if verbose: sys.stderr.write('geomods: chunking grid...')
        while True:
            if verbose: sys.stderr.write('.')
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
                if this_x_size == 0 or this_y_size == 0:
                    break
                
                srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)
                this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)
                dst_gt = [this_geo_x_origin, gt[1], 0.0, this_geo_y_origin, 0.0, gt[5]]
                
                #for band in bands:
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])
                if not np.all(band_data == band_data[0,:]):
                    o_chunk = '{}_chnk{}x{}.tif'.format(os.path.basename(src_fn).split('.')[0], x_i_chunk, i_chunk)
                    dst_fn = os.path.join(os.path.dirname(src_fn), o_chunk)
                    o_chunks.append(dst_fn)

                    dst_config = _cpy_infos(ds_config)
                    dst_config['nx'] = this_x_size
                    dst_config['ny'] = this_y_size
                    dst_config['geoT'] = dst_gt
                    
                    _write_gdal(band_data, dst_fn, dst_config)

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
                
        if verbose: sys.stderr.write('ok\n')
        src_ds = None
        return(o_chunks)

def crop(src_fn):
    '''Crop `src_fn` GDAL file by it's NoData value.
    Returns cropped array.'''
    
    src_ds = gdal.Open(src_fn)

    if srcds is None:
        return(None)

    ds_config = _gather_infos(src_ds)
    ds_arr = src_ds.GetRasterBand(1).ReadAsArray()

    src_ds = None
    
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

def dump(src_gdal, dst_xyz, dump_nodata = False, srcwin = None):
    '''Dump `src_gdal` GDAL file to ASCII XYZ'''

    status = 0
    band_nums = []
    delim = ' '
    skip = 1

    if band_nums == []: band_nums = [1]

    src_ds = gdal.Open(src_gdal)
    
    if src_ds is not None:
        bands = []
        for band_num in band_nums: 
            band = src_ds.GetRasterBand(band_num)
            if band is not None:
                bands.append(band)

        ds_config = _gather_infos(src_ds)
        gt = ds_config['geoT']

        if srcwin is None:
            srcwin = (0, 0, ds_config['nx'], ds_config['ny'])

        # fixme
        if dst_xyz is not None:
            try:
                dst_fh = open(dst_xyz, 'wt')
            except: dst_fh = dst_xyz
        else: dst_fh = sys.stdout

        band_format = (("%g" + delim) * len(bands)).rstrip(delim) + '\n'

        if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
           and abs(ds_config['nx'] * gt[1]) < 180 \
           and abs(ds_config['ny'] * gt[5]) < 180:
            format = '%.10g' + delim + '%.10g' + delim + '%s'
        else: format = '%.3f' + delim + '%.3f' + delim + '%s'

        for y in range(srcwin[1], srcwin[1] + srcwin[3], skip):
            nodata = ['-9999', 'nan']
            data = []
            for band in bands:
                if band.GetNoDataValue() is not None:
                    nodata.append(band_format % band.GetNoDataValue())
                band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)    
                band_data = np.reshape(band_data, (srcwin[2], ))
                data.append(band_data)

            for x_i in range(0, srcwin[2], skip):
                x = x_i + srcwin[0]

                geo_x = gt[0] + (x + 0.5) * gt[1] + (y + 0.5) * gt[2]
                geo_y = gt[3] + (x + 0.5) * gt[4] + (y + 0.5) * gt[5]

                x_i_data = []
                for i in range(len(bands)):
                    x_i_data.append(data[i][x_i])

                band_str = band_format % tuple(x_i_data)
                
                line = format % (float(geo_x), float(geo_y), band_str)

                if dump_nodata:
                    dst_fh.write(line)
                else:
                    if band_str not in nodata:
                        dst_fh.write(line)

        try:
            dst_fh.close()
        except: pass
        srcds = None

def dumpy(src_gdal, dump_nodata = False, srcwin = None):
    '''Dump `src_gdal` GDAL file to ASCII XYZ
    Function taken from GDAL's `gdal2xyz.py` script.'''

    delim = ' '
    skip = 1

    srcds = gdal.Open(src_gdal)

    if srcds is not None:
        band = srcds.GetRasterBand(1)
        gt = srcds.GetGeoTransform()
        if srcwin is None:
            srcwin = (0,0,srcds.RasterXSize,srcds.RasterYSize)

        band_format = (("%g" + delim) * len(bands)).rstrip(delim) + '\n'
        if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
           and abs(srcds.RasterXSize * gt[1]) < 180 \
           and abs(srcds.RasterYSize * gt[5]) < 180:
            format = '%.10g' + delim + '%.10g' + delim + '%s'
        else:
            format = '%.3f' + delim + '%.3f' + delim + '%s'

        for y in range(srcwin[1], srcwin[1] + srcwin[3], skip):

            nodata = ['-9999', 'nan']
            data = []

            if band.GetNoDataValue() is not None:
                nodata.append(band_format % band.GetNoDataValue())
            band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)    
            band_data = np.reshape(band_data, (srcwin[2], ))
            data.append(band_data)

            for x_i in range(0, srcwin[2], skip):
                x = x_i + srcwin[0]
                geo_x = gt[0] + (x + 0.5) * gt[1] + (y + 0.5) * gt[2]
                geo_y = gt[3] + (x + 0.5) * gt[4] + (y + 0.5) * gt[5]
                x_i_data = []
                for i in range(len(bands)):
                    x_i_data.append(data[i][x_i])

                band_str = band_format % tuple(x_i_data)
                line = format % (float(geo_x), float(geo_y), band_str)
                if dump_nodata:
                    yield(line)
                else:
                    if band_str not in nodata:
                        yield(line)

        srcds = None
    
def infos(src_fn):
    ds = gdal.Open(src_fn)
    if ds is not None:
        ds_config = _gather_infos(ds)
        print ds_config
        print _gt2extent(ds_config)
        
def null(dst_fn, extent, cellsize, nodata = -9999, outformat = 'GTiff'):
    '''generate a `null` grid with gdal'''

    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1
    
    nullArray = np.zeros((ycount, xcount))
    nullArray[nullArray == 0] = nodata

    dst_gt = (extent[0], cellsize, 0, extent[3], 0, (cellsize * -1.))
    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, sr_wkt(4326), gdal.GDT_Float32, -9999, outformat)
    _write_gdal(nullArray, dst_fn, ds_config)
    
def percentile(src_fn, perc = 95):

    ds = gdal.Open(src_fn)
    myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
    x_dim=myarray.shape[0]
    myarray[myarray==0]=np.nan

    myarray_flat=myarray.flatten()
    
    #p = np.nanpercentile(myarray_flat, perc)
    #p = np.percentile(myarray_flat, perc)
    myarray_flat = myarray_flat[~np.isnan(myarray_flat)]
    p = np.percentile(myarray_flat, perc)

    if p==p:
        percentile=p
        if percentile < 2:
            percentile = 2
    else: percentile = 1

    ds = myarray = None
    
    return(percentile)

def polygonize(src_gdal, dst_layer, verbose = False):
    '''run gdal.Polygonize on src_gdal and add polygon to dst_layer'''
    
    src_ds = gdal.Open(src_gdal)
    srcband = src_ds.GetRasterBand(1)

    if verbose: sys.stderr.write('geomods: polygonizing grid...')
    gdal.Polygonize(srcband, None, dst_layer, 0, [], None) #lambda x, y, z: sys.stderr.write('.'))
    if verbose: sys.stderr.write('ok\n')
    
    src_ds = srcband = None

def proximity(src_fn, dst_fn):
    '''Compute a proximity grid via GDAL'''

    prog_func = None

    src_ds = gdal.Open(src_fn)
    dst_ds = None
    
    if src_ds is not None:
        src_band = src_ds.GetRasterBand(1)
        ds_config = _gather_infos(src_ds)
        
        if dst_ds is None:
            drv = gdal.GetDriverByName('GTiff')
            dst_ds = drv.Create(dst_fn, ds_config['nx'], ds_config['ny'], 1, ds_config['dt'], [])

        dst_ds.SetGeoTransform(ds_config['geoT'])
        dst_ds.SetProjection(ds_config['proj'])

        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(ds_config['ndv'])

        gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)

        dst_band = src_band = dst_ds = src_ds = None

def query(src_xyz, src_grd, out_form):
    '''Query a gdal-compatible grid file with xyz data.'''

    xyzl = []
    out_array = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================

    ds = gdal.Open(src_grd)
    if ds is not None:
        ds_config = _gather_infos(ds)

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
            except: z = dsnodata

            if x > ds_gt[0] and y < float(ds_gt[3]):
                xpos, ypos = _geo2pixel(x, y, ds_gt)

                try: 
                    g = tgrid[ypos, xpos]
                except: g = ds_nd

                d = c = m = s = ds_nd
                if g != ds_nd:
                    d = z - g
                    m = z + g
                    c = _con_dec(math.fabs(float(d / (g + 0.00000001) * 100)), 2)
                    s = _con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    d = _con_dec(d, 4)

                    outs = []
                    for i in out_form:
                        outs.append(vars()[i])
                    xyzl.append(np.array(outs, dtype = ds_config['dtn']))
                    #xyzl.append(np.array(outs))
        dsband = ds = None
        out_array = np.array(xyzl, dtype = ds_config['dtn'])
        #out_array = np.array(xyzl)

    return(out_array)
    
def split(src_gdal, split_value = 0):
    '''split raster into two peices based on z value'''

    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_upper.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_lower.tif'.format(os.path.basename(src_gdal)[:-4]))
                                         
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = _gather_infos(src_ds)
        
        dst_config = _cpy_infos(src_config)
        dst_config['fmt'] = 'GTiff'

        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny']) 

        upper_array = ds_arr
        upper_array[upper_array <= split_value] = src_config['ndv'] 
        _write_gdal(upper_array, dst_upper, dst_config)
        upper_array = None

        lower_array = ds_arr
        lower_array[lower_array >= split_value] = src_config['ndv']
        _write_gdal(lower_array, dst_lower, dst_config)
        lower_array = src_ds = None

        return([dst_upper, dst_lower])

def sum(src_gdal):
    elev = gdal.Open(src_gdal)
    
    elev_array = elev.GetRasterBand(1).ReadAsArray() 
    sums = np.sum(elev_array)
    
    elev = elev_array = None
    return(sums)

def transform(src_gdal):
    ##Getting spatial reference of input raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)

    # WGS84 projection reference
    OSR_WGS84_REF = osr.SpatialReference()
    OSR_WGS84_REF.ImportFromEPSG(4326)
    
    # OSR transformation
    wgs84_to_image_trasformation = osr.CoordinateTransformation(OSR_WGS84_REF, srs)
    
    #Create geometry. Finally:
        
    wkt_geom.Transform(wgs84_to_image_trasformation)
        
def xyz2gdal(src_xyz, dst_gdal, extent, cellsize,
             dst_format='GTiff', zvalue='d', xloc=0, yloc=1, zloc=2, 
             delim=' ', verbose=False, overwrite=False):
    '''Create a GDAL supported grid from xyz data'''

    dst_nodata=-9999

    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1
    dst_gt = (extent[0], cellsize,0, extent[3], 0, (cellsize * -1.))
    
    if zvalue == 'z': sumArray = np.zeros((ycount, xcount))
    ptArray = np.zeros((ycount, xcount))
    
    pointnum = 0

    if verbose: sys.stderr.write('geomods: processing xyz data...')

    for this_xyz in xyz_parse(src_xyz):
        x = float(this_xyz[xloc])
        y = float(this_xyz[yloc])
        z = float(this_xyz[int(zloc)].strip())

        if x > extent[0] and x < extent[1]:
            if y > extent[2] and y < extent[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                ptArray[ypos, xpos] = 1

                if zvalue == 'z': sumArray[ypos, xpos] += float(z)
                if zvalue == 'd' or zvalue == 'z': 
                    ptArray[ypos, xpos] += 1
                else: ptArray[ypos, xpos] = 1
        pointnum += 1

    if zvalue == 'z':
        outarray = sumArray / ptArray
    elif zvalue == 'd': outarray = ptArray
    else: outarray = ptArray

    outarray[np.isnan(outarray)] = dst_nodata
    
    if verbose: sys.stderr.write('ok\n')

    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, sr_wkt(4326), gdal.GDT_Int32, dst_nodata, 'GTiff')
    _write_gdal(outarray, dst_gdal, ds_config)
    
    outarray = None 

def xyz2ogr(src_xyz, dst_ogr, xloc = 0, yloc = 1, zloc = 2, dst_fmt = 'ESRI Shapefile', overwrite = True, verbose = False):

    driver = ogr.GetDriverByName(dst_fmt)
    if os.path.exists(dst_ogr):
        driver.DeleteDataSource(dst_ogr)
    ds = driver.CreateDataSource(dst_ogr)

    layer = ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPoint25D)

    oft_title = "Longitude"
    oft_string = ogr.OFTReal
    fdw = 12
    fdp = 8
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)

    oft_title = "Latitude"
    oft_string = ogr.OFTReal
    fdw = 12
    fdp = 8
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)

    oft_title = "Elevation"
    oft_string = ogr.OFTReal
    fdw = 12
    fdp = 8
    fd = ogr.FieldDefn(oft_title, oft_string)
    fd.SetWidth(fdw)
    fd.SetPrecision(fdp)
    layer.CreateField(fd)
    
    f = ogr.Feature(feature_def = layer.GetLayerDefn())
    
    for this_xyz in xyz_parse(src_xyz):
        x = float(this_xyz[xloc])
        y = float(this_xyz[yloc])
        z = float(this_xyz[zloc])

        f.SetField(0, x)
        f.SetField(1, y)
        f.SetField(2, z)

        wkt = 'POINT(%f %f %f)' % (x,y,z)
        g = ogr.CreateGeometryFromWkt(wkt)
        f.SetGeometryDirectly(g)
        layer.CreateFeature(f)
        
## add option for grid/pixel node registration...currently pixel node only.
def xyz_mask(src_xyz, dst_gdal, extent, cellsize,
             dst_format='GTiff', xloc=0, yloc=1, zloc=2, 
             delim=' ', verbose=False):
    '''Create a num grid mask of xyz data. The output grid
    will contain 1 where data exists and 0 where no data exists.'''
    
    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1
    dst_gt = (extent[0], cellsize, 0, extent[3], 0, (cellsize * -1.))
    
    ptArray = np.zeros((ycount, xcount))

    if verbose: sys.stderr.write('geomods: processing xyz data...')

    for this_xyz in xyz_parse(src_xyz):
        x = float(this_xyz[xloc])
        y = float(this_xyz[yloc])

        if x > extent[0] and x < extent[1]:
            if y > extent[2] and y < extent[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                ptArray[ypos, xpos] = 1

    if verbose: sys.stderr.write('ok\n')
        
    ds_config = _set_infos(xcount, ycount, xcount * ycount, dst_gt, sr_wkt(4326), gdal.GDT_Int32, -9999, 'GTiff')
    _write_gdal(ptArray, dst_gdal, ds_config, verbose)
    
    ptArray = None

### End
