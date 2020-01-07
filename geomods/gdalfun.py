### gdalfun.py
##
## Copyright (c) 2012 - 2020 Matthew Love <matthew.love@colorado.edu>
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

import numpy as np
import osgeo.ogr as ogr
import osgeo.gdal as gdal
from gdalconst import *

GDAL_OPTS = ["COMPRESS=LZW", "INTERLEAVE=PIXEL", "TILED=YES",\
        "SPARSE_OK=TRUE", "BIGTIFF=YES" ]

def _ogr_get_fields(src_ogr):
    schema = []
    source = ogr.Open(src_ogr)
    layer = source.GetLayer()
    ldefn = layer.GetLayerDefn()

    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        schema.append(fdefn.name)
    
    return schema

def _gather_infos(src_ds):
    '''Gather information from `src_ds` GDAL dataset.'''

    ds_config = {}

    ds_config['nx'] = src_ds.RasterXSize
    ds_config['ny'] = src_ds.RasterYSize
    ds_config['nb'] = src_ds.RasterCount
    ds_config['geoT'] = src_ds.GetGeoTransform()
    ds_config['proj'] = src_ds.GetProjectionRef()
    ds_config['dt'] = src_ds.GetRasterBand(1).DataType
    ds_config['ndv'] = src_ds.GetRasterBand(1).GetNoDataValue()

    return ds_config

def _con_dec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''

    if x is None:
        return

    fstr = "%." + str( dec ) + "f"
    return fstr % x

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''Convert a geographic x,y value to a pixel location of geoTransform'''

    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else:
        pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt( geoTransform))

    return int(pixel_x), int(pixel_y)

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''Convert a pixel location to geographic coordinates given geoTransform'''

    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)

    return geo_x, geo_y

def _apply_gt(in_x, in_y, geoTransform):
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]

    return out_x, out_y

def _invert_gt(geoTransform):
    det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
    if abs( det ) < 0.000000000000001:
        return

    invDet = 1.0 / det

    # compute adjoint and divide by determinate
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet

    return outGeoTransform 

def chunks(src_fn, n_chunk = 10):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.'''

    band_nums = []
    status = 0
    
    if band_nums == []: band_nums = [1]

    src_ds = gdal.Open(src_fn)
    ds_config = _gather_infos(src_ds)

    bands = []
    for band_num in band_nums: 
        band = src_ds.GetRasterBand(band_num)
        if band is None:
            status = -1
            return status
        bands.append(band)

    gt = ds_config['geoT']

    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk

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

            srcwin = (this_x_origin, this_y_origin, this_x_size, this_y_size)

            # this_geo_x_origin = gt[0] + this_x_origin * gt[1] + this_y_origin * gt[2]
            # this_geo_y_origin = gt[3] + this_x_origin * gt[4] + this_y_origin * gt[5]

            this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)

            dst_gt = [this_geo_x_origin, gt[1], 0.0, this_geo_y_origin, 0.0, gt[5]]

            for band in bands:
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])    
                
                dst_fn = '%s_chnk%sx%s.tif' %(src_fn.split('.')[0], x_i_chunk, i_chunk)

                ods = gdal.GetDriverByName('GTiff').Create(dst_fn, this_x_size, this_y_size, 1, ds_config['dt'])
                ods.SetGeoTransform(dst_gt)
                ods.SetProjection(ds_config['proj'])
                ods.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
                ods.GetRasterBand(1).WriteArray(band_data)

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

def crop(src_fn):
    '''Crop `src_fn` GDAL file by it's NoData value.
    Returns cropped array.'''

    src_ds = gdal.Open(src_fn)

    if srcds is None:
        status = -1
        return status

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
    
    return out_array, ds_config

def dump(src_gdal, dst_xyz):
    '''Dump `src_gdal` GDAL file to ASCII XYZ
    Function taken from GDAL's `gdal2xyz.py` script.'''

    status = 0
    band_nums = []
    nodata = ['-9999', 'nan']
    srcwin = None
    delim = ' '
    skip = 1

    if band_nums == []: band_nums = [1]

    srcds = gdal.Open(src_gdal)

    if srcds is None:
        status = -1
        return status

    bands = []
    for band_num in band_nums: 
        band = srcds.GetRasterBand(band_num)
        if band is None:
            status = -1
            return status
        bands.append(band)

    gt = srcds.GetGeoTransform()
  
    if srcwin is None:
        srcwin = (0,0,srcds.RasterXSize,srcds.RasterYSize)

    if dst_xyz is not None:
        dst_fh = open(dst_xyz,'wt')
    else:
        dst_fh = sys.stdout

    band_format = (("%g" + delim) * len(bands)).rstrip(delim) + '\n'

    if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
       and abs(srcds.RasterXSize * gt[1]) < 180 \
       and abs(srcds.RasterYSize * gt[5]) < 180:
        format = '%.10g' + delim + '%.10g' + delim + '%s'
    else:
        format = '%.3f' + delim + '%.3f' + delim + '%s'

    for y in range(srcwin[1], srcwin[1] + srcwin[3], skip):

        data = []
        for band in bands:
            nodata.append(band.GetNoDataValue())
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

            if band_str not in nodata:
                dst_fh.write(line)

def proximity(self, src_fn, dst_fn):
    '''Compute a proximity grid via GDAL'''

    dst_ds = None
    prog_func = None
    dst_nodata = -9999

    src_ds = gdal.Open(src_fn)
    
    if src_ds is None:
        status = -1

    src_band = src_ds.GetRasterBand(1)

    if dst_ds is None:
        drv = gdal.GetDriverByName('GTiff')
        dst_ds = drv.Create(dst_fn, src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GetDataTypeByName('Float32'), [])

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())
    
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(dst_nodata)

    gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)

    dst_band = src_band = dst_ds = src_ds = None

def query(src_xyz, src_grd, out_form):
    '''Query a gdal-compatible grid file with xyz data.'''

    xyzl = []

    # Process the src grid file
    ds = gdal.Open(src_grd)
    dsgeot = ds.GetGeoTransform()
    dsband = ds.GetRasterBand(1)
    dsnodata = dsband.GetNoDataValue()

    dsdt = gdal.GetDataTypeName(dsband.DataType)

    # Load the src grid into a numpy array
    tgrid = dsband.ReadAsArray()

    cellsize = [float(dsgeot[1]), float(dsgeot[5])]
    xextent = float(dsgeot[0])
    yextent = float(dsgeot[3])
        
    # Process the src xyz data
    for xyz in src_xyz:
        x = xyz[0]
        y = xyz[1]
        try: 
            z = xyz[2]
        except: z = dsnodata

        # Continue if values are reasonable.
        if x > xextent and y < yextent:
            xpos,ypos = _geo2pixel(x, y, dsgeot)

            # Locate the grid cell and get it's value
            try: 
                g = tgrid[ypos,xpos]
            except: g = dsnodata

            d = c = m = s = dsnodata
                
            if g != dsnodata:
                d = z - g
                m = z + g
                c = _con_dec( math.fabs(float(d / (g + 0.00000001) * 100)), 2)
                s = _con_dec( math.fabs(d / (z + (g + 0.00000001))), 4)
                    
            d = _con_dec(d, 4)

            outs = []
            for i in out_form:
                outs.append(vars()[i])
            xyzl.append(np.array( outs, dtype = dsdt))

    dsband = ds = None
    return np.array(xyzl, dtype = dsdt)

### End
