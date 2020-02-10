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

def ogr_mask_union(src_layer, src_field, dst_defn = None):
    '''`union` a `src_layer`'s features based on `src_field` where
    `src_field` holds a value of 0 or 1. optionally, specify
    an output layer defn for the unioned feature.'''

    if dst_defn is None:
        dst_defn = src_layer.GetLayerDefn()

    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for f in src_layer:
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
    ds_config['fmt'] = src_ds.GetDriver().ShortName

    if ds_config['ndv'] is None:
        ds_config['ndv'] = -9999
    
    return(ds_config)

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

def _prj_file(dst_fn, epsg):
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(epsg)
    
    sr.MorphToESRI()
    sr_f = open(dst_fn, 'w')
    sr_f.write(sr.ExportToWkt())
    sr_f.close()

    sr = None

## =============================================================================
##
## GDAL Functions for external use
##
## TODO: Allow for specified output format
##       Allow for export to XYZ
##       Allow for return of chunks list
##
## =============================================================================

def chunks(src_fn, n_chunk = 10):
    '''split `src_fn` GDAL file into chunks with `n_chunk` cells squared.'''

    band_nums = []
    o_chunks = []
    status = 0
    
    if band_nums == []: band_nums = [1]

    src_ds = gdal.Open(src_fn)
    if src_ds is None:
        print('Error, unable to load file {}'.format(src_fn))
        return(o_chunks)
              
    ds_config = _gather_infos(src_ds)

    bands = []
    for band_num in band_nums: 
        band = src_ds.GetRasterBand(band_num)
        if band is None:
            status = -1
            return(status)
        bands.append(band)

    gt = ds_config['geoT']

    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk

    ## ==============================================
    ## parse through the grid and extract chunks
    ## ==============================================

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
            this_geo_x_origin, this_geo_y_origin = _pixel2geo(this_x_origin, this_y_origin, gt)
            dst_gt = [this_geo_x_origin, gt[1], 0.0, this_geo_y_origin, 0.0, gt[5]]

            for band in bands:
                band_data = band.ReadAsArray(srcwin[0], srcwin[1], srcwin[2], srcwin[3])    
                if not np.all(band_data == band_data[0,:]):
                    o_chunk = '{}_chnk{}x{}.tif'.format(os.path.basename(src_fn).split('.')[0], x_i_chunk, i_chunk)
                    dst_fn = os.path.join(os.path.dirname(src_fn), o_chunk)
                                          
                    o_chunks.append(dst_fn)

                    ods = gdal.GetDriverByName('GTiff').Create(dst_fn,
                                                               this_x_size,
                                                               this_y_size,
                                                               1,
                                                               ds_config['dt'])
                    ods.SetGeoTransform(dst_gt)
                    ods.SetProjection(ds_config['proj'])
                    ods.GetRasterBand(1).SetNoDataValue(ds_config['ndv'])
                    ods.GetRasterBand(1).WriteArray(band_data)

                    ods = None
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

def dump(src_gdal, dst_xyz, dump_nodata = False):
    '''Dump `src_gdal` GDAL file to ASCII XYZ
    Function taken from GDAL's `gdal2xyz.py` script.'''

    status = 0
    band_nums = []
    srcwin = None
    delim = ' '
    skip = 1

    if band_nums == []: band_nums = [1]

    srcds = gdal.Open(src_gdal)

    if srcds is None:
        status = -1
        return(status)

    bands = []
    for band_num in band_nums: 
        band = srcds.GetRasterBand(band_num)
        if band is None:
            status = -1
            return(status)
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

    srcds = None

def null(outfile, extent, cellsize, nodata = -9999, outformat = 'GTiff'):
    '''generate a `null` grid with gdal'''

    ## ==============================================    
    ## Set the rows and columns for the output grid
    ## ==============================================

    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1

    ## ==============================================
    ## Create the output GDAL Raster
    ## ==============================================

    if outformat == "AAIGrid":
        driver = gdal.GetDriverByName("MEM")
    else: driver = gdal.GetDriverByName(outformat)

    if os.path.exists(outfile):
        driver.Delete(outfile)

    dst_ds = driver.Create(outfile, xcount, ycount, 1, gdal.GDT_Float32)
    if dst_ds is None: sys.exit("Error: failed to open output file...%s" %(outfile))

    gt = (extent[0], cellsize, 0, extent[3], 0, (cellsize * -1.))
    dst_ds.SetGeoTransform(gt)
    dst_band = dst_ds.GetRasterBand(1)

    dst_band.SetNoDataValue(nodata)

    ## ==============================================
    ## Create a Numpy Array filled with 0's
    ## ==============================================

    nullArray = np.zeros((ycount, xcount))
    nullArray[nullArray == 0] = nodata

    ## ==============================================
    ## Write Numpy Array to the GDAL Raster Band
    ## ==============================================

    dst_band.WriteArray(nullArray)

    if outformat == "AAIGrid":
        driver = gdal.GetDriverByName(outformat)
        dst_ds_aai = driver.CreateCopy(outfile, dst_ds)
    
    ## ==============================================
    ## Clear the GDAL Raster(s)x
    ## ==============================================

    dst_ds = dst_ds_aai = None

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

def proximity(src_fn, dst_fn):
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
        dst_ds = drv.Create(dst_fn,
                            src_ds.RasterXSize,
                            src_ds.RasterYSize,
                            1,
                            gdal.GetDataTypeByName('Float32'),
                            [])

    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjectionRef())
    
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(dst_nodata)

    gdal.ComputeProximity(src_band, dst_band, ['DISTUNITS=PIXEL'], callback = prog_func)

    dst_band = src_band = dst_ds = src_ds = None

def query(src_xyz, src_grd, out_form):
    '''Query a gdal-compatible grid file with xyz data.'''

    xyzl = []

    ## ==============================================
    ## Process the src grid file
    ## ==============================================

    ds = gdal.Open(src_grd)
    dsgeot = ds.GetGeoTransform()
    dsband = ds.GetRasterBand(1)
    dsnodata = dsband.GetNoDataValue()

    dsdt = gdal.GetDataTypeName(dsband.DataType)

    ## ==============================================
    ## Load the src grid into a numpy array
    ## ==============================================

    tgrid = dsband.ReadAsArray()

    cellsize = [float(dsgeot[1]), float(dsgeot[5])]
    xextent = float(dsgeot[0])
    yextent = float(dsgeot[3])
     
    ## ==============================================   
    ## Process the src xyz data
    ## ==============================================

    for xyz in src_xyz:
        x = xyz[0]
        y = xyz[1]
        try: 
            z = xyz[2]
        except: z = dsnodata

        ## ==============================================
        ## Continue if values are reasonable.
        ## ==============================================

        if x > xextent and y < yextent:
            xpos,ypos = _geo2pixel(x, y, dsgeot)

            ## ==============================================
            ## Locate the grid cell and get it's value
            ## ==============================================

            try: 
                g = tgrid[ypos,xpos]
            except: g = dsnodata

            d = c = m = s = dsnodata
            if g != dsnodata:
                d = z - g
                m = z + g
                c = _con_dec(math.fabs(float(d / (g + 0.00000001) * 100)), 2)
                s = _con_dec(math.fabs(d / (z + (g + 0.00000001))), 4)
                    
                d = _con_dec(d, 4)

                outs = []
                for i in out_form:
                    outs.append(vars()[i])
                xyzl.append(np.array(outs, dtype = dsdt))

    dsband = ds = None

    return(np.array(xyzl, dtype = dsdt))

def scan(src_gdal, fail_if_nodata = False):
    '''Dump `src_gdal` GDAL file to ASCII XYZ
    Function taken from GDAL's `gdal2xyz.py` script.'''

    status = 0
    band_nums = []
    srcwin = None
    skip = 1
    
    if band_nums == []: band_nums = [1]

    srcds = gdal.Open(src_gdal)

    if srcds is None:
        status = -1
        return(status)

    bands = []
    for band_num in band_nums: 
        band = srcds.GetRasterBand(band_num)
        if band is None:
            status = -1
            return(status)
        bands.append(band)

    gt = srcds.GetGeoTransform()
  
    if srcwin is None:
        srcwin = (0,0,srcds.RasterXSize,srcds.RasterYSize)

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

            # geo_x = gt[0] + (x + 0.5) * gt[1] + (y + 0.5) * gt[2]
            # geo_y = gt[3] + (x + 0.5) * gt[4] + (y + 0.5) * gt[5]
            
            # x_i_data = []
            # for i in range(len(bands)):
            #     x_i_data.append(data[i][x_i])
            
    srcds = None

def _write_gdal(src_arr, dst_gdal, src_config, dst_fmt = 'GTiff'):

    UDataSet = gdal.GetDriverByName(dst_fmt).Create(dst_gdal, src_config['nx'], src_config['ny'], 1, src_config['dt'])
    UDataSet.SetGeoTransform(src_config['geoT'])
    UDataSet.SetProjection(src_config['proj'])
    UDataSet.GetRasterBand(1).SetNoDataValue(src_config['ndv'])
    UDataSet.GetRasterBand(1).WriteArray(src_arr)
    UDataset = None
    
def split(src_gdal, split_value = 0):
    '''split raster into two peices based on z value'''

    dst_upper = os.path.join(os.path.dirname(src_gdal), '{}_upper.tif'.format(os.path.basename(src_gdal)[:-4]))
    dst_lower = os.path.join(os.path.dirname(src_gdal), '{}_lower.tif'.format(os.path.basename(src_gdal)[:-4]))
                                         
    src_ds = gdal.Open(src_gdal)
    if src_ds is not None:
        src_config = _gather_infos(src_ds)

        ds_arr = src_ds.GetRasterBand(1).ReadAsArray(0, 0, src_config['nx'], src_config['ny']) 

        upper_array = ds_arr
        upper_array[upper_array <= split_value] = src_config['ndv'] 
        _write_gdal(upper_array, dst_upper, src_config)
        upper_array = None

        lower_array = ds_arr
        lower_array[lower_array >= split_value] = src_config['ndv']
        _write_gdal(lower_array, dst_lower, src_config)
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

    ## ==============================================    
    ## Set the rows and columns for the output grid
    ## ==============================================

    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1

    ## ==============================================
    ## Create the output GDAL Raster
    ## ==============================================

    if dst_format == "AAIGrid":
        driver = gdal.GetDriverByName("MEM")
    else: driver = gdal.GetDriverByName(dst_format)

    dst_ds = driver.Create(dst_gdal, xcount, ycount, 1, gdal.GDT_Float32)
    if dst_ds is None: sys.exit("geomods: failed to open output file...%s" %(outfile))

    dst_gt = (extent[0], cellsize,0, extent[3], 0, (cellsize * -1.))
    dst_ds.SetGeoTransform(dst_gt)

    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(dst_nodata)

    ## ==============================================
    ## Process the XYZ data
    ## ==============================================

    if zvalue == 'z': sumArray = np.zeros( (ycount, xcount) )
    ptArray = np.zeros( (ycount, xcount) )
    
    pointnum = 0

    if verbose:
        print('geomods: processing xyz data...')

    for this_xyz in xyz_parse(src_xyz):
        x = float(this_xyz[xloc])
        y = float(this_xyz[yloc])
        z = float(this_xyz[int(zloc)].strip())

        ## ==============================================
        ## Determine which cell to apply this count
        ## ==============================================

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

    if verbose:
        sys.stderr.write('ok\n')

    ## ==============================================
    ## write the output grid
    ## ==============================================

    outarray[np.isnan(outarray)] = dst_nodata

    dst_band.WriteArray(outarray)

    if dst_format == "AAIGrid":
        driver = gdal.GetDriverByName(dst_format)
        dst_ds_aai = driver.CreateCopy(dst_gdal, dst_ds)
    
    dst_ds = dst_ds_aai = None

def xyz_gmask(src_xyz, dst_gdal, extent, cellsize,
              dst_format='GTiff', xloc=0, yloc=1, zloc=2, 
              delim=' ', verbose=False):
    '''Create a num grid mask of xyz data'''
    
    ## ==============================================
    ## Set the rows and columns for the output grid
    ## ==============================================

    ysize = extent[3] - extent[2]
    xsize = extent[1] - extent[0]
    xcount = int(xsize / cellsize) + 1
    ycount = int(ysize / cellsize) + 1

    ## ==============================================
    ## Create the output GDAL Raster
    ## ==============================================

    dst_nodata=-9999

    if dst_format == "AAIGrid":
        driver = gdal.GetDriverByName("MEM")
    else: driver = gdal.GetDriverByName(dst_format)

    dst_ds = driver.Create(dst_gdal, xcount, ycount, 1, gdal.GDT_Int32)
    if dst_ds is None: sys.exit("geomods: failed to open output file...%s" %(outfile))

    dst_gt = (extent[0], cellsize,0, extent[3], 0, (cellsize * -1.))
    dst_ds.SetGeoTransform(dst_gt)

    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(dst_nodata)

    ## ==============================================
    ## Process the XYZ data
    ## ==============================================

    ptArray = np.zeros((ycount, xcount))

    if verbose:
        print('geomods: processing xyz data...')

    for this_xyz in xyz_parse(src_xyz):
        x = float(this_xyz[xloc])
        y = float(this_xyz[yloc])

        if x > extent[0] and x < extent[1]:
            if y > extent[2] and y < extent[3]:
                xpos, ypos = _geo2pixel(x, y, dst_gt)
                ptArray[ypos, xpos] = 1

    if verbose:
        sys.stderr.write('ok\n')

    ## ==============================================
    ## Write out the grid    
    ## ==============================================

    #ptArray[np.isnan(ptArray)] = dst_nodata

    dst_band.WriteArray(ptArray)

    if dst_format == "AAIGrid":
        driver = gdal.GetDriverByName(dst_format)
        dst_ds_aai = driver.CreateCopy(dst_gdal, dst_ds)
    
    dst_ds = dst_ds_aai = dst_band = ptArray = None

### End
