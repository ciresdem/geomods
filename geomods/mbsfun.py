### mbsfun.py
##
## Copyright (c) 2010 - 2020 CIRES Coastal DEM Team
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
### Code:

## ==============================================
## import gdal/numpy
## ==============================================
import gdal
import ogr
import numpy as np

## ==============================================
## import geomods
## ==============================================
from geomods import utils
from geomods import regions
from geomods import gdalfun

## ==============================================
## MB-System Wrapper Functions - mbsfun.py
##
## MS-System must be installed on the system to run
## these functions and commands.
## ==============================================
def mb_inf(src_xyz, src_fmt = 168):
    '''generate an info (.inf) file from a src_xyz file using MBSystem.

    return mb_inf_parse(inf_file)'''

    utils.run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz.name), verbose = False)
    return(mb_inf_parse('{}.inf'.format(src_xyz.name)))

def mb_inf_data_format(src_inf):
    with open(src_inf) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                if til[0] == 'MBIO':
                    return(til[4])
                    

def mb_inf_parse(src_inf):
    xyzi = {'name': src_inf, 'numpts': 0, 'minmax': [0,0,0,0,0,0], 'wkt': regions.region2wkt([0,0,0,0,0,0])}
    dims = []
    this_row = 0
    
    with open(src_inf) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                if til[0] == 'Swath':
                    if til[2] == 'File:':
                        xyzi['name'] = til[3]
                if til[0] == 'Number':
                    if til[2] == 'Records:':
                        xyzi['numpts'] = int(til[3])
                if til[0] == 'Minimum':
                    if til[1] == 'Longitude:':
                        xyzi['minmax'][0] = float(til[2])
                        xyzi['minmax'][1] = float(til[5])
                    elif til[1] == 'Latitude:':
                        xyzi['minmax'][2] = float(til[2])
                        xyzi['minmax'][3] = float(til[5])
                    elif til[1] == 'Depth:':
                        xyzi['minmax'][4] = float(til[5]) * -1
                        xyzi['minmax'][5] = float(til[2]) * -1
                if til[0] == 'CM':
                    if til[1] == 'dimensions:':
                        dims = [int(til[2]), int(til[3])]
                        cm_array = np.zeros((dims[0], dims[1]))
                if til[0] == 'CM:':
                    for j in range(0, dims[0]):
                        cm_array[this_row][j] = int(til[j+1])
                    this_row += 1

    xinc = (xyzi['minmax'][1] - xyzi['minmax'][0]) / dims[0]
    yinc = (xyzi['minmax'][2] - xyzi['minmax'][3]) / dims[1]
    xcount, ycount, dst_gt = regions.region2gt(xyzi['minmax'], xinc, y_inc = yinc)
    ds_config = gdalfun.gdal_set_infos(dims[0], dims[1], dims[1] * dims[0], dst_gt, gdalfun.gdal_sr_wkt(4326), gdal.GDT_Float32, 0, 'GTiff')

    driver = gdal.GetDriverByName('MEM')
    ds = driver.Create('tmp', ds_config['nx'], ds_config['ny'], 1, ds_config['dt'])
    ds.SetGeoTransform(ds_config['geoT'])
    ds.SetProjection(ds_config['proj'])
    ds_band = ds.GetRasterBand(1)
    ds_band.SetNoDataValue(ds_config['ndv'])
    ds_band.WriteArray(cm_array)
    
    tmp_ds = ogr.GetDriverByName('Memory').CreateDataSource('tmp_poly')
    tmp_layer = tmp_ds.CreateLayer('tmp_poly', None, ogr.wkbMultiPolygon)
    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
    
    gdal.Polygonize(ds_band, ds_band, tmp_layer, 0)

    ## TODO: scan all features
    feat = tmp_layer.GetFeature(0)
    geom = feat.GetGeometryRef()
    wkt = geom.ExportToWkt()
    tmp_ds = ds = None
    
    xyzi['wkt'] = wkt
    return(xyzi)
### End
