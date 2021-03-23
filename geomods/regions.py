### regions.py
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
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
## regions - regions are a bounding box list: [w, e, s, n]
##
### Code:

## imports
import os

## ==============================================
## import gdal/numpy
## ==============================================
import ogr
import osr
import numpy as np

## ==============================================
## import geomods
## ==============================================
from geomods import utils

def region_valid_p(region):
    """return True if `region` [xmin, xmax, ymin, ymax] appears to be valid

    Arguments:
      region(list): a region list
    """
    
    if region is not None:
        if region_wkt_p(region):
            return(False)
        
        try:
            if region[0] <= region[1] and region[2] <= region[3]:
                return(True)
            else:
                return(False)
        except:
            return(False)
        
    else:
        return(False)

def region_is_zeros(region):
    if region[0] == 0 and region[1] == 0 and region[2] == 0 and region[3] == 0:
        return(True)
    else: return(False)
                    
def region_wkt_p(region):
    
    try:
        if region.split()[0] == "POLYGON" or region.split()[0] == "MULTIPOLYGON":
            wkt_p = True
        else:
            wkt_p = False
    except:
        wkt_p = False
        
    return(wkt_p)
    
def str2region(region_str):
    """attempt to convert a string into a region

    Arguments:
      region(str): a region string

    Returns:
      region(list/wkt): a region list or a polygon wkt
    """
    
    try:
        this_region = [float(x) for x in region_str.split('/')]
    except ValueError:
        this_region = gdal_ogr_multi_poly(region_str)
    except Exception as e:
        utils.echo_error_msg('failed to parse region(s), {}'.format(e))
        this_region = None

    if this_region is None:
        if region_wkt_p(regions_str):
            this_region = region_str
            
    return(this_region)

def region_region(region):
    if region_valid_p(region):
        return(region)
    elif region_wkt_p(region):
        return(wkt2region(region))
    else:
        return(None)

def region_center(region):
    """find the center point [xc, yc] of the `region` [xmin, xmax, ymin, ymax]

    returns the center point [xc, yc]
    """
    
    xc = region[0] + (region[1] - region[0] / 2)
    yc = region[2] + (region[3] - region[2] / 2)
    return([xc, yc])

def region_pct(region, pctv):
    """calculate a percentage buffer for the `region` [xmin, xmax, ymin, ymax]

    returns the pctv buffer val of the region
    """
    
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv=0, pct=False):
    """return the region buffered by buffer-value `bv`
    if `pct` is True, attain the buffer-value via: region_pct(region, bv)

    returns the buffered region [xmin, xmax, ymin, ymax]
    """
    
    if pct:
        bv = region_pct(region, bv)
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def regions_reduce(region_a, region_b):
    """combine two regions and find their minimum combined region.
    if the regions don't overlap, will return an invalid region.
    check the result with region_valid_p()
    
    return the minimum region [xmin, xmax, ymin, ymax] when combining `region_a` and `region_b`
    """

    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] > region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] < region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] > region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] < region_b[3] else region_b[3]
    return(region_c)

def regions_merge(region_a, region_b):
    """combine two regions and find their maximum combined region.

    returns maximum region [xmin, xmax, ymin, ymax] when combining `region_a` `and region_b`
    """
    
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] < region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] > region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] < region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] > region_b[3] else region_b[3]
    if len(region_a) > 4 and len(region_b) > 4:
        region_c.append(region_a[4] if region_a[4] < region_b[4] else region_b[4])
        region_c.append(region_a[5] if region_a[5] > region_b[5] else region_b[5])
    return(region_c)

def regions_intersect_p(region_a, region_b):
    """check if two regions intersect.
    region_valid_p(regions_reduce(region_a, region_b))

    return True if `region_a` and `region_b` intersect else False
    """
    
    if region_a is not None and region_b is not None:
        return(region_valid_p(regions_reduce(region_a, region_b)))
    else: return(False)

def geoms_intersect_p(geom_a, geom_b):
    """check if OGR geometries `geom_a` and `geom_b` intersect

    returns True for intersection
    """
    
    if geom_a is not None and geom_b is not None:
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(True)
    
def regions_intersect_ogr_p(region_a, region_b):
    """check if two regions intersect.
    region_a_ogr_geom.Intersects(region_b_ogr_geom)
    
    return True if `region_a` and `region_b` intersect else False.
    """
    
    if region_a is not None and region_b is not None:
        geom_a = region2geom(region_a)
        geom_b = region2geom(region_b)
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(True)

def region_format(region, t='gmt'):
    """format region to string, defined by `t`

    t = 'str': xmin/xmax/ymin/ymax
    t = 'sstr': xmin xmax ymin ymax
    t = 'gmt': -Rxmin/xmax/ymin/ymax
    t = 'bbox': xmin,ymin,xmax,ymax
    t = 'osm_bbox': ymin,xmin,ymax,xmax
    t = 'te': xmin ymin xmax ymax
    t = 'ul_lr': xmin ymax xmax ymin
    t = 'fn': ymax_xmin

    returns the formatted region as str
    """
    
    if t == 'str': return('/'.join([str(x) for x in region[:4]]))
    elif t == 'sstr': return(' '.join([str(x) for x in region[:4]]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region[:4]]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'osm_bbox': return(','.join([str(region[2]), str(region[0]), str(region[3]), str(region[1])]))
    elif t == 'te': return(' '.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'ul_lr': return(' '.join([str(region[0]), str(region[3]), str(region[1]), str(region[2])]))
    elif t == 'fn':
        ns = 's' if region[3] < 0 else 'n'
        ew = 'e' if region[0] > 0 else 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(
            ns, abs(int(region[3])), abs(int(region[3] * 100)) % 100, 
            ew, abs(int(region[0])), abs(int(region[0] * 100)) % 100))
    elif t == 'inf': return(' '.join([str(x) for x in region]))

def region_warp(region, src_epsg=4326, dst_epsg=4326):
    
    if dst_epsg is None or src_epsg is None: return(region)
        
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(int(src_epsg))

    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(int(dst_epsg))
    try:
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except: pass
    dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)

    pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[0], region[2]))
    pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[1], region[3]))
    pointA.Transform(dst_trans)
    pointB.Transform(dst_trans)
    ds_region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(ds_region)
    
def region_chunk(region, inc, n_chunk=10, buff=0):
    """chunk the region [xmin, xmax, ymin, ymax] into 
    n_chunk by n_chunk cell regions, given inc.

    returns a list of chunked regions.
    """
    
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk
    o_chunks = []
    xcount, ycount, dst_gt = region2gt(region, inc)
    
    while True:
        y_chunk = n_chunk
        while True:
            this_x_origin = x_chunk - n_chunk
            this_y_origin = y_chunk - n_chunk
            this_x_size = x_chunk - this_x_origin
            this_y_size = y_chunk - this_y_origin
            
            geo_x_o = region[0] + this_x_origin * inc
            geo_x_t = geo_x_o + this_x_size * inc
            geo_y_o = region[2] + this_y_origin * inc
            geo_y_t = geo_y_o + this_y_size * inc

            if geo_y_t > region[3]: geo_y_t = region[3]
            if geo_y_o < region[2]: geo_y_o = region[2]
            if geo_x_t > region[1]: geo_x_t = region[1]
            if geo_x_o < region[0]: geo_x_o = region[0]
            chunkd_region = [geo_x_o, geo_x_t, geo_y_o, geo_y_t]
            region_buffer(chunkd_region, buff)
            o_chunks.append(chunkd_region)
        
            if y_chunk < ycount:
                y_chunk += n_chunk
                i_chunk += 1
            else: break
        if x_chunk < xcount:
            x_chunk += n_chunk
            x_i_chunk += 1
        else: break
    return(o_chunks)

def regions_sort(trainers, t_num=25, verbose=False):
    """sort regions by distance; regions is a list of regions [xmin, xmax, ymin, ymax].

    returns the sorted region-list
    """
    
    train_sorted = []
    for z, train in enumerate(trainers):
        train_d = []
        np.random.shuffle(train)
        train_total = len(train)
        while True:
            if verbose: utils.echo_msg_inline('sorting training tiles [{}]'.format(len(train)))
            if len(train) == 0: break
            this_center = region_center(train[0][0])
            train_d.append(train[0])
            train = train[1:]
            if len(train_d) > t_num or len(train) == 0: break
            dsts = [utils.euc_dst(this_center, region_center(x[0])) for x in train]
            min_dst = np.percentile(dsts, 50)
            d_t = lambda t: utils.euc_dst(this_center, region_center(t[0])) > min_dst
            np.random.shuffle(train)
            train.sort(reverse=True, key=d_t)
        if verbose: utils.echo_msg(' '.join([region_format(x[0], 'gmt') for x in train_d[:t_num]]))
        train_sorted.append(train_d)
    if verbose: utils.echo_msg_inline('sorting training tiles [OK]\n')
    return(train_sorted)

def z_region_valid_p(z_region):
    """return True if z_region appears to be valid"""

    if len(z_region) < 2: return(False)
    if z_region[0] > z_region[1]: return(False)
    return(True)

def z_region_pass(region, upper_limit=None, lower_limit=None):
    """return True if extended region [xmin, xmax, ymin, ymax, zmin, zmax] is 
    within upper and lower z limits
    """
    
    if region is not None:
        z_region = region[4:]
        if z_region is not None and len(z_region) >= 2:
            if upper_limit is not None:
                if z_region[0] >= upper_limit:
                    return(False)
            if lower_limit is not None:
                if z_region[1] <= lower_limit:
                    return(False)
    return(True)

def z_pass(z, upper_limit=None, lower_limit=None):
    """return True if z-value is between upper and lower z limits"""
    
    if z is None: return(True)
    if upper_limit is not None:
        if z >= upper_limit:
            return(False)
    if lower_limit is not None:
        if z <= lower_limit:
            return(False)
    return(True)


def create_wkt_polygon(coords, xpos=1, ypos=0):
    """convert coords to Wkt

    returns polygon as wkt
    """
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[xpos], coord[ypos])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def gt2region(ds_config):
    """convert a gdal geo-tranform to a region [xmin, xmax, ymin, ymax] via a data-source config dict.

    returns region of gdal data-source
    """
    
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + (geoT[1]*ds_config['nx']), geoT[3] + (geoT[5]*ds_config['ny']), geoT[3]])

def region2gt(region, inc, y_inc=None):
    """return a count info and a gdal geotransform based on extent and cellsize

    returns a list [xcount, ycount, geot]
    """
    
    if y_inc is None:
        y_inc = inc * -1.
        
    dst_gt = (region[0], inc, 0, region[3], 0, y_inc)    
    this_origin = utils._geo2pixel(region[0], region[3], dst_gt)
    this_end = utils._geo2pixel(region[1], region[2], dst_gt)
    this_size = (this_end[0] - this_origin[0], this_end[1] - this_origin[1])
    return(this_size[0], this_size[1], dst_gt)

def wkt2region(wkt):
    
    return(ogr.CreateGeometryFromWkt(wkt).GetEnvelope())
            
def wkt2geom(wkt):
    
    return(ogr.CreateGeometryFromWkt(wkt))

def region2wkt(region):
    
    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]    
    return(create_wkt_polygon(eg))
    
def region2geom(region):
    """convert an extent [west, east, south, north] to an OGR geometry

    Return: 
      ogr geometry
    """

    if region_wkt_p(region):
        geom = ogr.CreateGeometryFromWkt(region)

    else:
        eg = [[region[2], region[0]], [region[2], region[1]],
              [region[3], region[1]], [region[3], region[0]],
              [region[2], region[0]]]
        geom = ogr.CreateGeometryFromWkt(create_wkt_polygon(eg))
    return(geom)

def region2ogr(region, dst_ogr, append=False):
    """convert a region string to an OGR vector"""
    
    dst_wkt = region2wkt(region)
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

def gdal_ogr_regions(src_ds):
    """return the region(s) of the ogr dataset"""
    
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

def gdal_ogr_polys(src_ds):
    """return the region(s) of the ogr dataset"""
    
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                these_regions.append(pwkt)
                
        poly = None
    return(these_regions)

def gdal_ogr_multi_poly(src_ds):
    """return the region(s) of the ogr dataset"""
    
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for feat in p_layer:
                feat.geometry().CloseRings()
                wkt = feat.geometry().ExportToWkt()
                multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
            wkt = multi.ExportToWkt()
        poly = None
    return(wkt)


def gdal_region_warp(region, s_warp=4326, t_warp=4326):
    """warp region from source `s_warp` to target `t_warp`, using EPSG keys

    returns the warped region
    """
    
    if s_warp is None or s_warp == t_warp: return(region)
    
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(int(s_warp))

    try:
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except: pass
    
    if t_warp is not None:
        dst_srs = osr.SpatialReference()
        try:
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        dst_srs.ImportFromEPSG(int(t_warp))
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)        
        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[0], region[2]))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(region[1], region[3]))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        region = [pointA.GetX(), pointB.GetX(), pointA.GetY(), pointB.GetY()]
    return(region)

### End
