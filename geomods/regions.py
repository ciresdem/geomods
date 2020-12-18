### regions.py
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

import os

import ogr
import osr
import numpy as np
from geomods import gdalfun
from geomods import utils

## ==============================================
## regions - regions are a bounding box list:
## [w, e, s, n]
## -- regions.py --
## ==============================================
def region_valid_p(region):
    '''return True if `region` [xmin, xmax, ymin, ymax] appears to be valid'''
    
    if region is not None:
        if region[0] <= region[1] and region[2] <= region[3]: return(True)
        else: return(False)
    else: return(False)

def region_center(region):
    '''find the center point [xc, yc] of the `region` [xmin, xmax, ymin, ymax]

    returns the center point [xc, yc]'''
    
    xc = region[0] + (region[1] - region[0] / 2)
    yc = region[2] + (region[3] - region[2] / 2)
    return([xc, yc])

def region_pct(region, pctv):
    '''calculate a percentage buffer for the `region` [xmin, xmax, ymin, ymax]

    returns the pctv buffer val of the region'''
    
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv = 0, pct = False):
    '''return the region buffered by buffer-value `bv`
    if `pct` is True, attain the buffer-value via: region_pct(region, bv)

    returns the buffered region [xmin, xmax, ymin, ymax]'''
    
    if pct: bv = region_pct(region, bv)
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def regions_reduce(region_a, region_b):
    '''combine two regions and find their minimum combined region.
    if the regions don't overlap, will return an invalid region.
    check the result with region_valid_p()
    
    return the minimum region [xmin, xmax, ymin, ymax] when combining `region_a` and `region_b`'''
    
    region_c = [0, 0, 0, 0]
    region_c[0] = region_a[0] if region_a[0] > region_b[0] else region_b[0]
    region_c[1] = region_a[1] if region_a[1] < region_b[1] else region_b[1]
    region_c[2] = region_a[2] if region_a[2] > region_b[2] else region_b[2]
    region_c[3] = region_a[3] if region_a[3] < region_b[3] else region_b[3]
    return(region_c)

def regions_merge(region_a, region_b):
    '''combine two regions and find their maximum combined region.

    returns maximum region [xmin, xmax, ymin, ymax] when combining `region_a` `and region_b`'''
    
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
    '''check if two regions intersect.
    region_valid_p(regions_reduce(region_a, region_b))

    return True if `region_a` and `region_b` intersect else False'''
    
    if region_a is not None and region_b is not None:
        return(region_valid_p(regions_reduce(region_a, region_b)))
    else: return(False)
    
def regions_intersect_ogr_p(region_a, region_b):
    '''check if two regions intersect.
    region_a_ogr_geom.Intersects(region_b_ogr_geom)
    
    return True if `region_a` and `region_b` intersect else False.'''
    
    if region_a is not None and region_b is not None:
        geom_a = gdalfun.gdal_region2geom(region_a)
        geom_b = gdalfun.gdal_region2geom(region_b)
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(True)

def region_format(region, t = 'gmt'):
    '''format region to string, defined by `t`
    t = 'str': xmin/xmax/ymin/ymax
    t = 'sstr': xmin xmax ymin ymax
    t = 'gmt': -Rxmin/xmax/ymin/ymax
    t = 'bbox': xmin,ymin,xmax,ymax
    t = 'te': xmin ymin xmax ymax
    t = 'ul_lr': xmin ymax xmax ymin
    t = 'fn': ymax_xmin

    returns the formatted region as str'''

    if t == 'str': return('/'.join([str(x) for x in region[:4]]))
    elif t == 'sstr': return(' '.join([str(x) for x in region[:4]]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region[:4]]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'te': return(' '.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'ul_lr': return(' '.join([str(region[0]), str(region[3]), str(region[1]), str(region[2])]))
    elif t == 'fn':
        ns = 's' if region[3] < 0 else 'n'
        ew = 'e' if region[0] > 0 else 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100)) % 100, 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100)) % 100))
    elif t == 'inf': return(' '.join([str(x) for x in region]))

def region_chunk(region, inc, n_chunk = 10):
    '''chunk the region [xmin, xmax, ymin, ymax] into 
    n_chunk by n_chunk cell regions, given inc.

    returns a list of chunked regions.'''
    
    i_chunk = 0
    x_i_chunk = 0
    x_chunk = n_chunk
    o_chunks = []
    xcount, ycount, dst_gt = gdalfun.gdal_region2gt(region, inc)
    
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
            o_chunks.append([geo_x_o, geo_x_t, geo_y_o, geo_y_t])
        
            if y_chunk < ycount:
                y_chunk += n_chunk
                i_chunk += 1
            else: break
        if x_chunk < xcount:
            x_chunk += n_chunk
            x_i_chunk += 1
        else: break
    return(o_chunks)

def regions_sort(trainers):
    '''sort regions by distance; regions is a list of regions [xmin, xmax, ymin, ymax].

    returns the sorted region-list'''
    
    train_sorted = []
    for z, train in enumerate(trainers):
        train_d = []
        np.random.shuffle(train)
        while True:
            if len(train) == 0: break
            this_center = region_center(train[0][0])
            train_d.append(train[0])
            train = train[1:]
            if len(train) == 0: break
            dsts = [utils.hav_dst(this_center, region_center(x[0])) for x in train]
            min_dst = np.percentile(dsts, 50)
            d_t = lambda t: utils.hav_dst(this_center, region_center(t[0])) > min_dst
            np.random.shuffle(train)
            train.sort(reverse=True, key=d_t)
        #echo_msg(' '.join([region_format(x[0], 'gmt') for x in train_d[:25]]))
        train_sorted.append(train_d)
    return(train_sorted)

def region_warp(region, s_warp = 4326, t_warp = 4326):
    '''warp region from source `s_warp` to target `t_warp`, using EPSG keys

    returns the warped region'''
    
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

def z_region_valid_p(z_region):
    '''return True if z_region appears to be valid'''
    
    if len(z_region) < 2: return(False)
    if z_region[0] > z_region[1]: return(False)
    return(True)

def z_region_pass(region, upper_limit = None, lower_limit = None):
    '''return True if extended region [xmin, xmax, ymin, ymax, zmin, zmax] is 
    within upper and lower z limits'''
    
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

def z_pass(z, upper_limit = None, lower_limit = None):
    '''return True if z-value is between upper and lower z limits'''

    if z is None: return(True)
    if upper_limit is not None:
        if z >= upper_limit:
            return(False)
    if lower_limit is not None:
        if z <= lower_limit:
            return(False)
    return(True)

### End
