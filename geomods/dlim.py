### dlim.py - DataLists IMproved
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
### Code:

import os
import sys
import json
import glob
import hashlib

## ==============================================
## import gdal/numpy
## ==============================================
import gdal
import ogr
import osr
import numpy as np
from scipy import spatial

def dl_hash(fn, sha1 = False):

    BUF_SIZE = 65536  # lets read stuff in 64kbchunks!
    if sha1: this_hash = hashlib.sha1()
    else: this_hash = hashlib.md5()

    with open(fn, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            this_hash.update(data)
            
    return(this_hash.hexdigest())

def int_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(int(val))
    except: return(or_val)

def float_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(float(val))
    except: return(or_val)

def _geo2pixel(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = ((geo_x - geoTransform[0]) / geoTransform[1]) + .5
        pixel_y = ((geo_y - geoTransform[3]) / geoTransform[5]) + .5
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
    return(int(pixel_x), int(pixel_y))

def _geo2pixel_affine(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    note: use _geo2pixel instead

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    import affine
    forward_transform = affine.Affine.from_gdal(*geoTransform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x + 0.5), int(pixel_y + 0.5)
    return(pixel_x, pixel_y)

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    """convert a pixel location to geographic coordinates given geoTransform

    Args:
      pixel_x (int): the x pixel value
      pixel_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster
    
    Returns:
      list: [geographic-x, geographic-y]
    """
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    """apply geotransform to in_x,in_y
    
    Args:
      in_x (int): the x pixel value
      in_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: [geographic-x, geographic-y]
    """
    
    out_x = geoTransform[0] + int(in_x + 0.5) * geoTransform[1] + int(in_y + 0.5) * geoTransform[2]
    out_y = geoTransform[3] + int(in_x + 0.5) * geoTransform[4] + int(in_y + 0.5) * geoTransform[5]

    return(out_x, out_y)

def _invert_gt(geoTransform):
    """invert the geotransform
    
    Args:
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a geo-transform list describing a raster
    """
    
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

def create_wkt_polygon(coords, xpos = 1, ypos = 0):
    """convert coords to Wkt

    Args:
      coords (list): x/y geographic coords
      xpos (int): the position of the x value in coords
      ypos (int): the position of the y value in corrds

    Returns:
      wkt: polygon as wkt
    """

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[xpos], coord[ypos])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def remove_glob(*args):
    """glob `glob_str` and os.remove results, pass if error
    
    Args:
      *args (str): any number of pathname or dirname strings

    Returns:
      int: 0
    """
    
    for glob_str in args:
        try:
            globs = glob.glob(glob_str)
            for g in globs:
                if os.path.isdir(g):
                    remove_globs('{}/*'.format(g))
                    os.removedirs(g)
                else: os.remove(g)
        except: pass
    return(0)

def sr_wkt(epsg, esri = False):
    """convert an epsg code to wkt

    Returns:
      (str): wkt or None
    """
    
    try:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(epsg))
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)
    
## ==============================================
## verbosity functions
## ==============================================
def echo_warning_msg2(msg, prefix = 'waffles'):
    """echo warning msg to stderr using `prefix`

    >> echo_warning_msg2('message', 'test')
    test: warning, message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1mwarining\033[m, {}\n'.format(prefix, msg))

def echo_error_msg2(msg, prefix = 'waffles'):
    """echo error msg to stderr using `prefix`

    >> echo_error_msg2('message', 'test')
    test: error, message

    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles', nl = True):
    """echo `msg` to stderr using `prefix`

    >> echo_msg2('message', 'test')
    test: message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
      nl (bool): append a newline to the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
    sys.stderr.flush()
    
## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_msg_inline = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]), nl = False)

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix = os.path.basename(sys.argv[0]))

class _progress:
    """geomods minimal progress indicator"""

    def __init__(self, message = None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.opm = message
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
        
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def update_perc(self, p, msg = None):
        if len(p) == 2:
            self._clear_stderr()
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {:40}\r'.format(self.perc(p), msg if msg is not None else self.opm))
        else: self.update()
        
    def update(self, msg = None):

        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], msg if msg is not None else self.opm))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)

    def end(self, status, end_msg = None):
        self._clear_stderr()
        if end_msg is None: end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))

## ==============================================
## regions
## ==============================================
class region:
    """Representing a geographic region.
    
    Attributes:
      xmin (float): the minimum x(long) value
      xmax (float): the maximum x(long) value
      ymin (float): the minimum x(lat) value
      ymax (float): the maximum x(lat) value
      zmin (float): the minimum z(elev) value
      zmax (float): the maximum z(elev) value
      wmin (float): the minimum w(weight) value
      wmax (float): the maximum w(weight) value
      wkt (str): the wkt representation of the region
      epsg (int): the EPSG code representing the regions projection

    Methods:
      _valid_p(): True if the region is valid else False
      buffer(bv = 0): buffer the region by buffer value
      center(): find the center point of the region
    """

    full_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax, x.zmin, x.zmax, x.wmin, x.wmax]
    xy_region = lambda x: [x.xmin, x.xmax, x.ymin, x.ymax]
    z_region = lambda x: [x.zmin, x.zmax]
    w_region = lambda x: [x.wmin, x.wmax]
    
    def __init__(self, xmin = None, xmax = None, ymin = None, ymax = None,
                 zmin = None, zmax = None, wmin = None, wmax = None,
                 epsg = 4326, wkt = None):
        """
        Args:
          xmin (float): the minimum x(long) value
          xmax (float): the maximum x(long) value
          ymin (float): the minimum x(lat) value
          ymax (float): the maximum x(lat) value
          zmin (float): the minimum z(elev) value
          zmax (float): the maximum z(elev) value
          wmin (float): the minimum w(weight) value
          wmax (float): the maximum w(weight) value
          epsg (int): the EPSG code representing the regions projection
        """
        
        self.xmin = float_or(xmin)
        self.xmax = float_or(xmax)
        self.ymin = float_or(ymin)
        self.ymax = float_or(ymax)
        self.zmin = float_or(zmin)
        self.zmax = float_or(zmax)
        self.wmin = float_or(wmin)
        self.wmax = float_or(wmax)

        self.epsg = int_or(epsg)
        self.wkt = wkt
        
    def _valid_p(self, check_xy = False):
        """check the validity of a region
    
        Returns:
          bool: True if region  appears to be valid else False
        """
        
        if check_xy:
            if self.xmin is None: return(False)
            if self.xmax is None: return(False)
            if self.ymin is None: return(False)
            if self.ymax is None: return(False)

        if self.xmin is not None and self.xmax is not None:
            if self.xmin > self.xmax: return(False)
        if self.ymin is not None and self.ymax is not None:
            if self.ymin > self.ymax: return(False)
        if self.zmin is not None and self.zmax is not None:
            if self.zmin > self.zmax: return(False)
        if self.wmin is not None and self.wmax is not None:
            if self.wmin > self.wmax: return(False)

        if not any(self.full_region()): return(False)
        return(True)

    def from_list(self, region_list):
        """import a region from a region list 

        Args:
          region_list (list): [xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax]

        Returns:
          region-object: self
        """
        
        if len(region_list) >= 4:
            self.xmin = float_or(region_list[0])
            self.xmax = float_or(region_list[1])
            self.ymin = float_or(region_list[2])
            self.ymax = float_or(region_list[3])
            if len(region_list) >= 6:
                self.zmin = float_or(region_list[4])
                self.zmax = float_or(region_list[5])
                if len(region_list) >= 8:
                    self.wmin = float_or(region_list[6])
                    self.wmax = float_or(region_list[7])
            if self.wkt is None:
                if self._valid_p(check_xy = True):
                    self.wkt = self.export_as_wkt()
                else: self.wkt = None
        return(self)
        
    def from_string(self, region_str):
        """import a region from a region string 

        Args:
          region_str (str): <-R>xmin/xmax/ymin/ymax/zmin/zmax

        Returns:
          region-object: self
        """
        
        if region_str[:2] == '-R':
            region_str = region_str[2:]
            
        str_list = region_str.strip().split('/')
        if len(str_list) >= 4:
            r_list = [float_or(x) for x in str_list]
            self.from_list(r_list)
        elif region_str.split()[0] == "POLYGON":
            self.wkt = region_str
            self.from_list(ogr.CreateGeometryFromWkt(region_str).GetEnvelope())
        return(self)
    
    def from_geo_transform(self, geoT = None, x_count = None, y_count = None):
        """import a region from a region string 

        Args:
          geoT (list): a geo transform list
          x_count (int): length of x axis
          y_count (int): length of y axis

        Returns:
          region-object: self
        """
        
        if geoT is not None and x_count is not None and y_count is not None:
            self.xmin = geoT[0]
            self.xmax = geoT[0] + geoT[1] * x_count
            self.ymin = geoT[3] + geoT[5] * y_count
            self.ymax = geoT[3]
        if self.wkt is None:
            self.wkt = self.export_as_wkt()
        return(self)
    
    def format(self, t = 'gmt'):
        """format region to string, defined by `t`

        Args:
          t (str): output mode
            t = 'str': xmin/xmax/ymin/ymax
            t = 'sstr': xmin xmax ymin ymax
            t = 'gmt': -Rxmin/xmax/ymin/ymax
            t = 'bbox': xmin,ymin,xmax,ymax
            t = 'osm_bbox': ymin,xmin,ymax,xmax
            t = 'te': xmin ymin xmax ymax
            t = 'ul_lr': xmin ymax xmax ymin
            t = 'fn': ymax_xmin

        Returns:
          str: the formatted region as str or None if region is invalid
        """

        if self._valid_p():
            if t == 'str': return('/'.join([str(x) for x in self.region[:4]]))
            elif t == 'sstr': return(' '.join([str(x) for x in self.region[:4]]))
            elif t == 'gmt': return('-R' + '/'.join([str(self.xmin), str(self.xmax), str(self.ymin), str(self.ymax)]))
            elif t == 'bbox': return(','.join([str(self.xmin), str(self.ymin), str(self.xmax), str(self.ymax)]))
            elif t == 'osm_bbox': return(','.join([str(self.ymin), str(self.xmin), str(self.ymax), str(self.xmax)]))
            elif t == 'te': return(' '.join([str(self.xmin), str(self.ymin), str(self.xmax), str(self.ymax)]))
            elif t == 'ul_lr': return(' '.join([str(self.xmin), str(self.ymax), str(self.xmax), str(self.ymin)]))
            elif t == 'fn':
                ns = 's' if self.ymax < 0 else 'n'
                ew = 'e' if self.xmin > 0 else 'w'
                return('{}{:02d}x{:02d}_{}{:03d}x{:02d}\
                '.format(ns, abs(int(self.ymax)), abs(int(self.ymax * 100)) % 100, 
                         ew, abs(int(self.xmin)), abs(int(self.xmin * 100)) % 100))
            elif t == 'inf': return(' '.join([str(x) for x in self.region]))
            else: return('/'.join([str(x) for x in self.region[:4]]))
        else: return(None)
        
    def geo_transform(self, x_inc = 0, y_inc = None):
        """return a count info and a geotransform based on the region and a cellsize

        Args:
          x_inc (float): a x-axis gridding increment
          y_inc (float): a y-axis gridding increment

        Returns:
          list: [xcount, ycount, geot]
        """

        if y_inc is None: y_inc = x_inc * -1.
        dst_gt = (self.xmin, x_inc, 0, self.ymax, 0, y_inc)    
        this_origin = _geo2pixel(self.xmin, self.ymax, dst_gt)
        this_end = _geo2pixel(self.xmax, self.ymin, dst_gt)
        this_size = (this_end[0] - this_origin[0], this_end[1] - this_origin[1])
        return(this_end[0] - this_origin[0], this_end[1] - this_origin[1], dst_gt)

    def export_as_list(self, include_z = False, include_w = False):
        """export region as a list

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output

        Returns:
          list: the region values in a list
        """
        
        region_list = [self.xmin, self.xmax, self.ymin, self.ymax]
        if include_z:
            region_list.append(self.zmin)
            region_list.append(self.zmax)
        if include_w:
            region_list.append(self.wmin)
            region_list.append(self.wmax)
        return(region_list)
    
    def export_as_wkt(self):
        """transform a region to wkt

        Returns:
          wkt: a wkt polygon geometry
        """

        eg = [[self.ymin, self.xmin], [self.ymin, self.xmax],
              [self.ymax, self.xmax], [self.ymax, self.xmin],
              [self.ymin, self.xmin]]
        return(create_wkt_polygon(eg))

    def export_as_geom(self):
        """convert a region to an OGR geometry

        Returns:
          ogr-geom: an ogr polygon geometry
        """

        if self.wkt is not None:
            return(ogr.CreateGeometryFromWkt(self.wkt))
        else: return(None)

    def export_as_ogr(self, dst_ogr, dst_fmt = 'ESRI Shapefile', append = False):
        """convert a region string to an OGR vector

        Args:
          dst_ogr (str): destination ogr dataset pathname
          append (bool): Append to existing dataset
        """

        
        wkt = self.export_as_wkt()
        driver = ogr.GetDriverByName(dst_fmt)
        if os.path.exists(dst_ogr):
            driver.DeleteDataSource(dst_ogr)

        dst_ds = driver.CreateDataSource(dst_ogr)
        dst_lyr = dst_ds.CreateLayer(dst_ogr, geom_type = ogr.wkbPolygon)
        dst_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        dst_feat = ogr.Feature(dst_lyr.GetLayerDefn())
        dst_feat.SetGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
        dst_feat.SetField('id', 1)
        dst_lyr.CreateFeature(dst_feat)
        dst_feat = None
        dst_ds = None

    def srcwin(self, gt, x_count, y_count):
        """output the appropriate gdal srcwin for the region.

        Args:
          region (region): an input region

        returns the gdal srcwin
        """

        this_origin = [0 if x < 0 else x for x in _geo2pixel(self.xmin, self.ymax, gt)]
        this_end = [0 if x < 0 else x for x in _geo2pixel(self.xmax, self.ymin, gt)]
        this_size = [0 if x < 0 else x for x in ((this_end[0] - this_origin[0]), (this_end[1] - this_origin[1]))]
        if this_size[0] > x_count - this_origin[0]: this_size[0] = x_count - this_origin[0]
        if this_size[1] > y_count - this_origin[1]: this_size[1] = y_count - this_origin[1]
        return(this_origin[0], this_origin[1], this_size[0], this_size[1])
        
    def buffer(self, bv = 0, pct = None):
        """return the region buffered by buffer-value `bv`

        Args:
          bv (float): the buffer value
          pct (float): attain the buffer-value via percentage

        Returns:
          region-object: self
        """
        
        if self._valid:
            if pct is not None:
                ewp = (self.xmax - self.xmin) * (pctv * .01)
                nsp = (self.ymax - self.ymin) * (pctv * .01)
                bv = (ewp + nsp) / 2

            self.xmin -= bv
            self.xmax += bv
            self.ymin -= bv
            self.ymax += bv
        return(self)

    def center(self):
        """find the center point of the xy region

        Returns:
          list: the center point [xc, yc]
        """

        if self._valid_p():
            return([self.xmin + (self.xmax - self.xmin / 2),
                    self.ymax + (self.ymax - self.ymin / 2)])
        else: return(None)
        
    def chunk(self, inc, n_chunk = 10):
        """chunk the xy region [xmin, xmax, ymin, ymax] into 
        n_chunk by n_chunk cell regions, given inc.

        Args:
          inc (float): the chunking increment
          n_chunk (int): number of cells

        Returns:
          list: a list of chunked regions [<regions.region>, ...]
        """

        i_chunk = 0
        x_i_chunk = 0
        x_chunk = n_chunk
        o_chunks = []
        xcount, ycount, dst_gt = self.geo_transform(x_inc = inc)

        while True:
            y_chunk = n_chunk
            while True:
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = x_chunk - this_x_origin
                this_y_size = y_chunk - this_y_origin

                c_region = region()
                
                c_region.xmin = self.xmin + this_x_origin * inc
                c_region.xmax = c_region.xmin + this_x_size * inc
                c_region.ymin = self.ymin + this_y_origin * inc
                c_region.ymax = c_region.ymin + this_y_size * inc

                if c_region.ymax > self.ymax: c_region.ymax = self.ymax
                if c_region.ymin < self.ymin: c_region.ymin = self.ymin
                if c_region.xmax > self.xmax: c_region.xmax = self.xmax
                if c_region.xmin < self.xmin: c_region.xmin = self.xmin
                o_chunks.append(c_region)

                if y_chunk < ycount:
                    y_chunk += n_chunk
                    i_chunk += 1
                else: break
            if x_chunk < xcount:
                x_chunk += n_chunk
                x_i_chunk += 1
            else: break
        return(o_chunks)

    def warp(self, warp_epsg = 4326):
        """warp the region from to dst_epsg

        Args:
          warp_epsg (ing): the destination EPSG code

        Returns:
          region-object: warped self
        """

        warp_epsg = int_or(warp_epsg)
        if warp_epsg is None or self.epsg is None: return(self)

        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(self.epsg)

        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(warp_epsg)
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)

        pointA = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(self.xmin, self.ymin))
        pointB = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(self.xmax, self.ymax))
        pointA.Transform(dst_trans)
        pointB.Transform(dst_trans)
        self.xmin = pointA.GetX()
        self.xmax = pointB.GetX()
        self.ymin = pointA.GetY()
        self.ymax = pointB.GetY()
        self.epsg = warp_epsg
        return(self)

def ogr_wkts(src_ds):
    """return the wkt(s) of the ogr dataset
    """
    
    these_regions = []
    src_s = src_ds.split(':')
    if os.path.exists(src_s[0]):
        poly = ogr.Open(src_s[0])
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                r = region().from_string(pwkt)
                if len(src_s) > 1:
                    src_r = src_s[1].split('/')
                    if len(src_r) > 0: r.zmin = float_or(src_r[0])
                    if len(src_r) > 1: r.zmax = float_or(src_r[1])
                    if len(src_r) > 2: r.wmin = float_or(src_r[2])
                    if len(src_r) > 3:  r.wmax = float_or(src_r[3])
                these_regions.append(r)

        poly = None
    return(these_regions)
    
## ==============================================
## do things to and with regions 
## ==============================================
def regions_reduce(region_a, region_b):
    """combine two regions and find their minimum combined region.

    if the regions don't overlap, will return an invalid region.
    check the result with _valid_p()
    
    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      region: the minimum region when combining `region_a` and `region_b`
    """
    region_c = region()
    if region_a._valid_p() and region_b._valid_p():
        if region_a.xmin is not None and region_b.xmin is not None:
            region_c.xmin = region_a.xmin if region_a.xmin > region_b.xmin else region_b.xmin
        if region_a.xmax is not None and region_b.xmax is not None:
            region_c.xmax = region_a.xmax if region_a.xmax < region_b.xmax else region_b.xmax
        if region_a.ymin is not None and region_b.ymin is not None:
            region_c.ymin = region_a.ymin if region_a.ymin > region_b.ymin else region_b.ymin
        if region_a.ymax is not None and region_b.ymax is not None:
            region_c.ymax = region_a.ymax if region_a.ymax < region_b.ymax else region_b.ymax
        if region_a.zmin is not None and region_b.zmin is not None:
            region_c.zmin = region_a.zmin if region_a.zmin > region_b.zmin else region_b.zmin
        if region_a.zmax is not None and region_b.zmax is not None:
            region_c.zmax = region_a.zmax if region_a.zmax < region_b.zmax else region_b.zmax
        if region_a.wmin is not None and region_w.zmin is not None:
            region_c.wmin = region_a.wmin if region_a.wmin > region_b.wmin else region_b.wmin
        if region_a.wmax is not None and region_b.wmax is not None:
            region_c.wmax = region_a.wmax if region_a.wmax < region_b.wmax else region_b.wmax
    return(region_c)

def regions_merge(region_a, region_b):
    """combine two regions and find their maximum combined region.
    
    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      region: the maximum region [xmin, xmax, ymin, ymax] when combining `region_a` `and region_b`
    """
    
    region_c = region()
    if region_a._valid_p() and region_b._valid_p():
        region_c.xmin = region_a.xmin if region_a.xmin < region_b.xmin else region_b.xmin
        region_c.xmax = region_a.xmax if region_a.xmax > region_b.xmax else region_b.xmax
        region_c.ymin = region_a.ymin if region_a.ymin < region_b.ymin else region_b.ymin
        region_c.ymax = region_a.ymax if region_a.ymax > region_b.ymax else region_b.ymax
        if region_a.zmin is not None and region_b.zmin is not None:
            region_c.zmin = region_a.zmin if region_a.zmin < region_b.zmin else region_b.zmin
        if region_a.zmax is not None and region_b.zmax is not None:
            region_c.zmax = region_a.zmax if region_a.zmax < region_b.zmax else region_b.zmax
    return(region_c)

def regions_intersect_p(region_a, region_b):
    """check if two regions intersect.

    region_valid_p(regions_reduce(region_a, region_b))

    Args:
      region_a (region): a region object
      region_b (region): a region object

    Returns:
      bool: True if `region_a` and `region_b` intersect else False
    """
    
    if region_a._valid_p() and region_b._valid_p():
        return(regions_reduce(region_a, region_b)._valid_p())
    return(False)
    
def regions_intersect_ogr_p(region_a, region_b):
    """check if two regions intersect using ogr

    region_a_ogr_geom.Intersects(region_b_ogr_geom)
    
    Args:
      region_a (region): a region 
      region_b (region): a region 

    Returns:
      bool: True if `region_a` and `region_b` intersect else False.
    """

    if region_a._valid_p() and region_b._valid_p():
        if region_a.export_as_geom().Intersects(region_b.export_as_geom()):
            return(True)
    return(False)

def z_region_pass(region, upper_limit = None, lower_limit = None):
    """return True if extended region [xmin, xmax, ymin, ymax, zmin, zmax] is 
    within upper and lower z limits
    
    Args:
      region (list): a long-region [xmin, xmax, ymin, ymax, zmin, zmax]
      upper_limit (float): the z-max
      lower_limit (float): the z-min
    
    Returns:
      bool: True if region z values are within upper and lower limit
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

def xyz_in_region_p(xyz, this_region):
    """check if xyz point in inside the given region

    Args:
      xyz (xyz): an xyz data object
      this_region (region): a region object

    Returns:
      bool: True if xyz point inside region else False
    """

    pass_d = True
    xyz_wkt = xyz.export_as_wkt()
    p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
    r_geom = this_region.export_as_geom()
    if r_geom is not None:
        pass_d = p_geom.Within(r_geom)

    if pass_d:
        if this_region.zmin is not None:
            if xyz.z < this_region.zmin: pass_d = False
        if this_region.zmax is not None:
            if xyz.z > this_region.zmax: pass_d = False

    return(pass_d)

def gdal_ogr_regions(src_ds):
    """return the region(s) of the ogr dataset
    
    Args:
      src_ds (str): source ogr dataset pathname

    Returns:
      list: list of regions
    """
    
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                this_region = region()
                pgeom = pf.GetGeometryRef()
                pwkt = pgeom.ExportToWkt()
                this_region.from_list(ogr.CreateGeometryFromWkt(pwkt).GetEnvelope())
                these_regions.append(this_region)
        poly = None
    return(these_regions)

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

## ==============================================
## mbsystem parser
## ==============================================
class mbs_parser:
    """providing an mbsystem parser
    """

    def __init__(self, fn = None, epsg = None):
        self.fn = fn
        self.epsg = int_or(epsg)
        self.infos = {}
        
    def inf(self):
        pass

    def parse(self):
        pass

    def yield_xyz(self):
        pass
    
    def inf_parse(self):

        self.infos['name'] = self.fn
        self.infos['minmax'] = [0,0,0,0,0,0]
        self.infos['hash'] = None
        dims = []
        this_row = 0

        with open(self.fn) as iob:
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    if til[0] == 'Swath':
                        if til[2] == 'File:':
                            self.infos['name'] = til[3]
                    if til[0] == 'Number':
                        if til[2] == 'Records:':
                            self.infos['numpts'] = int_or(til[3])
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            self.infos['minmax'][0] = float_or(til[2])
                            self.infos['minmax'][1] = float_or(til[5])
                        elif til[1] == 'Latitude:':
                            self.infos['minmax'][2] = float_or(til[2])
                            self.infos['minmax'][3] = float_or(til[5])
                        elif til[1] == 'Depth:':
                            self.infos['minmax'][4] = float_or(til[5]) * -1
                            self.infos['minmax'][5] = float_or(til[2]) * -1
                    if til[0] == 'CM':
                        if til[1] == 'dimensions:':
                            dims = [int_or(til[2]), int_or(til[3])]
                            cm_array = np.zeros((dims[0], dims[1]))
                    if til[0] == 'CM:':
                        for j in range(0, dims[0]):
                            cm_array[this_row][j] = int_or(til[j+1])
                        this_row += 1

        mbs_region = region().from_list(self.infos['minmax'])
        xinc = (self.infos['minmax'][1] - self.infos['minmax'][0]) / dims[0]
        yinc = (self.infos['minmax'][2] - self.infos['minmax'][3]) / dims[1]

        if xinc >= 0 and yinc > 0:
            xcount, ycount, dst_gt = mbs_region.geo_transform(x_inc = xinc, y_inc = yinc)

            ds_config = {'nx': dims[0], 'ny': dims[1], 'nb': dims[1] * dims[0],
                         'geoT': dst_gt, 'proj': sr_wkt(self.epsg),
                         'dt': gdal.GDT_Float32, 'ndv': 0, 'fmt': 'GTiff'}

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
        else: wkt = mbs_region.export_as_wkt()

        self.infos['wkt'] = wkt
        return(self)    

## ==============================================
## xyz processing (datalists fmt:168)
## ==============================================
class xyz_point:
    """represnting an xyz data point
    
    Attributes:
      x (float): longitude/x
      y (float): latitude/y
      z (float): elevation/z
      w (float): weight/w
      epsg (int): the EPSG code representing the xy projection
      z_units (str): units of the z value
      z_datum (str): vertical datum of the z_value
    """
    
    def __init__(self, x = None, y = None, z = None, w = 1,
                 epsg = 4326, z_units = 'm', z_datum = 'msl'):
        """
        Args:
          x (float): longitude/x
          y (float): latitude/y
          z (float): elevation/z
          w (float): weight/w
          epsg (int): the EPSG code representing the xy projection
          z_units (str): units of the z value
          z_datum (str): vertical datum of the z_value
        """
        
        self.x = float_or(x)
        self.y = float_or(y)
        self.z = float_or(z)
        self.w = float_or(w)
        #self.epsg = int_or(epsg)
        #self.z_units = z_units
        #self.z_datum = z_datum

    def _valid_p(self):
        if self.x is None: return(False)
        if self.y is None: return(False)
        if self.z is None: return(False)
        return(True)
    
    def from_list(self, xyz_list, x_pos = 0, y_pos = 1, z_pos = 2, w_pos = 3):
        """load xyz data from a list

        Args:
          xyz_list (list): a list of xyz data [x,y,z,...]

        Returns:
          xyz: self
        """

        if len(xyz_list) > x_pos:
            self.x = float_or(xyz_list[x_pos])
        if len(xyz_list) > y_pos:
            self.y = float_or(xyz_list[y_pos])
        if len(xyz_list) > z_pos:
            self.z = float_or(xyz_list[z_pos])
        if len(xyz_list) > w_pos:
            self.w = float_or(xyz_list[w_pos])
        return(self)
    
    def from_string(self, xyz_str, x_pos = 0, y_pos = 1, z_pos = 2, w_pos = 3, delim = " "):
        """load xyz data from a string

        Args:
          xyz_str (str): a string representing delimited xyz data.

        Returns:
          xyz: self
        """
        
        this_line = xyz_str.strip()
        return(self.from_list(this_line.split(delim), x_pos, y_pos, z_pos, w_pos))
        
    def export_as_list(self, include_z = False, include_w = False):
        """export xyz as a list

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output

        Returns:
          list: the region values in a list
        """
        
        xyz_list = [self.x, self.y]
        if include_z:
            xyz_list.append(self.z)
        if include_w:
            xyz_list.append(self.w)
        return(xyz_list)

    def export_as_string(self, delim, include_z = False, include_w = False):
        """export xyz data as string

        Args:
          delim (str): the delimiter
        """

        l = self.export_as_list(include_z = include_z, include_w = include_w)
        return('{}\n'.format(delim.join([str(x) for x in l])))

    def export_as_wkt(self):
        return('POINT ({} {})'.format(self.x, self.y))
    
    def dump(self, delim = ' ', include_z = True, include_w = False, encode = False, dst_port = sys.stdout):
        """dump xyz as a string to dst_port

        Args:
          include_z: include the z-region in the output
          include_w: include the w-region in the output
          dst_port (port): an open destination port
          encode (bool): encode the output
            sys.stdout, etc prefer encoded data, while
            files usually like unencoded...
        """
    
        l = self.export_as_string(delim, include_z = include_z, include_w = include_w)
        if encode: l = l.encode('utf-8')
        dst_port.write(l)

    def transform(self, dst_trans):
        """transform the x/y using the dst_trans

        Args:
          dst_trans: an srs transformation object

        Returns:
          xyz: self
        """
        
        point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(self.x, self.y))
        point.Transform(dst_trans)
        self.x = point.GetX()
        self.y = point.GetY()

        return(self)
        
    def warp(self, warp_epsg = 4326):
        """transform the x/y using the warp_epsg code

        Args:
          warp_epsg (int): an epsg code to warp to

        Returns:
          xyz: self
        """

        warp_epsg = int_or(warp_epsg)
        if warp_epsg is None or self.epsg is None: return(self)

        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(self.epsg)

        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(warp_epsg)
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        
        if dst_trans is None: return(self)
        self.transform(dst_trans)
        
        return(self)

def xyz_line(xyz_line, dst_port = sys.stdout, encode = False):
    """write "xyz" `line` to `dst_port`
    `line` should be a list of xyz values [x, y, z, ...].
    
    Args:
      xyz_line (str): a string representing delimited data.
      dst_port (port): an open destination port
      encode (bool): encode the output
        sys.stdout, etc prefer encoded data, while
        files usually like unencoded...
    """
    
    delim = ' '
    
    l = '{}\n'.format(delim.join(['{:.7f}'.format(x) for x in xyz_line]))
    if encode: l = l.encode('utf-8')
    dst_port.write(l)
    
class xyz_parser:
    """representing an xyz dataset stream
    """

    def __init__(self, ds = None, name = "<xyz-dataset>", delim = None, xpos = 0, ypos = 1, zpos = 2,
                 skip = 0, x_scale = 1, y_scale = 1, z_scale = 1, x_offset = 0, y_offset = 0, epsg = 4326,
                 src_region = None, warp = None, verbose = False):
        self.ds = ds
        self.fn = ds.fn
        
        self.name = name
        self.delim = delim
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.skip = skip
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.z_scale = z_scale
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.epsg = int_or(epsg)
        self.warp = int_or(warp)
        self.region = src_region
        self.verbose = verbose
        self.infos = {}
        self.data_entries = []
        
        self._known_delims = [',', '/', ':']
        self._known_fmts = ['xyz', 'csv', 'dat', 'ascii']

    def inf(self):
        """generate a infos dictionary from the xyz dataset

        Returns:
          dict: a data-entry infos dictionary
        """
                
        pts = []
        self.infos['name'] = self.fn
        self.infos['hash'] = dl_hash(self.fn)
        self.infos['numpts'] = 0
        this_region = region()

        for i, l in enumerate(self.yield_xyz()):
            if i == 0:
                this_region.from_list([l.x, l.x, l.y, l.y, l.z, l.z])
            else:
                try:
                    if l.x < this_region.xmin: this_region.xmin = l.x
                    elif l.x > this_region.xmax:  this_region.xmax = l.x
                    if l.y < this_region.ymin: this_region.ymin = l.y
                    elif l.y > this_region.ymax: this_region.ymax = l.y
                    if l.z < this_region.zmin: this_region.zmin = l.z
                    elif l.z > this_region.zmax: this_region.zmax = l.z
                except: pass
            pts.append(l.export_as_list(include_z = True))
            self.infos['numpts'] = i

        self.infos['minmax'] = this_region.export_as_list(include_z = True)
        if self.infos['numpts'] > 0:
            try:
                out_hull = [pts[i] for i in spatial.ConvexHull(pts, qhull_options='Qt').vertices]
                out_hull.append(out_hull[0])
                self.infos['wkt'] = regions.create_wkt_polygon(out_hull, xpos = 0, ypos = 1)

            except: self.infos['wkt'] = this_region.export_as_wkt()
        return(self.infos)
        
    def line_delim(self, xyz_line):
        """guess a line delimiter
        Args:
          xyz_line (str): a string representing delimited data.

        Returns:
          str: delimiter (or None)
        """

        for delim in self._known_delims:
            this_xyz = xyz_line.split(delim)
            if len(this_xyz) > 1:
                self.delim = delim

    def parse(self):
        if self.region is not None:
            self.ds.inf()
            inf_region = region().from_list(self.ds.infos['minmax'])
            if regions_intersect_p(inf_region, self.region):
                self.data_entries = [self.ds]
        else: self.data_entries = [self.ds]
        return(self)
                
    def yield_xyz(self):
        """xyz file parsing generator

        Yields:
          xyz: xyz data
        """
        
        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                self.src_data = open(self.fn, "r")
            else: self.src_data = sys.stdin
        else: self.src_data = sys.stdin
        
        sk = self.skip
        this_xyz = xyz_point()
        if self.region is not None:
            if self.region.epsg != self.epsg:
                if self.warp is not None:
                    if self.region.epsg != self.warp:
                        self.region.warp(warp_epsg = self.epsg)
                else: self.region.warp(warp_epsg = self.epsg)

        warp_epsg = int_or(self.warp)
        if warp_epsg is not  None and self.epsg is not None:
            src_srs = osr.SpatialReference()
            src_srs.ImportFromEPSG(self.epsg)

            dst_srs = osr.SpatialReference()
            dst_srs.ImportFromEPSG(warp_epsg)
            try:
                src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            except: pass
            dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        else: dst_trans = None
        ln = 0
        for xyz_line in self.src_data:
            if ln >= sk:
                if self.delim is None: self.line_delim(xyz_line)
                this_xyz.from_string(xyz_line, delim = self.delim, x_pos = self.xpos, y_pos = self.ypos)
                if this_xyz._valid_p():
                    this_xyz.x = (this_xyz.x + self.x_offset) * self.x_scale
                    this_xyz.y = (this_xyz.y + self.y_offset) * self.y_scale
                    this_xyz.z *= self.z_scale
                    if self.region is not None and self.region._valid_p():
                        if xyz_in_region_p(this_xyz, self.region):
                            if dst_trans is not None:
                                this_xyz.transform(dst_trans)
                            ln += 1
                            yield(this_xyz)
                    else:
                        if dst_trans is not None:
                            this_xyz.transform(dst_trans)
                        ln +=1
                        yield(this_xyz)        
            else: sk -= 1
        if self.verbose: echo_msg('parsed {} data records from {}'.format(ln, self.fn))
        self.src_data.close()

## ==============================================
## raster processing (fmt: 200)
## ==============================================
class raster_parser:
    """providing a raster parser
    """

    def __init__(self, ds = None, src_region = None, mask = None, warp = None, verbose = False, step = 1, epsg = None):
        self.ds = ds
        self.fn = ds.fn
        self.mask = mask
        self.verbose = verbose
        self.step = 1
        self.epsg = int_or(epsg)
        self.warp = int_or(warp)
        self.region = src_region
        self.infos = {}
        self.data_entries = []
        
        self.src_ds = None
        self.ds_config = None
        self.ds_open_p = True
        
    def _open_ds(self):
        """open the gdal datasource and gather infos 

        Returns:
          raster_parser: self
        """
        
        if self.fn is not None:
            if os.path.exists(str(self.fn)):
                try:
                    self.src_ds = gdal.Open(self.fn)
                except: self.src_ds = None
            else: self.src_ds = None
        else: self.src_ds = None

        if self.src_ds is not None:
            self.gather_infos()
            self.ds_open_p = True
        else: self.ds_open_p = False
        return(self)

    def _close_ds(self):
        """close the gdal datasource

        Returns:
          raster_parser: self
        """
        
        self.src_ds = None
        self.ds_config = None
        self.ds_open_p = False
        return(self)

    def inf(self):
        """generate a infos dictionary from the raster dataset

        Returns:
          dict: a data-entry infos dictionary
        """

        self.infos['name'] = self.fn
        self.infos['hash'] = dl_hash(self.fn)
        self._open_ds()
        if self.ds_open_p:
            
            this_region = region().from_geo_transform(geoT = self.gt, x_count = self.x_count, y_count = self.y_count)
            try: zr = self.src_ds.GetRasterBand(1).ComputeRasterMinMax()
            except: zr = [None, None]
            this_region.zmin = zr[0]
            this_region.zmax = zr[1]
            self.infos['minmax'] = this_region.export_as_list(include_z = True)
            self.infos['numpts'] = self.x_count * self.y_count
            self.infos['wkt'] = this_region.export_as_wkt()
        self._close_ds()
        return(self.infos)
    
    def gather_infos(self, scan = False):
        """gather information from `src_ds` GDAL dataset

        Returns:
          raster_parser: self
        """

        if self.ds_open_p:
            gt = self.src_ds.GetGeoTransform()
            if self.region is not None and self.region._valid_p(check_xy = True):
                self.srcwin = self.region.srcwin(gt, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
            else: self.srcwin = (0, 0, self.src_ds.RasterXSize, self.src_ds.RasterYSize)
            src_band = self.src_ds.GetRasterBand(1)
            self.gt = (gt[0] + (self.srcwin[0] * gt[1]), gt[1], 0., gt[3] + (self.srcwin[1] * gt[5]), 0., gt[5])

            proj = self.src_ds.GetProjectionRef()
            src_srs = osr.SpatialReference()
            src_srs.ImportFromWkt(proj)
            src_srs.AutoIdentifyEPSG()
            srs_auth = src_srs.GetAuthorityCode(None)

            self.epsg = int_or(srs_auth)
            self.x_count = self.srcwin[2]
            self.y_count = self.srcwin[3]
            self.dt = src_band.DataType
            self.dtn = gdal.GetDataTypeName(src_band.DataType)
            self.ndv = src_band.GetNoDataValue()
            if self.ndv is None: self.ndv = -9999
            self.fmt = self.src_ds.GetDriver().ShortName
            self.zr = None

            if scan:
                src_arr = src_band.ReadAsArray(srcwin[0], self.srcwin[1], self.srcwin[2], self.srcwin[3])
                self.zr = (np.amin(src_arr), np.amax(src_arr))
                src_arr = None
        return(self)

    def parse(self):
        if self.region is not None:
            self.ds.inf()
            inf_region = region().from_list(self.ds.infos['minmax'])
            if regions_intersect_p(inf_region, self.region):
                self.data_entries = [self.ds]
        else: self.data_entries = [self.ds]
        return(self)
    
    def yield_xyz(self):
        """parse the data from gdal dataset src_ds (first band only)

        Yields:
          xyz: the parsed xyz data
        """

        self._open_ds()
        out_xyz = xyz_point()
        if self.src_ds is not None:
            ln = 0
            band = self.src_ds.GetRasterBand(1)
            gt = self.gt
            warp_epsg = self.warp
            
            if warp_epsg is not  None and self.epsg is not None:
                src_srs = osr.SpatialReference()
                src_srs.ImportFromEPSG(self.epsg)

                dst_srs = osr.SpatialReference()
                dst_srs.ImportFromEPSG(warp_epsg)
                try:
                    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                    dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                except: pass
                dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
            else: dst_trans = None
                        
            msk_band = None
            if self.mask is not None:
                src_mask = gdal.Open(self.mask)
                msk_band = src_mask.GetRasterBand(1)

            nodata = ['{:g}'.format(-9999), 'nan', float('nan')]
            if self.ndv is not None: nodata.append('{:g}'.format(self.ndv))
            for y in range(self.srcwin[1], self.srcwin[1] + self.srcwin[3], 1):
                band_data = band.ReadAsArray(self.srcwin[0], y, self.srcwin[2], 1)
                if self.region is not None:
                    z_region = self.region.z_region()
                    if z_region[0] is not None:
                        band_data[band_data < z_region[0]] = -9999
                    if z_region[1] is not None:
                        band_data[band_data > z_region[1]] = -9999
                if msk_band is not None:
                   msk_data = msk_band.ReadAsArray(self.srcwin[0], y, self.srcwin[2], 1)
                   band_data[msk_data==0]=-9999
                band_data = np.reshape(band_data, (self.srcwin[2], ))
                for x_i in range(0, self.srcwin[2], 1):
                    x = x_i + self.srcwin[0]
                    z = band_data[x_i]
                    if '{:g}'.format(z) not in nodata:
                        ln += 1
                        out_xyz.x, out_xyz.y = _pixel2geo(x, y, gt)
                        out_xyz.z = z
                        if dst_trans is not None: out_xyz.transform(dst_trans)
                        yield(out_xyz)
            band = None
            src_mask = None
            msk_band = None
            if self.verbose: echo_msg('parsed {} data records from {}'.format(ln, self.fn))
        self._close_ds()

## ==============================================
## datalist processing (fmt: -1)
## ==============================================        
class datalist_parser:
    """representing a datalist parser
    
    A datalist is an MB-System style datalist, or
    a archive (zip/tar.gz/gz) or just a datfile.

    """

    def __init__(self, ds = None, src_region = None, verbose = False, epsg = None, warp = None):
        self.ds = ds
        self.fn = ds.fn
        self.name = '.'.join(self.fn.split('.')[:-1])
        self.weight = ds.weight
        self.region = src_region
        self.data_entries = []
        self.infos = {}
        self.epsg = int_or(epsg)
        self.warp = int_or(warp)
        self.verbose = verbose
        
    def parse(self):
        """import a datalist entry from a string
    
        Returns:
          datalist_parser: self
        """
        
        _prog = _progress("parsing datalist {}".format(self.fn))
        these_entries = []
        this_dir = os.path.dirname(self.fn)
        if os.path.exists(self.fn):
            with open(self.fn, 'r') as op:
                for this_line in op:
                    _prog.update()
                    if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                        data_set = xyz_dataset(parent = self.ds, weight = self.weight).from_string(this_line, this_dir)
                        if data_set.valid_p():
                            dls = data_set.data_types[data_set.data_format]['parser'](
                                data_set,
                                {'src_region': self.region, 'verbose': self.verbose, 'epsg': self.epsg, 'warp': self.warp}
                            )                           
                            for entry in dls.data_entries:
                                these_entries.append(entry)
                        #else: echo_warning_msg('invalid dataset: `{}`'.format(data_set.fn))
        else: echo_warning_msg('could not open datalist/entry {}'.format(self.fn))
        self.data_entries = these_entries
        _prog.end(0, "parsed datalist {}".format(self.fn))
        return(self)
    
    def yield_xyz(self):
        """parse the data from the datalist

        Yields:
          xyz: the parsed xyz data
        """

        if self.verbose: _prog = _progress('parsing data from {}'.format(self.fn))
        for i, this_entry in enumerate(self.data_entries):
            if self.verbose: _prog.update_perc((i, len(self.data_entries)))
            for xyz in this_entry.yield_xyz(src_region = self.region, verbose = self.verbose, epsg = self.epsg, warp = self.warp):
                yield(xyz)
        if self.verbose: _prog.end(0, 'parsed data from {}'.format(self.fn))
        
    def inf(self):
        """return the region of the datalist and generate
        an associated `.inf` file if `inf_file` is True.

        Args:
          dl (str): a datalist pathname
          inf_file (bool): generate an inf file
          epsg (int): EPSG code
          overwrite (bool): overwrite a possibly existing inf_file

        Returns:
          list: the region [xmin, xmax, ymin, ymax]
        """

        out_regions = []
        self.infos['name'] = self.fn
        self.infos['numpts'] = 0
        self.infos['hash'] = dl_hash(self.fn)
        for entry in self.data_entries:
            entry.inf()
            out_regions.append(entry.infos['minmax'])
            self.infos['numpts'] += entry.infos['numpts']

        l = 0
        for this_region in out_regions:
            if l == 0:
                tmp_region = region().from_list(this_region)
                if tmp_region._valid_p():
                    out_region = region().from_list(this_region)
                    l += 1
            else:
                tmp_region = region().from_list(this_region)
                if tmp_region._valid_p():
                    out_region = regions_merge(out_region, tmp_region)
        self.infos['minmax'] = out_region.export_as_list(include_z = True)
        self.infos['wkt'] = out_region.export_as_wkt()
        return(self.infos)

## ==============================================
## datalist dataset
## ==============================================
class xyz_dataset:
    """representing an xyz-able parser or a data-list entry
    """

    data_types = {
        -1: {'name': 'datalist',
            'fmts': ['datalist', 'mb-1'],
             'parser': lambda s, x: datalist_parser(ds = s, **x).parse(),
             },
        167: {'name': 'yxz',
              'fmts': ['yxz'],
              'parser': lambda s, x: xyz_parser(ds = s, xpos = 1, ypos = 0, **x).parse(),
              },
        168: {'name': 'xyz',
              'fmts': ['xyz', 'csv', 'dat', 'ascii', 'txt'],
              'parser': lambda s, x: xyz_parser(ds = s, **x).parse(),
              },
        200: {'name': 'raster',
            'fmts': ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'],
              'parser': lambda s, x: raster_parser(ds = s, **x).parse(),
              },
    }
    
    def __init__(self, fn = None, data_format = None, weight = 1, epsg = 4326, name = "<xyz-dataset>", parent = 'top-level', verbose = False):

        self.fn = fn
        self.name = name
        self.data_format = data_format
        self.epsg = epsg
        self.infos = {}
        self.parent = parent
        self.weight = weight
        self.verbose = verbose

    def guess_data_format(self):
        if self.fn is not None:
            for key in self.data_types.keys():
                if self.fn.split('.')[-1] in self.data_types[key]['fmts']:
                    self.data_format = key
                    break
                
    def valid_p(self):
        """check if self is a valid datalist entry

        Returns:
          bools: True if valid else False
        """
        
        if self.fn is None: return(False)
        if self.data_format is None: return(False)
        if self.fn is not None:
            if not os.path.exists(self.fn): return (False)
            if os.stat(self.fn).st_size == 0: return(False)
        return(True)

    def _hash(self, sha1 = False):
        """generate a hash of the xyz-dataset source file

        Returns:
          str: hexdigest
        """
        
        BUF_SIZE = 65536  # lets read stuff in 64kbchunks!
        if sha1: this_hash = hashlib.sha1()
        else: this_hash = hashlib.md5()
        
        with open(self.fn, 'rb') as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data:
                    break
                this_hash.update(data)
            
        return(this_hash.hexdigest())
    
    def echo(self, **kwargs):
        """print self as a datalist entry string.
        """

        dls = self.data_types[self.data_format]['parser'](self, kwargs)
        for entry in dls.data_entries:
            print('{} {} {}'.format(entry.fn, entry.data_format, entry.weight))
                    
    def from_string(self, dl_e, directory = None):
        """import a datalist entry string 'fn fmt wt' to an xyz-dataset object
        minimally, this is a string of a filename/fetch module.

        Returns:
          xyz_dataset: self
        """

        this_entry = dl_e.rstrip().split()
        try:
            entry = [x if n == 0 else int_or(x) if n < 2 else float_or(x) if n < 3 else x for n, x in enumerate(this_entry)]
        except Exception as e:
            echo_error_msg('could not parse entry {}'.format(dl_e))
            return(self)
        if len(entry) < 2:
            for key in self.data_types.keys():
                se = entry[0].split('.')
                see = se[-1] if len(se) > 1 else entry[0].split(":")[0]
                if see in self.data_types[key]['fmts']:
                    entry.append(int(key))
                    break
            if len(entry) < 2:
                echo_error_msg('could not parse entry {}'.format(dl_e))
                return(self)
        if len(entry) < 3: entry.append(self.weight)

        self.fn = entry[0] if directory is None else os.path.join(directory, entry[0])
        self.data_format = entry[1]
        if self.data_format is None:
            self.guess_data_format()
        
        if self.weight is not None:
            self.weight *= entry[2]
        return(self)
        
    def inf(self, **kwargs):
        """read/write an inf file

        Args:
          kwargs (dict): any arguments to pass to the dataset parser
        
        Returns:
          dict: xyz-dataset info dictionary
        """
        
        inf_path = '{}.inf'.format(self.fn)
        if os.path.exists(inf_path):
            try:
                with open(inf_path) as i_ob:
                    self.infos = json.load(i_ob)
            except ValueError:
                try:
                    self.infos = mbs_parser(fn = inf_path, epsg = self.epsg).inf_parse().infos
                except:
                    echo_error_msg('failed to parse inf {}'.format(inf_path))
            except: echo_error_msg('failed to parse inf {}'.format(inf_path))
        else: self.infos = {}
        #if 'hash' not in self.infos.keys() or self._hash() != self.infos['hash']:
        if 'hash' not in self.infos.keys():
            echo_msg("generating inf for {}".format(self.fn))
            self.infos = self.data_types[self.data_format]['parser'](self, kwargs).inf()
            self.infos['format'] = self.data_format
            if self.infos['minmax'] is not None:
                with open('{}.inf'.format(self.fn), 'w') as inf:
                    inf.write(json.dumps(self.infos))
            if self.parent != 'top-level':
                remove_glob('{}.inf'.format(self.parent.fn))
                self.parent.infos = {}
            self.infos['epsg'] = self.epsg
        self.infos['format'] = self.data_format
        return(self.infos)
        
    def yield_xyz(self, **kwargs):
        """parse the xyz data from the xyz-dataset, yield results

        Args:
          kwargs (dict): any arguments to pass to the dataset parser

        Yields:
          xyz: xyz data
        """
        
        for xyz in self.data_types[self.data_format]['parser'](self, kwargs).yield_xyz():
            xyz.w = self.weight
            yield(xyz)

    def dump(self, dst_port = sys.stdout, encode = False, **kwargs):
        """dump the xyz-dataset to dst_port
        
        Args:
          dst_port (port): an open file port to write to
          encode (bool): True to encode the output
          kwargs (dict): any arguments to pass to the dataset parser
        """
        
        for this_xyz in self.yield_xyz(**kwargs):
            this_xyz.dump(include_w = True if self.weight is not None else False, dst_port = dst_port, encode = False)
            
_datalist_fmts_short_desc = lambda: '\n  '.join(['{}\t{}'.format(key, xyz_dataset().data_types[key]['fmts']) for key in xyz_dataset().data_types])
    
## ==============================================
## dadtalists cli
## ==============================================    
datalists_version = '0.1.0'
datalists_usage = '''{cmd} ({dl_version}): DataLists IMproved; Process and generate datalists

usage: {cmd} [ -ghiqwPRW [ args ] ] DATALIST ...

Options:
  -R, --region\t\tRestrict processing to the desired REGION 
\t\t\tWhere a REGION is [ xmin/xmax/ymin/ymax/[ zmin/zmax/[ wmin/wmax ] ] ]
\t\t\tUse '-' to indicate no bounding range; e.g. -R -/-/-/-/-10/10/1/-
\t\t\tOR an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied, will use each region found therein.
\t\t\tAppend :zmin/zmax/[ wmin/wmax ] to the file path to extended region.
  -P, --s_epsg\t\tSet the projection EPSG code of the datalist
  -W, --t_epsg\t\tSet the output warp projection EPSG code

  --glob\t\tGlob the datasets in the current directory to stdout
  --info\t\tGenerate and return an info dictionary of the dataset
  --weights\t\tOutput weight values along with xyz
  --quiet\t\tLower the verbosity to a quiet

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Supported datalist formats: 
  {dl_formats}

Examples:
  % {cmd} my_data.datalist -R -90/-89/30/31
  % {cmd} -R-90/-89/30/31/-100/100 *.tif -l -w > tifs_in_region.datalist
  % {cmd} -R my_region.shp my_data.xyz -w -s_epsg 4326 -t_epsg 3565 > my_data_3565.xyz

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format(cmd =  os.path.basename(sys.argv[0]), 
           dl_version = datalists_version,
           dl_formats =  _datalist_fmts_short_desc())

def datalists_cli(argv = sys.argv):
    """run datalists from command-line

    See `datalists_cli_usage` for full cli options.
    """

    dls = []
    epsg = None
    warp = None
    i_regions = []
    these_regions = []
    want_weights = False
    want_inf = False
    want_list = False
    want_glob = False
    want_verbose = True
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            i_regions.append(str(argv[i + 1]))
            i = i + 1
        elif arg[:2] == '-R':
            i_regions.append(str(arg[2:]))
        elif arg == '-s_epsg' or arg == '--s_epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg == '-t_epsg' or arg == '--t_epsg' or arg == '-W':
            warp = argv[i + 1]
            i = i + 1
        elif arg == '--weights' or arg == '-w':
            want_weights = True
        elif arg == '--info' or arg == '-i':
            want_inf = True
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--glob' or arg == '-g':
            want_glob = True
        elif arg == '--quiet' or arg == '-q':
            want_verbose = False
        elif arg == '--help' or arg == '-h':
            print(datalists_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), datalists_version))
            sys.exit(1)
        elif arg[0] == '-':
            print(datalists_usage)
            sys.exit(0)
        else: dls.append(arg)
        i = i + 1

    if want_glob:
        for key in xyz_dataset().data_types.keys():
            if key != -1:
                for f in xyz_dataset().data_types[key]['fmts']:
                    globs = glob.glob('*.{}'.format(f))
                    [sys.stdout.write('{}\n'.format(' '.join([x, str(key), '1']))) for x in globs]
        sys.exit(0)

    for i_region in i_regions:
        tmp_region = region().from_string(i_region)
        if tmp_region._valid_p():
            these_regions.append(tmp_region)
        else:
            tmp_region = ogr_wkts(i_region)
            for i in tmp_region:
                if i._valid_p():
                    these_regions.append(i)
                    
    if len(these_regions) == 0:
        these_regions = [None]
    else:
        if want_verbose: echo_msg('parsed {} region(s)'.format(len(these_regions)))
    
    for rn, this_region in enumerate(these_regions):
        if len(dls) == 0:
            print(datalists_usage)
            echo_error_msg("you must specify some type of data")
        xdls = [xyz_dataset().from_string(" ".join(["-" if x == "" else x for x in dl.split(":")])) for dl in dls]
        for xdl in xdls:
            if xdl.valid_p():
                if not want_weights: xdl.weight = None
                if want_inf: print(xdl.inf())
                elif want_list: xdl.echo(src_region = this_region, verbose = want_verbose)
                else: xdl.dump(src_region = this_region, verbose = want_verbose, epsg = epsg, warp = warp)
### End
