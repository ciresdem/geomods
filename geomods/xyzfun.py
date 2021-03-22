### xyzfun.py
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

import sys
import copy
import json

## ==============================================
## import gdal/numpy
## ==============================================
import gdal
import ogr
import osr
import numpy as np
from scipy import spatial

## ==============================================
## import geomods
## ==============================================
from geomods import utils
from geomods import regions
from geomods import mbsfun

## ==============================================
## xyz processing (datalists fmt:168)
## ==============================================
_xyz_config = {
    'delim': None,
    'xpos': 0,
    'ypos': 1,
    'zpos': 2,
    'skip': 0,
    'z-scale': 1,
    'x-off': 0,
    'y-off': 0,
    'name': '<xyz-data-stream>',
    'upper_limit': None,
    'lower_limit': None,
    'epsg': 4326,
    'warp': None,
    'verbose': False,
}

#_known_delims = [',', ' ', '\t', '/', ':']
_known_delims = [',', '/', ':']
_known_xyz_fmts = ['xyz', 'csv', 'dat', 'ascii']

def xyz_line_delim(xyz_line):
    """guess a line delimiter

    Args:
      xyz_line (str): a string representing delimited data.
    
    Returns:
      str: delimiter (or None)
    """
    
    for delim in _known_delims:
        this_xyz = xyz_line.split(delim)
        if len(this_xyz) > 1: return(delim)
        
    return(None)

def xyz_warp(xyz, dst_trans):
    """transform the x/y using the dst_trans

    Args:
      xyz (list): xyz data
      dst_trans (srs-trans): srs transformation object

    Returns:
      list: xyz data; [x, y, z]
    """
    
    if dst_trans is None: return(xyz)
    point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(xyz[0], xyz[1]))
    point.Transform(dst_trans)
    return([point.GetX(), point.GetY(), xyz[2]])

def xyz_parse_line(xyz_line, xyz_c = _xyz_config):
    """parse an xyz line-string, using _xyz_config

    Args:
      xyz_line (str): a string representing delimited data.
      xyz_c (dict): xyz config dictionary

    Returns:
      list: xyz data; [x, y, z]
    """
    
    try:
        this_line = xyz_line.strip()
    except AttributeError as e:
        utils.echo_error_msg('input is list, should be xyz-line: {}'.format(e))
        return(None)
    except Exception as e:
        utils.echo_error_msg(e)
        return(None)
    if xyz_c['delim'] is None:
        xyz_c['delim'] = xyz_line_delim(this_line)
    this_xyz = this_line.split(xyz_c['delim'])
    try:
        o_xyz = [float(this_xyz[xyz_c['xpos']]) + float(xyz_c['x-off']),
                 float(this_xyz[xyz_c['ypos']]),
                 float(this_xyz[xyz_c['zpos']]) * float(xyz_c['z-scale'])]
    except IndexError as e:
        if xyz_c['verbose']: utils.echo_error_msg(e)
        return(None)
    except Exception as e:
        if xyz_c['verbose']: utils.echo_error_msg(e)
        return(None)
    return(o_xyz)
   
def xyz_parse(src_xyz, xyz_c = _xyz_config, region = None, verbose = False):
    """xyz file parsing generator

    Args:
      src_xyz (generataor): list/generator of xyz data
      xyz_c (dict): xyz config dictionary
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      verbose (bool): increase verbosity

    Yields:
      list: xyz data [x, y, z, ...]
    """
    
    ln = 0
    pass_d = True
    skip = int(xyz_c['skip'])
    
    if xyz_c['epsg'] == xyz_c['warp'] or xyz_c['epsg'] is None:
        xyz_c['warp'] = None
    
    if xyz_c['warp'] is not None:
        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(int(xyz_c['epsg']))
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(int(xyz_c['warp']))
    
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except: pass
        
        dst_trans = osr.CoordinateTransformation(src_srs, dst_srs)
        
    else: src_srs = dst_srs = dst_trans = None
    
    for xyz in src_xyz:
        pass_d = True
        
        if ln >= skip:
            this_xyz = xyz_parse_line(xyz, xyz_c)
            
            if this_xyz is not None:
                if xyz_c['warp'] is not None:
                    this_xyz = xyz_warp(this_xyz, dst_trans)
                    
                if region is not None:
                    if regions.region_valid_p:
                        if not xyz_in_region_p(this_xyz, region):
                            pass_d = False
                    
                if xyz_c['upper_limit'] is not None or xyz_c['lower_limit'] is not None:
                    if not regions.z_pass(this_xyz[2], upper_limit = xyz_c['upper_limit'], lower_limit = xyz_c['lower_limit']):
                        pass_d = False
                        
            else: pass_d = False
            
            if pass_d:
                ln += 1
                yield(this_xyz)
                
        else: skip -= 1
        
    if verbose:
        if ln == 0:
            status = -1            
        else:
            status = 0
            
        utils.echo_msg('parsed {} data records from {}'.format(ln, xyz_c['name']))

def xyz_dump(src_xyz, xyz_c = _xyz_config, region = None, verbose = False, dst_port = sys.stdout):
    """dump the xyz data from the xyz datalist entry to dst_port

    Args:
      src_xyz (generataor): list/generator of xyz data
      xyz_c (dict): xyz config dictionary
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      verbose (bool): increase verbosity    
      dst_port (port): an open destination port
    """
    
    for xyz in xyz_parse(src_xyz, xyz_c, region, verbose):
        xyz_line(xyz, dst_port, True)
        
def xyz2py(src_xyz):
    """return src_xyz as a python list

    Args:
      src_xyz (generataor): list/generator of xyz data

    Returns:
      list: list of xyz data
    """
    
    xyzpy = []
    return([xyzpy.append(xyz) for xyz in xyz_parse(src_xyz)])

def xyz_block(src_xyz, region, inc, weights = False, verbose = False):
    """block the src_xyz data to the mean block value

    Args:
      src_xyz (generataor): list/generator of xyz data
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      inc (float): blocking increment, in native units
      weights (bool): block using weights
      verbose (bool): increase verbosity    

    Yields:
      list: xyz data for each block with data
    """
    
    xcount, ycount, dst_gt = regions.region2gt(region, inc)
    sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    ptArray = np.zeros((ycount, xcount))
    if weights: wtArray = np.zeros((ycount, xcount))
    if verbose: utils.echo_msg('blocking data to {}/{} grid'.format(ycount, xcount))
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if weights: z * this_xyz[3]
        #w = this_xyz[3]
        #z = z * w
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = utils._geo2pixel(x, y, dst_gt)
                try:
                    sumArray[ypos, xpos] += z
                    ptArray[ypos, xpos] += 1
                    if weights: wtArray[ypos, xpos] += this_xyz[3]
                except: pass
    ptArray[ptArray == 0] = np.nan
    if weights:
        wtArray[wtArray == 0] = 1
        outarray = (sumArray / wtArray) / ptArray
    else: outarray = sumArray / ptArray

    sumArray = ptArray = None
    if weights: wtArray = None

    outarray[np.isnan(outarray)] = -9999
    
    for y in range(0, ycount):
        for x in range(0, xcount):
            geo_x, geo_y = utils._pixel2geo(x, y, dst_gt)
            z = outarray[y,x]
            if z != -9999:
                yield([geo_x, geo_y, z])
    
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
    
    delim = _xyz_config['delim'] if _xyz_config['delim'] is not None else ' '
    
    l = '{}\n'.format(delim.join([str(x) for x in xyz_line]))
    if encode: l = l.encode('utf-8')
    dst_port.write(l)

def xyz2wkt(xyz):
    return('POINT ({} {})'.format(xyz[0], xyz[1]))
    
def xyz_in_region_p(xyz, region):
    """check if xyz point in inside the given region

    Args:
      xyz (list): an xyz data list
      region (list): a `region` list [xmin, xmax, ymin, ymax]

    Returns:
      bool: True if xyz point inside region else False
    """
    
    if regions.region_wkt_p(region):
        xyz_wkt = xyz2wkt(xyz)
        p_geom = ogr.CreateGeometryFromWkt(xyz_wkt)
        r_geom = ogr.CreateGeometryFromWkt(region)

        return(p_geom.Within(r_geom))
    
    if xyz[0] < region[0]: return(False)
    elif xyz[0] > region[1]: return(False)
    elif xyz[1] < region[2]: return(False)
    elif xyz[1] > region[3]: return(False)
    else: return(True)

def xyz_inf(src_xyz):
    """generate and return or read and return an xyz inf file.

    Args:
      src_xyz (generataor): list/generator of xyz data

    Returns:
      dict: an xyz infos dictionary
    """
    
    pts = []
    xyzi = {}
    xyzi['name'] = src_xyz.name
    xyzi['numpts'] = 0
    xyzi['minmax'] = [0, 0, 0, 0, 0, 0]
    
    for i, l in enumerate(xyz_parse(src_xyz)):
        if i == 0:
            xyzi['minmax'] = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            try:
                if l[0] < xyzi['minmax'][0]: xyzi['minmax'][0] = l[0]
                elif l[0] > xyzi['minmax'][1]: xyzi['minmax'][1] = l[0]
                if l[1] < xyzi['minmax'][2]: xyzi['minmax'][2] = l[1]
                elif l[1] > xyzi['minmax'][3]: xyzi['minmax'][3] = l[1]
                if l[2] < xyzi['minmax'][4]: xyzi['minmax'][4] = l[2]
                elif l[2] > xyzi['minmax'][5]: xyzi['minmax'][5] = l[2]
            except: pass
        pts.append(l)
        xyzi['numpts'] = i

    if xyzi['numpts'] > 0:
        try:
            out_hull = [pts[i] for i in spatial.ConvexHull(pts, qhull_options='Qt').vertices]
            out_hull.append(out_hull[0])
            xyzi['wkt'] = regions.create_wkt_polygon(out_hull, xpos = 0, ypos = 1)

            with open('{}.inf'.format(src_xyz.name), 'w') as inf:
                inf.write(json.dumps(xyzi))
        except: xyzi['wkt'] = regions.region2wkt(xyzi['minmax'])
    return(xyzi)
            
def xyz_inf_entry(entry):
    """find the region of the xyz datalist entry.
    
    Args:
      entry (list): a datalist entry [path, fmt, weight, ...]

    Returns:
      list: the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry
    """
    
    with open(entry[0]) as infile:
        try:
            minmax = mbsfun.mb_inf(infile)
        except: minmax = xyz_inf(infile)
    return(minmax)        

def xyz_yield_entry(entry, region = None, verbose = False, z_region = None, epsg = None):
    """yield the xyz data from the xyz datalist entry

    Args:
      entry (list): a datalist entry [path, fmt, weight, ...]
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      verbose (bool): increase verbosity
      z_region (str): 'min-z/max-z'
      epsg (int): an output EPSG for xyz data

    Yields:
      list: xyz data list [x, y, z, <w, ...>]
    """
    
    xyzc = copy.deepcopy(_xyz_config)
    xyzc['name'] = entry[0]
    if z_region is not None and len(z_region) >= 2:
        xyzc['lower_limit'] = z_region[0]
        xyzc['upper_limit'] = z_region[1]
    
    with open(entry[0]) as infile:
        for line in xyz_parse(infile, xyz_c = xyzc, region = region, verbose = verbose):
            yield(line + [entry[2]] if entry[2] is not None else line)
    
def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, z_region = None, epsg = None):
    """dump the xyz data from the xyz datalist entry to dst_port

    Args:
      entry (list): a datalist entry [path, fmt, weight, ...]
      region (list): a `region` list [xmin, xmax, ymin, ymax]
      verbose (bool): increase verbosity
      z_region (str): 'min-z/max-z'
      epsg (int): an output EPSG for xyz data
    """
    
    for xyz in xyz_yield_entry(entry, region, verbose, z_region, epsg):
        xyz_line(xyz, dst_port, True, None)        
### End
