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
import numpy as np
from scipy import spatial
import gdal
import ogr
import osr

from geomods import utils
from geomods import regions
from geomods import gdalfun
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
    'verbose': False,
    'epsg': 4326,
    'warp': None,
}

#_known_delims = [',', ' ', '\t', '/', ':']
_known_delims = [',', '/', ':']

def xyz_line_delim(xyz):
    for delim in _known_delims:
        this_xyz = xyz.split(delim)
        if len(this_xyz) > 1: return(delim)
    return(None)

def xyz_warp(xyz, dst_trans):
    '''transform the x/y using the dst_trans'''
    
    if dst_trans is None: return(xyz)
    #if xyz is None: return(xyz)
    point = ogr.CreateGeometryFromWkt('POINT ({} {})'.format(xyz[0], xyz[1]))
    point.Transform(dst_trans)
    return([point.GetX(), point.GetY(), xyz[2]])

def xyz_parse_line(xyz, xyz_c = _xyz_config):
    '''parse an xyz line-string, using _xyz_config

    returns [x, y, z]'''

    this_line = xyz.strip()
    if xyz_c['delim'] is None:
        xyz_c['delim'] = xyz_line_delim(this_line)
    this_xyz = this_line.split(xyz_c['delim'])
    try:
        o_xyz = [float(this_xyz[xyz_c['xpos']]) + float(xyz_c['x-off']), float(this_xyz[xyz_c['ypos']]), float(this_xyz[xyz_c['zpos']]) * float(xyz_c['z-scale'])]
    except IndexError as e:
        if xyz_c['verbose']: utils.echo_error_msg(e)
        return(None)
    except Exception as e:
        if xyz_c['verbose']: utils.echo_error_msg(e)
        return(None)
    return(o_xyz)
   
def xyz_parse(src_xyz, xyz_c = _xyz_config, region = None, verbose = False):
    '''xyz file parsing generator
    `src_xyz` is a file object or list of xyz data.

    yields each xyz line as a list [x, y, z, ...]'''
    
    ln = 0
    skip = int(xyz_c['skip'])
    xpos = xyz_c['xpos']
    ypos = xyz_c['ypos']
    zpos = xyz_c['zpos']
    verbose = xyz_c['verbose']
    pass_d = True

    if xyz_c['epsg'] == xyz_c['warp'] or xyz_c['epsg'] is None: xyz_c['warp'] = None
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
    
    #if verbose: utils.echo_msg('parsing xyz data from {}...'.format(xyz_c['name']))
    for xyz in src_xyz:
        pass_d = True
        if ln >= skip:
            this_xyz = xyz_parse_line(xyz, xyz_c)
            if this_xyz is not None:
                if xyz_c['warp'] is not None:
                    this_xyz = xyz_warp(this_xyz, dst_trans)
                if region is not None:
                    if not xyz_in_region_p(this_xyz, region): pass_d = False
                if xyz_c['upper_limit'] is not None or xyz_c['lower_limit'] is not None:
                    if not regions.z_pass(this_xyz[2], upper_limit = xyz_c['upper_limit'], lower_limit = xyz_c['lower_limit']): pass_d = False
            else: pass_d = False
            if pass_d:
                ln += 1
                yield(this_xyz)
        else: skip -= 1
    if verbose: utils.echo_msg('parsed {} data records from {}'.format(ln, xyz_c['name']))

def xyz2py(src_xyz):
    '''return src_xyz as a python list'''
    
    xyzpy = []
    return([xyzpy.append(xyz) for xyz in xyz_parse(src_xyz)])

def xyz_block(src_xyz, region, inc, dst_xyz = sys.stdout, weights = False, verbose = False):
    '''block the src_xyz data to the mean block value

    yields the xyz value for each block with data'''
    
    xcount, ycount, dst_gt = gdalfun.gdal_region2gt(region, inc)
    sumArray = np.zeros((ycount, xcount))
    gdt = gdal.GDT_Float32
    ptArray = np.zeros((ycount, xcount))
    if weights: wtArray = np.zeros((ycount, xcount))
    if verbose: utils.echo_msg('blocking data to {}/{} grid'.format(ycount, xcount))
    for this_xyz in src_xyz:
        x = this_xyz[0]
        y = this_xyz[1]
        z = this_xyz[2]
        if weights:
            w = this_xyz[3]
            z = z * w
        if x > region[0] and x < region[1]:
            if y > region[2] and y < region[3]:
                xpos, ypos = gdalfun._geo2pixel(x, y, dst_gt)
                try:
                    sumArray[ypos, xpos] += z
                    ptArray[ypos, xpos] += 1
                    if weights: wtArray[ypos, xpos] += w
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
            geo_x, geo_y = gdalfun._pixel2geo(x, y, dst_gt)
            z = outarray[y,x]
            if z != -9999:
                yield([geo_x, geo_y, z])
    
def xyz_line(line, dst_port = sys.stdout, encode = False):
    '''write "xyz" `line` to `dst_port`
    `line` should be a list of xyz values [x, y, z, ...].'''
    delim = _xyz_config['delim'] if _xyz_config['delim'] is not None else ' '
    
    l = '{}\n'.format(delim.join([str(x) for x in line]))
    if encode: l = l.encode('utf-8')
    dst_port.write(l)

def xyz_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n], else False.'''
    
    if src_xy[0] < src_region[0]: return(False)
    elif src_xy[0] > src_region[1]: return(False)
    elif src_xy[1] < src_region[2]: return(False)
    elif src_xy[1] > src_region[3]: return(False)
    else: return(True)

def xyz_inf(src_xyz):

    pts = []
    xyzi = {}
    xyzi['name'] = src_xyz.name
    xyzi['numpts'] = 0
    xyzi['minmax'] = [0, 0, 0, 0, 0, 0]
    
    #utils.echo_msg('generating inf file for {}'.format(src_xyz.name))
    
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
            xyzi['wkt'] = gdalfun.gdal_create_polygon(out_hull, xpos = 0, ypos = 1)

            with open('{}.inf'.format(src_xyz.name), 'w') as inf:
                inf.write(json.dumps(xyzi))
        except: xyzi['wkt'] = gdalfun.gdal_region2wkt(xyzi['minmax'])
    return(xyzi)
            
def xyz_inf_entry(entry):
    '''find the region of the xyz datalist entry
    
    returns the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry'''
    
    with open(entry[0]) as infile:
        #try:
        minmax = mbsfun.mb_inf(infile)
        #except: minmax = xyz_inf(infile)
    return(minmax)        

def xyz_yield_entry(entry, region = None, verbose = False, z_region = None, epsg = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''

    xyzc = copy.deepcopy(_xyz_config)
    xyzc['name'] = entry[0]
    if z_region is not None and len(z_region) >= 2:
        xyzc['lower_limit'] = z_region[0]
        xyzc['upper_limit'] = z_region[1]
    
    with open(entry[0]) as infile:
        for line in xyz_parse(infile, xyz_c = xyzc, region = region, verbose = verbose):
            yield(line + [entry[2]] if entry[2] is not None else line)
    
def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, z_region = None):
    '''dump the xyz data from the xyz datalist entry to dst_port'''
    
    for xyz in xyz_yield_entry(entry, region, verbose, z_region):
        xyz_line(xyz, dst_port, True, None)

def gdal_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, epsg = None, z_region = None):
    '''dump the xyz data from the gdal entry to dst_port'''
    
    for xyz in gdalfun.gdal_yield_entry(entry, region, verbose, epsg, z_region):
        xyz_line(xyz, dst_port, True)
        
### End
