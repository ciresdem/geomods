### xyzfun.py
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

import sys
import copy

from geomods import utils

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
    'name': '<xyz-data-stream>',
    'upper_limit': None,
    'lower_limit': None}

_known_delims = [',', ' ', '\t', '/', ':']

def xyz_line_delim(xyz):
    for delim in _known_delims:
        this_xyz = xyz.split(delim)
        if len(this_xyz) > 1: return(delim)
    return(None)

def xyz_parse_line(xyz, xyz_c = _xyz_config):
    '''parse an xyz line-string, using _xyz_config

    returns [x, y, z]'''

    this_line = xyz.strip()
    if xyz_c['delim'] is None:
        xyz_c['delim'] = xyz_line_delim(this_line)
    this_xyz = this_line.split(xyz_c['delim'])
    try:
        o_xyz = [float(this_xyz[xyz_c['xpos']]), float(this_xyz[xyz_c['ypos']]), float(this_xyz[xyz_c['zpos']]) * float(xyz_c['z-scale'])]
    except IndexError as e:
        print(this_xyz)
        utils.echo_error_msg(e)
        return(None)
    except Exception as e:
        print(this_xyz)
        utils.echo_error_msg(e)
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
    pass_d = True
    #if verbose: utils.echo_msg('parsing xyz data from {}...'.format(xyz_c['name']))
    for xyz in src_xyz:
        pass_d = True
        if ln >= skip:
            this_xyz = xyz_parse_line(xyz, xyz_c)
            if this_xyz is not None:
                if region is not None:
                    if not xyz_in_region_p(this_xyz, region): pass_d = False
                if xyz_c['upper_limit'] is not None or xyz_c['lower_limit'] is not None:
                    if not z_pass(this_xyz[2], upper_limit = xyz_c['upper_limit'], lower_limit = xyz_c['lower_limit']): pass_d = False
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
    
    xcount, ycount, dst_gt = gdal_region2gt(region, inc)
    #xcount += 1
    #ycount += 1
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
                xpos, ypos = _geo2pixel(x, y, dst_gt)
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
    '''scan an xyz file and find it's min/max values and
    write an associated inf file for the src_xyz file.

    returns region [xmin, xmax, ymin, ymax] of the src_xyz file.'''
    
    minmax = []
    for i,l in enumerate(xyz_parse(src_xyz)):
        if i == 0:
            minmax = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            try:
                if l[0] < minmax[0]: minmax[0] = l[0]
                elif l[0] > minmax[1]: minmax[1] = l[0]
                if l[1] < minmax[2]: minmax[2] = l[1]
                elif l[1] > minmax[3]: minmax[3] = l[1]
                if l[2] < minmax[4]: minmax[4] = l[2]
                elif l[2] > minmax[5]: minmax[5] = l[2]
            except: pass
    if len(minmax) == 6:
        with open('{}.inf'.format(src_xyz.name), 'w') as inf:
            #utils.echo_msg('generating inf file for {}'.format(src_xyz.name))
            inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
        return(minmax)
    else: return(0,0,0,0,0,0)

def xyz_inf_entry(entry):
    '''find the region of the xyz datalist entry
    
    returns the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry'''
    
    with open(entry[0]) as infile:
        try:
            minmax = mbsfun.mb_inf(infile)
        except: minmax = xyz_inf(infile)
    return(minmax)        

def xyz_yield_entry(entry, region = None, verbose = False, z_region = None):
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
        
### End
