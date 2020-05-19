### datalists.py
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
### Commentary:
##
## DATALIST and REGION functions - datalists.py
##
## MBSystem/Waffles style datalists.
## Recurse through a datalist file and process the results.
##
## a datalist '*.datalist' file should be formatted as in MBSystem:
## ~path ~format ~weight
##
## with the additional columns of ~"meta,data" ~etc
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
## a format of 200 represents GDAL data
## a format of 300 represents LAS/LAZ data
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
## or via _datalist_inf_hooks
##
## if 'region' is specified, will only process data that falls within
## the given region
##
## Adjust the _datalist*hooks dictionary to add custom processing/testing
##
### Code:

import os
import sys
import gdal
import ogr
import numpy as np

import utils

_version = '0.3.0'

_known_dl_delims = [':', ' ']
_known_delims = [',', ' ', '\t', '/', ':']
_known_datalist_fmts = {-1: ['datalist', 'mb-1'], 168: ['xyz', 'dat'], 200: ['tif', 'img', 'grd', 'nc']}

## ==============================================
## datalist processing and testing hooks
## ==============================================

_datalist_pass_hooks = {-1: lambda e: os.path.exists(e[0]), 168: lambda e: os.path.exists(e[0]), 200: lambda e: os.path.exists(e[0])}
_datalist_inf_hooks = {-1: lambda e: datalist_inf(*e), 168: lambda e: xyz_inf_entry(e), 200: lambda e: gdal_inf_entry(e)}
_datalist_hooks = {-1: lambda e, v: datalist(*e, verbose = v), 168: lambda e, v: xyz_dump_entry(e, verbose = v), 200: lambda e, v: gdal_dump_entry(e, verbose = v)}

## ==============================================
## stderr messaging
## ==============================================

echo_msg = lambda m: utils.echo_msg(m, prefix = 'datalists')
echo_error_msg = lambda m: utils.echo_msg(m, prefix = 'datalists')

## ==============================================
## regions - regions are a bounding box list:
## [w, e, s, n]
## -- regions.py --
## ==============================================

def region_valid_p(region):
    '''return True if region appears valid'''
    if region[0] < region[1] and region[2] < region[3]: return(True)
    else: return(False)

def region_center(region):
    '''return the center point of the region'''
    xc = region[0] + (region[1] - region[0] / 2)
    yc = region[2] + (region[3] - region[2] / 2)
    return([xc, yc])

def region_pct(self, pctv):
    '''return the pctv val of the region'''
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv = 0):
    '''return the region buffered by bv'''
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def regions_reduce(region_a, region_b):
    '''return the minimum region when combining region_a and region_b'''
    region_c = [0, 0, 0, 0]
    if region_a[0] > region_b[0]: region_c[0] = region_a[0]
    else: region_c[0] = region_b[0]
    if region_a[1] < region_b[1]: region_c[1] = region_a[1]
    else: region_c[1] = region_b[1]
    if region_a.south > region_b[2]: region_c[2] = region_a[2]
    else: region_c[2] = region_b[2]
    if region_a[3] < region_b[2]: region_c[3] = region_a[3]
    else: region_c[3] = region_b[3]
    return(region_c)

def regions_merge(region_a, region_b):
    '''merge two regions into a single region'''
    region_c = [0, 0, 0, 0]
    if region_a[0] < region_b[0]: region_c[0] = region_a[0]
    else: region_c[0] = region_b[0]
    if region_a[1] > region_b[1]: region_c[1] = region_a[1]
    else: region_c[1] = region_b[1]
    if region_a[2] < region_b[2]: region_c[2] = region_a[2]
    else: region_c[2] = region_b[2]
    if region_a[3] > region_b[3]: region_c[3] = region_a[3]
    else: region_c[3] = region_b[3]
    return(region_c)

def regions_intersect_p_depr(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''    
    return(region_valid_p(regions_reduce(region_a, region_b)))
    
def regions_intersect_p(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''        
    geom_a = gdal_region2geom(region_a)
    geom_b = gdal_region2geom(region_b)
    if geom_a.Intersects(geom_b):
        return(True)
    else: return(False)

def region_format(region, t = 'gmt'):
    '''format region to string'''
    if t == 'gmt': return('-R' + '/'.join([str(x) for x in region]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'fn':
        if region[3] < 0: ns = 's'
        else: ns = 'n'
        if region[0] > 0: ew = 'e'
        else: ew = 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100) % 100), 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100) % 100)))

## ==============================================
## inf files (data info) inf.py
## mbsystem/gmt/waffles infos
## ==============================================

def inf_generate(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    return(inf_entry([data_path, data_fmt]))

def inf_parse(src_inf):
    '''parse an inf file (mbsystem or gmt) and return minmax'''
    minmax = [0, 0, 0, 0]
    with open(src_inf) as iob:
        for il in iob:
            til = il.split()
            if len(til) > 1:
                try: 
                    minmax = [float(x) for x in til]
                except:
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            minmax[0] = til[2]
                            minmax[1] = til[5]
                        elif til[1] == 'Latitude:':
                            minmax[2] = til[2]
                            minmax[3] = til[5]
    return([float(x) for x in minmax])

def inf_entry(src_entry):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file
    or the result of `gmt gmtinfo file.xyz -C`, which is
    a 6 column line with minmax info, etc.
    returns the region of the inf file.'''
    minmax = None
    if os.path.exists(src_entry[0]):
        path_i = src_entry[0] + '.inf'
        if not os.path.exists(path_i):
            _datalist_inf_hooks[src_entry[1]](src_entry)
        if os.path.exists(path_i): minmax = inf_parse(path_i)[:4]
    if region_valid_p(minmax):
        return(minmax)
    else: return(None)
    
## ==============================================
## xyz processing (fmt:168)
## ==============================================

def xyz_parse(src_xyz):
    '''xyz file parsing generator'''
    for xyz in src_xyz:
        this_line = xyz.strip()
        for delim in _known_delims:
            this_xyz = this_line.split(delim)
            if len(this_xyz) > 1:
                this_delim = delim
                break
        yield [float(x) for x in this_xyz]

def xyz_line(line, dst_port = sys.stdout, delimiter = ' ', weight = None):
    '''write XYZ `line` to `dst_port` using `delimiter` and `weight`'''
    if weight is None:
        w_string = ''
    else: w_string = '{}{}'.format(delimiter, weight)
    l = '{}{}\n'.format(delimiter.join([str(x) for x in line]), w_string)#.encode(dst_port.encoding)
    if dst_port != sys.stdout: l = l.encode('utf-8')
    dst_port.write(l)

def xyz_inf(src_xyz):
    '''return minmax info from a src_xyz file.'''
    minmax = []
    for i,l in enumerate(xyz_parse(src_xyz)):
        if i == 0:
            minmax = [l[0], l[0], l[1], l[1], l[2], l[2]]
        else:
            if l[0] < minmax[0]: minmax[0] = l[0]
            elif l[0] > minmax[1]: minmax[1] = l[0]
            if l[1] < minmax[2]: minmax[2] = l[1]
            elif l[1] > minmax[3]: minmax[3] = l[1]
            if l[1] < minmax[2]: minmax[4] = l[2]
            elif l[1] > minmax[3]: minmax[5] = l[2]
    with open('{}.inf'.format(src_xyz.name), 'w') as inf:
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    return(minmax)

def xyz_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n], else False.'''
    if src_xy[0] < src_region[0]: return(False)
    elif src_xy[0] > src_region[1]: return(False)
    elif src_xy[1] < src_region[2]: return(False)
    elif src_xy[1] > src_region[3]: return(False)
    else: return(True)

def xyz_inf_entry(entry):
    with open(entry[0]) as infile:
        minmax = xyz_inf(infile)
            
def xyz_dump_entry(entry, dst_port = sys.stdout, region = None, delimiter = ' ', verbose = False):
    '''dump ascii xyz data to dst_port'''
    with open(entry[0]) as infile:
        for line in xyz_parse(infile):
            if region is not None:
                if xyz_in_region_p(line, region):
                    xyz_line(line, dst_port, delimiter, entry[2])
            else: xyz_line(line, dst_port, delimiter, entry[2])

## ==============================================
## gdal processing (fmt:200)
## ==============================================

def gdal_infos(src_gdal):
    if os.path.exists(src_gdal):
        ds = gdal.Open(src_gdal)
        dsc = None
        if ds is not None:
            dsc = gdal_gather_infos(ds)
        ds = None
        return(dsc)
    else: return(None)

def gdal_gather_infos(src_ds):
    '''gather information from `src_ds` GDAL dataset'''
    ds_config = {
        'nx': src_ds.RasterXSize,
        'ny': src_ds.RasterYSize,
        'nb':src_ds.RasterCount,
        'geoT': src_ds.GetGeoTransform(),
        'proj': src_ds.GetProjectionRef(),
        'dt': src_ds.GetRasterBand(1).DataType,
        'dtn': gdal.GetDataTypeName(src_ds.GetRasterBand(1).DataType),
        'ndv': src_ds.GetRasterBand(1).GetNoDataValue(),
        'fmt': src_ds.GetDriver().ShortName,
    }
    if ds_config['ndv'] is None: ds_config['ndv'] = -9999
    return(ds_config)

def gdal_set_infos(nx, ny, nb, geoT, proj, dt, ndv, fmt):
    '''Set a datasource config dictionary'''
    return({'nx': nx, 'ny': ny, 'nb': nb, 'geoT': geoT, 'proj': proj, 'dt': dt, 'ndv': ndv, 'fmt': fmt})

def gdal_gt2region(ds_config):
    '''convert a gdal geo-tranform to an extent [w, e, s, n]'''
    geoT = ds_config['geoT']
    return([geoT[0], geoT[0] + geoT[1] * ds_config['nx'], geoT[3] + geoT[5] * ds_config['ny'], geoT[3]])

def gdal_create_polygon(coords):
    '''convert coords to Wkt'''
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[1], coord[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def gdal_region2geom(region):
    '''convert an extent [west, east, south, north] to an OGR geometry'''
    eg = [[region[2], region[0]], [region[2], region[1]],
          [region[3], region[1]], [region[3], region[0]],
          [region[2], region[0]]]
    geom = ogr.CreateGeometryFromWkt(gdal_create_polygon(eg))
    return(geom)

def gdal_region(src_ds):
    '''return the extent of the src_fn gdal file.'''
    ds_config = gdal_gather_infos(src_ds)
    return(gdal_gt2region(ds_config))

def gdal_inf(src_ds):
    '''generate an info (.inf) file from a src_gdal file using gdal'''
    minmax = [str(x) for x in gdal_region(src_ds)]
    with open('{}.inf'.format(src_ds.GetDescription()), 'w') as inf:
        inf.write('{}\n'.format(' '.join(minmax)))
    return(minmax)

def _geo2pixel(geo_x, geo_y, geoTransform):
    '''convert a geographic x,y value to a pixel location of geoTransform'''
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = (geo_x - geoTransform[0]) / geoTransform[1]
        pixel_y = (geo_y - geoTransform[3]) / geoTransform[5]
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt( geoTransform))
    return(int(round(pixel_x)), int(round(pixel_y)))

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    '''convert a pixel location to geographic coordinates given geoTransform'''
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    '''apply geotransform to in_x,in_y'''
    out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
    out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]
    return(out_x, out_y)

def _invert_gt(geoTransform):
    '''invert the geotransform'''
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

def gdal_srcwin(src_ds, region):
    '''given a gdal file src_fn and a region [w, e, s, n],
    output the appropriate gdal srcwin.'''    
    ds_config = gdal_gather_infos(src_ds)
    this_origin = _geo2pixel(region[0], region[3], ds_config['geoT'])
    this_end = _geo2pixel(region[1], region[2], ds_config['geoT'])
    this_size = (int(this_end[0] - this_origin[0]), int(this_end[1] - this_origin[1]))
    this_origin = [0 if x < 0 else x for x in this_origin]
    this_size = [0 if x < 0 else x for x in this_size]
    if this_size[0] > ds_config['nx'] - this_origin[0]: this_size[0] = ds_config['nx'] - this_origin[0]
    if this_size[1] > ds_config['ny'] - this_origin[1]: this_size[1] = ds_config['ny'] - this_origin[1]
    return(this_origin[0], this_origin[1], this_size[0], this_size[1])

def gdal_parse(src_ds, dst_xyz = sys.stdout, delim = ' ', weight = None, dump_nodata = False, srcwin = None):
    '''send the data from gdal file src_gdal to dst_xyz port'''
    if weight is None:
        w_string = ''
    else: w_string = '{}{}'.format(delim, weight)

    band = src_ds.GetRasterBand(1)
    ds_config = gdal_gather_infos(src_ds)
    gt = ds_config['geoT']
        
    if srcwin is None: srcwin = (0, 0, ds_config['nx'], ds_config['ny'])

    for y in range(srcwin[1], srcwin[1] + srcwin[3], 1):
        nodata = ['{:.10f}',format(-9999), 'nan']
        data = []

        if band.GetNoDataValue() is not None: nodata.append('{:.10f}'.format(band.GetNoDataValue()))
        
        band_data = band.ReadAsArray(srcwin[0], y, srcwin[2], 1)
        band_data = np.reshape(band_data, (srcwin[2], ))
        data.append(band_data)

        for x_i in range(0, srcwin[2], 1):
            x = x_i + srcwin[0]
            geo_x, geo_y = _pixel2geo(x, y, gt)
            x_i_data = [data[0][x_i]]
            z = x_i_data[0]
            line = [geo_x, geo_y, z, weight]
            
            if dump_nodata:
                xyz_line(line, dst_xyz)
            else:
                if '{:.10f}'.format(z) not in nodata:
                    xyz_line(line, dst_xyz)

def gdal_inf_entry(entry):
    ds = gdal.Open(entry[0])
    minmax = gdal_inf(ds)
    ds = None
                    
def gdal_dump_entry(entry, dst_port = sys.stdout, region = None, delimiter = ' ', verbose = False):
    ds = gdal.Open(entry[0])
    if region is not None:
        srcwin = gdal_srcwin(ds, region)
    else: srcwin = None
    gdal_parse(ds, dst_xyz = dst_port, delim = delimiter, weight = entry[2], dump_nodata = False, srcwin = srcwin)
    ds = None
            
## ==============================================
## datalist processing (fmt:-1)
## ==============================================

def datalist_inf(dl, fmt = -1, wt = 1):
    '''return the region of the datalist.'''
    out_regions = []
    _datalist_hooks[168] = lambda e, v: out_regions.append(inf_entry(e))
    _datalist_hooks[200] = lambda e, v: out_regions.append(inf_entry(e))
    datalist(dl)
    out_regions = [x for x in out_regions if x is not None]
    out_region = out_regions[0]
    for i in out_regions[1:]:
        out_region = regions_merge(out_region, i)
    with open('{}.inf'.format(dl), 'w') as inf:
        inf.write('{}\n'.format(' '.join([str(x) for x in out_region])))
    return(out_region)

def datalist_dump(dl, dst_port = sys.stdout, region = None, verbose = False):
    '''dump the data from datalist to dst_port'''
    _datalist_hooks[168] = lambda e, v: xyz_dump_entry(e, region = region, dst_port = dst_port, verbose = v)
    _datalist_hooks[200] = lambda e, v: gdal_dump_entry(e, region = region, dst_port = dst_port, verbose = v)
    datalist(dl, verbose = verbose)
    
def entry2py(dle):
    '''convert a datalist entry to python'''
    for delim in _known_dl_delims:
        this_entry = dle.rstrip().split(delim)
        if len(this_entry) > 1: break
    try:
        entry = [x if n == 0 else float(x) if n < 3 else x for n, x in enumerate(this_entry)]
    except ValueError as e: return(None)
    return(entry)
    
def datalist2py(dl):
    '''convert a datalist to python data'''
    these_entries = []
    this_entry = entry2py(dl)
    if len(this_entry) < 2: this_entry.append(-1)
    if len(this_entry) < 3: this_entry.append(1)
    if this_entry[1] == -1:
        with open(this_entry[0], 'r') as op:
            while True:
                this_line = op.readline().rstrip()
                if not this_line: break
                if this_line[0] != '#' and this_line[0] != '\n':
                    these_entries.append(entry2py(this_line))
    else: these_entries.append(this_entry)
    return(these_entries)

def datalist(dl, fmt = -1, wt = 1, verbose = False):
    '''recurse a datalist/entry'''
    this_dir = os.path.dirname(dl)
    these_entries = datalist2py(dl)
    if len(these_entries) == 0: these_entries = [entry2py(dl)]
    for this_entry in these_entries:
        this_entry[0] = os.path.join(this_dir, this_entry[0])
        this_entry[2] = wt * this_entry[2]
        if _datalist_pass_hooks[this_entry[1]](this_entry):
            if verbose is True: echo_msg('using datafile {}'.format(this_entry[0]))
            _datalist_hooks[this_entry[1]](this_entry, verbose)
                
## ==============================================
## mainline - datalists
## ==============================================

datalists_cli_usage = '''datalists [-wR] <datalist/entry>

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''

def datalists_cli():
    argv = sys.argv
    dls = []
    region = None
    want_weight = False
    want_verbose = False
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '--region' or arg == '-R':
            region = str(argv[i + 1])
            i += 1
        elif arg[:2] == '-R': region = str(arg[2:])
        elif arg == 'w': want_weight = True
        elif arg == '--verbose' or arg == '-V': want_verbose = True
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(datalists_cli_usage)
            sys.exit(0)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}'.format(_version))
            sys.exit(0)
        else: dls.append(arg)
        i += 1
    if len(dls) == 0:
        sys.stderr.write(datalists_cli_usage)
        echo_error_msg('''must specify a datalist/entry file''')
        sys.exit(-1)
    if region is not None:
        region = [float(x) for x in region.split('/')]
        _datalist_pass_hooks[168] = lambda e: regions_intersect_p(region, inf_entry(e))
        _datalist_pass_hooks[200] = lambda e: regions_intersect_p(region, inf_entry(e))
        _datalist_hooks[168] = lambda e, v: xyz_dump_entry(e, region = region, verbose = want_verbose)
        _datalist_hooks[200] = lambda e, v: gdal_dump_entry(e, region = region, verbose = want_verbose)
    ## ==============================================
    ## recurse the datalist
    ## ==============================================
    master = '.master.datalist'
    with open(master, 'w') as md:
        for dl in dls:
            if len(dl.split(':')) == 1:
                for key in _known_datalist_fmts.keys():
                    if dl.split(':')[0].split('.')[-1] in _known_datalist_fmts[key]:
                        md.write('{} {} 1\n'.format(dl, key))
    try:
        datalist(master, verbose = want_verbose)
    except KeyboardInterrupt as e:
        echo_error_msg('user killed process')
        utils.remove_glob(master)
        sys.exit(-1)
    utils.remove_glob(master)

if __name__ == '__main__': datalists_cli()
### End
