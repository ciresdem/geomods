### datalists.py
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
## datalists and entries - datalists.py
##
## a datalists is a space-delineated file containing:
## /path/to/data data-fmt data-weight data,meta,data
##
## datalist processing (datalists fmt:-1)
## entry processing fmt:*
## TODO -> pass_h to list for all entry passing
##
### Code:

import os
import sys
import glob
import json

from geomods import utils
from geomods import regions
from geomods import gdalfun
from geomods import mbsfun
from geomods import gmtfun
from geomods import xyzfun
from geomods import lasfun
from geomods import fetches

## ==============================================
## Datalist formats and lambdas
## ==============================================
_known_dl_delims = [' ']
_dl_dl_h = {
    -1: {
        'fmts': ['datalist', 'mb-1'],
        'inf': lambda e, p: datalist_inf_entry(e),
        'yield': lambda e, r, v, z, w, p: datalist_yield_entry(e, r, v, z, w, p),
    },
    168: {
        'fmts': ['xyz', 'csv', 'dat', 'ascii'],
        'inf': lambda e, p: xyzfun.xyz_inf_entry(e),
        'yield': lambda e, r, v, z, w, p: xyzfun.xyz_yield_entry(e, region = r, verbose = v, z_region = z),
    },
    200: {
        'fmts': ['tif', 'img', 'grd', 'nc', 'vrt', 'bag'],
        'inf': lambda e, p: gdalfun.gdal_inf_entry(e, p),
        'yield': lambda e, r, v, z, w, p: gdalfun.gdal_yield_entry(e, region = r, verbose = v, z_region = z, epsg = p),
    },
    300: {
        'fmts': ['las', 'laz'],
        'inf': lambda e, p: lasfun.las_inf_entry(e),
        'yield': lambda e, r, v, z, w, p: lasfun.las_yield_entry(e, region = r, verbose = v, z_region = z),
    },
    400: {
        'fmts': ['nos', 'dc', 'gmrt', 'srtm_cgiar', 'srtm_plus', 'mar_grav', 'charts', 'mb', 'tnm', 'emodnet', 'chs', 'hrdem', 'cudem'],
        'inf': lambda e, p: fetches.fetch_inf_entry(e, p),
        'yield': lambda e, r, v, z, w, p: fetches.fetch_yield_entry(e, region = r, warp = p, verbose = v),
    },
}

_datalist_fmts_short_desc = lambda: '\n  '.join(['{}\t{}'.format(key, _dl_dl_h[key]['fmts']) for key in _dl_dl_h])

def path_exists_or_url(src_str):
    '''return True if src_str is a valid path, url or fetches module'''
    
    if os.path.exists(src_str): return(True)
    if src_str[:4] == 'http': return(True)
    if src_str.split(':')[0] in _dl_dl_h[400]['fmts']: return(True)
    utils.echo_error_msg('invalid datafile/datalist: {}'.format(src_str))
    return(False)

def intersect_p(r, e, p = None):
    '''return True if region r intersects with the wkt found in the 
    inf file for datalist entry e'''
    
    dl_i = inf_entry(e, epsg = p)
    r_geom = gdalfun.gdal_region2geom(r)
    e_geom = gdalfun.gdal_wkt2geom(dl_i['wkt']) if 'wkt' in dl_i.keys() else r_geom
    return(regions.geoms_intersect_p(r_geom, e_geom))

def datalist_default_hooks():
    return([lambda e: path_exists_or_url(e[0])])

_dl_pass_h = [lambda e: path_exists_or_url(e[0])]

## ==============================================
## inf files (data info) inf.py
## mbsystem/gmt/waffles infos
## ==============================================
def inf_generate(data_path, data_fmt = 168):
    '''generate an info (.inf) file from the data_path'''
    
    return(inf_entry([data_path, data_fmt]), True)

def inf_parse(src_inf):
    '''parse an inf file (mbsystem or gmt)
    
    returns region: [xmin, xmax, ymin, ymax, zmin, zmax]'''

    xyzi = {'name': src_inf, 'numpts': 0, 'minmax': [0,0,0,0,0,0], 'wkt': gdalfun.gdal_region2wkt([0,0,0,0,0,0])}
    try:
        with open(src_inf) as iob:
            xyzi = json.load(iob)
    except:
        try:
            xyzi = gmtfun.gmt_inf_parse(src_inf)
        except:
            try:
                xyzi = mbsfun.mb_inf_parse(src_inf)
            except: pass

    return(xyzi)

def inf_entry(src_entry, overwrite = False, epsg = None):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file or the 
    result of `gmt gmtinfo file.xyz -C`, which is a 6 column line 
    with minmax info, etc.

    returns the region of the inf file.'''
    
    ei = {'name': src_entry, 'numpts': 0, 'minmax': [0,0,0,0,0,0], 'wkt': gdalfun.gdal_region2wkt([0,0,0,0,0,0])}
    if os.path.exists(src_entry[0]) or src_entry[1] == 400:
        path_i = src_entry[0] + '.inf'
        if not os.path.exists(path_i) or overwrite:
            ei = _dl_dl_h[src_entry[1]]['inf'](src_entry, epsg)
        else: ei = inf_parse(path_i)
        if not regions.region_valid_p(ei['minmax']): return({})
    return(ei)
    
## ==============================================
## datalists
## ==============================================
def datalist_inf(dl, inf_file = True, epsg = None, overwrite = False):
    '''return the region of the datalist and generate
    an associated `.inf` file if `inf_file` is True.'''
    
    out_regions = []
    dl_i = {'name': dl, 'minmax': None, 'numpts': 0, 'wkt': None}    
    utils.echo_msg('generating inf for datalist {}'.format(dl))
    
    for entry in datalist(dl, pass_h = _dl_pass_h):
        #if entry[1] == 200:
        entry_inf = inf_entry(entry, epsg = epsg, overwrite = overwrite)
        #else: entry_inf = inf_entry(entry, overwrite = overwrite)
        if entry_inf is not None:
            out_regions.append(entry_inf['minmax'][:6])
            dl_i['numpts'] += entry_inf['numpts']
    
    out_regions = [x for x in out_regions if x is not None]
    if len(out_regions) == 0:
        dl_i['minmax'] = None
    elif len(out_regions) == 1:
        dl_i['minmax'] = out_regions[0]
    else:
        out_region = out_regions[0]
        for x in out_regions[1:]:
            out_region = regions.regions_merge(out_region, x)
        dl_i['minmax'] = out_region

    dl_i['wkt'] = gdalfun.gdal_region2wkt(dl_i['minmax'])
    
    if dl_i['minmax'] is not None and inf_file:
        with open('{}.inf'.format(dl), 'w') as inf:
            inf.write(json.dumps(dl_i))
    return(dl_i)

def datalist_inf_entry(e, inf_file = True, epsg = None, overwrite = False):
    '''write an inf file for datalist entry e
    
    return the region [xmin, xmax, ymin, ymax]'''
    if e[0] is not None:
        return(datalist_inf(e[0], epsg = epsg, overwrite = overwrite))

def datalist_append_entry(entry, datalist):
    '''append entry to datalist file `datalist`'''
    
    with open(datalist, 'a') as outfile:
        outfile.write('{}\n'.format(' '.join([str(x) for x in entry])))
        
def datalist_echo(entry):
    '''echo datalist entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))
    datalist(entry[0])

def datafile_echo(entry):
    '''echo datafile entry to stderr'''
    
    sys.stderr.write('{}\n'.format([str(x) for x in entry]))

def datalist_major(dls, major = '.mjr.datalist', region = None):
    '''set the major datalist
    `dls` is a list of datalist entries, minimally: ['datafile.xyz']

    returns the major datalist filename'''
    
    with open(major, 'w') as md:        
        for dl in dls:
            entries = datalist2py(dl, region)
            for entry in entries:
                md.write('{}\n'.format(' '.join([str(e) for e in entry])))
    if os.stat(major).st_size == 0:
        utils.remove_glob(major)
        utils.echo_error_msg('bad datalist/entry, {}'.format(dls))
        return(None)    
    return(major)

def entry2py(dl_e):
    '''convert a datalist entry to python

    return the entry as a list [fn, fmt, wt, ...]'''
    
    this_entry = dl_e.rstrip().split()
    try:
        entry = [x if n == 0 else int(x) if n < 2 else float(x) if n < 3 else x for n, x in enumerate(this_entry)]
    except Exception as e:
        utils.echo_error_msg('could not parse entry {}'.format(dl_e))
        return(None)
    if len(entry) < 2:
        for key in _dl_dl_h.keys():
            se = entry[0].split('.')
            if len(se) == 1:
                see = entry[0].split(':')[0]
            else: see = se[-1]
            if see in _dl_dl_h[key]['fmts']:
                entry.append(int(key))
    if len(entry) < 3: entry.append(1)
    return(entry)

def datalist2py(dl, region = None):
    '''convert a datalist to python data
    
    returns a list of datalist entries.'''
    
    these_entries = []
    this_entry = entry2py(dl)
    this_dir = os.path.dirname(this_entry[0])
    if this_entry[1] == -1:
        if os.path.exists(this_entry[0]):
            with open(this_entry[0], 'r') as op:
                for this_line in op:
                    if this_line[0] != '#' and this_line[0] != '\n' and this_line[0].rstrip() != '':
                        these_entries.append([os.path.join(this_dir, x) if n == 0 else x for n,x in enumerate(entry2py(this_line.rstrip()))])
        else: utils.echo_error_msg('could not open datalist/entry {}'.format(this_entry[0]))
                
    else: these_entries.append(this_entry)
    return(these_entries)
            
def datalist(dl, wt = None, pass_h = _dl_pass_h,
             yield_dl_entry = False, verbose = False):
    '''recurse a datalist/entry
    for entry in datalist(dl): do_something_with entry

    yields entry [path, fmt, wt, ...]'''

    these_entries = datalist2py(dl)
    if len(these_entries) == 0: these_entries = [entry2py(dl)]
    for this_entry in these_entries:
        if this_entry is not None:
            this_entry[2] = wt if wt is None else wt * this_entry[2] 
            this_entry = this_entry[:3] + [' '.join(this_entry[3:]).split(',')] + [os.path.basename(dl).split('.')[0]]
            if not False in [x(this_entry) for x in pass_h]:
                if verbose and this_entry[1] == -1: utils.echo_msg('parsing datalist ({}) {}'.format(this_entry[2], this_entry[0]))
                if this_entry[1] == -1:
                    if yield_dl_entry: yield(this_entry)
                    for entry in datalist(this_entry[0], wt = this_entry[2], pass_h = pass_h, yield_dl_entry = yield_dl_entry, verbose = verbose):
                        yield(entry)
                else: yield(this_entry)

def datalist_archive_yield_entry(entry, dirname = 'archive', region = None, inc = 1, weight = None, verbose = None, z_region = None, epsg = None):
    '''archive a datalist entry.
    a datalist entry is [path, format, weight, ...]

    yield the xyz line data'''
    
    if region is None:
        a_name = entry[-1]
    else: a_name = '{}_{}_{}'.format(entry[-1], regions.region_format(region, 'fn'), this_year())
    
    i_dir = os.path.dirname(entry[0])
    i_xyz = os.path.basename(entry[0]).split('.')[0]
    i_xyz = ''.join(x for x in i_xyz if x.isalnum())
    a_dir = os.path.join(dirname, a_name, 'data', entry[-1])
    a_xyz_dir = os.path.join(a_dir, 'xyz')
    a_xyz = os.path.join(a_xyz_dir, i_xyz + '.xyz')
    a_dl = os.path.join(a_xyz_dir, '{}.datalist'.format(entry[-1]))
    
    if not os.path.exists(a_dir): os.makedirs(a_dir)
    if not os.path.exists(a_xyz_dir): os.makedirs(a_xyz_dir)

    with open(a_xyz, 'w') as fob:
        for xyz in datalist_yield_entry(entry, region = region, verbose = verbose, z_region = z_region, epsg = epsg):
            xyzfun.xyz_line(xyz, fob)
            yield(xyz)
            
    inf = mbsfun.mb_inf(a_xyz)
    datalist_append_entry([i_xyz + '.xyz', 168, entry[2] if entry[2] is not None else 1], a_dl)
    
def datalist_yield_entry(this_entry, region = None, verbose = False, z_region = None, w_region = None, epsg = None):
    '''yield the xyz data from the datalist entry [entry/path, entry-format, entry-weight]

    yields xyz line data [x, y, z, ...]'''

    if this_entry[1] != -1:
        for xyz in _dl_dl_h[this_entry[1]]['yield'](this_entry, region, verbose, z_region, w_region, epsg):
            yield(xyz)
        
def datalist_yield_xyz(dl, wt = None, pass_h = _dl_pass_h,
                       region = None, archive = False,
                       mask = False, verbose = False, z_region = None, epsg = None):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''
    
    for this_entry in datalist(dl, wt = wt, pass_h = pass_h, verbose = verbose):
        dly = datalist_yield_entry(this_entry, region, verbose = verbose, z_region = z_region, epsg = epsg)
        if archive: dly = datalist_archive_yield_entry(this_entry, dirname = 'archive', region = region, weight = wt, verbose = verbose, z_region = z_region, epsg = None)
        for xyz in dly:
            yield(xyz)

def datalist_dump_xyz(dl, wt = None,  pass_h = _dl_pass_h,
                      region = None, archive = False, mask = False,
                      verbose = False, dst_port = sys.stdout, z_region = None, epsg = None):
    '''parse out the xyz data from the datalist
    for xyz in datalist_yield_xyz(dl): xyz_line(xyz)

    yields xyz line data [x, y, z, ...]'''

    for xyz in datalist_yield_xyz(dl, wt, pass_h, region, archive, mask, verbose, z_region, epsg):
        xyzfun.xyz_line(xyz, dst_port, verbose)

## ==============================================
## dadtalists cli
## ==============================================    
datalists_version = '0.0.1'
datalists_usage = '''{} ({}): Process and generate datalists

usage: {} [ -aghirvFPR [ args ] ] DATALIST ...

Options:
  -R, --region\t\tSpecifies the desired REGION;
  -P, --epsg\t\tSpecify the EPSG of the DATALIST.
  -F, --format\t\tOnly process FORMAT data type.

  -a, --archive\t\tARCHIVE the data from the DATALIST.
  -g, --glob\t\tGLOB FORMAT data into the DATALIST.
  -i, --info-file\tGenerate INF files for the data in the DATALIST
  -l, --list\t\tLIST the datafiles from the DATALIST.
  -m, --mask\t\tMASK the datafiles from the DATALIST.
  -r, --region-info\tReturn the full REGION of the DATALIST.
  -w, --weights\t\toutput weights along with each datalist.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

Supported datalist formats: 
  {}

 Examples:
 % {} my_data.datalist -R -90/-89/30/31 -g -i

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            datalists_version, 
            os.path.basename(sys.argv[0]),
            _datalist_fmts_short_desc(),
            os.path.basename(sys.argv[0]))

def datalists_cli(argv = sys.argv):

    status = 0
    dls = []
    i_region = None
    i_inc = 0.000277777
    o_bn = None
    epsg = None
    these_regions = []
    want_verbose = False
    want_inf = False
    want_region = False
    want_sm = False
    want_glob = False
    want_list = False
    want_archive = False
    want_mask = False
    want_weights = False
    z_region = None
    w_region = None
    dl_fmt = None
    #dl_fmts = []
    
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================

    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            i_region = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            i_region = str(arg[2:])

        elif arg == '--output-name' or arg == '-O':
            o_bn = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-O':
            o_bn = str(arg[2:])

        elif arg == '--increment' or arg == '-E':
            try:
                i_inc = float(argv[i + 1])
            except:
                sys.stederr.write('error, -E should be a float value\n')
                sys.exit(1)
            i = i + 1
        elif arg[:2] == '-E':
            try:
                i_inc = float(arg[2:])
            except:
                sys.stederr.write('error, -E should be a float value\n')
                sys.exit(1)

        elif arg == '--epsg' or arg == '-P':
            epsg = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-P':
            epsg = arg[2:]

        elif arg == '--format' or arg == '-F':
            dl_fmt = int(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-F':
            dl_fmt = int(arg[2:])

        elif arg == '--glob' or arg == '-g':
            want_glob = True

        elif arg == '--spatial-md' or arg == '-s':
            want_sm = True

        elif arg == '--list' or arg == '-l':
            want_list = True

        elif arg == '--mask' or arg == '-m':
            want_mask = True

        elif arg == '--archive' or arg == '-a':
            want_archive = True

        elif arg == '--info-file' or arg == '-i':
            want_inf = True

        elif arg == '--region-info' or arg == '-r':
            want_region = True

        elif arg == '--weights' or arg == '-w':
            want_weights = True

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}'.format(os.path.basename(sys.argv[0]), datalists_version))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        else:
            dls.append(arg)

        i = i + 1

    if want_glob:
        if dl_fmt is None:
            dl_fmts = list(_dl_dl_h.keys())
            dl_fmts.remove(-1)
        else:
            dl_fmts = [dl_fmt]
        # for f in _known_datalist_fmts[-1]:
        #     globs = glob.glob('*.{}'.format(f))
        #     [sys.stdout.write('{}\n'.format(' '.join([x, str(-1), '1']))) for x in globs]
        for key in dl_fmts:
            for f in _dl_dl_h[key]['fmts']:
                globs = glob.glob('*.{}'.format(f))
                [sys.stdout.write('{}\n'.format(' '.join([x, str(key), '1']))) for x in globs]
        sys.exit(0)
        
    if len(dls) == 0:
        print(datalists_usage)
        sys.exit(1)
        
    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if i_region is not None:
        try:
            these_regions = [[float(x) for x in i_region.split('/')]]
        except ValueError: these_regions = gdalfun.gdal_ogr_regions(i_region)
        except Exception as e:
            utils.echo_error_msg('failed to parse region(s), {}'.format(e))
        #else: these_regions = [[-180,180,-90,90]]
    else: these_regions = [None]
    if len(these_regions) == 0: utils.echo_error_msg('failed to parse region(s), {}'.format(i_region))

    for rn, this_region in enumerate(these_regions):
        
        ## ==============================================
        ## Load the input datalist
        ## ==============================================
        dl_m = datalist_major(dls, region = this_region)
        utils.echo_msg('processed datalist')
        dlp_hooks = datalist_default_hooks()
        if regions.region_valid_p(this_region):
            dlp_hooks.append(lambda e: intersect_p(this_region, e, epsg))
        if z_region is not None:
            dlp_hooks.append(lambda e: regions.z_region_pass(inf_entry(e)['minmax'], upper_limit = z_region[1], lower_limit = z_region[0]))
        if w_region is not None:
            dlp_hooks.append(lambda e: regions.z_pass(e[2], upper_limit = w_region[1], lower_limit = w_region[0]))
        
        if want_inf:
            datalist_inf_entry([dl_m, -1, 1], epsg = epsg, overwrite = True)
        elif want_region:
            print(datalist_inf(dl_m, inf_file = False, epsg = epsg, overwrite = False))
        elif want_list:
            for this_entry in datalist(dl_m, wt = 1, pass_h = dlp_hooks):
                print(' '.join([','.join(x) if i == 3 else os.path.abspath(str(x)) if i == 0 else str(x) for i,x in enumerate(this_entry[:-1])]))
        elif want_mask: pass
        else:
            datalist_dump_xyz(dl_m, wt = 1 if want_weights else None, pass_h = dlp_hooks, region = this_region, epsg = epsg)
        utils.remove_glob(dl_m)
        
### End
