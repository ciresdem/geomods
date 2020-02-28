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
### Code:

import sys
import os
import glob
import ogr
import regions
import gdalfun
import utils

_version = '0.1.2'

## =============================================================================
##
## Datalist Class
##
## MBSystem style datalists.
## Recurse through a datalist file and process the results.
##
## a datalist '*.datalits' file should be formatted as in MBSystem:
## ~path ~format ~weight
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
##
## if 'i_region' is specified, will only process data that falls within
## the given region
##
## =============================================================================

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

def xyz_region(src_xyz):
    out, status = utils.run_cmd('gmt gmtinfo {} -I-'.format(src_xyz), False, True)
    o_region = regions.region(out[2:])
    return(o_region)

def xyz_inf(src_xyz):
    minmax = []
    with open(src_xyz, 'r') as in_file:
        for i,l in enumerate(xyz_parse(in_file)):
            if i == 0:
                minmax[0] = l[0], minmax[1] = l[0]
                minmax[2] = l[1], minmax[3] = l[1]
                minmax[4] = l[2], minmax[5] = l[2]
            else:
                if l[0] < minmax[0]: minmax[0] = l[0]
                elif l[0] > minmax[1]: minmax[1] = l[0]
                if l[1] < minmax[2]: minmax[2] = l[1]
                elif l[1] > minmax[3]: minmax[3] = l[1]
                if l[1] < minmax[2]: minmax[4] = l[2]
                elif l[1] > minmax[3]: minmax[5] = l[2]
                
    return(minmax)

def xyz_inf_gmt(src_xyz):
    out, status = utils.run_cmd('gmt gmtinfo {} -C > {}.inf'.format(src_xyz, src_xyz), False, True)
    return(out, status)

def xyz_inf_mb(src_xyz):
    out, status = utils.run_cmd('mbdatalist -O -F168 -I{}'.format(src_xyz), False, True)
    return(out, status)

def xy_in_region_p(src_xy, src_region):
    '''return True if point [x, y] is inside region [w, e, s, n]'''

    x = src_xyz[0]
    y = src_xyz[1]

    if x < src_region[0]:
        return False
    elif x > src_region[1]:
        return False
    elif y < src_region[2]:
        return False
    elif y > src_region[3]:
        return False
    else: return True

def gdal_inf_gmt(src_gdal):
    out, status = utils.run_cmd('gmt grdinfo {} -C > {}.inf'.format(src_gdal, src_gdal), False, True)
    return(out, status)

def generate_inf(data_path, data_fmt = 168):    
    this_region = datafile_region(data_path, data_fmt)
    if this_region is None:
        if data_fmt == 168:
            xyz_inf_mb(data_path)
        elif data_fmt == 200:
            gdal_inf_gmt(data_path)
    
def datafile_region(data_path, data_fmt = 168):
    '''Read .inf file and extract minmax info.
    the .inf file can either be an MBSystem style inf file
    or the result of `gmt gmtinfo file.xyz -C`, which is
    a 6 column line with minmax info, etc.'''

    if data_fmt == 200:
        minmax = (map(str, gdalfun._extent(data_path)))
    elif data_fmt == 168:    
        path_i = data_path + '.inf'

        minmax = [0, 0, 0, 0]
        if not os.path.exists(path_i):
            xyz_inf_mb(data_path)

        if os.path.exists(path_i):
            iob = open(path_i)
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    ## GMT inf
                    try:
                        minmax[0] = float(til[0])
                        minmax[1] = float(til[1])
                        minmax[2] = float(til[2])
                        minmax[3] = float(til[3])
                    ## mbsystem inf
                    except:
                        if til[0] == 'Minimum':
                            if til[1] == 'Longitude:':
                                minmax[0] = til[2]
                                minmax[1] = til[5]
                            elif til[1] == 'Latitude:':
                                minmax[2] = til[2]
                                minmax[3] = til[5]

    try: 
        o_region = regions.region('/'.join(minmax))
    except: o_region = None

    return(o_region)

def datalist_set_weight(data_e, weight = None):
    if weight is not None:
        try:
            dweight = int(weight)
        except: dweight = 1
    else:
        if len(data_e) > 2:
            try:
                dweight = int(data_e[2])
            except: dweight = 1
        else: dweight = 1
    return dweight


def datalist_bounds(dl, layer, region, inc):
    '''gather geometries from datalist `dl` and append
    results to ogr `layer`. Load the datalist, generate
    a NUM-MSK grid, polygonize said NUM-MSK then union
    the polygon and add it to the output layer.
    `dl` is a datalist line [dl_path, dl_fmt, dl_weight, metadata..., ...]
    `layer` is an OGR layer object.'''

    status = 0
    defn = layer.GetLayerDefn()
    this_datalist = datalists.datalist(dl[0], region)
    this_datalist._load_data()

    try:
        o_v_fields = [dl[3], dl[4], dl[5], dl[6], dl[7], dl[8], dl[9], dl[10].strip()]
    except: o_v_fields = [this_datalist._name, 'Unknown', 0, 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

    pb = utils._progress('gathering geometries from \033[1m{}\033[m...'.format(this_datalist._path_basename))
    this_mask = this_datalist.mask(inc = inc)
    if not os.path.exists(this_mask):
        status = -1
    else:
        utils.remove_glob('{}_poly.*'.format(this_o_name))

        tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_datalist._name))
        tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_datalist._name), None, ogr.wkbMultiPolygon)
        tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))
        
        pb1 = utils._progress('polygonizing \033[1m{}\033[m...'.format(this_datalist._path_basename))
        gdalfun.polygonize(this_mask, tmp_layer, verbose = True)
        pb1.opm = 'polygonized \033[1m{}\033[m.'.format(this_datalist._path_basename)
        pb1.end(status)

        if len(tmp_layer) > 1:
            out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn)
            for i, f in enumerate(gdalfun._ogr_get_layer_fields(layer)):
                out_feat.SetField(f, o_v_fields[i])

            layer.CreateFeature(out_feat)
            
        tmp_ds = tmp_layer = out_feat = None
        utils.remove_glob('{}_poly.*'.format(this_datalist._name))
        os.remove(this_mask)

    pb.opm = 'gathered geometries from \033[1m{}\033[m.'.format(this_datalist._path_basename)
    pb.end(status)

## {fmt_id, [fmt_exts, ...]}
known_fmts = {168: ['xyz', 'dat'],
              200: ['tif', 'img']}
    
class datalist:
    '''MBSystem style datalists for elevation data.
    Supports XYZ and GDAL data.'''

    def __init__(self, i_datalist, i_region = None, i_fmt = None, verbose = False):
        if not os.path.exists(i_datalist): 
            open(i_datalist, 'a').close()

        self.region = i_region
        
        self._path = i_datalist
        self._path_dirname = os.path.dirname(self._path)
        self._path_basename = os.path.basename(self._path)
        self._name = os.path.basename(self._path).split('.')[0]
        self.i_fmt = i_fmt
        self.datalist = []
        self.datafiles = []
        self.verbose = verbose
        
        self._load()

    def _reset(self):
        '''reload the datalist'''

        self.datalist = []
        self.datafiles = []
        self._load()

    def _valid_p(self):
        '''validate the datalist'''

        ## an empty datalist should technically be valid
        ## also check if lines have at least 2 columns
        if len(self.datalist) > 0:
            return(True)
        else: return(False)

    def _valid_d(self, data_l):
        '''validate a dataline'''

        if data_l[0] == '#' or data_l[0] == '\n' or data_l[0] == '': return(False)
        dl_cols = [x.strip() for x in data_l.split(' ')]
        if len(dl_cols) < 2: return(False)
        path_d = os.path.join(self._path_dirname, dl_cols[0])
        if not os.path.exists(path_d): return(False)
        try:
            int(dl_cols[1])
        except: return(False)
        return(True)
        
    def _load(self):
        '''read a datalist and gather all the datalists found therein.'''

        with open(self._path, 'r') as fob:
            for dl in fob:
                if self._valid_d(dl):
                    dl_cols = [x.strip() for x in dl.split(' ')]
                    self.datalist.append(dl_cols)

    def _load_data(self):
        '''load a datalist and gather all datafiles found therein.'''

        if self.verbose:
            status = 0
            pb = utils._progress('loading datalist `\033[1m{}\033[m`...'.format(self._name))
            
        self._proc_data(lambda x, t, w: self.datafiles.append([x,t,w]))

        if self.verbose:
            if len(self.datafiles) == 0: status = -1
            pb.end(status, 'loaded datalist `\033[1m{}\033[m`.'.format(self._name))
        
    def _proc_data(self, proc = lambda x, t, w: None, weight = None):
        '''Recurse through the datalist and run proc on each data file.'''

        for data_e in self.datalist:
            dpath = data_e[0]
            try:
                dformat = int(data_e[1])
            except: dformat = 168
            dweight = datalist_set_weight(data_e, weight)
            
            ## ==============================================
            ## Datalist format -1
            ## Open as a datalist and recurse
            ## ==============================================

            if dformat == -1:
                d = datalist(dpath, self.region, self.i_fmt, self.verbose)
                if self.i_fmt == -1:
                    proc(dpath, dformat, dweight)
                d._proc_data(proc, dweight)

            ## ==============================================
            ## Data formats 0+
            ## check if passes tests and append to datafiles if so.
            ## ==============================================

            else:
                usable = True
                path_d = os.path.join(self._path_dirname, dpath)
                
                if not os.path.exists(path_d): usable = False
                if self.i_fmt is not None and self.i_fmt != dformat: usable = False 

                if usable:
                    #dinf_region = None
                    dinf_region = datafile_region(path_d, data_fmt = dformat)
                    # if dformat == 168:
                    #     dinf_region = datafile_region(path_d)
                    # elif dformat == 200:
                    #     dinf_region = datafile_region(path_d, data_fmt = 200)

                    if self.region is not None and dinf_region is not None:
                        if not regions.regions_intersect_p(self.region, dinf_region):
                            usable = False

                if usable: proc(path_d, dformat, dweight)

    def _gen_inf(self):
        '''load a dalist and process datafiles'''

        self._proc_data(lambda x, t, w: generate_inf(x, t))
                
    def _append_datafile(self, dfile, dformat = 168, dweight = 1):
        '''append a datalist entry to a datalist, given filename, format and weight'''

        with open(self._path, 'a') as outfile:
            outfile.write('{} {} {}\n'.format(dfile, dformat, dweight))

        self._reset()

    def _echo_datafiles(self, osep='/n'):
        '''return a string of datafiles in the datalist'''

        data_paths = [x[0] for x in self.datafiles]
        return(osep.join(data_paths))

    def _cat_port(self, dst_port):
        '''Catenate the xyz data from the datalist'''
        try:
            for df in self.datafiles:
                if df[1] == 168: # xyz
                    with open(df[0]) as infile:
                        for line in infile:
                            dst_port.write(line)
                elif df[1] == 200: #gdal
                    gdalfun.dump(df[0], dst_port)
        except: sys.stderr.write('geomods: error, pipe broken.\n')

    def _caty(self):
        '''Catenate the xyz data from the datalist to a generator
        usage: for line in self._caty(): proc(line)'''

        for df in self.datafiles:
            if df[1] == 168: # xyz
                with open(df[0]) as infile:
                    for line in infile:
                        yield(line)
            elif df[1] == 200: #gdal
                gdalfun.dumpy(df[0])
        
    def _caty_xyz(self): #depr
        '''Catenate the xyz data from the datalist to a generator
        usage: for line in self._caty_xyz(): proc(line)'''
        
        data_paths = [x[0] for x in self.datafiles]
        for fn in data_paths:
            with open(fn) as infile:
                for line in infile:
                    yield(line)

    def _join(self, osep='\n'):
        df = []
        data_paths = [x[0] for x in self.datafiles]
        for fn in data_paths:
            with open(fn) as infile:
                for line in infile:
                    df.append(line)

        return osep.join(df)

    def gather_region(self):

        out_regions = []
        self._proc_data(lambda x, t, w: out_regions.append(datafile_region(x, t)))
        out_regions = [x for x in out_regions if x is not None]
        out_region = out_regions[0]
        for i in out_regions[1:]:
            out_region = regions.regions_merge(out_region, i)
                
        self.region = out_region
        
    def mask(self, inc = .000277777):
        '''Generate a num-msk with GDAL from a datalist.'''

        o_mask = '{}_num_msk.tif'.format(self._name)
        gdalfun.xyz_mask(self._caty(), o_mask, self.region.region, inc, verbose = False)
        
        return(o_mask)    
    
## =============================================================================
##
## Mainline - run datalists from console.
##
## =============================================================================

_usage = '''{} ({}): Process and generate datalists

usage: {} [ -hisvEOPR [ args ] ] datalist ...

Options:
  -R, --region\t\tSpecifies the desired region;
  -E, --increment\tCell size in native units
  -O, --output-name\tThe output file name.
  -P, --epsg\t\tSpecify the projection of the datalist.
  -F, --format\t\tSpecify the input datalist format.

Modules:

  -i, --info-file\tGenerate inf files for the data in the datalist
  -r, --region-info\tReturn the full region of the datalist
  -s, --spatial-md\tGenerate spatial metadata from the datalist

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} my_data.datalist -R -90/-89/30/31
 % {} my_data.datalist -R -90/-89/30/31 -E.000925 -s

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def main():

    status = 0
    i_datalist = None
    i_region = None
    i_inc = 0.000277777
    o_bn = None
    epsg = '4326'
    these_regions = []
    want_verbose = False
    want_inf = False
    want_region = False
    want_sm = False
    want_glob = False
    dl_fmt = 168
    
    argv = sys.argv
        
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

        elif arg == '--info-file' or arg == '-i':
            want_inf = True

        elif arg == '--region-info' or arg == '-r':
            want_region = True

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        else: 
            i_datalist = arg

        i = i + 1

    if i_datalist is None:
        print (_usage)
        sys.exit(1)
        
    utils.check_config()

    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if want_verbose: pb = utils._progress('loading region(s)...')
    if i_region is None:
        these_regions = [None]
    else:
        try: 
            these_regions = [regions.region(i_region)]
        except: these_regions = [regions.region('/'.join(map(str, x))) for x in gdalfun._ogr_extents(i_region)]

    if len(these_regions) == 0:
        these_regions = [None]

    for this_region in these_regions:
        if this_region is not None:
            if not this_region._valid:
                status = -1
    if want_verbose: pb.end(status, 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions)))
    
    if status == -1:
        utils._error_msg('failed to load region(s)')
        print(_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):
        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        this_datalist = datalist(i_datalist, this_region, dl_fmt, verbose = want_verbose)

        if want_glob:
            for f in known_fmts[dl_fmt]:
                globs = glob.glob('*.{}'.format(f))
                [this_datalist._append_datafile(x, 168, 1) for x in globs]
        
        if this_datalist._valid_p():

            ## ==============================================
            ## Generate spatial metadata for the datalist
            ## ==============================================

            if want_sm:
                import metadata

                if this_datalist.region is None:
                    this_datalist.gather_region()

                print this_datalist.region.region
                this_datalist.i_fmt = -1
                this_datalist._load_data()

                print this_datalist.datafiles
                sm = metadata.spatial_metadata(this_datalist, this_datalist.region, i_inc, o_bn, lambda: False, verbose = want_verbose)
                sm.run()

            ## ==============================================
            ## Generate inf files for the datalist
            ## ==============================================

            elif want_inf:
                this_datalist._gen_inf()

            ## ==============================================
            ## Print thre maximum datalist region
            ## ==============================================

            elif want_region:
                this_datalist.gather_region()
                sys.stdout.write('{}\n'.format(this_datalist.region.gmt))

            ## ==============================================
            ## load the data from the datalist and
            ## print to stdout
            ## ==============================================

            else:
                #print this_datalist.datalists
                #if this_datalist.region is None:
                #    this_datalist.gather_region()
                #this_datalist._load_data()
                #print this_datalist.datalist
                #print this_datalist.datalists
                print this_datalist.datafiles
                #this_datalist._cat_port(sys.stdout)

if __name__ == '__main__':
    main()    

### End
