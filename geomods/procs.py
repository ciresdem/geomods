### procs.py
##
## Copyright (c) 2020 CIRES Coastal DEM Team
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

import os
import sys
import time
import zipfile
import gzip
import threading
import Queue as queue
import ogr

import regions
import datalists
import gdalfun
import utils
import vdatum

import json

_version = '0.1.1'

## =============================================================================
##
## Data Processing Classes and functions.
## 
## =============================================================================

proc_infos = { 
    'gdal':[lambda x, a: x.proc_gdal(*a), 'process gdal data to chunked XYZ [:split_value:increment:input_vdatum]'],
    'lidar':[lambda x, a: x.proc_las(*a), 'process lidar las/laz data to XYZ [:increment:input_vdatum]'],
    'ascii':[lambda x, a: x.proc_ascii(*a), 'process ascii xyz data to XYZ [:delim:x_loc,y_loc,z_loc:skip_lines:increment:input_vdatum]'],
    'ogr':[lambda x, a: x.proc_ogr(*a), 'process OGR data. [layer:height]'],
    'mb':[lambda x, a: x.proc_mb(*a), 'process Multibeam data. []'],
    'ngs':[lambda x, a: x.proc_ngs(*a), 'process NGS json data. []'],
}

def _proc_queue(q):
    '''process queue `q` of fetched data'''

    while True:
        work = q.get()
        if not work[0][-1]():
            proc_args = tuple(work[0])
            proc_mod_args = tuple(work[1])
            if proc_args[1] is not None:
                #print proc_args
                proc_d = procs(*proc_args)
                #try:
                proc_d.run(proc_mod_args)
                #except:
                #   sys.stderr.write('\x1b[2K\r')
                #  sys.stderr.write('procs: error, unknown error processing {}\n'.format(work[0][0]))
                # sys.stderr.flush()
       
        q.task_done()

class procs_from_queue(threading.Thread):
    
    def __init__(self, proc_mod, proc_mod_opts, these_regions, data_files = [], callback = lambda: False):
        threading.Thread.__init__(self)
        self.proc_q = queue.Queue()
        self.data_files = data_files
        self.these_regions = these_regions
        self.proc_mod = proc_mod
        self.proc_mod_opts = proc_mod_opts
        self.stop = callback

    def run(self):
        for _ in range(3):
            t = threading.Thread(target = _proc_queue, args = (self.proc_q,))
            t.daemon = True
            t.start()
                
        for df in self.data_files:
            self.proc_q.put([[df, self.proc_mod, self.these_regions, self.stop], self.proc_mod_opts])

        self.proc_q.join()
        
class procs:
    '''process data to XYZ.'''

    def __init__(self, src_file, proc_mod, regions = [], callback = lambda: False):

        if callback is None:
            self.stop = lambda: False
        else: self.stop = callback
        self.status = 0
        self.xyzs = []

        ## Set src and dst variables
        self.src_file = os.path.abspath(src_file)
        self.src_dir = os.path.dirname(src_file)
        self.src_bn = os.path.basename(src_file)
        self.src_root = self.src_bn.split('.')[0]
        self.proc_mod = proc_mod
        self.dst_regions = regions
        self.proc_dir = os.path.join(self.src_dir, self.proc_mod)
        self.xyz_dir = os.path.join(self.proc_dir, 'xyz')
                                    
        ## Initialize VDatum
        self.this_vd = vdatum.vdatum()

        ## make the output xyz directory
        if not os.path.exists(self.xyz_dir):
            try:
                os.makedirs(self.xyz_dir)
            except: self.status = -1

        self.gdal_exts = ['.tif', '.img', '.asc']
        self.lidar_exts = ['.las', '.laz']
        self.ascii_exts = ['.xyz', '.ascii', '.csv', '.dat']
        self.ogr_exts = ['.000', '.shp', '.gmt']

    def run(self, args):
        '''Run the elevation data processor.
        Process data to WGS84/NAVD88 XYZ and append
        to a DATALIST.'''

        pb = utils._progress('processing local file: \033[1m{}\033[m...'.format(self.src_bn))
        if self.proc_mod in proc_infos.keys():
            #try:
            proc_infos[self.proc_mod][0](self, args)
            #except: self.status = -1
        else: self.status = -1

        if self.status == 0:
            self.add_to_datalist()

        pb.opm = 'processed local file: \033[1m{}\033[m.'.format(self.src_bn)
        pb.end(self.status)

    def _clean_zips(self, zip_files):
        '''remove all files in `zip_files`'''
        
        for i in zip_files:
            i_file = os.path.join(self.s_dir, self.s_t, i)
            if os.path.isfile(i_file):
                os.remove(i_file)
                zip_files = [x for x in zip_files if x != i]

        if len(zip_files) > 0:
            for i in zip_files:
                i_file = os.path.join(self.s_dir, self.s_t, i)
                if os.path.isdir(i_file):
                    try:
                        os.removedirs(i_file)
                    except: pass

    def unzip(self, zip_file):
        '''unzip `zip_file` and return a list of extracted file names.'''

        zip_ref = zipfile.ZipFile(zip_file)
        zip_files = zip_ref.namelist()
        zip_ref.extractall(self.proc_dir)
        zip_ref.close()

        return(zip_files)

    def gunzip(self, gz_file):
        '''gunzip `gz_file` and return the extracted file name.'''
        
        gz_split = gz_file.split('.')[:-1]
        guz_file = '{}.{}'.format(gz_split[0], gz_split[1])
        with gzip.open(self.src_file, 'rb') as in_gz, \
             open(guz_file, 'w') as f:
            while True:
                block = in_gz.read(65536)
                if not block:
                    break
                else: f.write(block)
        return(guz_file)
        
    def add_to_datalist(self):
        '''Add xyz files in self.xyzs to the datalist.'''

        for o_xyz in self.xyzs:
            if os.stat(o_xyz).st_size != 0:
                dl_f = '{}.datalist'.format(os.path.abspath(self.src_dir).split('/')[-1])
                s_datalist = datalists.datalist(os.path.join(self.xyz_dir, dl_f))
                s_datalist._append_datafile(['{}'.format(os.path.basename(o_xyz)), 168, 1])
                s_datalist._reset()

                out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(self.xyz_dir, dl_f)), False, False)
            else: os.remove(o_xyz)

    def run_vdatum(self, src_xyz, dst_xyz, i_vert = 'MLLW', o_vert = 'NAVD88'):
        '''run vdatum wrapper'''

        if self.this_vd.vdatum_path is not None:
            self.this_vd.ivert = i_vert
            self.this_vd.overt = o_vert
            self.this_vd.ds_dir = os.path.relpath(os.path.join(os.path.dirname(src_xyz), 'result'))

            tmp_xyz = os.path.basename(src_xyz.lower())
            os.rename(src_xyz, tmp_xyz)
            self.this_vd.run_vdatum(os.path.relpath(tmp_xyz))
            os.rename(os.path.join(self.this_vd.ds_dir.lower(), os.path.basename(tmp_xyz)), dst_xyz)
            os.remove(tmp_xyz)
        else: self.status = -1
        
    def xyz_region(self, src_xyz):
        out, self.status = utils.run_cmd('gmt gmtinfo {} -I-'.format(src_xyz), False, False)
        o_region = regions.region(out[2:])
        return(o_region)
        
    def proc_xyz(self, src_xyz, block = None, i_vert = None):
        '''process geographic XYZ data to NAVD88 XYZ via blockmedian and vdatum.'''
        
        if os.stat(src_xyz).st_size != 0:
            data_region = self.xyz_region(src_xyz)

            for rn, this_region in enumerate(self.dst_regions):
                self.status = 0

                if regions.regions_intersect_p(this_region, data_region):
                    tmp_xyz = os.path.join(self.xyz_dir, '{}_tmp.xyz'.format(os.path.basename(src_xyz).split('.')[0]))
                    dst_xyz = os.path.join(self.xyz_dir, '{}_{}.xyz'.format(os.path.basename(src_xyz).split('.')[0], this_region.fn))

                    out, self.status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR=space', False, False)
                    if block is not None:
                        try:
                            inc = float(block)
                        except: inc = .000030864
                        out, self.status = utils.run_cmd('gmt blockmedian {} -I{} {} -r -V > {}\
                        '.format(src_xyz, inc, this_region.gmt, tmp_xyz), False, False)
                    else:
                        out, self.status = utils.run_cmd('gmt gmtselect {} {} -V -o0,1,2 > {}\
                        '.format(src_xyz, this_region.gmt, tmp_xyz), False, False)

                    if os.stat(tmp_xyz).st_size != 0:
                        if i_vert is not None:
                            self.run_vdatum(tmp_xyz, dst_xyz, i_vert = i_vert)
                        else: os.rename(tmp_xyz, dst_xyz)
                        
                    if os.path.exists(tmp_xyz):
                        os.remove(tmp_xyz)

                    if not os.path.exists(dst_xyz):
                        self.status = -1
                    elif os.stat(dst_xyz).st_size == 0:
                        self.status = -1
                        os.remove(dst_xyz)
                    else: self.xyzs.append(dst_xyz)

    def proc_ngs(self):
        import csv
        
        with open(self.src_file, 'r') as json_file:
            r = json.load(json_file)
        
        if len(r) > 0:
            for rn, this_region in enumerate(self.dst_regions):
                outfile = open(os.path.join(self.src_dir, 'ngs_results_{}.csv'.format(this_region.fn)), 'w')
                outcsv = csv.writer(outfile)
                outcsv.writerow(r[0].keys())
                [outcsv.writerow(row.values()) for row in r]
                outfile.close()
                    
    def proc_las(self, block = None, i_vert = None, classes = '2 29'):
        '''process las/laz lidar data to XYZ.'''

        out, self.status = utils.run_cmd('las2txt -verbose -parse xyz -keep_class {} -i {}'.format(classes, self.src_file), True, True)

        if self.status == 0:
            o_fn_txt = os.path.join(self.src_dir, '{}.txt'.format(self.src_file.split('.')[0]))

            self.proc_xyz(o_fn_txt, block, i_vert)
            os.remove(o_fn_txt)

    def proc_ascii(self, delim = ' ', xyz_cols = '0,1,2', skip = 0, block = None, i_vert = None):

        if self.src_file.split('.')[-1] == 'zip':
            zips = self.unzip(self.src_file)

            for ext in self.ascii_exts:
                for zf in zips:
                    if ext in zf:
                        src_ascii = os.path.join(self.proc_dir, zf)
                        print src_ascii
                        break
        elif self.src_file.split('.')[-1] == 'gz':
            tmp_ascii = self.gunzip(self.src_file)
            src_ascii = os.path.join(self.proc_dir, os.path.basename(tmp_ascii))
            os.rename(tmp_ascii, src_ascii)
        else: src_ascii = self.src_file

        xyz_loc = map(int, xyz_cols.split(','))
        tmp_ascii = os.path.join(self.proc_dir, '{}_tmp.ascii'.format(os.path.basename(src_ascii).split('.')[0]))
        
        with open(src_ascii, 'r') as in_ascii,\
             open(tmp_ascii, 'w') as out_ascii:
            l_n = 1
            for line in in_ascii:
                if l_n > int(skip):
                    row = line.rstrip().split(delim)
                    if len(row) > 2:
                        out_row = '{} {} {}\n'.format(row[xyz_loc[0]], row[xyz_loc[1]], row[xyz_loc[2]])
                        out_ascii.write(out_row)
                l_n += 1

        self.proc_xyz(tmp_ascii, block, i_vert)
        os.remove(tmp_ascii)

    def proc_mb(self, region = None):
        '''process mb data'''

        ## mblist to xyz
        #print self.dst_regions[0].gmt
        src_xyz = os.path.basename(self.src_file).split('.')[0] + '.xyz'
        mbl_cmd = 'mblist -MX20 -OXYZ -I{}  > {}'.format(self.src_file, src_xyz)
        out, self.status = utils.run_cmd(mbl_cmd, True, True)
        
        self.proc_xyz(src_xyz, block = None, i_vert = 'lmsl')
        
    def proc_ogr(self, layer = None, z = 'Height', i_vert = None):
        '''process OGR data to XYZ via OGR.'''

        if os.path.basename(self.src_file).split('.')[-1].lower() == 'zip':
            zips = self.unzip(self.src_file)
            for ext in self.ogr_exts:
                for zf in zips:
                    if ext in zf.lower():
                        src_ogr = os.path.join(self.proc_dir, zf)
                        break
        else: src_ogr = self.src_file
        
        dst_xyz = os.path.join(self.proc_dir, os.path.basename(src_ogr).split('.')[0] + '.xyz')
        
        ds_ogr = ogr.Open(src_ogr)
        if layer is not None:
            layer_s = ds_ogr.GetLayerByName(layer) #'SOUNDG'
        else: layer_s = ds_ogr.GetLayer()
        
        if layer_s is not None:
            o_xyz = open(dst_xyz, 'w')
            for f in layer_s:
                g = json.loads(f.GetGeometryRef().ExportToJson())
                for xyz in g['coordinates']:
                    if z != 'Height':
                        xyz_l = '{} {} {}\n'.format(xyz[0], xyz[1], xyz[2])
                    else: xyz_l = '{} {} {}\n'.format(xyz[0], xyz[1], xyz[2]*-1)
                    o_xyz.write(xyz_l)
            o_xyz.close()
        else: self.status = -1

        ds_ogr = layer_s = None

        if os.path.exists(dst_xyz):
            self.proc_xyz(dst_xyz, block = None, i_vert = i_vert)
    
    def proc_gdal(self, split = None, block = None, i_vert = None):
        '''process gdal grid data to XYZ.'''
        
        if self.src_file.split('.')[-1].lower() == 'zip':
            zips = self.unzip(self.src_file)

            for ext in self.gdal_exts:
                for zf in zips:
                    if ext in zf.lower():
                        src_gdal = os.path.join(self.proc_dir, zf)
                        break
        elif self.src_file.split('.')[-1].lower() == 'gz':
            tmp_gdal = self.gunzip(self.src_file)
            src_gdal = os.path.join(self.proc_dir, os.path.basename(tmp_gdal))
            os.rename(tmp_gdal, src_gdal)
        else: src_gdal = self.src_file

        out_chunks = gdalfun.chunks(src_gdal, 1000)
        for chunk in out_chunks:
            if not self.stop():
                if split is not None:
                    split_chunk = gdalfun.split(chunk, int(split))
                    chunk_gdal = split_chunk[0]
                else: chunk_gdal = chunk
                
                #if reproject:
                tmp_gdal = os.path.join(self.proc_dir, '{}_tmp.tif'.format(os.path.basename(chunk_gdal).split('.')[0]))
                out, self.status = utils.run_cmd('gdalwarp {} {} -t_srs EPSG:4326\
                '.format(chunk_gdal, tmp_gdal), False, False)
                #else: tmp_gdal = chunk_gdal

                tmp_xyz = os.path.join(self.xyz_dir, '{}.tmp'.format(os.path.basename(chunk_gdal).split('.')[0]))
                tmp_xyz_ob = open(tmp_xyz, 'wt')
                gdalfun.dump(tmp_gdal, tmp_xyz_ob)
                tmp_xyz_ob.close()
                self.proc_xyz(tmp_xyz, block, i_vert)

                os.remove(tmp_xyz)
                os.remove(tmp_gdal)
                #if not reproject:
                os.remove(chunk)
                if split is not None:
                    os.remove(split_chunk[0])
                    os.remove(split_chunk[1])

    ## NHDS Shapefiles [ TNM:2:Shapefile ]
    ## unzip only nhd_waterbodies
    ## clip to region and remove null feature name
    ## merge and union
    ## unzip only nhd_area
    ## merge and uin
    ## merge and union nhd_area and nhd_waterbody
    ## edit via gis to cover region in water.
    ## make a waffles 'coastline' module using tnm and/or gmt gsshg
                    
## =============================================================================
##
## Mainline - run  from console.
##
## =============================================================================

def procs_desc(x):
    fd = []
    for key in x: 
        fd.append('{:10}\t\t{}'.format(key, x[key][1]))
    return fd

_usage = '''{} ({}): Process and generate datalists

usage: {} [ -hvMR [ args ] ] datalist ...

Modules:
  {}

Options:
  -R, --region\t\tSpecifies the desired region.
  -M, --module\t\tData procs module.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} -R -90/-89/30/31 -M gdal::.0000925 dc_raster.img
 % {} -R tiles.shp -M gdal:0::mllw ned13.zip
 % {} -M ascii:,:2,1,3:1::mllw *.xyz.gz
 % {} -M lidar:.000030864 *.laz

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]),
            '\n  '.join(procs_desc(proc_infos)),
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]),
            os.path.basename(sys.argv[0]))

def main():

    status = 0
    i_data = []
    i_region = None
    proc_opts = {}
    want_verbose = False
    stop_threads = False
    these_regions = []
    
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

        elif arg == '--module' or arg == '-M':
            opts = argv[i + 1].split(':')
            proc_opts[opts[0]] = list(opts[1:])
            i = i + 1
        elif arg[:2] == '-M':
            opts = arg[2:].split(':')
            proc_opts[opts[0]] = list(opts[1:])

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
            i_data.append(arg)

        i = i + 1

    if len(i_data) == 0 or len(proc_opts.keys()) == 0:
        print(_usage)
        sys.exit(0)
        
    for key in proc_opts.keys():
        proc_opts[key] = [None if x == '' else x for x in proc_opts[key]]
        
    if i_region is None:
        i_region = '-180/180/0/90'

        #this_region = regions.region(i_region)

    utils.check_config()
    
    pb = utils._progress('loading region(s)...')
    try: 
        these_regions = [regions.region(i_region)]
    except: these_regions = [regions.region('/'.join(map(str, x))) for x in gdalfun._ogr_extents(i_region)]

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
          status = -1
    pb.opm = 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions))
    pb.end(status)
            
    if status == -1:
        print(_usage)
        sys.exit(1)

    for proc_mod in proc_opts:
        if proc_mod in proc_infos.keys():
            pb = utils._progress('running procs module \033[1m{}\033[m on \033[1m{}\033[m data file(s) in \033[1m{}\033[m region(s)...\
            '.format(proc_mod, len(i_data),len(these_regions)))
            args = proc_opts[proc_mod]
            df = procs_from_queue(proc_mod, args, these_regions, i_data, lambda: stop_threads)

            try:
                df.start()
                while True:
                    time.sleep(1)
                    pb.update()
                    if not df.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): 
                status = -1
                stop_threads = True

            pb.opm = 'ran procs module \033[1m{}\033[m on \033[1m{}\033[m data file(s) in \033[1m{}\033[m region(s).\
            '.format(proc_mod, len(i_data), len(these_regions))
            pb.end(status)

            
if __name__ == '__main__':
    main()

        
## End
