### fetches.py
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

import os
import sys
import time
import math

import requests
import lxml.html as lh
import lxml.etree

import zipfile
import gzip
import csv
import json
import threading
import Queue as queue

import numpy as np

try:
    import osgeo.ogr as ogr
    import osgeo.gdal as gdal
except ImportError:
    try:
        import ogr
        import gdal
    except ImportError:
        sys.exit(-1)

import regions
import datalists
import gdalfun
import utils

_version = '0.3.0'

gdal.PushErrorHandler('CPLQuietErrorHandler')

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## =============================================================================

r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(_version) }

namespaces = {'gmd': 'http://www.isotc211.org/2005/gmd', 
              'gmi': 'http://www.isotc211.org/2005/gmi', 
              'gco': 'http://www.isotc211.org/2005/gco',
              'gml': 'http://www.isotc211.org/2005/gml'}

def fetch_queue(q, p = None):
    '''fetch queue `q` of fetch results'''

    while True:
        fetch_args = q.get()
        this_region = fetch_args[2]
        this_dt = fetch_args[4]
        fetch_args[2] = None

        if not fetch_args[3]():
            fetch_file(*tuple(fetch_args))
            if p is not None:
                outf = fetch_args[1]
                surv_dir = os.path.dirname(outf)
                surv_fn = os.path.basename(outf)
                #surv_t = fetch_args[0].split('/')[-2]

                if os.path.exists(outf):
                    p.put([[surv_dir, surv_fn, outf, this_dt, this_region], fetch_args[3]])

        q.task_done()

def fetch_file(src_url, dst_fn, params = None, callback = None, datatype = None):
    '''fetch src_url and save to dst_fn'''
    
    status = 0
    req = None
    halt = callback
    pb = utils._progress('fetching remote file: \033[1m{}\033[m...'.format(os.path.basename(src_url)))

    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 

    try:
        req = requests.get(src_url, stream = True, params = params, headers = r_headers)
    except requestion.exceptions.RequestExceptions as e:
        print 'Error: {}'.format(e)
        status = -1

    if req is not None:
        with open(dst_fn, 'wb') as local_file:
            for chunk in req.iter_content(chunk_size = 50000):
                if chunk:
                    if halt(): 
                        status = -1
                        break
                    local_file.write(chunk)

    pb.opm = 'fetched remote file: \033[1m{}\033[m.'.format(os.path.basename(src_url))
    pb.end(status)
    return(status)

def fetch_req(src_url, params = None, tries = 5, timeout = 2):
    '''fetch src_url and return the requests object'''

    if tries <= 0: return(None)
    try:
        return(requests.get(src_url, stream = True, params = params, timeout = timeout, headers = r_headers))
    except: return(fetch_req(src_url, params = params, tries = tries - 1, timeout = timeout + 1))

def fetch_nos_xml(src_url):
    '''fetch src_url and return it as an XML object'''

    results = lxml.etree.fromstring('<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8'))

    try:
        req = fetch_req(src_url)
        results = lxml.etree.fromstring(req.text.encode('utf-8'))
    except: pass

    return(results)
        
def fetch_html(src_url):
    '''fetch src_url and return it as an HTML object'''

    req = fetch_req(src_url)
    if req:
        return(lh.document_fromstring(req.text))
    else: return(None)

def fetch_csv(src_url):
    '''fetch src_url and return it as a CSV object'''

    req = fetch_req(src_url)
    if req:
        return(csv.reader(req.text.split('\n'), delimiter = ','))
    else: return(None)

class fetch_results(threading.Thread):
    '''fetch results'''

    def __init__(self, results, region, out_dir, want_proc = False, stop_threads = lambda: False):
        threading.Thread.__init__(self)

        self.fetch_q = queue.Queue()

        self.results = results
        self.region = region
        self._outdir = out_dir
        self.want_proc = want_proc
        self.stop_threads = stop_threads
        
    def run(self):
        if self.want_proc:
            proc_q = queue.Queue()
            for _ in range(3):
                p = threading.Thread(target = proc_queue, args = (proc_q,))
                p.daemon = True
                p.start()
        else: proc_q = None

        for _ in range(3):
            t = threading.Thread(target = fetch_queue, args = (self.fetch_q, proc_q))
            t.daemon = True
            t.start()

        for row in self.results:
            self.fetch_q.put([row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2]])

        self.fetch_q.join()
        if self.want_proc:
            proc_q.join()

## =============================================================================
##
## Processing Classes and functions.
## 
## =============================================================================

proc_infos = { 
    'GEODAS_XYZ':lambda x: x.proc_geodas(),
    'GRID_BAG':lambda x: x.proc_bag(),
    'ENC':lambda x: x.proc_enc(),
    'lidar':lambda x: x.proc_dc_las(),
    'raster':lambda x: x.proc_dc_raster(),
    'srtm':lambda x: x.proc_dc_raster(),
}

def proc_queue(q):
    '''process queue `q` of fetched data'''

    while True:
        work = q.get()

        proc_d = proc(work[0], work[1])
        if not work[1]():
            proc_d.run()
        
        q.task_done()

class proc:
    '''Process fetched data to XYZ'''

    def __init__(self, s_file_info, callback = None):

        self.stop = callback
        self.status = 0
        self.xyzs = []

        self.s_dir = s_file_info[0]
        self.s_fn = s_file_info[1]
        self.o_fn = s_file_info[2]
        self.s_t = s_file_info[3]
        self.s_region = s_file_info[4]

        self.this_vd = utils.vdatum()

        self.xyz_dir = os.path.join(self.s_dir, self.s_t.lower(), 'xyz')

        if not os.path.exists(self.xyz_dir):
            try:
                os.makedirs(self.xyz_dir)
            except: self.status = -1

        self.o_fn_bn = os.path.basename(self.o_fn).split('.')[0]
        self.o_fn_tmp = os.path.join(self.s_dir, '{}_tmp.xyz'.format(self.o_fn_bn.lower()))
        self.o_fn_p_xyz = os.path.join(self.s_dir, '{}.xyz'.format(self.o_fn_bn))
        self.o_fn_xyz = os.path.join(self.xyz_dir, '{}_{}.xyz'.format(self.o_fn_bn, self.s_region.fn))
        self.o_fn_p_ras = os.path.join(self.s_dir, '{}.tif'.format(self.o_fn_bn))
    
    def run(self):
        '''Run the elevation data processor.
        Process data to WGS84/NAVD88 XYZ and append
        to a DATALIST.'''

        pb = utils._progress('processing local file: \033[1m{}\033[m'.format(self.s_fn))
        if self.s_t in proc_infos.keys():
            proc_infos[self.s_t](self)

            if self.status == 0:
                self._add_to_datalist()
                if os.path.exists(self.o_fn):
                    os.remove(self.o_fn)
        else: self.status = -1

        pb.opm = 'processed local file: \033[1m{}\033[m'.format(self.s_fn)
        pb.end(self.status)

    def _add_to_datalist(self):
        '''Add xyz files in self.xyzs to the datalist.'''

        for o_xyz in self.xyzs:
            if os.stat(o_xyz).st_size != 0:
                sdatalist = datalists.datalist(os.path.join(self.xyz_dir, '{}.datalist'.format(os.path.abspath(self.s_dir).split('/')[-1])))
                sdatalist._append_datafile('{}'.format(os.path.basename(o_xyz)), 168, 1)
                sdatalist._reset()

                out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(self.xyz_dir, '{}.datalist'.format(os.path.basename(self.s_dir)))), False, None)
            else: os.remove(o_xyz)

    def proc_gmrt(self):
        '''Process GMRT geotiff to XYZ.'''

        gdalfun.dump(self.o_fn, self.o_fn_xyz)
        if os.stat(self.o_fn_xyz).st_size != 0:
            self.xyzs.append(self.o_fn_xyz)
        else: 
            self.status = -1
            os.remove(self.o_fn_xyz)
        
    def proc_geodas(self):
        '''Process NOAA NOS data to WGS84/NAVD88 XYZ'''

        try:
            in_f = gzip.open(os.path.join(self.s_dir, self.s_fn), 'rb')
            s = in_f.read()
            s_csv = csv.reader(s.split('\n'), delimiter = ',')
            next(s_csv, None)
        except: self.status = -1

        if self.status == 0 and not self.stop():
            with open(os.path.join(self.o_fn_tmp), 'w') as out_xyz:
                d_csv = csv.writer(out_xyz, delimiter = ' ')

                for row in s_csv:
                    if len(row) > 2:
                        this_xyz = [float(row[2]), float(row[1]), float(row[3]) * -1]
                        d_csv.writerow(this_xyz)
            in_f.close()

            if self.this_vd.vdatum_path is not None:
                self.this_vd.ivert = 'mllw'
                self.this_vd.overt = 'navd88'
                self.this_vd.ds_dir = os.path.relpath(os.path.join(self.xyz_dir, 'result'))

                self.status = self.this_vd.run_vdatum(os.path.relpath(self.o_fn_tmp))

                if self.status == 0 and os.path.exists(os.path.join(self.xyz_dir, 'result', os.path.basename(self.o_fn_tmp))): 
                    os.rename(os.path.join(self.xyz_dir, 'result', os.path.basename(self.o_fn_tmp)), self.o_fn_xyz)
            else: os.rename(self.o_fn_tmp, self.o_fn_xyz)

            if os.path.exists(self.o_fn_tmp):
                os.remove(self.o_fn_tmp)

            if os.path.exists(self.o_fn_xyz):
                self.xyzs.append(self.o_fn_xyz)

    def proc_bag(self):
        '''Process NOAA NOS BAG data to WGS84/NAVD88 XYZ.'''

        s_gz = os.path.join(self.s_dir, self.s_fn)

        if s_gz.split('.')[-1] == 'gz':
            with gzip.open(s_gz, 'rb') as in_bag:
                s = in_bag.read()
                with open(s_gz[:-3], 'w') as f:
                    f.write(s)

            s_bag = s_gz[:-3]
        elif s_gz.split('.')[-1] == 'bag':
            s_bag = s_gz
        else: self.status = -1

        if s_bag.split('.')[-1] != 'bag':
            self.status = -1

        s_tif = os.path.join(self.s_dir, self.s_fn.split('.')[0].lower() + '.tif')            
        s_xyz = s_gz.split('.')[0] + '.xyz'

        if self.status == 0 and not self.stop():
            out_chunks = gdalfun.chunks(s_bag, 5000)

            if s_gz.split('.')[-1] == 'gz':
                os.remove(s_bag)

            for chunk in out_chunks:
                if self.stop():
                    break

                i_xyz = '{}_{}.xyz'.format(chunk.split('.')[0], self.s_region.fn)
                i_tif = chunk.split('.')[0] + '_wgs.tif'
                o_xyz = os.path.join(self.xyz_dir, os.path.basename(i_xyz).lower())

                out, self.status = utils.run_cmd('gdalwarp {} {} -t_srs EPSG:4326'.format(chunk, i_tif), False, None)
                
                if self.status != 0 or self.stop():
                    os.remove(chunk)
                    break

                os.remove(chunk)
                gdalfun.dump(i_tif, i_xyz)
                os.remove(i_tif)

                if os.stat(i_xyz).st_size == 0:
                    self.status = -1
                    os.remove(i_xyz)
                    break

                if self.this_vd.vdatum_path is not None:
                    self.this_vd.ivert = 'mllw'
                    self.this_vd.overt = 'navd88'
                    self.this_vd.ds_dir = os.path.relpath(os.path.join(self.xyz_dir, 'result'))
                    
                    self.this_vd.run_vdatum(os.path.relpath(i_xyz))
                    
                    os.rename(os.path.join(self.xyz_dir, 'result', os.path.basename(i_xyz)), o_xyz)
                    os.remove(i_xyz)
                else: os.rename(i_xyz, o_xyz)

                self.xyzs.append(o_xyz)

    def proc_dc_las(self):
        '''Process NOAA Digital Coast lidar data (las/laz) to Ground XYZ.'''

        out, self.status = utils.run_cmd('las2txt -verbose -parse xyz -keep_class 2 29 -i {}'.format(self.o_fn), False, None)
        if self.status == 0:
            self.o_fn_bn = os.path.basename(self.o_fn).split('.')[0]
            self.o_fn_txt = os.path.join(self.s_dir, '{}.txt'.format(self.o_fn_bn))

            out, self.status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR=space', False, None)
            out, self.status = utils.run_cmd('gmt blockmedian {} -I.1111111s {} -r -V > {}'.format(self.o_fn_txt, self.s_region.gmt, self.o_fn_xyz), False, None)

            os.remove(self.o_fn_txt)

            if self.status == 0:
                self.xyzs.append(self.o_fn_xyz)

    def proc_dc_raster(self):
        '''Process NOAA Digital Coast raster data to WGS84 XYZ chunks.'''

        out, self.status = utils.run_cmd('gdalwarp {} {} -t_srs EPSG:4326'.format(self.o_fn, self.o_fn_p_ras), False, None)

        out_chunks = gdalfun.chunks(self.o_fn_p_ras, 1000)
        os.remove(self.o_fn_p_ras)

        for chunk in out_chunks:

            i_xyz = chunk.split('.')[0] + '.xyz'
            o_xyz = os.path.join(self.xyz_dir, os.path.basename(i_xyz))

            gdalfun.dump(chunk, o_xyz)

            if os.stat(o_xyz).st_size != 0:
                self.xyzs.append(o_xyz)
            else:
                os.remove(o_xyz)

            #os.remove(i_xyz)

    def proc_enc(self):
        '''Proces Electronic Nautical Charts to NAVD88 XYZ.'''

        zip_ref = zipfile.ZipFile(self.o_fn)
        zip_ref.extractall(os.path.join(self.s_dir, 'enc'))
        zip_ref.close()

        s_fn_000 = os.path.join(self.s_dir, 'enc/ENC_ROOT/', self.o_fn_bn, '{}.000'.format(self.o_fn_bn))

        ds_000 = ogr.Open(s_fn_000)
        layer_s = ds_000.GetLayerByName('SOUNDG')
        if layer_s is not None:
            o_xyz = open(self.o_fn_p_xyz, 'w')
            for f in layer_s:
                g = json.loads(f.GetGeometryRef().ExportToJson())
                for xyz in g['coordinates']:
                    xyz_l = '{} {} {}\n'.format(xyz[0], xyz[1], xyz[2]*-1)
                    o_xyz.write(xyz_l)
            o_xyz.close()
        else: self.status = -1

        if os.path.exists(self.o_fn_p_xyz) and os.stat(self.o_fn_p_xyz).st_size != 0:

            out, self.status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR=space', False, None)
            out, self.status = utils.run_cmd('gmt gmtselect {} {} -V > {}'.format(self.o_fn_p_xyz, self.s_region.gmt, self.o_fn_tmp), False, None)

            if os.stat(self.o_fn_tmp).st_size != 0:

                if self.this_vd.vdatum_path is not None:
                    self.this_vd.ivert = 'mllw'
                    self.this_vd.overt = 'navd88'
                    self.this_vd.ds_dir = os.path.relpath(os.path.join(self.xyz_dir, 'result'))

                    self.this_vd.run_vdatum(os.path.relpath(self.o_fn_tmp))

                    os.rename(os.path.join(self.xyz_dir, 'result', os.path.basename(self.o_fn_tmp)), self.o_fn_xyz)
                else: os.rename(self.o_fn_tmp, self.o_fn_xyz)

                os.remove(self.o_fn_p_xyz)

                if os.stat(self.o_fn_xyz).st_size == 0:
                    os.remove(self.o_fn_xyz)
                    self.status = -1
                else:
                    self.xyzs.append(self.o_fn_xyz)

        if os.path.exists(self.o_fn_tmp):
            os.remove(self.o_fn_tmp)

## =============================================================================
##
## Reference Vector
##
## the reference vector location and related functions
##
## =============================================================================

this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

def create_polygon(coords):
    '''convert coords to Wkt'''

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[1], coord[0])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
        
    return(poly.ExportToWkt())

def bounds2geom(bounds):
    '''convert a bounds [west, east, south, north] to an 
    OGR geometry'''

    b1 = [[bounds[2], bounds[0]],
          [bounds[2], bounds[1]],
          [bounds[3], bounds[1]],
          [bounds[3], bounds[0]],
          [bounds[2], bounds[0]]]

    return(ogr.CreateGeometryFromWkt(create_polygon(b1)))

def addf_ref_vector(ogr_layer, survey):
    '''add a survey to the reference vector layer'''

    layer_defn = ogr_layer.GetLayerDefn()
    
    feat = ogr.Feature(layer_defn)
    feat.SetGeometry(survey[0])
    feat.SetField('Name', str(survey[1]))
    feat.SetField('ID', str(survey[2]))
    feat.SetField('Date', int(survey[3]))
    feat.SetField('XML', str(survey[4]))
    feat.SetField('Data', str(survey[5]))
    feat.SetField('Datatype', str(survey[6]))
        
    ogr_layer.CreateFeature(feat)
    
    feat = geom = None

def update_ref_vector(src_vec, surveys, update=True):
    '''update or create a reference vector'''

    layer = None

    if update:
        ds = ogr.GetDriverByName('GMT').Open(src_vec, 1)
        if ds is not None:
            layer = ds.GetLayer()
    else:
        ds = ogr.GetDriverByName('GMT').CreateDataSource(src_vec)
        layer = ds.CreateLayer('%s' %(os.path.basename(src_vec).split('.')[0]), None, ogr.wkbPolygon)
        layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
        layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('XML', ogr.OFTString))
        data_field = ogr.FieldDefn('Data', ogr.OFTString)
        data_field.SetWidth(1000)
        layer.CreateField(data_field)
        layer.CreateField(ogr.FieldDefn('Datatype', ogr.OFTString))
  
    if layer is not None:
        for feature in layer:
            layer.SetFeature(feature)

        for i in surveys:
            addf_ref_vector(layer, i)

    ds = layer = None

## =============================================================================
##
## DC Fetch - NOAA Digital Coast
##
## -- ftp.coast.noaa.gov/pub/DigitalCoast --
## -- https://coast.noaa.gov/digitalcoast --
##
## Raster and Lidar data from the Digital Coast
## Most raster data is Tiff and most lidar data is LAZ
## Will need some laz processing utilities to process the lidar data,
## such as lastools or liblas, etc.
##
## =============================================================================

class dc:
    '''Fetch elevation data from the Digital Coast'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading Digital Coast fetch module...')
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dav_id = 'https://coast.noaa.gov/dataviewer/#/lidar/search/where:ID='
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'raster2']

        self._ref_vector = os.path.join(fetchdata, 'dc.gmt')
        self._outdir = os.path.join(os.getcwd(), 'dc')

        self._status = 0
        self._surveys = []
        self._results = []

        self._index = False

        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)
        else: self._status = -1

        pb.opm = 'loaded Digital Coast fetch module.'
        pb.end(self._status)
        
    def run(self):
        '''Run the Digital Coast fetching module'''

        if self._want_update:
            self._update()
        else:
            self.search_gmt()            
            return(self._results)

    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================

    def _update(self):
        '''Update the DC reference vector after scanning
        the relevant metadata from Digital Coast.'''

        for ld in self._dc_dirs:
            if self._has_vector:
                gmt2 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
                layer = gmt2.GetLayer()
            else: layer = []

            page = fetch_html(self._dc_htdata_url + ld)

            tables = page.xpath('//table')
            tr = tables[0].xpath('.//tr')

            if len(tr) <= 0:
                break

            cols = []
            for i in tr[0]: cols.append(i.text_content())

            ## ==============================================
            ## Collect relavant info from metadata
            ## ==============================================

            for i in range(1, len(tr)):
                if not self.stop():
                    cells = tr[i].getchildren()
                    dc = {}
                    for j, cell in enumerate(cells):
                        cl = cell.xpath('a')
                        if len(cl) > 0:
                            if cols[j] == 'Dataset Name':
                                dc[cols[j]] = cell.text_content()
                                dc['Metadata'] = cl[0].get('href')
                            else: dc[cols[j]] = cl[0].get('href')
                        else: dc[cols[j]] = cell.text_content()

                    if self._has_vector: 
                        layer.SetAttributeFilter('ID = "%s"' %(dc['ID #']))

                    if len(layer) == 0:
                        if 'Metadata' in dc.keys():
                            xml_doc = fetch_nos_xml(dc['Metadata'])

                            wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal',
                                              namespaces = namespaces)
                            el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal',
                                              namespaces = namespaces)
                            sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal',
                                              namespaces = namespaces)
                            nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal',
                                              namespaces = namespaces)
                            if wl is not None:
                                obbox = bounds2geom([float(wl.text), float(el.text),
                                                     float(sl.text), float(nl.text)])

                                try: 
                                    odate = int(dc['Year'][:4])
                                except: odate = 1900

                                ## ==============================================
                                ## Append survey to surveys list
                                ## ==============================================

                                out_s = [obbox, 
                                         dc['Dataset Name'], 
                                         dc['ID #'], 
                                         odate, 
                                         dc['Metadata'], 
                                         dc['https'], 
                                         ld.split("_")[0]]

                                self._surveys.append(out_s)

            ## ==============================================
            ## Add matching surveys to reference vector
            ## ==============================================

            if len(self._surveys) > 0:
                update_ref_vector(self._ref_vector, self._surveys, self._has_vector)

            if os.path.exists(self._ref_vector): 
                self._has_vector = True
            else: self._has_vector = False

            self._surveys = []
            gmt2 = layer = None

    ## ==============================================
    ## Filter for results
    ## ==============================================

    def search_gmt(self):
        '''Search for data in the reference vector file'''

        tw = utils._progress('filtering Digital Coast reference vector...')
        gmt1 = ogr.Open(self._ref_vector)
        layer = gmt1.GetLayer(0)

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        ## ==============================================
        ## Find surveys that fit all filters and within
        ## the region of interest.
        ## ==============================================

        for feature1 in layer:
            if self.stop() or self._status != 0:
                break

            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                surv_url = feature1.GetField('Data')
                surv_id = surv_url.split('/')[-2]
                surv_dt = feature1.GetField('Datatype')
                suh = fetch_html(surv_url)
                if suh is None: 
                    self._status = -1
                    break

                #if self._status == 0:
                if 'lidar' in surv_dt:

                    ## ==============================================
                    ## Lidar data has a minmax.csv file to get extent
                    ## for each survey file.
                    ## ==============================================

                    scsv = suh.xpath('//a[contains(@href, ".csv")]/@href')[0]
                    dc_csv = fetch_csv(surv_url + scsv)
                    next(dc_csv, None)

                    for tile in dc_csv:
                        if len(tile) > 3:
                            tb = [float(tile[1]), float(tile[2]),
                                  float(tile[3]), float(tile[4])]
                            tile_geom = bounds2geom(tb)
                            if tile_geom.Intersects(self._boundsGeom):
                                self._results.append([os.path.join(surv_url, tile[0]), '{}/{}'.format(surv_id, tile[0]), 'lidar'])

                elif 'raster' in surv_dt:

                    ## ==============================================
                    ## Raster data has a tileindex shapefile 
                    ## to get extent
                    ## ==============================================

                    sshpz = suh.xpath('//a[contains(@href, ".zip")]/@href')[0]
                    fetch_file(surv_url + sshpz, os.path.join('.', sshpz),
                               callback = self.stop)

                    zip_ref = zipfile.ZipFile(sshpz)
                    zip_ref.extractall('dc_tile_index')
                    zip_ref.close()

                    if os.path.exists('dc_tile_index'):
                        ti = os.listdir('dc_tile_index')

                    ts = None
                    for i in ti:
                        if ".shp" in i:
                            ts = os.path.join('dc_tile_index/', i)

                    if ts is not None:
                        shp1 = ogr.Open(ts)
                        slay1 = shp1.GetLayer(0)

                        for sf1 in slay1:
                            geom = sf1.GetGeometryRef()

                            if geom.Intersects(self._boundsGeom):
                                tile_url = sf1.GetField('URL').strip()
                                self._results.append([tile_url, '{}/{}'.format(surv_id, tile_url.split('/')[-1]), 'raster'])

                        shp1 = slay1 = None

                    for i in ti:
                        ts = os.remove(os.path.join('dc_tile_index/', i))
                    os.removedirs(os.path.join('.', 'dc_tile_index'))
                    os.remove(os.path.join('.', sshpz))
        if len(self._results) == 0: self._status = -1
        gmt1 = layer = None

        tw.opm = 'filtered \033[1m{}\033[m data files from Digital Coast reference vector'.format(len(self._results))
        tw.end(self._status)

## =============================================================================
##
## NOS Fetch
##
## fetch NOS BAG and XYZ sounding data from NOAA
## BAG data is in projected units and MLLW (height)
## XYZ data is CSV in MLLW (Sounding)
##
## =============================================================================

class nos:
    '''Fetch NOS BAG and XYZ sounding data from NOAA'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading NOS fetch module...')
        self._nos_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso_u/xml/' %(nd)
        self._nos_directories = ["B00001-B02000/", "D00001-D02000/", "F00001-F02000/", \
                                 "H00001-H02000/", "H02001-H04000/", "H04001-H06000/", \
                                 "H06001-H08000/", "H08001-H10000/", "H10001-H12000/", \
                                 "H12001-H14000/", "L00001-L02000/", "L02001-L04000/", \
                                 "T00001-T02000/", "W00001-W02000/"]

        self._outdir = os.path.join(os.getcwd(), 'nos')
        self._ref_vector = os.path.join(fetchdata, 'nos.gmt')
        self._local_ref_vector = 'nos.gmt'

        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self._status = 0
        self._surveys = []
        self._results = []

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None: 
            self._bounds = bounds2geom(extent.region)
        else: self._status = -1

        self.fetch_q = queue.Queue()
        pb.opm = 'loaded NOS fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the NOS fetching module.'''

        if self._want_update:
            self._update()
        else:
            self.search_gmt()
            #self._results = list(set(self._results))
            t = np.asarray(self._results)
            p = t != ''
            r = t[p.all(1)]
            
            return(r.tolist())

    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================

    def _parse_nos_xml(self, xml_url, sid):
        '''pare the NOS XML file and extract relavant infos.'''

        xml_doc = fetch_nos_xml(xml_url)
        title = 'Unknown'

        title_ = xml_doc.find('.//gmd:title/gco:CharacterString', namespaces = namespaces)
        if title_ is not None:
            title = title_.text

        obbox = None
        wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
        el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
        sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
        nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
        if wl is not None and el is not None and sl is not None and nl is not None:
            obbox = bounds2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)])

        odt = '0000'        
        dt = xml_doc.find('.//gmd:date/gco:Date', namespaces = namespaces)
        if dt is not None:
            odt = dt.text[:4]

        dts = ['GEODAS_XYZ', 'BAG', 'GRID_BAG']
        xml_dsf = 'None'
        xml_dsu = []

        ## ==============================================
        ## Get all the data URLs
        ## ==============================================

        dfs = xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString',
                              namespaces = namespaces)
        dus = xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL',
                              namespaces = namespaces)
        if dus is not None:
            for i,j in enumerate(dus):
                if dfs is not None:
                    for dt in dts:
                        if dfs[i].text in dts:
                            xml_dsu.append(j.text)
                            xml_dsf = dfs[i].text

            return([obbox, title, sid, odt, xml_url, ','.join(list(set(xml_dsu))), xml_dsf])
        return([None])

    def _scan_directory(self, nosdir):
        '''Scan an NOS directory and parse the XML for each survey.'''

        if self._has_vector:
            gmt1 = ogr.GetDriverByName('GMT').Open(self._local_ref_vector, 0)
            layer = gmt1.GetLayer()
        else: layer = []

        xml_catalog = self._nos_xml_url(nosdir)
        page = fetch_html(xml_catalog)
        rows = page.xpath('//a[contains(@href, ".xml")]/@href')

        ## ==============================================
        ## Parse each survey found in the directory
        ## and append it to the surveys list
        ## ==============================================

        for survey in rows:
            if self.stop():
                break
            sid = survey[:-4]

            if self._has_vector:
                layer.SetAttributeFilter('ID = "{}"'.format(sid))

            if len(layer) == 0:
                xml_url = xml_catalog + survey
                s_entry = self._parse_nos_xml(xml_url, sid)

                if s_entry[0]:
                    self._surveys.append(s_entry)
        gmt1 = layer = None

    def _update(self):
        '''Crawl the NOS database and update the NOS reference vector.'''

        for j in self._nos_directories:
            tw = utils._progress('scanning {}...'.format(j))
            if self.stop():
                break

            if os.path.exists(self._local_ref_vector): 
                self._has_vector = True
            else: self._has_vector = False

            self._scan_directory(j)
            update_ref_vector(self._local_ref_vector, self._surveys, self._has_vector)

            self._surveys = []
            tw.opm = 'scanned {}...'.format(j)
            tw.end(self._status)

    ## ==============================================
    ## Filter for results
    ## ==============================================

    def search_gmt(self):
        '''Search the NOS reference vector and append the results
        to the results list.'''

        tw = utils._progress('filtering NOS reference vector...')
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            if not self.stop():
                geom = feature.GetGeometryRef()
                if geom.Intersects(self._bounds):
                    fldata = feature.GetField('Data').split(',')

                    for i in fldata:
                        self._results.append([i, i.split('/')[-1], feature.GetField('Datatype')])

        gmt1 = layer = None
        tw.opm = 'filtered \033[1m{}\033[m data files from NOS reference vector'.format(len(self._results))
        tw.end(self._status)

## =============================================================================
##
## Chart Fetch - ENC & RNC
##
## Fetch digital charts from NOAA, including ENC and RNC
##
## =============================================================================

class charts():
    '''Fetch digital chart data from NOAA'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading CHARTS fetch module...')
        self.fetch_q = queue.Queue()

        self._enc_data_catalog = 'http://www.charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'http://www.charts.noaa.gov/RNCs/RNCProdCat_19115.xml'

        self._outdir = os.path.join(os.getcwd(), 'charts')
        self._ref_vector = os.path.join(fetchdata, 'charts.gmt')

        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self._dt_xml = { 'ENC':self._enc_data_catalog,
                         'RNC':self._rnc_data_catalog }
        self._checks = self._dt_xml.keys()[0]

        self._status = 0
        self._results = []
        self._chart_feats = []

        self.stop = callback
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)
        else: self._status = -1

        pb.opm = 'loaded CHARTS fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the charts fetching module.'''

        if self._want_update:
            self._update()
        else:
            self.search_gmt()

            return self._results
      
    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================

    def _parse_charts_xml(self, update = True):
        '''parse the charts XML and extract the survey results'''

        if update:
            ds = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
            layer = ds.GetLayer()
        else: layer = []

        namespaces = {'gmd': 'http://www.isotc211.org/2005/gmd', 
                      'gco': 'http://www.isotc211.org/2005/gco',
                      'gml': 'http://www.isotc211.org/2005/gml' }
        
        charts = self.chart_xml.findall('.//{*}has', namespaces = namespaces)

        for chart in charts:
            opoly = []
            titles = chart.findall('.//gmd:title', namespaces = namespaces)
            title = titles[0].find('gco:CharacterString', namespaces = namespaces).text
            id = titles[1].find('gco:CharacterString', namespaces = namespaces).text

            cd = chart.find('.//gmd:date/gco:Date', namespaces = namespaces)
            if cd is not None:
                cd = cd.text

            if update:
                layer.SetAttributeFilter('Name = "{}"'.format(title))

            if len(layer) == 0:
                polygon = chart.find('.//{*}Polygon', namespaces = namespaces)
                if polygon is not None:
                    nodes = polygon.findall('.//{*}pos', namespaces = namespaces)

                    for node in nodes:
                        opoly.append(map(float, node.text.split(' ')))

                    linkage = chart.find('.//{*}linkage/{*}URL', namespaces = namespaces)
                    if linkage is not None:
                        linkage = linkage.text

                    if self._checks == 'RNC': opoly.append(opoly[0])
                    geom = ogr.CreateGeometryFromWkt(create_polygon(opoly))
                    self._chart_feats.append([geom, title, id, cd[:4], self._dt_xml[self._checks], linkage, self._checks])

        ds = layer = None

    def _update(self):
        '''Update or create the reference vector file'''

        for dt in self._dt_xml.keys():
            self._checks = dt

            self.chart_xml = fetch_nos_xml(self._dt_xml[self._checks])
            self._parse_charts_xml(self._has_vector)

            if len(self._chart_feats) > 0:
                update_ref_vector(self._ref_vector, self._chart_feats, self._has_vector)

            if os.path.exists(self._ref_vector): 
                self._has_vector = True
            else: self._has_vector = False

            self._chart_feats = []

    ## ==============================================
    ## Filter for results
    ## ==============================================

    def search_gmt(self):
        '''Search for data in the reference vector file'''

        tw = utils._progress('filtering CHARTS reference vector...')
        ds = ogr.Open(self._ref_vector)
        layer = ds.GetLayer(0)

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                self._results.append([feature1.GetField('Data'), feature1.GetField('Data').split('/')[-1], feature1.GetField('Datatype')])

        ds = layer = None
        tw.opm = 'filtered \033[1m{}\033[m data files from CHARTS reference vector.'.format(len(self._results))
        tw.end(self._status)

## =============================================================================
##
## SRTM Fetch (cgiar)
##
## Fetch srtm tiles from cgiar.
##
## =============================================================================

class srtm_cgiar:
    '''Fetch SRTM data from CGIAR'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading SRMT (CGIAR) fetch module...')
        self._srtm_url = 'http://srtm.csi.cgiar.org'
        self._srtm_dl_url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/'

        self._status = 0
        self._results = []
        self._outdir = os.path.join(os.getcwd(), 'srtm')
        self._ref_vector = os.path.join(fetchdata, 'srtm.gmt')
        if not os.path.exists(self._ref_vector):
            self._status = -1

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self._boundsGeom = None
        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)

        pb.opm = 'loaded SRTM (CGIAR) fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the SRTM fetching module.'''

        self.search_gmt()

        return(self._results)

        if self._want_list:
            self.print_results()
        else: self.fetch_results()

    ## ==============================================
    ## Filter for results
    ## ==============================================

    def search_gmt(self):
        tw = utils._progress('filtering SRTM reference vector...')
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            geom = feature.GetGeometryRef()

            if geom.Intersects(self._boundsGeom):
                geo_env = geom.GetEnvelope()
                srtm_lon = int(math.ceil(abs((-180 - geo_env[1]) / 5)))
                srtm_lat = int(math.ceil(abs((60 - geo_env[3]) / 5)))
                out_srtm = 'srtm_{:02}_{:02}.zip'.format(srtm_lon, srtm_lat)
                self._results.append(['{}{}'.format(self._srtm_dl_url, out_srtm), out_srtm, 'srtm'])

        tw.opm = 'filtered \033[1m{}\033[m data files from SRTM reference vector.'.format(len(self._results))
        tw.end(self._status)

## =============================================================================
##
## The National Map (TNM) - USGS
##
## Fetch elevation data from The National Map
## NED, 3DEP, Etc.
##
## TODO: Break up search regions on large regions 
##
## =============================================================================

class tnm:
    '''Fetch elevation data from The National Map'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading The National Map fetch module...')
        self._tnm_api_url = "http://viewer.nationalmap.gov/tnmaccess/"
        self._tnm_dataset_url = "http://viewer.nationalmap.gov/tnmaccess/api/datasets?"
        self._tnm_product_url = "http://viewer.nationalmap.gov/tnmaccess/api/products?"

        self._status = 0
        self._req = None
        self._results = []
        self._dataset_results = {}

        self._outdir = os.path.join(os.getcwd(), 'tnm')
        self._ref_vector = None

        self._tnm_ds = [1, 2]
        self._tnm_df = ['IMG']

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        if extent is not None: 
            self._req = fetch_req(self._tnm_dataset_url)
            self._datasets = self._req.json()

            sbDTags = []
            for ds in self._tnm_ds:
                sbDTags.append(self._datasets[ds]['sbDatasetTag'])
        
            self.data = { 'datasets':sbDTags,
                          'bbox':extent.bbox }
        
            self._req = fetch_req(self._tnm_product_url, params = self.data)
            if self._req is not None:
                self._dataset_results = self._req.json()
            else: self._status = -1

        pb.opm = 'loaded The National Map fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the TNM (National Map) fetching module.'''

        if self._status == 0:
            for i in self._dataset_results['items']:
                f_url = i['downloadURL']
                self._results.append([f_url, f_url.split('/')[-1], 'tnm'])

        return(self._results)

## =============================================================================
##
## MB Fetch
##
## Fetch Multibeam bathymetric surveys from NOAA
## MBSystem is required to process the resulting data
##
## =============================================================================

class mb:
    '''Fetch multibeam bathymetry from NOAA'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading Multibeam fetch module...')
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._ref_vector = None
        self._outdir = os.path.join(os.getcwd(), 'mb')

        self._status = 0
        self._req = None
        self._results = []
        self._surveys = []
        self._survey_list = []

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None:
            self.data = { 'geometry':extent.bbox }
            self._req = fetch_req(self._mb_search_url, params = self.data)
            if self._req is not None:
                self._survey_list = self._req.content.split('\n')[:-1]
            else: self._status = -1
        else: self._status = -1

        pb.opm = 'loaded Multibeam fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the MB (multibeam) fetching module.'''
        if self._status == 0:
            for r in self._survey_list:
                dst_pfn = r.split(' ')[0]
                dst_fn = r.split(' ')[0].split('/')[-1:][0]
                survey = r.split(' ')[0].split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                self._results.append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])

        return(self._results)
        
    def proc_results(self, local):
        for res in self._survey_list:
            survey = res.split(' ')[0].split('/')[6]
            if survey not in self._surveys:
                self._surveys.append(survey)

            if local: 
                sf = open(survey + '.mb-1', 'a')
                sf.write('.' + res)
                sf.write('\n')
                sf.close()
            else:
                su = open(survey + '.url', 'a')
                data_url = self._mb_data_url + '/'.join(res.split('/')[3:])
                su.write(data_url.split(' ')[0])
                su.write('\n')
                su.close()

## =============================================================================
##
## USACE Fetch
##
## Fetch USACE bathymetric surveys
##
## =============================================================================

class usace:
    '''Fetch USACE bathymetric surveys'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading USACE fetch module...')
        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self._ref_vector = None

        self._status = 0
        self._req = None
        self._results = []
        self._survey_list = []

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None:
            self.data = { 'geometry':extent.bbox,
                          'inSR':4326,
                          'f':'pjson' }
 
            self._req = fetch_req(self._usace_gs_api_url, params = self.data)
            if self._req is not None:
                self._survey_list = self._req.json()
            else: self._status = -1
        else: self._status = -1

        pb.opm = 'loaded USACE fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the USACE fetching module'''

        if self._status == 0:
            for feature in self._survey_list['features']:
                fetch_fn = feature['attributes']['SOURCEDATALOCATION']
                self._results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])

        return(self._results)

## =============================================================================
##
## GMRT Fetch
##
## fetch extracts of the GMRT.
##
## =============================================================================

class gmrt:
    '''Fetch raster data from the GMRT'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading GMRT fetch module...')
        self._gmrt_grid_url = "https://www.gmrt.org/services/GridServer?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')

        self._status = 0
        self._results = []
        self._ref_vector = None

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None: 
            self.data = { 'north':extent.north,
                          'west':extent.west,
                          'south':extent.south,
                          'east':extent.east,
                          'layer':'topo',
                          'resolution':'max',
                          'format':'geotiff' }
        
            self._req = fetch_req(self._gmrt_grid_url, params = self.data, tries = 2)

            gmrt_fn = self._req.headers['content-disposition'].split('=')[1].strip()
            outf = os.path.join(self._outdir,
                                '{}_{}.{}'.format(gmrt_fn.split('.')[0],
                                                  self.region.fn,
                                                  gmrt_fn.split('.')[1]))

            self._results = [self._req.url, outf, 'gmrt']
            
        pb.opm = 'loaded GMRT fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the GMRT fetching module'''
        
        if self._req is not None:
            return([self._results])
        else: return([])            

## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS
##
## =============================================================================

class ngs:
    '''Fetch NGS monuments from NGS'''

    def __init__(self, extent = None, filters = [], callback = None):
        pb = utils._progress('loading NGS Monument fetch module...')
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')
        self._ref_vector = None

        self._status = 0
        self._req = None
        self._results = []

        self.stop = callback
        self._want_proc = True
        self._want_update = False
        self._filters = filters

        self.region = extent
        if extent is not None:
            self.data = { 'maxlon':extent.east,
                          'minlon':extent.west,
                          'maxlat':extent.north,
                          'minlat':extent.south }

            self._req = fetch_req(self._ngs_search_url, params = self.data)

        pb.opm = 'loaded NGS Monument fetch module.'
        pb.end(self._status)

    def run(self):
        '''Run the NGS (monuments) fetching module.'''

        self._results.append([self._req.url, 'ngs_results_{}.txt'.format(self.region.fn), 'ngs'])
        return(self._results)
        
    def proc_results(self):
        try:
            r = self._results.json()

            if len(r) > 0:
                if not os.path.exists(self._outdir):
                    os.makedirs(self._outdir)

                outfile = open(os.path.join(self._outdir,
                                            'ngs_results_{}.csv'.format(self.region.fn)), 'w')
            
                outcsv = csv.writer(outfile)
                outcsv.writerow(r[0].keys())
                for row in r:
                    outcsv.writerow(row.values())

                outfile.close()
        except: self._status = -1

## =============================================================================
##
## Run fetch from command-line
##
## =============================================================================

fetch_infos = { 
    'dc':[lambda x, f, c: dc(x, f, callback = c), 'digital coast'],
    'nos':[lambda x, f, c: nos(x, f, callback = c), 'noaa nos bathymetry'],
    'charts':[lambda x, f, c: charts(x, f, callback = c), 'noaa nautical charts'],
    'srtm':[lambda x, f, c: srtm_cgiar(x, f, callback = c), 'srtm from cgiar'],
    'tnm':[lambda x, f, c: tnm(x, f, callback = c), 'the national map'],
    'mb':[lambda x, f, c: mb(x, f, callback = c), 'noaa multibeam'],
    'gmrt':[lambda x, f, c: gmrt(x, f, callback = c), 'gmrt'],
    'usace':[lambda x, f, c: usace(x, f, callback = c), 'usace bathymetry'],
    'ngs':[lambda x, f, c: ngs(x, f, callback = c), 'ngs monuments'],
}

def fetch_desc(x):
    fd = []
    for key in x: 
        fd.append('{:10}\t\t{}'.format(key, x[key][1]))
    return fd

_usage = '''{} ({}): Fetch geographic elevation data.

usage: {} [ -fhlpRuv [ args ] ] module(s) ...

Modules:
  {}

Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -l, --list-only\tOnly fetch a list of surveys in the given region.
  -f, --filter\t\tSQL style attribute filter; use ? to see available field names.
  -p, --process\t\tProcess fetched data to ASCII XYZ format.

  --update\t\tUpdate the stored list of surveys.
  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Examples:
 % sudo {} --update
 % {} nos charts -R -90.75/-88.1/28.7/31.25 -f "Date > 2000"
 % {} dc -R tiles.shp -p
 % {} dc -R tiles.shp -f "Datatype LIKE 'lidar%'" -l > dc_lidar.urls

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format(os.path.basename(sys.argv[0]),
           _version, 
           os.path.basename(sys.argv[0]), 
           '\n  '.join(fetch_desc(fetch_infos)),
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]))

def main():
    status = 0
    extent = None
    want_list = False
    want_update = False
    want_proc = False
    stop_threads = False

    f = []
    fetch_class = []
    these_regions = []

    ## ==============================================
    ## process Command-Line
    ## ==============================================

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '--region' or arg == '-R':
            extent = str(sys.argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            iregion = str(arg[2:])

        elif arg == '--list-only' or arg == '-l':
            want_list = True

        elif arg == '--filter' or arg == '-f':
            f.append(sys.argv[i + 1])
            i = i + 1

        elif arg == '--process' or arg == '-p':
            want_proc = True

        elif arg == '--update' or arg == '-u':
            want_update = True

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}\
            '.format(os.path.basename(sys.argv[0]), 
                     _version, 
                     utils._license))
            sys.exit(1)

        else: fetch_class.append(sys.argv[i])

        i = i + 1

    if extent is None and want_update is False:
        print(_usage)
        sys.exit(1)

    if len(fetch_class) == 0:
        fetch_class = fetch_infos.keys()

    ## ==============================================
    ## check platform and installed software
    ## ==============================================

    tw = utils._progress('checking system status...')
    platform = sys.platform
    tw.opm = 'platform is {}.'.format(platform)
    tw.end(status)

    if want_proc:
        tw.opm = 'checking for GMT'
        if utils.cmd_exists('gmt'): 
            gmt_vers, status = utils.run_cmd('gmt --version')
        else: status = -1
        tw.end(status)

        tw.opm = 'checking for MBSystem'
        if utils.cmd_exists('mbgrid'): 
            mbs_vers, status = utils.run_cmd('mbgrid -version')
        else: status = -1
        tw.end(status)

        tw.opm = 'checking for LASTools'
        if utils.cmd_exists('las2txt'): 
            status = 0
        else: status = -1
        tw.end(status)

        tw.opm = 'checking for GDAL command-line'
        if utils.cmd_exists('gdal-config'): 
            gdal_vers, status = utils.run_cmd('gdal-config --version')
        else: status = -1
        tw.end(status)

    tw.opm = 'checked system status.'
    tw.end(status)

    ## ==============================================
    ## process input region(s)
    ## ==============================================

    if extent is not None:
        tw = utils._progress('processing region(s)...')
        try: 
            these_regions = [regions.region(extent)]
        except:
            if os.path.exists(extent):
                _poly = ogr.Open(extent)
                _player = _poly.GetLayer(0)
                for pf in _player:
                    _pgeom = pf.GetGeometryRef()
                    these_regions.append(regions.region("/".join(map(str, _pgeom.GetEnvelope()))))

        if len(these_regions) == 0:
            status = -1

        for this_region in these_regions:
            if not this_region._valid: 
                status = -1

        tw.opm = 'processed {} region(s)'.format(len(these_regions))
        tw.end(status)
    else: these_regions = [regions.region('-180/180/0/90')]

    ## ==============================================
    ## fetch some data in each of the input regions
    ## ==============================================

    for rn, this_region in enumerate(these_regions):
        if stop_threads:
            return
        for fc in fetch_class:

            ## ==============================================
            ## Run the Fetch Module
            ## ==============================================

            status = 0
            pb = utils._progress('running fetch module \033[1m{}\033[m on region \033[1m{}\033[m ({}/{})...\
            '.format(fc, this_region.region_string, rn+1, len(these_regions)))
            fl = fetch_infos[fc][0](this_region.buffer(5, percentage = True), f, lambda: stop_threads)
            fl._want_update = want_update
            r = fl.run()

            if len(r) == 0:
                status = -1
            else:
                if want_list:
                    for result in r:
                        print(result[0])
                else:
                    fr = fetch_results(r, this_region, fl._outdir, want_proc, lambda: stop_threads)

                    try:
                        fr.start()
                        while True:
                            time.sleep(1)
                            pb.update()
                            if not fr.is_alive():
                                break
                    except (KeyboardInterrupt, SystemExit): 
                        fr._status = -1
                        stop_threads = True

                    fr.join()
            pb.opm = 'ran fetch module \033[1m{}\033[m on region \033[1m{}\033[m ({}/{})...\
            '.format(fc, this_region.region_string, rn+1, len(these_regions))
            pb.end(status)

if __name__ == '__main__':
    main()

### End
