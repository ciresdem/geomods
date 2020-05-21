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
### Commentary:
##
## current fetch modules: dc, nos, mb, charts, usace, srtm, tnm, gmrt
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
import csv
import json
import threading

import numpy as np
import ogr
import gdal

#import procs

try:
    import Queue as queue
except: import queue as queue

_version = '0.3.2'

## =============================================================================
## functions from waffles:
## =============================================================================

def region_format(region, t = 'gmt'):
    '''format region to string'''
    if t == 'str': return('/'.join([str(x) for x in region]))
    elif t == 'gmt': return('-R' + '/'.join([str(x) for x in region]))
    elif t == 'bbox': return(','.join([str(region[0]), str(region[2]), str(region[1]), str(region[3])]))
    elif t == 'fn':
        if region[3] < 0: ns = 's'
        else: ns = 'n'
        if region[0] > 0: ew = 'e'
        else: ew = 'w'
        return('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(region[3])), abs(int(region[3] * 100) % 100), 
                                                        ew, abs(int(region[0])), abs(int(region[0] * 100) % 100)))

def region_valid_p(region):
    '''return True if region appears valid'''
    if region is not None:
        if region[0] < region[1] and region[2] < region[3]: return(True)
        else: return(False)
    else: return(False)

def region_pct(region, pctv):
    '''return the pctv val of the region'''
    ewp = (region[1] - region[0]) * (pctv * .01)
    nsp = (region[3] - region[2]) * (pctv * .01)
    return((ewp + nsp) / 2)

def region_buffer(region, bv = 0, pct = False):
    '''return the region buffered by bv'''
    if pct: bv = region_pct(region, bv)
    return([region[0] - bv, region[1] + bv, region[2] - bv, region[3] + bv])

def gdal_ogr_regions(src_ds):
    '''return the region(s) of the ogr dataset'''
    these_regions = []
    if os.path.exists(src_ds):
        poly = ogr.Open(src_ds)
        if poly is not None:
            p_layer = poly.GetLayer(0)
            for pf in p_layer:
                pgeom = pf.GetGeometryRef()
                these_regions.append(pgeom.GetEnvelope())
        poly = None
    return(these_regions)

def gdal_fext(src_drv_name):
    '''return the common file extention given a GDAL driver name'''
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER): fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None: return(fexts.split()[0])
        else: return(None)
    except:
        if src_drv_name == 'GTiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        return(fext)

def echo_error_msg2(msg, prefix = 'waffles'):
    '''echo error msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: error, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix = 'waffles'):
    '''echo msg to stderr'''
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.flush()
    sys.stderr.write('{}: {}\n'.format(prefix, msg))

echo_msg = lambda m: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
echo_error_msg = lambda m: echo_error_msg2(m, prefix = os.path.basename(sys.argv[0]))

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
        this_dt = fetch_args[4].lower()
        fetch_args[2] = None

        if not fetch_args[3]():
            if fetch_args[0].split(':')[0] == 'ftp':
                fetch_ftp_file(*tuple(fetch_args))
            else: fetch_file(*tuple(fetch_args))
            
            ## ==============================================
            ## initiate the processing module if desired
            ## ==============================================
            if p is not None:
                if this_dt == 'tnm':
                    proc_mod = 'gdal'
                    proc_opts = [0, True, None]
                elif  this_dt == 'grid_bag':
                    proc_mod = 'gdal'
                    proc_opts = [None, True, 'mllw']
                elif this_dt == 'gmrt' or this_dt == 'raster' or this_dt == 'srtm':
                    proc_mod = 'gdal'
                    proc_opts = [None, None, None]
                elif this_dt == 'geodas_xyz':
                    proc_mod = 'ascii'
                    proc_opts = [',', '2,1,3', 1, True, 'mllw']
                elif this_dt == 'lidar':
                    proc_mod = 'lidar'
                    proc_opts = [None, None]
                elif this_dt == 'enc':
                    proc_mod = 'ogr'
                    proc_opts = ['SOUNDG', 'Depth', 'mllw']
                elif this_dt == 'mb':
                    proc_mod = 'mb'
                    proc_opts = []
                elif this_dt == 'ngs':
                    proc_mod = 'ngs'
                    proc_opts = []
                else:
                    proc_mod = None
                    proc_opts = []
                
                if os.path.exists(fetch_args[1]):
                    p.put([[fetch_args[1], proc_mod, [this_region], fetch_args[3]], proc_opts])
        q.task_done()

def fetch_ftp_file(src_url, dst_fn, params = None, callback = None, datatype = None):
    '''fetch an ftp file via urllib2'''
    import urllib2
    status = 0
    f = None
    halt = callback
    echo_msg('fetching remote ftp file: {}...'.format(os.path.basename(src_url)))

    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    try:
        f = urllib2.urlopen(src_url)
    except: status - 1
    
    if f is not None:
        with open(dst_fn, 'wb') as local_file:
            local_file.write(f.read())
    echo_msg('fetched remote ftp file: {}.'.format(os.path.basename(src_url)))
    return(status)

def fetch_file(src_url, dst_fn, params = None, callback = None, datatype = None):
    '''fetch src_url and save to dst_fn'''
    status = 0
    req = None
    halt = callback
    echo_msg('fetching remote file: {}...'.format(os.path.basename(src_url)))

    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 

    if not os.path.exists(dst_fn):
        try:
            req = requests.get(src_url, stream = True, params = params, headers = r_headers)
        except requests.ConnectionError as e:
            echo_error_msg('Error: {}'.format(e))
            status = -1

        if req is not None:
            with open(dst_fn, 'wb') as local_file:
                for chunk in req.iter_content(chunk_size = 50000):
                    if chunk:
                        if halt(): 
                            status = -1
                            break
                        local_file.write(chunk)
    else: status = -1
    echo_msg('fetched remote file: {}.'.format(os.path.basename(src_url)))
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
    def __init__(self, results, region, out_dir, want_proc = False, callback = lambda: False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.results = results
        self.region = region
        self._outdir = out_dir
        self.want_proc = want_proc
        self.stop_threads = callback
        
    def run(self):
        if self.want_proc:
            proc_q = queue.Queue()
            for _ in range(3):
                p = threading.Thread(target = procs._proc_queue, args = (proc_q,))
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
## Reference Vector
##
## the reference vector location and related functions
## dc, nos, charts and srtm use reference vectors
##
## =============================================================================
this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

def _ogr_create_polygon(coords):
    '''convert coords to Wkt'''

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords:
        ring.AddPoint(coord[1], coord[0])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    poly_wkt = poly.ExportToWkt()
    poly = None
    
    return(poly_wkt)

def bounds2geom(bounds):
    '''convert a bounds [west, east, south, north] to an 
    OGR geometry'''
    b1 = [[bounds[2], bounds[0]],
          [bounds[2], bounds[1]],
          [bounds[3], bounds[1]],
          [bounds[3], bounds[0]],
          [bounds[2], bounds[0]]]
    return(ogr.CreateGeometryFromWkt(_ogr_create_polygon(b1)))

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
        echo_msg('loading Digital Coast fetch module...')
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dav_id = 'https://coast.noaa.gov/dataviewer/#/lidar/search/where:ID='
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'raster2']

        self._ref_vector = os.path.join(fetchdata, 'dc.gmt')
        self._outdir = os.path.join(os.getcwd(), 'dc')

        self._status = 0
        self._surveys = []
        self._results = []

        self._index = False
        self._filters = filters
        
        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self.stop = callback

        self.region = extent
        if extent is not None: 
            self._boundsGeom = bounds2geom(self.region)
        else: self._status = -1
        
    def run(self, datatype = None, index = False, update = False):
        '''Run the Digital Coast fetching module'''

        if index:
            if index.lower() == 'false':
                index = False
        
        self._index = index
        
        if update:
            if update.lower() == 'false':
                update = False

        self._want_update = update
                
        if self._want_update:
            self._update()
            return([])

        if datatype is not None:
            if datatype.lower() == 'none':
                datatype = None
            else: self.datatype = datatype
            
        self.datatype = datatype
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
                                print(out_s)
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

        echo_msg('filtering Digital Coast reference vector...')
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

                if self.datatype is not None:
                    if self.datatype.lower() not in surv_dt:
                        next(layer, None)
                        continue
            
                if self._index:
                    print("%s (%s): %s (%s) - %s" %(feature1.GetField("ID"),feature1.GetField("Datatype"),feature1.GetField("Name"),feature1.GetField("Date"),feature1.GetField("Data")))
                else:
                
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
        echo_msg('filtered \033[1m{}\033[m data files from Digital Coast reference vector'.format(len(self._results)))
        
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
        echo_msg('loading NOS fetch module...')
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
        self._filters = filters
        
        self.region = extent
        self._bounds = None

    def run(self, datatype = None, update = False):
        '''Run the NOS fetching module.'''
        if update:
            if update.lower() == 'false':
                update = False
        
        self._want_update = update
                
        if self._want_update:
            self._update()
            return([])

        if self.region is None:
            return([])
        
        self._bounds = bounds2geom(self.region)

        if datatype is not None:
            if datatype.lower() == 'none':
                datatype = None
            else: self.datatype = datatype
            
        self.datatype = datatype

        self.search_gmt()

        if len(self._results) > 0:
            t = np.asarray(self._results)
            p = t != ''
            r = t[p.all(1)]
            
            return(r.tolist())
        else: return([])

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
            echo_msg('scanning {}...'.format(j))
            if self.stop():
                break

            if os.path.exists(self._local_ref_vector): 
                self._has_vector = True
            else: self._has_vector = False

            self._scan_directory(j)
            update_ref_vector(self._local_ref_vector, self._surveys, self._has_vector)

            self._surveys = []
            echo_msg('scanned {}...'.format(j))

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def search_gmt(self):
        '''Search the NOS reference vector and append the results
        to the results list.'''

        echo_msg('filtering NOS reference vector...')
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            if not self.stop():
                geom = feature.GetGeometryRef()
                if geom.Intersects(self._bounds):
                    fldata = feature.GetField('Data').split(',')
                    if feature.GetField('Datatype').lower() != 'none':
                        if self.datatype is not None:
                            if self.datatype.lower() in feature.GetField('Datatype').lower():
                                for i in fldata:
                                    self._results.append([i, i.split('/')[-1], feature.GetField('Datatype')])
                        else:
                            for i in fldata:
                                self._results.append([i, i.split('/')[-1], feature.GetField('Datatype')])

        gmt1 = layer = None
        tw.end(self._status, 'filtered \033[1m{}\033[m data files from NOS reference vector'.format(len(self._results)))

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
        echo_msg('loading CHARTS fetch module...')
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
        self._filters = filters
        
        self.region = extent
        self._boundsGeom = None

    def run(self, datatype = None, update = False):
        '''Run the charts fetching module.'''
        if update:
            if update.lower() == 'false':
                update = False
        
        self._want_update = update
        
        if self._want_update:
            self._update()
            return([])

        if self.region is None:
            return([])

        if datatype is not None:
            if datatype.lower() == 'none':
                datatype = None
        
        self.dt = datatype
        self._boundsGeom = bounds2geom(self.region)
        self.search_gmt()
        return(self._results)
      
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
        echo_msg('updating reference vector: {}'.format(self._ref_vector))
        
        for dt in self._dt_xml.keys():
            echo_msg('updating {}'.format(dt))
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

        echo_msg('filtering CHARTS reference vector...')
        ds = ogr.Open(self._ref_vector)
        layer = ds.GetLayer(0)

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                if self.dt is not None:
                    if feature1.GetField('Datatype').lower() == self.dt.lower():
                        self._results.append([feature1.GetField('Data'), feature1.GetField('Data').split('/')[-1], feature1.GetField('Datatype')])
                else: self._results.append([feature1.GetField('Data'), feature1.GetField('Data').split('/')[-1], feature1.GetField('Datatype')])

        ds = layer = None
        echo_msg('filtered \033[1m{}\033[m data files from CHARTS reference vector.'.format(len(self._results)))

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
        echo_msg('loading SRMT (CGIAR) fetch module...')
        self._srtm_url = 'http://srtm.csi.cgiar.org'
        self._srtm_dl_url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/'

        self._status = 0
        self._outdir = os.path.join(os.getcwd(), 'srtm')
        self._ref_vector = os.path.join(fetchdata, 'srtm.gmt')
        if not os.path.exists(self._ref_vector):
            self._status = -1

        self._results = []
        self._filters = filters

        self._boundsGeom = None
        self.region = extent

    def run(self):
        '''Run the SRTM fetching module.'''

        if self.region is None: return([])
        self._boundsGeom = bounds2geom(self.region)
        self.search_gmt()
        return(self._results)

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def search_gmt(self):
        echo_msg('filtering SRTM reference vector...')
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

        echo_msg('filtered \033[1m{}\033[m data files from SRTM reference vector.'.format(len(self._results)))

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
        echo_msg('loading The National Map fetch module...')
        self._tnm_api_url = "http://viewer.nationalmap.gov/tnmaccess/"
        self._tnm_dataset_url = "http://viewer.nationalmap.gov/tnmaccess/api/datasets?"
        self._tnm_product_url = "http://viewer.nationalmap.gov/tnmaccess/api/products?"
        self._outdir = os.path.join(os.getcwd(), 'tnm')
        
        self._status = 0
        self._req = None
        self._results = []
        self._dataset_results = {}

        self._tnm_ds = []
        self._tnm_df = []
        
        self._req = None        
        self.region = extent

    def run(self, index = False, ds = 3, sub_ds = None, formats = None, extent = None):
        '''Run the TNM (National Map) fetching module.'''
        if self.region is None:
            return([])
        
        self._req = fetch_req(self._tnm_dataset_url)

        if self._req is not None:
            try:
                self._datasets = self._req.json()
            except:
                echo_error_msg('error, try again')
                self._status = -1
        else: self._status = -1

        if index:
            if index.lower() == 'false':
                index = False

        if self._status == 0:
            if index:
                self.print_datasets()
                sys.exit()

            this_ds = [int(ds)]
            if sub_ds is not None:
                this_ds.append(int(sub_ds))

            self._tnm_ds = [this_ds]

            if formats is not None:
                self._tnm_df = formats.split(',')
            else: self._tnm_df = []

            if extent is not None:
                self._extents = [str(extent)]
            else: self._extents = []
            
            self.filter_datasets()
            
        return(self._results)

    def filter_datasets(self):
        '''Filter TNM datasets.'''
        echo_msg('filtering TNM dataset results...')
        sbDTags = []
        
        for ds in self._tnm_ds:
            dtags = self._datasets[ds[0]]['tags']
            
            if len(ds) > 1:
                if len(dtags) == 0:
                    sbDTags.append( self._datasets[ds[0]]['sbDatasetTag'])
                else:
                    dtag = dtags.keys()[ds[1]]
                    sbDTag = self._datasets[ds[0]]['tags'][dtag]['sbDatasetTag']
                    
                    if len(self._tnm_df) == 0:
                        formats = self._datasets[ds[0]]['tags'][dtag]['formats']
                        self._tnm_df = formats
                    sbDTags.append(sbDTag)
            else:
                if len(self._tnm_df) == 0:
                    all_formats = True
                else: all_formats = False

                if len(dtags) == 0:
                    sbDTags.append( self._datasets[ds[0]]['sbDatasetTag'])
                else:
                    for dtag in dtags.keys():
                        sbDTag = self._datasets[ds[0]]['tags'][dtag]['sbDatasetTag']

                        if len(self._tnm_df) == 0:
                            formats = self._datasets[ds[0]]['tags'][dtag]['formats']

                            for ff in formats:
                                self._tnm_df.append(ff)
                                
                        sbDTags.append(sbDTag)

        self.data = { 'datasets': sbDTags,
                      'bbox': region_format(self.region, 'bbox') }

        if len(self._tnm_df) > 0:
            self.data['prodFormats'] = ','.join(self._tnm_df)

        req = fetch_req(self._tnm_product_url, params = self.data)

        if req is not None:
            try:
                self._dataset_results = req.json()
            except ValueError:
                echo_error_msg("tnm server error, try again")
                
        else: self._status = -1

        if len(self._dataset_results) > 0:
            for item in self._dataset_results['items']:
                if len(self._extents) > 0:
                    for extent in self._extents:
                        if item['extent'] == extent:
                            try:
                                f_url = item['downloadURL']
                                self._results.append([f_url, f_url.split('/')[-1], 'tnm'])
                            except: pass
                else:
                    try:
                        f_url = item['downloadURL']
                        self._results.append([f_url, f_url.split('/')[-1], 'tnm'])
                    except: pass

        echo_msg('filtered \033[1m{}\033[m data files from TNM dataset results.'.format(len(self._results)))

    def print_datasets(self):
        for i,j in enumerate(self._datasets):
            print('%s: %s [ %s ]' %(i, j['title'], ", ".join(j['formats'])))
            for m,n in enumerate(j['tags']):
                print('\t%s: %s [ %s ]' %(m, n, ", ".join(j['tags'][n]['formats'])))

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
        echo_msg('loading Multibeam fetch module...')
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._outdir = os.path.join(os.getcwd(), 'mb')
        self._status = 0
        self._req = None
        self._results = []
        self.region = extent

    def run(self):
        '''Run the MB (multibeam) fetching module.'''
        if self.region is None: return([])
        
        self.data = { 'geometry': region_format(self.region, 'bbox') }
        self._req = fetch_req(self._mb_search_url, params = self.data)
        if self._req is not None:
            survey_list = self._req.content.split('\n')[:-1]
            for r in survey_list:
                dst_pfn = r.split(' ')[0]
                dst_fn = dst_pfn.split('/')[-1:][0]
                survey = dst_pfn.split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                self._results.append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])
        return(self._results)
        
## =============================================================================
##
## USACE Fetch
##
## Fetch USACE bathymetric surveys via eHydro
##
## =============================================================================
class usace:
    '''Fetch USACE bathymetric surveys'''

    def __init__(self, extent = None, filters = [], callback = None):
        echo_msg('loading USACE fetch module...')
        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self._status = 0
        self._req = None
        self._results = []
        self.region = extent

    def run(self, stype = None):
        '''Run the USACE fetching module'''

        if self.region is None:
            return([])
        
        self.data = { 'geometry': region_format(self.region, 'bbox'),
                      'inSR':4326,
                      'outSR':4326,
                      'f':'pjson' }
 
        self._req = fetch_req(self._usace_gs_api_url, params = self.data)

        if self._req is not None:
            survey_list = self._req.json()
            for feature in survey_list['features']:
                fetch_fn = feature['attributes']['SOURCEDATALOCATION']
                if stype is not None:
                    if feature['attributes']['SURVEYTYPE'].lower() == stype.lower():
                        self._results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                else: self._results.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
        return(self._results)

## =============================================================================
##
## GMRT Fetch
##
## fetch extracts of the GMRT.
## https://www.gmrt.org/index.php
##
## =============================================================================
class gmrt:
    '''Fetch raster data from the GMRT'''
    def __init__(self, extent = None, filters = [], callback = None):
        echo_msg('loading GMRT fetch module...')
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')
        self._status = 0
        self._results = []
        self.region = extent

    def run(self, res = 'max', fmt = 'geotiff'):
        '''Run the GMRT fetching module'''
        if self.region is None: return([])
        self.data = { 'north':self.region[3],
                      'west':self.region[0],
                      'south':self.region[2],
                      'east':self.region[1],
                      'mformat':'json',
                      'resolution':res,
                      'format':fmt}
        
        self._req = fetch_req(self._gmrt_grid_urls_url, params = self.data, tries = 10, timeout = 2)
        if self._req is not None:
            gmrt_urls = self._req.json()
            for url in gmrt_urls:
                opts = {}
                for url_opt in url.split('?')[1].split('&'):
                    opt_kp = url_opt.split('=')
                    opts[opt_kp[0]] = opt_kp[1]

                url_region = [float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])]
                outf = 'gmrt_{}_{}.{}'.format(opts['layer'], region_format(url_region, 'fn'), gdal_fext(opts['format']))
                self._results.append([url, outf, 'gmrt'])
        return(self._results)
    
## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS
##
## =============================================================================
class ngs:
    '''Fetch NGS monuments from NOAA'''
    def __init__(self, extent = None, filters = [], callback = None):
        echo_msg('loading NGS Monument fetch module...')
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')
        self._status = 0
        self._req = None
        self._results = []
        self.region = extent

    def run(self, csv = False):
        '''Run the NGS (monuments) fetching module.'''
        if self.region is None: return([])
        self.data = { 'maxlon':self.region[0],
                      'minlon':self.region[1],
                      'maxlat':self.region[3],
                      'minlat':self.region[2] }

        self._req = fetch_req(self._ngs_search_url, params = self.data)
        if self._req is not None: self._results.append([self._req.url, 'ngs_results_{}.json'.format(region_format(self.region, 'fn')), 'ngs'])
        return(self._results)
        
## =============================================================================
##
## Run fetches from command-line
##
## =============================================================================
fetch_infos = { 
    'dc':[lambda x, f, c: dc(x, f, c), '''NOAA Digital Coast access
    \t\t\t< dc:datatype=None:index=False:update=False >
    \t\t\t:datatype=[lidar/raster] - Only fetch lidar or raster data.
    \t\t\t:index=[True/False] - True to display indexed results.
    \t\t\t:update=[True/False] - True to update stored reference vector.'''],
    'nos':[lambda x, f, c: nos(x, f, c), '''NOAA NOS bathymetry surveys and data (hydro & BAG)
    \t\t\t< nos:datatype=None:update=False >
    \t\t\t:datatype=[bag/xyz] - Only fetch BAG or Hydro-XYZ data.
    \t\t\t:update=[True/False] - True to update stored reference vector.'''],
    'charts':[lambda x, f, c: charts(x, f, c), '''NOAA Nautical CHARTS (RNC & ENC)
    \t\t\t< charts:datatype=None:update=False >
    \t\t\t:dataype=[ENC/RNC] - Only fetch either ENC or RNC data.
    \t\t\t:update=[True/False] - True to update stored reference vector.'''],
    'srtm':[lambda x, f, c: srtm_cgiar(x, f, c), '''SRTM elevation data from CGIAR'''],
    'tnm':[lambda x, f, c: tnm(x, f, c), '''The National Map (TNM) from USGS
    \t\t\t< tnm:ds=1:sub_ds=None:formats=None:index=False >
    \t\t\t:index=[True/False] - True to display an index of available datasets.
    \t\t\t:ds=[dataset index value (0-15)] - see :index=True for dataset index values.
    \t\t\t:sub_ds=[sub-dataset index value (0-x)] - see :index=True for sub-dataset index values.
    \t\t\t:formats=[data-set format] - see :index=True for dataset format options.'''],
    'mb':[lambda x, f, c: mb(x, f, c), '''NOAA MULTIBEAM survey data'''],
    'gmrt':[lambda x, f, c: gmrt(x, f, c), '''The Global Multi-Reosolution Topography Data Synthesis (GMRT) 
    \t\t\t< gmrt:res=max:fmt=geotiff >
    \t\t\t:res=[an Integer power of 2 zoom level (<=1024)]
    \t\t\t:fmt=[netcdf/geotiff/esriascii/coards]'''],
    'usace':[lambda x, f, c: usace(x, f, c), '''USACE bathymetry surveys via eHydro'''],
    'ngs':[lambda x, f, c: ngs(x, f, c), '''NOAA NGS monuments''']
}

def fetch_desc(x):
    fd = []
    for key in x: 
        fd.append('{:18}{}'.format(key, x[key][-1]))
    return fd

_usage = '''{} ({}): Fetch geographic elevation data.

usage: {} [ -hlpRv [ args ] ] module[:parameter=value]* ...

General Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -f, --filter\t\tSQL style filters for certain datasets.
\t\t\tFields to filter include: 'Name', 'Date' and 'Datatype'
  -p, --process\t\tProcess fetched data to ASCII XYZ format. <beta>
\t\t\tRequires external data processing programs (GMT/GDAL/MBSYSTEM/LASTOOLS/VDATUM)

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Modules and their options:
  {}

Examples:
 % {} -R -90.75/-88.1/28.7/31.25 nos -f "Date > 2000"
 % {} -R region.shp -p dc nos:datatype=bag charts:datatype=enc
 % {} -R region.shp dc:datatype=lidar -l > dc_lidar.urls
 % {} -R -89.75/-89.5/30.25/30.5 tnm:ds=4:formats=IMG gmrt:res=max:fmt=geotiff

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(os.path.basename(sys.argv[0]),
           _version, 
           os.path.basename(sys.argv[0]),
           '\n  '.join(fetch_desc(fetch_infos)),
           os.path.basename(sys.argv[0]),
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
    mod_opts = {}

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
            extent = str(arg[2:])
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
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format( _version))
            sys.exit(0)
        else: 
            opts = arg.split(':')
            if opts[0] in fetch_infos.keys():
                mod_opts[opts[0]] = list(opts[1:])
            else: echo_error_msg('invalid module name `{}`'.format(opts[0]))
        i = i + 1

    if len(mod_opts) == 0:
        echo_error_msg('you must select a fetch module')
        sys.stderr.write(_usage)
        sys.exit(1)
        
    for key in mod_opts.keys():
        mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]
        
    ## ==============================================
    ## process input region(s)
    ## ==============================================
    if extent is not None:
        try:
            these_regions = [[float(x) for x in extent.split('/')]]
        except ValueError: these_regions = gdal_ogr_regions(extent)
        if len(these_regions) == 0: status = -1
        for this_region in these_regions:
            if not region_valid_p(this_region): status = -1
        echo_msg('loaded {} region(s)'.format(len(these_regions)))
    else: these_regions = [[-180, 180, 0, 90]]

    ## ==============================================
    ## fetch some data in each of the input regions
    ## ==============================================
    for rn, this_region in enumerate(these_regions):
        if stop_threads: return

        ## ==============================================
        ## Run the Fetch Module(s)
        ## ==============================================
        for fetch_mod in mod_opts.keys():
            status = 0
            args = tuple(mod_opts[fetch_mod])
            echo_msg('running fetch module {} on region {} ({}/{})...\
            '.format(fetch_mod, region_format(this_region, 'str'), rn+1, len(these_regions)))
            fl = fetch_infos[fetch_mod][0](region_buffer(this_region, 5, pct = True), f, lambda: stop_threads)

            args_d = {}
            for arg in args:
                p_arg = arg.split('=')
                if len(p_arg) > 1:
                    args_d[p_arg[0]] = p_arg[1]

            try:
                r = fl.run(**args_d)
            except ValueError as e:
                echo_error_msg('something went wrong, {}'.format(e))
                sys.exit(-1)
            echo_msg('found {} data files.'.format(len(r)))
            
            if want_list:
                for result in r:
                    print(result[0])
            else:
                fr = fetch_results(r, this_region, fl._outdir, want_proc, lambda: stop_threads)
                try:
                    fr.start()
                    while True:
                        time.sleep(5)
                        sys.stderr.write('.')
                        if not fr.is_alive():
                            break
                except (KeyboardInterrupt, SystemExit):
                    echo_error_msg('user breakage')
                    fl._status = -1
                    stop_threads = True
                fr.join()
            echo_msg('ran fetch module {} on region {} ({}/{})...\
            '.format(fetch_mod, region_format(this_region, 'str'), rn+1, len(these_regions)))

if __name__ == '__main__':  main()
### End
