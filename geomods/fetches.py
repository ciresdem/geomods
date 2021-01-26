### fetches.py
##
## Copyright (c) 2012 - 2021 CIRES Coastal DEM Team
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
## current fetch modules: dc, nos, mb, charts, usace, srtm_cgiar, srtm_plus, tnm, gmrt, emodnet, mar_grav
##
## Possible BAG errors with GDAL >= 3
##
### Code:

import os
import sys
import time
import math
import requests
import ftplib
import lxml.html as lh
import lxml.etree
import zipfile
import csv
import json
import threading
import copy
import numpy as np
import ogr
import gdal
try:
    import Queue as queue
except: import queue as queue
from xml.dom import minidom

from geomods import utils
from geomods import regions
from geomods import vdatumfun
from geomods import gdalfun
from geomods import xyzfun

_version = '0.5.0'

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## =============================================================================
r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(_version) }
namespaces = {
    'gmd': 'http://www.isotc211.org/2005/gmd', 
    'gmi': 'http://www.isotc211.org/2005/gmi', 
    'gco': 'http://www.isotc211.org/2005/gco',
    'gml': 'http://www.isotc211.org/2005/gml',
}
thredds_namespaces = {
    'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
}

def fetch_queue(q, p):
    '''fetch queue `q` of fetch results'''
    while True:
        fetch_args = q.get()
        this_region = fetch_args[2]
        this_dt = fetch_args[4].lower()
        fetch_args[2] = None
        #utils.echo_msg(fetch_args[-1])
        if not fetch_args[3]():
            if p is None:
                if fetch_args[0].split(':')[0] == 'ftp':
                    fetch_ftp_file(*tuple(fetch_args))
                else: fetch_file(*tuple(fetch_args))
            else:
                if not os.path.exists(os.path.dirname(fetch_args[1])):
                    try:
                        os.makedirs(os.path.dirname(fetch_args[1]))
                    except: pass
                #utils.echo_msg(fetch_args[1])
                o_x_fn = '.'.join(fetch_args[1].split('.')[:-1]) + '.xyz'
                if not os.path.exists(o_x_fn):
                    with open(o_x_fn, 'w') as out_xyz:
                        p._dump_xyz([fetch_args[0], fetch_args[1], fetch_args[-1]], dst_port = out_xyz)
                    
        q.task_done()

def fetch_ftp_file(src_url, dst_fn, params = None, callback = None, datatype = None):
    '''fetch an ftp file via urllib2'''
    import urllib2
    status = 0
    f = None
    halt = callback
    utils.echo_msg('fetching remote ftp file: {}...'.format(os.path.basename(src_url)))
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
    #utils.echo_msg('fetched remote ftp file: {}.'.format(os.path.basename(src_url)))
    return(status)

def fetch_file(src_url, dst_fn, params = None, callback = lambda: False, datatype = None, overwrite = False, verbose = False):
    '''fetch src_url and save to dst_fn'''
    status = 0
    req = None
    halt = callback

    if verbose: utils.echo_msg('fetching remote file: {}...'.format(os.path.basename(src_url)))

    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 

    if not os.path.exists(dst_fn) or overwrite:
        try:
            with requests.get(src_url, stream = True, params = params, headers = r_headers, timeout=(45,320)) as req:
                req_h = req.headers
                if req.status_code == 200:
                    curr_chunk = 0
                    with open(dst_fn, 'wb') as local_file:
                        for chunk in req.iter_content(chunk_size = 8196):
                            if chunk:
                                if halt(): 
                                    status = -1
                                    break
                                local_file.write(chunk)
                else: utils.echo_error_msg('server returned: {}'.format(req.status_code))
                
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
    else:
        if os.path.exists(dst_fn): return(status)
        status = -1
    if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0: status = -1
    if verbose: utils.echo_msg('fetched remote file: {}.'.format(os.path.basename(dst_fn)))
    return(status)

def fetch_req(src_url, params = None, tries = 5, timeout = 2, read_timeout = 10):
    '''fetch src_url and return the requests object'''
    if tries <= 0:
        utils.echo_error_msg('max-tries exhausted')
        return(None)
    try:
        return(requests.get(src_url, stream = True, params = params, timeout = (timeout,read_timeout), headers = r_headers))
    except: return(fetch_req(src_url, params = params, tries = tries - 1, timeout = timeout + 1, read_timeout = read_timeout + 10))

def fetch_nos_xml(src_url, verbose = False):
    '''fetch src_url and return it as an XML object'''
    results = lxml.etree.fromstring('<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8'))
    try:
        req = fetch_req(src_url, timeout = .25)
        results = lxml.etree.fromstring(req.text.encode('utf-8'))
    except:
        if verbose: utils.echo_error_msg('could not access {}'.format(src_url))
    return(results)
        
def fetch_html(src_url):
    '''fetch src_url and return it as an HTML object'''
    req = fetch_req(src_url, timeout = 2)
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
    '''fetch results gathered from a fetch module.
    results is a list of URLs with data type'''
    def __init__(self, results, region, out_dir, proc = None, vdc = None, callback = lambda: False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.results = results
        self.region = region
        self._outdir = out_dir
        self.stop_threads = callback
        self.proc = proc
        
    def run(self):
        for _ in range(3):
            t = threading.Thread(target = fetch_queue, args = (self.fetch_q, self.proc))
            t.daemon = True
            t.start()
            
        fq = [[row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2]] for row in self.results]
        [self.fetch_q.put([row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2]]) for row in self.results]
        self.fetch_q.join()
    
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
    '''convert coords to Wkt

    returns WKT polygon'''
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[1], coord[0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def bounds2geom(bounds):
    '''convert a bounds [west, east, south, north] to an 
    OGR geometry
    
    returns OGR geometry'''
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
        if ds is not None: layer = ds.GetLayer()
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
        [layer.SetFeature(f) for f in layer]
        [addf_ref_vector(layer, i) for i in surveys]
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
class dc_ftp:
    def __init__(self):
        self._dc_ftp_url = 'ftp.coast.noaa.gov'
        self.ftp = ftplib.FTP(self._dc_ftp_url)
        self.ftp.login()
        self.ftp.cwd("/pub/DigitalCoast")
        self.parent = self.ftp.pwd()

    def _home(self):
        self.ftp.cwd(self.parent)

    def _list(self):
        self.ftp.retrlines("LIST")

    def fetch_file(self, filename, full=True):
        if full:
            dirname = "." + self.ftp.pwd()
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        else: dirname = "."
            
        with open(dirname + "/" + os.path.basename(filename), "wb") as local_file:
            self.ftp.retrbinary('RETR {}'.format(filename), local_file.write)

    def fetch_res(self, fn):
        response = []
        self.ftp.retrbinary('RETR {}'.format(fn), lambda d: response.append(d))
        return(''.join(response))
    
class dc:
    '''Fetch elevation data from the Digital Coast'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading Digital Coast fetch module...')

        ## ==============================================
        ## URLs and directories
        ## ==============================================
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_ftp_url = 'ftp://ftp.coast.noaa.gov/pub/DigitalCoast'
        self._dc_ftp_url_full = "ftp://ftp.coast.noaa.gov"
        self._dc_dav_id = 'https://coast.noaa.gov/dataviewer/#/lidar/search/where:ID='
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'lidar3_z', 'lidar4_z', 'raster1', 'raster2', 'raster5']
        #self._dc_dirs = ['raster1', 'raster2', 'raster5']
        self._dc_subdirs = ['geoid12a', 'geoid12b', 'geoid18', 'elevation']

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._local_ref_vector = 'dc.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'dc.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'nos.gmt'
        
        self._outdir = os.path.join(os.getcwd(), 'dc')

        ## ==============================================
        ## _surveys holds the survey info for reference vector updates
        ## _results holds the survey info for reference vector filtering
        ## ==============================================
        self._surveys = []
        self._results = []
        
        self._index = False
        self._filters = filters
        #self._has_vector = False

        ## ==============================================
        ## region of insterest.
        ## ==============================================
        self.region = extent
        if self.region is not None: 
            self._boundsGeom = bounds2geom(self.region)
        else: self._status = -1

        self._status = 0
        self.stop = callback
        self._verbose = False

        self.dcftp = dc_ftp()
        
    def run(self, datatype = None, index = False, update = False):
        '''Run the Digital Coast fetching module'''
        self._index = index
        self._want_update = update
        if self._want_update:
            #self._ref_vector = self._local_ref_vector
            #self._has_vector = False
            self._update_from_ftp()
            
        utils.echo_msg('using reference vector: {}'.format(self._ref_vector))
            
        self.datatype = datatype
        self._filter_reference_vector()
        return(self._results)
        
    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================

    def _fetch_res(self, url):
        import urllib2
        
        response = urllib2.urlopen(url)
        results = response.read()
        response.close()
        return(results)

    def _fetch_xml(self, url):
        try:
            return lxml.etree.fromstring(self._fetch_res(url).encode('utf-8'))
        except: return(None)
        
    def _update_directory(self, s_dtype):

        # if self._has_vector:
        #     gmt2 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        #     layer = gmt2.GetLayer()
        #     layerDef = layer.GetLayerDefn()
        # else: layer = []
        layer = []
        
        pdir = self.dcftp.ftp.pwd()
        if s_dtype == 'lidar': self.dcftp.ftp.cwd("data")
        data_dir = self.dcftp.ftp.pwd()
        data_list = self.dcftp.ftp.nlst()
        utils.echo_msg(data_dir)
        utils.echo_msg_inline('scanning {} surveys in {} [    ].'.format(len(data_list), data_dir))
        
        for it, dataset in enumerate(data_list):
            try: did = int(dataset.split("_")[-1])
            except: did = None
            perc = int((float(it)/len(data_list)) * 100)
            utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(data_list), data_dir, perc, did))
            if did is not None:
                #if self._has_vector: layer.SetAttributeFilter('ID = "{}"'.format(did))
                #if self._has_vector: layer.SetAttributeFilter('ID = "%s"' %(did))
                if len(layer) == 0:
                    self.dcftp.ftp.cwd(dataset)
                    dc_files = self.dcftp.ftp.nlst()
                    for dc_file in dc_files:
                        if "metadata.xml" in dc_file or 'met.xml' in dc_file:
                            xml_url = self._dc_ftp_url + "/" + dataset + "/" + dc_file
                            ds_url = self._dc_ftp_url + "/" + dataset + "/"
                            xml_res = self.dcftp.fetch_res(dc_file)
                            if len(xml_res) == 0:
                                utils.echo_error_msg('xml_doc is none')
                                break
                            
                            xml_doc = lxml.etree.fromstring(xml_res)

                            wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
                            el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
                            sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
                            nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
                            try:
                                this_region = [float(x) for x in [wl.text, el.text, sl.text, nl.text]]
                            except: this_region = None

                            if regions.region_valid_p(this_region):
                                obbox = bounds2geom(this_region)
                                try:
                                    odate = int(xml_doc.find('.//gco:Date', namespaces = namespaces).text[:4])
                                except: odate = 1900

                                try:
                                    otitle = xml_doc.find('.//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString', namespaces = namespaces).text
                                except: otitle = 'UNKNOWN DATASET'
                                
                                ds_urls = xml_doc.findall('.//gmd:CI_OnlineResource', namespaces=namespaces)
                                for i in ds_urls:
                                    i_name = i.find('.//gmd:name/gco:CharacterString',namespaces=namespaces)
                                    if i_name is not None:
                                        if i_name.text == 'Bulk Download':
                                            ds_url = i.find('.//gmd:linkage/gmd:URL', namespaces=namespaces).text
                                        
                                utils.echo_msg(ds_url)
                                ## ==============================================
                                ## Append survey to surveys list
                                ## ==============================================
                                out_s = [obbox, otitle, did, odate, xml_url, ds_url, s_dtype]
                                self._surveys.append(out_s)

            self.dcftp.ftp.cwd(data_dir)
        self.dcftp.ftp.cwd(pdir)
        utils.echo_msg_inline('scanning {} surveys in {} [ OK ].\n'.format(len(data_list), data_dir))
        gmt2 = layer = None

    def _update_from_ftp(self):
        file_list = self.dcftp.ftp.nlst()
        for this_dir in self._dc_dirs:
            utils.echo_msg('scanning {}'.format(this_dir))
            s_dtype = 'lidar' if 'lidar' in this_dir else 'raster' if 'raster' in this_dir else None
            self.dcftp.ftp.cwd(this_dir)
            sub_dir = self.dcftp.ftp.pwd()
            for this_sub_dir in self._dc_subdirs:
                try:
                    self.dcftp.ftp.cwd(this_sub_dir)
                    #utils.echo_msg('scanning {}'.format(self.dcftp.ftp.pwd()))
                    self._update_directory(s_dtype)
                    #except error_perm: pass
                except Exception as e:
                    pass
                self.dcftp.ftp.cwd(sub_dir)
            self.dcftp._home()
                
            ## ==============================================
            ## Add matching surveys to reference vector
            ## ==============================================
            if len(self._surveys) > 0: update_ref_vector(self._local_vector, self._surveys, False)
            self._has_vector = True if os.path.exists(self._local_vector) else False
            self._surveys = []
            gmt2 = layer = None
    
    def _update(self):
        '''Update the DC reference vector after scanning
        the relevant metadata from Digital Coast.'''
        utils.echo_msg('updating Digital Coast fetch module')
        for ld in self._dc_dirs:
            utils.echo_msg('scanning {}'.format(ld))
            #if self._has_vector:
            #    gmt2 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 1)
            #    layer = gmt2.GetLayer()
            #    layerDef = layer.GetLayerDefn()
            #    #print([layerDef.GetFieldDefn(i).GetName() for i in range(layerDef.GetFieldCount())])
            #else: layer = []
            layer = []
            page = fetch_html(self._dc_htdata_url + ld)
            tables = page.xpath('//table')
            tr = tables[0].xpath('.//tr')
            if len(tr) <= 0: break
            cols = []
            [cols.append(i.text_content()) for i in tr[0]]
            #for i in tr[0]: cols.append(i.text_content())
            #sys.stderr.write('fetches: scanning {} surveys in {} [    ].'.format(len(tr), ld))
            utils.echo_msg_inline('scanning {} surveys in {} [    ].'.format(len(tr), ld))
            
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
                    perc = int((float(i)/len(tr)) * 100)
                    #sys.stderr.write('\x1b[2K\rfetches: scanning {} surveys in {} [{:3}%] - {}.'.format(len(tr), ld, perc, dc['ID #']))
                    utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(tr), ld, perc, dc['ID #']))
                    #print(dc['ID #'])
                    if self._has_vector: layer.SetAttributeFilter('ID = "%s"' %(dc['ID #']))
                    #if self._has_vector:
                        #print('ID = {}'.format(dc['ID #']))
                    #    print(layer.SetAttributeFilter("ID = '{}'".format(dc['ID #'])))
                        
                    #for feature in layer:
                     #   print('dcID - ' + feature.GetField('ID') + '-' + dc['ID #'])
                    #print(len(layer))
                    if len(layer) == 0:
                        if 'Metadata' in dc.keys():
                            xml_doc = fetch_nos_xml(dc['Metadata'], verbose = False)
                            wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
                            el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
                            sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
                            nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
                            try:
                                this_region = [float(x) for x in [wl.text, el.text, sl.text, nl.text]]
                            except: this_region = None
                            if regions.region_valid_p(this_region):
                                obbox = bounds2geom(this_region)
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
                                #print(out_s)
                                self._surveys.append(out_s)

            ## ==============================================
            ## Add matching surveys to reference vector
            ## ==============================================
            #print(self._surveys)
            utils.echo_msg_inline('scanning {} surveys in {} [ OK ].\n'.format(len(tr), ld))
            gmt2 = layer = None
        if len(self._surveys) > 0: update_ref_vector(self._ref_vector, self._surveys, self._has_vector)
        #self._has_vector = True if os.path.exists(self._ref_vector) else False
        #print(self._has_vector)
        #self._surveys = []

    ## ==============================================
    ## Filter reference vector for results
    ## ==============================================
    def _filter_reference_vector(self):
        '''Search for data in the reference vector file'''
        utils.echo_msg('filtering Digital Coast reference vector...')
        gmt1 = ogr.Open(self._ref_vector)
        layer = gmt1.GetLayer(0)
        [layer.SetAttributeFilter('{}'.format(filt)) for filt in self._filters]

        ## ==============================================
        ## Find surveys that fit all filters and within
        ## the region of interest.
        ## ==============================================
        for feature1 in layer:
            if self.stop() or self._status != 0: break
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
                    utils.echo_msg('{} ({}): {} ({}) - {}\
                    '.format(feature1.GetField("ID"),feature1.GetField("Datatype"),
                             feature1.GetField("Name"),feature1.GetField("Date"),
                             feature1.GetField("Data")))
                else:
                    suh = fetch_html(surv_url)
                    if suh is None: 
                        self._status = -1
                        break
                    if 'lidar' in surv_dt:
                        ## ==============================================
                        ## Lidar data has a minmax.csv file to get extent
                        ## for each survey file.
                        ## ==============================================
                        scsv = suh.xpath('//a[contains(@href, ".csv")]/@href')[0]
                        dc_csv = fetch_csv(surv_url + "/" + scsv)
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
                        fetch_file(surv_url + sshpz, os.path.join('.', sshpz), callback = self.stop)
                        zip_ref = zipfile.ZipFile(sshpz)
                        zip_ref.extractall('dc_tile_index')
                        zip_ref.close()
                        if os.path.exists('dc_tile_index'): ti = os.listdir('dc_tile_index')

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
                        for i in ti: ts = os.remove(os.path.join('dc_tile_index/', i))
                        os.removedirs(os.path.join('.', 'dc_tile_index'))
                        os.remove(os.path.join('.', sshpz))
                        
        if len(self._results) == 0: self._status = -1
        gmt1 = layer = None
        utils.echo_msg('filtered \033[1m{}\033[m data files from Digital Coast reference vector'.format(len(self._results)))

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry, datatype = None):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        if src_ext == 'laz' or src_ext == 'las': dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img': dt = 'raster'
        else: dt = None
        if dt == 'lidar':
            if fetch_file(entry[0], src_dc, callback = lambda: False, verbose = self._verbose) == 0:
                xyz_dat = utils.yield_cmd('las2txt -verbose -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                '.format(regions.region_format(self.region, 'te'), '2 29', src_dc), verbose = False)
                xyzc = copy.deepcopy(xyzfun._xyz_config)
                xyzc['name'] = src_dc
                for xyz in xyzfun.xyz_parse(xyz_dat, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        elif dt == 'raster':
            try:
                src_ds = gdal.Open(entry[0])
            except Exception as e:
                fetch_file(entry[0], src_dc, callback = lambda: False, verbose = self._verbose)
                try:
                    src_ds = gdal.Open(src_dc)
                except Exception as e:
                    utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                    src_ds = None
            except Exception as e:
                utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                src_ds = None
            
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.region_warp(self.region, s_warp = 4326, t_warp = gdalfun.gdal_getEPSG(src_ds)))
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = 4326, verbose = self._verbose):
                    yield(xyz)
            src_ds = None
        utils.remove_glob(src_dc)

    def _dump_xyz(self, entry, datatype = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, datatype):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    def _yield_results_to_xyz(self, datatype = None):
        if len(self._results) == 0: self.run(datatype)
        for entry in self._results:
            for xyz in self._yield_xyz(entry, datatype):
                yield(xyz)
                
    def _dump_results_to_xyz(self, datatype = None, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz(datatype):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
## =============================================================================
##
## NOS Fetch
##
## fetch NOS BAG and XYZ sounding data from NOAA
## BAG data is in projected units and MLLW (height)
## XYZ data is CSV in MLLW (Sounding)
##
## https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html
##
## =============================================================================
class nos:
    '''Fetch NOS BAG and XYZ sounding data from NOAA'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading NOS fetch module...')

        ## ==============================================
        ## NOS urls and directories
        ## ==============================================
        self._nos_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso_u/xml/' %(nd)
        self._nos_directories = [
            "B00001-B02000/", "D00001-D02000/", "F00001-F02000/", \
            "H00001-H02000/", "H02001-H04000/", "H04001-H06000/", \
            "H06001-H08000/", "H08001-H10000/", "H10001-H12000/", \
            "H12001-H14000/", "L00001-L02000/", "L02001-L04000/", \
            "T00001-T02000/", "W00001-W02000/" \
        ]

        ## ==============================================
        ## NOS data comes as xyz or bag, sometimes gzipped
        ## ==============================================
        self._nos_fmts = ['.xyz.gz', '.bag.gz', '.bag']

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._local_ref_vector = 'nos.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'nos.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'nos.gmt'
                
        self._outdir = os.path.join(os.getcwd(), 'nos')

        ## ==============================================
        ## misc. variables
        ## ==============================================
        self._status = 0
        self._surveys = []
        self._results = []
        self.stop = callback
        self._want_proc = True
        self._filters = filters
        self.region = extent
        self._bounds = None
        self._datalists_code = 401
        self._verbose = False

    def run(self, datatype = None, update = False):
        '''Run the NOS fetching module.'''

        ## ==============================================
        ## update/generate the reference vector
        ## ==============================================
        if update:
            if str(update).lower() == 'false':
                update = False
                
        self._want_update = update                
        if self._want_update:
            self._update()
        
        ## ==============================================
        ## filter the reference vector and return a
        ## list of survey results
        ## ==============================================
        utils.echo_msg('using reference vector: {}'.format(self._ref_vector))
        if self.region is None: return([])
        self._bounds = bounds2geom(self.region)
        if datatype is not None:
            if datatype.lower() == 'none':
                datatype = None
            else: self.datatype = datatype
        self.datatype = datatype
        self._filter_reference_vector()
        if len(self._results) > 0:
            t = np.asarray(self._results)
            p = t != ''
            r = t[p.all(1)]
            return(r.tolist())
        else: return([])

    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================
    def _parse_nos_xml(self, xml_url, sid, nosdir):
        '''pare the NOS XML file and extract relavant infos.'''
        xml_doc = fetch_nos_xml(xml_url)
        title_ = xml_doc.find('.//gmd:title/gco:CharacterString', namespaces = namespaces)
        title = 'Unknown' if title_ is None else title_.text

        wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
        el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
        sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
        nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
        if wl is not None and el is not None and sl is not None and nl is not None:
            obbox = bounds2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)])
        else: obbox = None
        dt = xml_doc.find('.//gmd:date/gco:Date', namespaces = namespaces)
        odt = dt.text[:4] if dt is not None else '0000'
        dts = ['GEODAS_XYZ', 'BAG', 'GRID_BAG']
        xml_dsf = 'None'
        xml_dsu = []
        
        ## ==============================================
        ## Get all the data URLs
        ## ==============================================        
        dfs = xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces = namespaces)
        dus = xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces = namespaces)
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
        ## ==============================================
        ## load the reference vector if it exists and
        ## scan the remote nos directories
        ## ==============================================
        #if self._has_vector:
        #    gmt1 = ogr.GetDriverByName('GMT').Open(self._local_ref_vector, 0)
        #    layer = gmt1.GetLayer()
        #else: layer = []
        layer = []
        
        xml_catalog = self._nos_xml_url(nosdir)
        page = fetch_html(xml_catalog)
        rows = page.xpath('//a[contains(@href, ".xml")]/@href')
        utils.echo_msg_inline('scanning {} surveys in {} [    ].'.format(len(rows), nosdir))
        
        ## ==============================================
        ## Parse each survey found in the directory
        ## and append it to the surveys list
        ## ==============================================
        for i, survey in enumerate(rows):
            if self.stop(): break
            sid = survey[:-4]
            perc = int((float(i)/len(rows)) * 100)
            utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(rows), nosdir, perc, sid))

            #if self._has_vector: layer.SetAttributeFilter('ID = "{}"'.format(sid))
            if len(layer) == 0:
                xml_url = xml_catalog + survey
                s_entry = self._parse_nos_xml(xml_url, sid, nosdir)
                if s_entry[0]: self._surveys.append(s_entry)
        gmt1 = layer = None
        utils.echo_msg_inline('scanning {} surveys in {} [ OK ].\n'.format(len(rows), nosdir))

    def _update(self):
        '''Crawl the NOS database and update/generate the NOS reference vector.'''
        for j in self._nos_directories:
            if self.stop(): break
            #self._has_vector = True if os.path.exists(self._local_ref_vector) else self._has_vector
            self._scan_directory(j)
        update_ref_vector(self._local_ref_vector, self._surveys, False)
        #self._surveys = []

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def _filter_reference_vector(self):
        '''Search the NOS reference vector and append the results
        to the results list.'''
        utils.echo_msg('filtering NOS reference vector...')
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()
        for filt in self._filters: layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            if not self.stop():
                geom = feature.GetGeometryRef()
                if geom.Intersects(self._bounds):
                    fldata = feature.GetField('Data').split(',')
                    if feature.GetField('Datatype').lower() != 'none':
                        if self.datatype is not None:
                            if self.datatype.lower() in feature.GetField('Datatype').lower():
                                [self._results.append([i, i.split('/')[-1], feature.GetField('Datatype')]) for i in fldata]
                        else:
                            [self._results.append([i, i.split('/')[-1], feature.GetField('Datatype')]) for i in fldata]
        gmt1 = layer = None
        utils.echo_msg('filtered \033[1m{}\033[m data files from NOS reference vector'.format(len(self._results)))

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry, datatype = None, vdc = None, xyzc = None):
        if vdc is None: vdc = copy.deepcopy(vdatumfun._vd_config)
        if xyzc is None: xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_nos = os.path.basename(entry[1])
        dt = None
        if fetch_file(entry[0], src_nos, callback = lambda: False, verbose = self._verbose) == 0:
            src_ext = src_nos.split('.')
            if len(src_ext) > 2:
                if src_ext[-2] == 'bag': dt = 'grid_bag'
                elif src_ext[-2] == 'xyz': dt = 'geodas_xyz'
                else: dt = None
            elif len(src_ext) == 2:
                if src_ext[-1] == 'bag': dt = 'grid_bag'
                elif src_ext[-1] == 'xyz': dt = 'geodas_xyz'
                else: dt = None
            else: dt = None
            if dt == 'geodas_xyz':
                nos_f, nos_zips = utils.procs_unzip(src_nos, ['xyz', 'dat'])
                vdc['ivert'] = 'mllw:m:sounding'
                vdc['overt'] = 'lmsl:m:height'
                vdc['delim'] = 'comma'
                vdc['xyzl'] = '2,1,3'
                vdc['skip'] = '1'
                vdc['region'] = '5'
                #out, status = vdatumfun.run_vdatum(nos_f, vdc)
                #nos_f_r = os.path.join('result', os.path.basename(nos_f))
                nos_f_r = nos_f

                xyzc['delim'] = ','
                xyzc['skip'] = 1
                xyzc['xpos'] = 2
                xyzc['ypos'] = 1
                xyzc['zpos'] = 3
                xyzc['z-scale'] = -1
                xyzc['name'] = nos_f_r
                if os.path.exists(nos_f_r):
                    with open(nos_f_r, 'r') as in_n:
                        for xyz in xyzfun.xyz_parse(in_n, xyz_c = xyzc, region = self.region, verbose = self._verbose):
                            yield(xyz)
                    
            elif dt == 'grid_bag':
                src_bag, nos_zips = utils.procs_unzip(src_nos, ['gdal', 'tif', 'img', 'bag', 'asc'])
                nos_f = '{}.tmp'.format(os.path.basename(src_bag).split('.')[0])
                vdc['ivert'] = 'mllw:m:height'
                vdc['overt'] = 'lmsl:m:height'
                vdc['region'] = '5'
                vdc['delim'] = 'space'
                vdc['skip'] = '0'
                vdc['xyzl'] = '0,1,2'
                xyzc['delim'] = ' '
                xyzc['skip'] = 0
                xyzc['xpos'] = 0
                xyzc['ypos'] = 1
                xyzc['zpos'] = 2
                xyzc['z-scale'] = 1
                xyzc['name'] = src_bag
                #, gdalfun.gdal_region(src_ds, warp = 4326):
                src_ds = gdal.Open(src_bag)
                if src_ds is not None:
                    srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.region_warp(self.region, s_warp = 4326, t_warp = gdalfun.gdal_getEPSG(src_ds)))
                    with open(nos_f, 'w') as cx:
                        for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = 4326):
                            xyzfun.xyz_line(xyz, cx)
                    src_ds = None
                    if os.stat(nos_f).st_size != 0:
                        #out, status = vdatumfun.run_vdatum(nos_f, vdc)
                        #src_r_bag = os.path.join('result', os.path.basename(nos_f))
                        src_r_bag = nos_f
                        if os.path.exists(src_r_bag):
                            with open(src_r_bag, 'r') as in_b:
                                for xyz in xyzfun.xyz_parse(in_b, xyz_c = xyzc, verbose = self._verbose):
                                    yield(xyz)
                utils.remove_glob(src_bag)
            utils.remove_glob(nos_f)
        utils.remove_glob(src_nos)
        #vdatumfun.vdatum_clean_result()
        
    def _dump_xyz(self, entry, dst_port = sys.stdout):
         for xyz in self._yield_xyz(entry):
             xyzfun.xyz_line(xyz, dst_port, self._verbose)

    def _dump_xyz_entry(self, entry, dst_port = sys.stdout):
             
        src_nos = os.path.basename(entry[1])
        dt = None
        if fetch_file(entry[0], src_nos, callback = lambda: False, verbose = self._verbose) == 0:
            src_ext = src_nos.split('.')
            if len(src_ext) > 2:
                if src_ext[-2] == 'bag': dt = 'grid_bag'
                elif src_ext[-2] == 'xyz': dt = 'geodas_xyz'
                else: dt = None
            elif len(src_ext) == 2:
                if src_ext[-1] == 'bag': dt = 'grid_bag'
                elif src_ext[-1] == 'xyz': dt = 'geodas_xyz'
                else: dt = None
            else: dt = None

            if dt == 'grid_bag':
                gdalfun.gdal2xyz_chunks(src_nos, 2000, None, 4269, 'mllw,navd88', 'bags.datalist', True)
            elif dt == 'geodas_xyz':
                xyzc = copy.deepcopy(xyzfun._xyz_config)
                xyzc['delim'] = ','
                xyzc['skip'] = 1
                xyzc['xpos'] = 2
                xyzc['ypos'] = 1
                xyzc['zpos'] = 3
                gdalfun.xyz_transform(src_nos, xyzc, None, 'mllw:m:sounding,navd88:m:height', self.region, 'hydronos.datalist', True)
            utils.remove_glob(src_nos)
            
    def _yield_results_to_xyz(self, datatype = None):
        self.run(datatype)
        vdc = copy.deepcopy(vdatumfun._vd_config)
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        if vdc['jar'] is None: vdc['jar'] = vdatumfun.vdatum_locate_jar()[0]
        
        for entry in self._results:
            for xyz in self._yield_xyz(entry, datatype, vdc, xyzc):
                yield(xyz)

    def _dump_results_to_xyz(self, datatype = None, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz(datatype):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                
class cudem:
    def __init__(self, extent = None, filters = [], callback = None):

        self._cudem_catalog = "https://www.ngdc.noaa.gov/thredds/demCatalog.xml"
        self._ngdc_url = "https://www.ngdc.noaa.gov"

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._local_ref_vector = 'cudem.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'cudem.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'cudem.gmt'

        self._outdir = os.path.join(os.getcwd(), 'cudem')

        self._surveys = []
        self._results = []
        self._wcs_results = []
        
        self._status = 0
        self._filters = filters
        self.region = extent
        self._boundsGeom = None
        self._verbose = False
        self.stop = callback

    def run(self, datatype = None, update = False):
        '''Run the cudem fetching module.'''
        self._want_update = update
        if self._want_update:
            self._update()
            
        utils.echo_msg('using reference vector: {}'.format(self._ref_vector))

        if self.region is None: return([])
        self.dt = datatype
        self._boundsGeom = bounds2geom(self.region)
        self._filter_reference_vector()
        return(self._results)

    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================

    def _parse_ds_catalog(self, ds_xml, ds_url, want_update = False):
        utils.echo_msg('scanning {}'.format(ds_url))
        #if want_update:
        #    gmt1 = ogr.GetDriverByName('GMT').Open(fetchdata + "cudem.gmt", 0)
        #    layer = gmt1.GetLayer()
        #else:
        layer = []
        ds = ds_xml.findall('.//th:dataset', namespaces = thredds_namespaces)
        ds_services = ds_xml.findall('.//th:service', namespaces = thredds_namespaces)

        for node in ds:
            ds_name = node.attrib['name']
            ds_id = node.attrib['ID']
            sub_catalogRefs = node.findall('.//th:catalogRef', namespaces = thredds_namespaces)
            if len(sub_catalogRefs) > 0:
                self._parse_catalog(node, ds_url)
                break
            #if self.want_update:
            #    layer.SetAttributeFilter("Name = '%s'" %(str(ds_name)))
            if len(layer) == 0:
                try:
                    ds_path = node.attrib['urlPath']
                except: continue
                if ds_path:
                    iso_url = False
                    wcs_url = False
                    http_url = False
                    for service in ds_services:
                        service_name = service.attrib['name']
                        if service_name == 'iso':
                            iso_url = self._ngdc_url + service.attrib['base'] + ds_path
                        if service_name == 'wcs':
                            wcs_url = self._ngdc_url + service.attrib['base'] + ds_path
                        if service_name == 'http':
                            http_url = self._ngdc_url + service.attrib['base'] + ds_path
                    if iso_url and http_url and wcs_url:
                        iso_xml = fetch_nos_xml(iso_url)
                        #iso_xml = self._fetch_xml(iso_url)

                        wl = iso_xml.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
                        el = iso_xml.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
                        sl = iso_xml.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
                        nl = iso_xml.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
                        if wl is not None and el is not None and sl is not None and nl is not None:
                            obbox = bounds2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)])
                        else: obbox = None
                        dt = iso_xml.find('.//gmd:date/gco:Date', namespaces = namespaces)
                        odt = dt.text[:4] if dt is not None else '0000'

                        zv = iso_xml.find('.//gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString', namespaces = namespaces)
                        if zv is not None:
                            zvar = zv.text
                        else: zvar = 'z'

                        survey = [obbox, ds_name, ds_id, odt, iso_url, wcs_url, zvar]
                        #print(survey)
                        self._surveys.append(survey)
        utils.echo_msg('found {} dems in {}'.format(len(self._surveys), ds_url))
    
    def _parse_catalog(self, catalog_xml, catalog_url, want_update = False):

        catalogRefs = catalog_xml.findall('.//th:catalogRef', namespaces = thredds_namespaces)
        for catalog in catalogRefs:
            cat_href = catalog.attrib['{http://www.w3.org/1999/xlink}href']
            if cat_href[0] == "/":
                cat_url = self._ngdc_url + cat_href
            else: cat_url = os.path.dirname(catalog_url) + "/" + cat_href
            xmldoc = fetch_nos_xml(cat_url)
            self._parse_ds_catalog(xmldoc, cat_url, want_update)
    
    def _update(self):
        cudem_cat_xml = fetch_nos_xml(self._cudem_catalog)
        self._parse_catalog(cudem_cat_xml, self._cudem_catalog, True)
        update_ref_vector(self._local_ref_vector, self._surveys, False)

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def _filter_reference_vector(self):
        '''Search for data in the reference vector file'''

        utils.echo_msg('filtering CUDEM reference vector...')
        ds = ogr.Open(self._ref_vector)
        layer = ds.GetLayer(0)

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                wcs_url_service = "{}?request=GetCoverage&version=1.0.0&service=WCS&coverage={}&bbox={}&format=NetCDF3"\
                                               .format(feature1.GetField("Data"), feature1.GetField("Datatype"), regions.region_format(self.region, 'bbox'))
                if self.dt is not None:
                    if feature1.GetField('Datatype').lower() == self.dt.lower():
                        self._results.append([wcs_url_service, feature1.GetField('Data').split('/')[-1], feature1.GetField('Datatype')])
                else: self._results.append([wcs_url_service, feature1.GetField('Data').split('/')[-1], feature1.GetField('Datatype')])
        ds = layer = None
        utils.echo_msg('filtered \033[1m{}\033[m data files from CUDEM reference vector.'.format(len(self._results)))

## =============================================================================
##
## HRDEM Fetch - Canada High Resolution DEM dataset
##
## Fetch Canadian HRDEM data.
## https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6
##
## =============================================================================
class hrdem():
    '''Fetch HRDEM data from Canada'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading HRDEM fetch module...')

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._hrdem_footprints_url = 'ftp://ftp.maps.canada.ca/pub/elevation/dem_mne/highresolution_hauteresolution/Datasets_Footprints.zip'
        self._local_ref_vector = 'hrdem.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'hrdem.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'hrdem.gmt'

        self._outdir = os.path.join(os.getcwd(), 'hrdem')

        self._results = []
        
        self._filters = filters
        self.region = extent
        self._boundsGeom = None
        self._verbose = False
        self.stop = callback
        
    def run(self, update = False):
        '''Run the hrdem fetching module.'''

        if self.region is None: return([])
        self._boundsGeom = bounds2geom(self.region)
        if update:
            utils.echo_msg('updating')
            self._update()
            #self._ref_vector = self._local_ref_vector
            
        utils.echo_msg('using reference vector: {}'.format(self._ref_vector))
        self._filter_reference_vector()
        return(self._results)

    def _update(self):
        v_zip = os.path.basename(self._hrdem_footprints_url)
        fetch_ftp_file(self._hrdem_footprints_url, v_zip)
        v_files = utils.unzip(v_zip)
        v_shp = [s for s  in v_files if '.shp' in s]
        utils.run_cmd('ogr2ogr {} {} -f GMT'.format(self._local_ref_vector, v_shp[0], verbose = True))
        utils._clean_zips(v_files)
        utils.remove_glob(v_zip)
    
    ## ==============================================
    ## Filter for results
    ## ==============================================
    def _filter_reference_vector(self):
        '''Search for data in the reference vector file'''

        utils.echo_msg('filtering HRDEM reference vector...')
        ds = ogr.Open(self._ref_vector)
        layer = ds.GetLayer(0)

        for filt in self._filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                ftp_url = feature1.GetField('Ftp_dtm').replace('http', 'ftp')
                self._results.append([ftp_url, ftp_url.split('/')[-1], 'hrdem'])

        ds = layer = None
        utils.echo_msg('filtered \033[1m{}\033[m data files from HRDEM reference vector.'.format(len(self._results)))

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        try:
            src_ds = gdal.Open(entry[0])
        except Exception as e:
            fetch_file(entry[0], src_dc, callback = lambda: False, verbose = True)
            try:
                src_ds = gdal.Open(src_dc)
            except Exception as e:
                utils.echo_error_msg('could not read hrdem raster file: {}, {}'.format(entry[0], e))
                src_ds = None
        except Exception as e:
            utils.echo_error_msg('could not read hrdem raster file: {}, {}'.format(entry[0], e))
            src_ds = None

        if src_ds is not None:
            srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.region_warp(self.region, s_warp = 4326, t_warp = gdalfun.gdal_getEPSG(src_ds)))
            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = 4326, verbose = True):
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_dc)

    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    def _yield_results_to_xyz(self, datatype = None):
        if len(self._results) == 0: self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)
                
    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
        
        
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
        utils.echo_msg('loading CHARTS fetch module...')
        self._enc_data_catalog = 'http://www.charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'http://www.charts.noaa.gov/RNCs/RNCProdCat_19115.xml'

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._local_ref_vector = 'charts.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'charts.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'charts.gmt'

        self._outdir = os.path.join(os.getcwd(), 'charts')
        
        self._dt_xml = { 'ENC':self._enc_data_catalog,
                         'RNC':self._rnc_data_catalog }
        self._checks = self._dt_xml.keys()[0]
        
        self._results = []
        self._chart_feats = []
        
        self._filters = filters
        self.region = extent
        self._boundsGeom = None
        self._verbose = False
        self._status = 0
        self.stop = callback

    def run(self, datatype = None, update = False):
        '''Run the charts fetching module.'''
        self._want_update = update
        if self._want_update:
            self._update()
            
        utils.echo_msg('using reference vector: {}'.format(self._ref_vector))
        
        if self.region is None: return([])
        self.dt = datatype
        self._boundsGeom = bounds2geom(self.region)
        self._filter_reference_vector()
        return(self._results)
      
    ## ==============================================
    ## Reference Vector Generation
    ## ==============================================
    def _parse_charts_xml(self, update=True):
        '''parse the charts XML and extract the survey results'''
        # if update:
        #     ds = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        #     layer = ds.GetLayer()
        # else: layer = []
        layer = []
        namespaces = {
            'gmd': 'http://www.isotc211.org/2005/gmd', 
            'gco': 'http://www.isotc211.org/2005/gco',
            'gml': 'http://www.isotc211.org/2005/gml',
        }
        charts = self.chart_xml.findall('.//{*}has', namespaces = namespaces)

        for chart in charts:
            opoly = []
            titles = chart.findall('.//gmd:title', namespaces = namespaces)
            title = titles[0].find('gco:CharacterString', namespaces = namespaces).text
            id = titles[1].find('gco:CharacterString', namespaces = namespaces).text
            cd = chart.find('.//gmd:date/gco:Date', namespaces = namespaces)
            if cd is not None: cd = cd.text
            #if update: layer.SetAttributeFilter('Name = "{}"'.format(title))
            if len(layer) == 0:
                polygon = chart.find('.//{*}Polygon', namespaces = namespaces)
                if polygon is not None:
                    nodes = polygon.findall('.//{*}pos', namespaces = namespaces)
                    for node in nodes:  opoly.append(map(float, node.text.split(' ')))
                    linkage = chart.find('.//{*}linkage/{*}URL', namespaces = namespaces)
                    if linkage is not None: linkage = linkage.text
                    if self._checks == 'RNC': opoly.append(opoly[0])
                    geom = ogr.CreateGeometryFromWkt(_ogr_create_polygon(opoly))
                    self._chart_feats.append([geom, title, id, cd[:4], self._dt_xml[self._checks], linkage, self._checks])
        ds = layer = None

    def _update(self):
        '''Update or create the reference vector file'''
        utils.echo_msg('updating reference vector: {}'.format(self._local_ref_vector))
        for dt in self._dt_xml.keys():
            utils.echo_msg('updating {}'.format(dt))
            self._checks = dt
            self.chart_xml = fetch_nos_xml(self._dt_xml[self._checks])
            self._parse_charts_xml(self._has_vector)
        if len(self._chart_feats) > 0: update_ref_vector(self._local_ref_vector, self._chart_feats, False)
        #self._has_vector = True if os.path.exists(self._local_ref_vector) else False
        #self._chart_feats = []

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def _filter_reference_vector(self):
        '''Search for data in the reference vector file'''

        utils.echo_msg('filtering CHARTS reference vector...')
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
        utils.echo_msg('filtered \033[1m{}\033[m data files from CHARTS reference vector.'.format(len(self._results)))

    ## ==============================================
    ## Process results to xyz
    ## ==============================================
    def _yield_xyz(self, entry, vdc = None, xyzc = None):
        if vdc is None: vdc = vdatumfun._vd_config
        if xyzc is None: xyzc = xyzfun._xyz_config
        xyzc['z-scale'] = -1
        src_zip = os.path.basename(entry[1])
        
        if fetch_file(entry[0], src_zip, callback = lambda: False) == 0:
            if entry[-1].lower() == 'enc':
                src_ch, src_zips = utils.procs_unzip(src_zip, ['.000'])
                dst_xyz = src_ch.split('.')[0] + '.xyz'
                ds_ogr = ogr.Open(src_ch)
                layer_s = ds_ogr.GetLayerByName('SOUNDG')
                if layer_s is not None:
                    with open(dst_xyz, 'w') as o_xyz:
                        for f in layer_s:
                            g = json.loads(f.GetGeometryRef().ExportToJson())
                            for xyz in g['coordinates']:
                                xyzfun.xyz_line([float(x) for x in xyz], o_xyz, self._verbose)

                ds_ogr = layer_s = None

                vdc['ivert'] = 'mllw:m:sounding'
                vdc['overt'] = 'lmsl:m:height'
                vdc['delim'] = 'space'
                vdc['xyzl'] = '0,1,2'
                vdc['skip'] = '0'
                vdc['region'] = '5'
                #out, status = vdatumfun.run_vdatum(dst_xyz, vdc)
                #ch_f_r = os.path.join('result', os.path.basename(dst_xyz))
                ch_f_r = dst_xyz

                if os.path.exists(ch_f_r):
                    with open(ch_f_r, 'r') as in_c:
                        for xyz in xyzfun.xyz_parse(in_c, xyz_c = xyzc, verbose = self._verbose):
                            yield(xyz)

                utils.remove_glob(src_ch)
                utils.remove_glob(dst_xyz)
                utils._clean_zips(src_zips)
                vdatumfun.vdatum_clean_result()
        utils.remove_glob(src_zip)
        
    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
        
    def _yield_results_to_xyz(self, datatype = None):
        if len(self._results) == 0: self.run(datatype)
        vdc = copy.deepcopy(vdatumfun._vd_config)            
        if vdc['jar'] is None: vdc['jar'] = vdatumfun.vdatum_locate_jar()[0]
        for entry in self._results:
            for xyz in self._yield_xyz(entry, datatype):
                yield(xyz)
                
    def _dump_results_to_xyz(self, datatype = None, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz(datatype):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)

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
        utils.echo_msg('loading SRTM (CGIAR) fetch module...')
        self._srtm_url = 'http://srtm.csi.cgiar.org'
        self._srtm_dl_url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/'

        ## ==============================================
        ## Reference vector
        ## ==============================================
        self._local_ref_vector = 'srtm.gmt'
        if not os.path.exists(self._local_ref_vector):
            self._ref_vector = os.path.join(fetchdata, 'srtm.gmt')
        else: self._ref_vector = self._local_ref_vector

        self._has_vector = True if os.path.exists(self._ref_vector) else False
        if not self._has_vector: self._ref_vector = 'srtm.gmt'

        self._outdir = os.path.join(os.getcwd(), 'srtm_cgiar')
        
        self._results = []
        self._filters = filters
        self._boundsGeom = None
        self.region = extent
        self._verbose = False

    def run(self):
        '''Run the SRTM fetching module.'''
        if self.region is None: return([])
        self._boundsGeom = bounds2geom(self.region)
        self._filter_reference_vector()
        return(self._results)

    ## ==============================================
    ## Filter for results
    ## ==============================================
    def _filter_reference_vector(self):
        utils.echo_msg('filtering SRTM reference vector...')
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()
        for filt in self._filters: layer.SetAttributeFilter('{}'.format(filt))
        for feature in layer:
            geom = feature.GetGeometryRef()
            if geom.Intersects(self._boundsGeom):
                geo_env = geom.GetEnvelope()
                srtm_lon = int(math.ceil(abs((-180 - geo_env[1]) / 5)))
                srtm_lat = int(math.ceil(abs((60 - geo_env[3]) / 5)))
                out_srtm = 'srtm_{:02}_{:02}.zip'.format(srtm_lon, srtm_lat)
                self._results.append(['{}{}'.format(self._srtm_dl_url, out_srtm), out_srtm, 'srtm'])
        utils.echo_msg('filtered \033[1m{}\033[m data files from SRTM reference vector.'.format(len(self._results)))

    ## ==============================================
    ## Process results to xyz
    ## ==============================================
    def _yield_xyz(self, entry):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = lambda: False, verbose = self._verbose) == 0:
            src_srtm, src_zips = utils.procs_unzip(os.path.basename(entry[1]), ['tif', 'img', 'gdal', 'asc', 'bag'])
            try:
                src_ds = gdal.Open(src_srtm)
            except:
                utils.echo_error_msg('could not read srtm data: {}'.format(src_srtm))
                src_ds = None

            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
                    yield(xyz)
                src_ds = None

            utils.remove_glob(src_srtm)
            utils._clean_zips(src_zips)
        utils.remove_glob(entry[1])
    
    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                    
    def _yield_results_to_xyz(self):
        if len(self._results) == 0: self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)            

    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
## =============================================================================
##
## mar_grav - Sattelite Altimetry Topography from Scripps
##
## https://topex.ucsd.edu/WWW_html/mar_grav.html
## ftp://topex.ucsd.edu/pub/global_grav_1min/
## https://topex.ucsd.edu/marine_grav/explore_grav.html
## https://topex.ucsd.edu/marine_grav/white_paper.pdf
##
## =============================================================================
class mar_grav:
    '''Fetch mar_grav sattelite altimetry topography'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading mar_grav fetch module...')
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'
        self._outdir = os.path.join(os.getcwd(), 'mar_grav')
        self._results = []
        self.region = extent
        self._datalists_code = 414
        self._verbose = False

    def run(self):
        '''Run the mar_grav fetching module.'''
        if self.region is None: return([])

        self.data = {
            'north':self.region[3],
            'west':self.region[0],
            'south':self.region[2],
            'east':self.region[1],
            'mag':1,
        }

        self._req = fetch_req(self._mar_grav_url, params = self.data, tries = 10, timeout = 2)
        if self._req is not None:
            url = self._req.url
            outf = 'mar_grav_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            self._results.append([url, outf, 'mar_grav'])
        return(self._results)
        
    ## ==============================================
    ## Process results to xyz
    ## ==============================================
    def _yield_xyz(self, entry, xyzc = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = lambda: False, verbose = self._verbose) == 0:
            if xyzc is None: xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            xyzc['x-off'] = -360
            xyzc['verbose'] = True
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
    
    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                    
    def _yield_results_to_xyz(self):
        if len(self._results) == 0: self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)            

    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
## =============================================================================
##
## SRTM Plus
##
## Fetch srtm+ data
## https://topex.ucsd.edu/WWW_html/srtm15_plus.html
## http://topex.ucsd.edu/sandwell/publications/180_Tozer_SRTM15+.pdf
##
## =============================================================================
class srtm_plus:
    '''Fetch SRTM+ data'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading SRTM+ fetch module...')
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'
        self._outdir = os.path.join(os.getcwd(), 'srtm_plus')
        self._results = []
        self.region = extent
        self._datalists_code = 413
        self._verbose = False

    def run(self):
        '''Run the SRTM fetching module.'''
        if self.region is None: return([])

        self.data = {
            'north':self.region[3],
            'west':self.region[0],
            'south':self.region[2],
            'east':self.region[1],
        }

        self._req = fetch_req(self._srtm_url, params = self.data, tries = 10, timeout = 2)
        if self._req is not None:
            url = self._req.url
            outf = 'srtm_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            self._results.append([url, outf, 'srtm'])
        return(self._results)
        
    ## ==============================================
    ## Process results to xyz
    ## ==============================================
    def _yield_xyz(self, entry, xyzc = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = lambda: False, verbose = self._verbose) == 0:
            if xyzc is None: xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
    
    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                    
    def _yield_results_to_xyz(self):
        if len(self._results) == 0: self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)            

    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                        
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
        self._datalists_code = 405
        self._verbose = False

    def run(self, index = False, ds = 3, sub_ds = None, formats = None, extent = None):
        '''Run the TNM (National Map) fetching module.'''
        utils.echo_msg('loading The National Map fetch module...')
        if self.region is None: return([])
        self._req = fetch_req(self._tnm_dataset_url)
        if self._req is not None:
            try:
                self._datasets = self._req.json()
            except Exception as e:
                utils.echo_error_msg('try again, {}'.format(e))
                self._status = -1
        else: self._status = -1

        if self._status == 0:
            if index:
                self.print_dataset_index()
                sys.exit()

            this_ds = [int(ds)]
            if sub_ds is not None: this_ds.append(int(sub_ds))
            self._tnm_ds = [this_ds]
            self._tnm_df = [] if formats is None else formats.split(',')
            self._extents = [] if extent is None else [str(extent)]            
            self.filter_datasets()
        return(self._results)

    def filter_datasets(self):
        '''Filter TNM datasets.'''
        utils.echo_msg('filtering TNM dataset results...')
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
                all_formats = False if len(self._tnm_df) == 0 else True
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

        self.data = {
            'datasets': sbDTags,
            'bbox': regions.region_format(self.region, 'bbox'),
        }
        if len(self._tnm_df) > 0: self.data['prodFormats'] = ','.join(self._tnm_df)
        req = fetch_req(self._tnm_product_url, params = self.data)

        if req is not None:
            try:
                self._dataset_results = req.json()
            except ValueError:
                utils.echo_error_msg('tnm server error, try again')
            except Exception as e:
                utils.echo_error_msg('error, {}'.format(e))
                
        else: self._status = -1

        if len(self._dataset_results) > 0:
            for item in self._dataset_results['items']:
                #print(item)
                if len(self._extents) > 0:
                    for extent in self._extents:
                        if item['extent'] == extent:
                            try:
                                f_url = item['downloadURL']
                                #self._results.append([f_url, f_url.split('/')[-1], 'tnm'])
                                if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                                    tnm_ds = 'ned'
                                elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                                    tnm_ds = 'lidar'
                                else: tnm_ds = 'tnm'
                                self._results.append([f_url, os.path.join(item['datasets'][0].replace('/', '').replace(' ', '_'), f_url.split('/')[-1]), tnm_ds])
                            except: pass
                else:
                    try:
                        f_url = item['downloadURL']
                        if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                            tnm_ds = 'ned'
                        elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                            tnm_ds = 'lidar'
                        else: tnm_ds = 'tnm'
                        self._results.append([f_url, os.path.join(item['datasets'][0].replace('/', '').replace(' ', '_'), f_url.split('/')[-1]), tnm_ds])
                    except: pass
        utils.echo_msg('filtered \033[1m{}\033[m data files from TNM dataset results.'.format(len(self._results)))

    def print_dataset_index(self):
        for i,j in enumerate(self._datasets):
            print('%s: %s [ %s ]' %(i, j['title'], ", ".join(j['formats'])))
            for m,n in enumerate(j['tags']):
                print('\t%s: %s [ %s ]' %(m, n, ", ".join(j['tags'][n]['formats'])))

    def _yield_xyz(self, entry):
        '''yield the xyz data from the tnm fetch module'''
        
        if fetch_file(entry[0], entry[1], callback = lambda: False, verbose = self._verbose) == 0:
            datatype = entry[-1]
            if datatype == 'ned':
                src_tnm, src_zips = utils.procs_unzip(entry[1], ['tif', 'img', 'gdal', 'asc', 'bag'])
                try:
                    src_ds = gdal.Open(src_tnm)
                except:
                    utils.echo_error_msg('could not read tnm data: {}'.format(src_tnm))
                    src_ds = None

                if src_ds is not None:
                    srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                    for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
                        if xyz[2] != 0:
                            yield(xyz)
                    src_ds = None

                utils.remove_glob(src_tnm)
                utils._clean_zips(src_zips)
        utils.remove_glob(entry[1])

    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
    
    def _yield_results_to_xyz(self, index = False, ds = 3, sub_ds = None, formats = None, extent = None):
        if len(self._results) == 0: self.run(index, ds, sub_ds, formats, extent)
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)

    def _dump_results_to_xyz(self, index = False, ds = 3, sub_ds = None, formats = None, extent = None, dst_port = sys.stdout):
        if len(self._results) == 0: self.run(index, ds, sub_ds, formats, extent)
        for xyz in self._yield_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
                
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
        utils.echo_msg('loading Multibeam fetch module...')
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._outdir = os.path.join(os.getcwd(), 'mb')
        self._req = None
        self._results = []
        self.region = extent
        self._datalists_code = 406
        self._verbose = False
        self._stop = callback

    def run(self):
        '''Run the MB (multibeam) fetching module.'''
        if self.region is None: return([])        
        self.data = { 'geometry': regions.region_format(self.region, 'bbox') }
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

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry, vdc = None, xyzc = None):
        if vdc is None: vdc = vdatumfun._vd_config
        if xyzc is None: xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_mb = os.path.basename(entry[1])

        if fetch_file(entry[0], src_mb, callback = self._stop, verbose = self._verbose) == 0:
            src_xyz = os.path.basename(src_mb).split('.')[0] + '.xyz'
            out, status = utils.run_cmd('mblist -MX20 -OXYZ -I{}  > {}'.format(src_mb, src_xyz), verbose = False)
            vdc['ivert'] = 'lmsl:m:height'
            vdc['overt'] = 'navd88:m:height'
            vdc['delim'] = 'tab'
            vdc['xyzl'] = '0,1,2'
            vdc['skip'] = '0'
            xyzc['delim'] = '\t'
            xyzc['z-scale'] = 1
            #out, status = vdatumfun.run_vdatum(src_xyz, vdc)
            #mb_r = os.path.join('result', os.path.basename(src_xyz))
            mb_r = src_xyz

            with open(mb_r, 'r') as in_m:
                for xyz in xyzfun.xyz_parse(in_m, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
            utils.remove_glob(src_xyz)
            #vdatumfun.vdatum_clean_result()
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_mb))
        utils.remove_glob(src_mb)

    def _dump_xyz(self, entry, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
    
    def _yield_results_to_xyz(self):
        if len(self._results) == 0: self.run()
        vdc = copy.deepcopy(vdatumfun._vd_config)
        xyzc = copy.deepcopy(xyzfun._xyz_config)        
        if vdc['jar'] is None: vdc['jar'] = vdatumfun.vdatum_locate_jar()[0]
        
        for entry in self._results:
            for xyz in self._yield_xyz(vdc, xyzc):
                yield(xyz)

    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
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
        utils.echo_msg('loading USACE fetch module...')
        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self._status = 0
        self._req = None
        self._results = []
        self._datalists_code = 407
        self._verbose = False
        self.region = extent

    def run(self, stype = None):
        '''Run the USACE fetching module'''
        if self.region is None: return([])
        self.data = {
            'geometry': regions.region_format(self.region, 'bbox'),
            'inSR':4326,
            'outSR':4326,
            'f':'pjson',
        }
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
## fetch extracts of the GMRT. - Global Extents
## https://www.gmrt.org/index.php
##
## =============================================================================
class gmrt:
    '''Fetch raster data from the GMRT'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading GMRT fetch module...')
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')
        self._status = 0
        self._results = []
        self.region = extent
        self._datalists_code = 408
        self._verbose = False

    def run(self, res = 'max', fmt = 'geotiff'):
        '''Run the GMRT fetching module'''
        if self.region is None: return([])
        self.data = {
            'north':self.region[3],
            'west':self.region[0],
            'south':self.region[2],
            'east':self.region[1],
            'mformat':'json',
            'resolution':res,
            'format':fmt,
        }
        self._req = fetch_req(self._gmrt_grid_urls_url, params = self.data, tries = 10, timeout = 2)
        if self._req is not None:
            gmrt_urls = self._req.json()
            for url in gmrt_urls:
                opts = {}
                for url_opt in url.split('?')[1].split('&'):
                    opt_kp = url_opt.split('=')
                    opts[opt_kp[0]] = opt_kp[1]
                url_region = [float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])]
                outf = 'gmrt_{}_{}.{}'.format(opts['layer'], regions.region_format(url_region, 'fn'), gdalfun.gdal_fext(opts['format']))
                self._results.append([url, outf, 'gmrt'])
        return(self._results)

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry, res = 'max', fmt = 'geotiff'):
        src_gmrt = 'gmrt_tmp.{}'.format(gdalfun.gdal_fext(fmt))
        if fetch_file(entry[0], src_gmrt, callback = lambda: False, verbose = self._verbose) == 0:
            #try:
            src_ds = gdal.Open(src_gmrt)
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
                    yield(xyz)
            src_ds = None
            #except: utils.echo_error_msg('could not read gmrt data: {}'.format(src_gmrt))
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_gmrt))
        utils.remove_glob(src_gmrt)

    def _dump_xyz(self, src_gmrt, res = 'max', fmt = 'geotiff', dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_gmrt, res, fmt):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    def _yield_results_to_xyz(self, res = 'max', fmt = 'geotiff'):
        self.run(res, fmt)
        for entry in self._results:
            for xyz in self._yield_xyz(entry, res, fmt):
                yield(xyz)
                
    def _dump_results_to_xyz(self, res = 'max', fmt = 'geotiff', dst_port = sys.stdout):
        for xyz in self._yield_to_xyz(res, fmt):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)


## =============================================================================
##
## EMODNET Fetch
##
## fetch extracts of the EMOD DTM - Mostly European extents
## https://portal.emodnet-bathymetry.eu/
##
## =============================================================================
class emodnet:
    '''Fetch raster data from the EMODNET DTM'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading EMODNET fetch module...')
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._outdir = os.path.join(os.getcwd(), 'emodnet')
        self._results = []
        self.region = extent
        self._datalists_code = 415
        self._verbose = True

    def run(self):
        '''Run the EMODNET fetching module'''
        if self.region is None: return([])

        desc_data = {
            'request': 'DescribeCoverage',
            'version': '2.0.1',
            'CoverageID': 'emodnet:mean',
            'service': 'WCS',
            }

        desc_req = fetch_req(self._emodnet_grid_url, params = desc_data)
        desc_results = lxml.etree.fromstring(desc_req.text.encode('utf-8'))
        g_env = desc_results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces = namespaces)[0]
        hl = map(float, g_env.find('{http://www.opengis.net/gml/3.2}high').text.split())

        g_bbox = desc_results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
        lc = map(float, g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split())
        uc = map(float, g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split())

        ds_region = [lc[1], uc[1], lc[0], uc[0]]
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]

        if regions.regions_intersect_ogr_p(self.region, ds_region):
            emodnet_wcs = '{}service=WCS&request=GetCoverage&version=1.0.0&Identifier=emodnet:mean&coverage=emodnet:mean&format=GeoTIFF&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
                                      .format(self._emodnet_grid_url, regions.region_format(self.region, 'bbox'), resx, resy)
            outf = 'emodnet_{}.tif'.format(regions.region_format(self.region, 'fn'))
            self._results.append([emodnet_wcs, outf, 'emodnet'])
            
        return(self._results)

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry):
        src_emodnet = 'emodnet_tmp.tif'
        if fetch_file(entry[0], src_emodnet, callback = lambda: False, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_emodnet)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
        utils.remove_glob(src_emodnet)

    def _dump_xyz(self, src_emodnet, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_emodnet):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    def _yield_results_to_xyz(self):
        self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)
                
    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)

## =============================================================================
##
## CHS Fetch
##
## fetch bathymetric soundings from the Canadian Hydrographic Service (CHS) - Canada Only
## https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d
##
## =============================================================================
class chs:
    '''Fetch bathymetric soundings from the CHS'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading CHS fetch module...')
        #self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._outdir = os.path.join(os.getcwd(), 'chs')
        self._results = []
        self.region = extent
        self._datalists_code = 418
        self._verbose = True

    def run(self):
        '''Run the CHS fetching module'''
        if self.region is None: return([])

        desc_data = {
            'request': 'DescribeCoverage',
            'version': '2.0.1',
            'CoverageID': 'caris:NONNA 100',
            'service': 'WCS',
            }
        desc_req = fetch_req(self._chs_url, params = desc_data)

        desc_results = lxml.etree.fromstring(desc_req.text.encode('utf-8'))
        g_env = desc_results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope', namespaces = namespaces)[0]
        hl = map(float, g_env.find('{http://www.opengis.net/gml/3.2}high').text.split())

        g_bbox = desc_results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
        lc = map(float, g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split())
        uc = map(float, g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split())

        ds_region = [lc[1], uc[1], lc[0], uc[0]]
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]

        if regions.regions_intersect_ogr_p(self.region, ds_region):
            chs_wcs = '{}service=WCS&request=GetCoverage&version=1.0.0&Identifier=caris:NONNA+100&coverage=caris:NONNA+100&format=GeoTIFF&bbox={}&resx={}&resy={}&crs=EPSG:4326'\
                                  .format(self._chs_url, regions.region_format(self.region, 'bbox'), resx, resy)
            outf = 'chs_{}.tif'.format(regions.region_format(self.region, 'fn'))
            self._results.append([chs_wcs, outf, 'chs'])
            
        return(self._results)

    ## ==============================================
    ## Process results to xyz
    ## ==============================================    
    def _yield_xyz(self, entry):
        src_chs = 'chs_tmp.tif'
        if fetch_file(entry[0], src_chs, callback = lambda: False, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_chs)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_chs))
        utils.remove_glob(src_chs)

    def _dump_xyz(self, src_chs, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_chs):
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    def _yield_results_to_xyz(self):
        self.run()
        for entry in self._results:
            for xyz in self._yield_xyz(entry):
                yield(xyz)
                
    def _dump_results_to_xyz(self, dst_port = sys.stdout):
        for xyz in self._yield_results_to_xyz():
            xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
## =============================================================================
##
## GEBCO Fetch
##
## fetch extracts of the GEBCO Bathymetry Grid. - Global Extents
## https://www.gebco.net/
##
## =============================================================================
class gebco:
    '''Fetch raster data from the GEBCO Bathymetric Grid'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading GEBCO fetch module...')
        self._gebco_grid_url = "https://www.gebco.net/data_and_products/gebco_web_services/2019/mapserv?"
        self._outdir = os.path.join(os.getcwd(), 'gebco')
        self._status = 0
        self._results = []
        self.region = extent
        self._datalists_code = 412
        self._verbose = False

    def run(self, fmt = 'image/tiff'):
        '''Run the GEBCO fetching module'''
        if self.region is None: return([])

        self.data = {
            'request': 'getmap',
            'service': 'wms',
            'crs': 'EPSG:4326',
            'format': fmt,
            'layers': 'gebco_2019_grid',
            'version': '1.3.0',
            'BBOX': regions.region_format(self.region, 'bbox'),
        }

        self._req = fetch_req(self._gebco_grid_url, params = self.data, tries = 10, timeout = 2)
        if self._req is not None:
            #gebco_results = self._req.json()
            url = self._req.url
            outf = 'gebco_{}.tif'.format(regions.region_format(self.region, 'fn'))
            self._results.append([url, outf, 'gebco'])
        return(self._results)

    # ## ==============================================
    # ## Process results to xyz
    # ## ==============================================    
    # def _yield_xyz(self, entry, res = 'max', fmt = 'geotiff'):
    #     src_gmrt = 'gmrt_tmp.{}'.format(gdalfun.gdal_fext(fmt))
    #     print(entry)
    #     if fetch_file(entry[0], src_gmrt, callback = lambda: False, verbose = self._verbose) == 0:
    #         #try:
    #         src_ds = gdal.Open(src_gmrt)
    #         if src_ds is not None:
    #             srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
    #             for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose):
    #                 yield(xyz)
    #         src_ds = None
    #         #except: utils.echo_error_msg('could not read gmrt data: {}'.format(src_gmrt))
    #     else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_gmrt))
    #     utils.remove_glob(src_gmrt)

    # def _dump_xyz(self, src_gmrt, res = 'max', fmt = 'geotiff', dst_port = sys.stdout):
    #     for xyz in self._yield_xyz(src_gmrt, res, fmt):
    #         xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
    # def _yield_results_to_xyz(self, res = 'max', fmt = 'geotiff'):
    #     self.run(res, fmt)
    #     for entry in self._results:
    #         for xyz in self._yield_xyz(entry, res, fmt):
    #             yield(xyz)
                
    # def _dump_results_to_xyz(self, res = 'max', fmt = 'geotiff', dst_port = sys.stdout):
    #     for xyz in self._yield_to_xyz(res, fmt):
    #         xyzfun.xyz_line(xyz, dst_port, self._verbose)
            
## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS - US Only
##
## =============================================================================
class ngs:
    '''Fetch NGS monuments from NOAA'''
    def __init__(self, extent = None, filters = [], callback = None):
        utils.echo_msg('loading NGS Monument fetch module...')
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')
        self._status = 0
        self._req = None
        self._results = []
        self.region = extent
        self._datalists_code = 409
        self._verbose = False

        self._desc = '''
        National Geodetic Service (NGS) Monuments 

        options:
         <csv=False>
        '''

    def run(self, csv = False):
        '''Run the NGS (monuments) fetching module.'''
        if self.region is None: return([])
        self.data = { 'maxlon':self.region[0],
                      'minlon':self.region[1],
                      'maxlat':self.region[3],
                      'minlat':self.region[2] }

        self._req = fetch_req(self._ngs_search_url, params = self.data)
        if self._req is not None: self._results.append([self._req.url, 'ngs_results_{}.json'.format(regions.region_format(self.region, 'fn')), 'ngs'])
        return(self._results)

    def proc_ngs(self):
        '''process ngs monuments'''
        
        with open(self.src_file, 'r') as json_file: r = json.load(json_file)
        if len(r) > 0:
            for rn, this_region in enumerate(self.dst_regions):
                outfile = open(os.path.join(self.src_dir, 'ngs_results_{}.csv'.format(this_region.fn)), 'w')
                outcsv = csv.writer(outfile)
                outcsv.writerow(r[0].keys())
                [outcsv.writerow(row.values()) for row in r]
                outfile.close()

## ==============================================
## fetches processing (datalists fmt:400 - 499)
## ==============================================
def fetch_inf_entry(entry = []):
    return([-180,180,-90,90])

def fetch_yield_entry(entry = ['nos:datatype=xyz'], region = None, verbose = False):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''
    
    fetch_mod = entry[0].split(':')[0]
    fetch_args = entry[0].split(':')[1:]
    
    fl = fetch_infos[fetch_mod][0](regions.region_buffer(region, 5, pct = True), [], lambda: False)
    args_d = utils.args2dict(fetch_args, {})
    fl._verbose = verbose

    for xyz in fl._yield_results_to_xyz(**args_d):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dump_entry(entry = ['nos:datatype=nos'], dst_port = sys.stdout, region = None, verbose = False):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    
    for xyz in fetch_yield_entry(entry, region, verbose):
        xyz_line(xyz, dst_port, True)
        
def fetch_module_yield_entry(entry, region = None, verbose = False, module = 'dc'):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''

    fl = fetch_infos[module][0](regions.region_buffer(region, 5, pct = True), [], lambda: False)
    fl._verbose = verbose
    fetch_entry = [entry[0], entry[0].split('/')[-1], module]
    
    for xyz in fl._yield_xyz(fetch_entry):
        yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dc_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, module = 'dc'):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    
    for xyz in fetch_module_yield_entry(entry, region, verbose, module):
        xyz_line(xyz, dst_port, True)        
                
## =============================================================================
##
## Run fetches from command-line
##
## =============================================================================
fetch_infos = { 
    'dc':[lambda x, f, c: dc(x, f, c), '''NOAA Digital Coast
    Lidar and Raster data from NOAA's Digital Coast

    < dc:datatype=None:index=False:update=False >
     :datatype=[lidar/raster] - Only fetch lidar or raster data.
     :index=[True/False] - True to display indexed results.
     :update=[True/False] - True to update stored reference vector.'''],
     'nos':[lambda x, f, c: nos(x, f, c), '''NOAA NOS Bathymetric Data
    Bathymetry surveys and data (xyz & BAG)

    < nos:datatype=None:update=False >
     :datatype=[bag/xyz] - Only fetch BAG or Hydro-XYZ data.
     :update=[True/False] - True to update stored reference vector.'''],
     'charts':[lambda x, f, c: charts(x, f, c), '''NOAA Nautical CHARTS (RNC & ENC)
    Raster and Vector U.S. Nautical Charts

    < charts:datatype=None:update=False >
     :dataype=[ENC/RNC] - Only fetch either ENC or RNC data.
     :update=[True/False] - True to update stored reference vector.'''],
    'srtm_cgiar':[lambda x, f, c: srtm_cgiar(x, f, c), '''SRTM elevation data (CGIAR)
    Global Topography at 90m DEMs

    < srtm_cgiar >'''],
    'srtm_plus':[lambda x, f, c: srtm_plus(x, f, c), '''SRTM15+ elevation data (Scripps)
    Global Bathymetry and Topography at 15 Arc Sec:SRTM15+
    
    < srtm_plus >'''],
    'tnm':[lambda x, f, c: tnm(x, f, c), '''The National Map (TNM) from USGS
    Various datasets from USGS's National Map.
    The National Map is a collaborative effort among the USGS and other Federal, State, and local partners to improve and deliver topographic information for the Nation.

    https://www.usgs.gov/core-science-systems/national-geospatial-program/national-map
    https://viewer.nationalmap.gov/advanced-viewer/

    < tnm:ds=1:sub_ds=None:formats=None:index=False >
     :index=[True/False] - True to display an index of available datasets.
     :ds=[dataset index value (0-15)] - see :index=True for dataset index values.
     :sub_ds=[sub-dataset index value (0-x)] - see :index=True for sub-dataset index values.
     :formats=[data-set format] - see :index=True for dataset format options.'''],
    'mb':[lambda x, f, c: mb(x, f, c), '''NOAA MULTIBEAM survey data
    Multibeam Bathymetric Data from NOAA.
    In addition to deepwater data, the Multibeam Bathymetry Database (MBBDB) includes hydrographic multibeam survey data from NOAA's National Ocean Service (NOS).

    https://www.ngdc.noaa.gov/mgg/bathymetry/multibeam.html

    < mb >'''],
    'gmrt':[lambda x, f, c: gmrt(x, f, c), '''The Global Multi-Reosolution Topography Data Synthesis (GMRT) 
    Global Elevation raster dataset via GMRT.

    < gmrt:res=max:fmt=geotiff >
     :res=[an Integer power of 2 zoom level (<=1024)]
     :fmt=[netcdf/geotiff/esriascii/coards]'''],
    'cudem':[lambda x, f, c: cudem(x, f, c), '''ncei cudem thredds catalog
    NOAA NCEI THREDDS DEM Catalog Access

    < cudem >'''],
    'usace':[lambda x, f, c: usace(x, f, c), '''USACE bathymetry surveys via eHydro
    Bathymetric Channel surveys from USACE - U.S. only.
    The hydrographic surveys provided by this application are to be used for informational purposes only and should not be used as a navigational aid. 
    Channel conditions can change rapidly and the surveys may or may not be accurate.

    https://navigation.usace.army.mil/Survey/Hydro

    < usace >'''],
    'ngs':[lambda x, f, c: ngs(x, f, c), '''NOAA NGS Monument Data
    Monument data from NOAA's Nagional Geodetic Survey (NGS) monument dataset.

    < ngs >'''],
    'mar_grav':[lambda x, f, c: mar_grav(x, f, c), '''Marine Gravity from Sattelite Altimetry topographic data.
    Elevation data from Scripps Marine Gravity dataset.

    < mar_grav >'''],
    'emodnet':[lambda x, f, c: emodnet(x, f, c), '''EMODNET Elevation Data
    European Bathymetry/Topographic data from EMODNET

    https://portal.emodnet-bathymetry.eu/help/help.html

    < emodnet >'''],
    'chs':[lambda x, f, c: chs(x, f, c), '''CHS Bathymetry
    CHS NONNA 100m Bathymetric Survey Grids
    Non-Navigational gridded bathymetric data based on charts and soundings.

    https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d

    < chs >'''],
    'hrdem':[lambda x, f, c: hrdem(x, f, c), '''Canada HRDEM Elevation Data
    Canada High-Resolution Digital Elevation Model (HRDEM).
    Collection of lidar-derived DTMs across Canada.

    https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995

    < hredem >'''],
}

def fetch_desc(x):
    fd = []
    for key in x: 
        fd.append('{:18}{}'.format(key, x[key][-1]))
    return fd

_fetch_long_desc = lambda x: 'fetches modules:\n% fetches ... <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['{:14}{}\n'.format(key, x[key][-1]) for key in x]) + '\n'
_fetch_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

_usage = '''{} [OPTIONS] <module[:parameter=value]* ...>

Fetch geographic elevation data.

General Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -f, --filter\t\tSQL style filters for certain datasets.
\t\t\tFields to filter include: 'Name', 'Date' and 'Datatype'
  -p, --process\t\tProcess fetched data to ASCII XYZ format in WGS84. <beta>

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Modules (see fetches --modules <module-name> for more info):
  {}

Examples:
 % {} -R -90.75/-88.1/28.7/31.25 nos -f "Date > 2000"
 % {} -R region.shp -p dc nos:datatype=bag charts:datatype=enc
 % {} -R region.shp dc:datatype=lidar -l > dc_lidar.urls
 % {} -R -89.75/-89.5/30.25/30.5 tnm:ds=4:formats=IMG gmrt:res=max:fmt=geotiff

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(os.path.basename(sys.argv[0]),
           _fetch_short_desc(fetch_infos),
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]))
           
#'\n  '.join(fetch_desc(fetch_infos)),
def fetches_cli(argv = sys.argv):
    status = 0
    extent = None
    want_list = False
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
    while i < len(argv):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            extent = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            extent = str(arg[2:])
        elif arg == '--list-only' or arg == '-l':
            want_list = True
        elif arg == '--filter' or arg == '-f':
            f.append(argv[i + 1])
            i = i + 1
        elif arg == '--process' or arg == '-p':
            want_proc = True
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format( _version))
            sys.exit(0)
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in fetch_infos.keys():
                    sys.stderr.write(_fetch_long_desc({k: fetch_infos[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_fetch_long_desc(fetch_infos))
            except: sys.stderr.write(_fetch_long_desc(fetch_infos))
            sys.exit(0)
        else: 
            opts = arg.split(':')
            if opts[0] in fetch_infos.keys():
                mod_opts[opts[0]] = list(opts[1:])
            else: utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
        i = i + 1

    if len(mod_opts) == 0:
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must select a fetch module')
        sys.exit(1)
        
    for key in mod_opts.keys():
        mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]
        
    ## ==============================================
    ## process input region(s)
    ## ==============================================
    if extent is not None:
        try:
            these_regions = [[float(x) for x in extent.split('/')]]
        except ValueError: these_regions = gdalfun.gdal_ogr_regions(extent)
        if len(these_regions) == 0: status = -1
        for this_region in these_regions:
            if not regions.region_valid_p(this_region): status = -1
        utils.echo_msg('loaded {} region(s)'.format(len(these_regions)))
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
            utils.echo_msg('running fetch module {} on region {} ({}/{})...\
            '.format(fetch_mod, regions.region_format(this_region, 'str'), rn+1, len(these_regions)))
            fl = fetch_infos[fetch_mod][0](regions.region_buffer(this_region, 5, pct = True), f, lambda: stop_threads)
            #fl._verbose = True
            args_d = utils.args2dict(args)
            try:
                r = fl.run(**args_d)
            except ValueError as e:
                utils.echo_error_msg('something went wrong, {}'.format(e))
                sys.exit(-1)
            except Exception as e:
                utils.echo_error_msg('{}'.format(e))
                sys.exit(-1)
            utils.echo_msg('found {} data files.'.format(len(r)))
            
            if want_list:
                for result in r:
                    print(result[0])
            else:
                fr = fetch_results(r, this_region, fl._outdir, fl if want_proc else None, lambda: stop_threads)
                fr.daemon = True
                try:
                    fr.start()
                    while True:
                        time.sleep(2)
                        sys.stderr.write('\x1b[2K\r')
                        perc = float((len(r) - fr.fetch_q.qsize())) / len(r) * 100 if len(r) > 0 else 1
                        sys.stderr.write('fetches: fetching remote data files [{}%]'.format(perc))
                        sys.stderr.flush()
                        if not fr.is_alive():
                            break
                except (KeyboardInterrupt, SystemExit):
                    utils.echo_error_msg('user breakage...please wait for while fetches exits.')
                    fl._status = -1
                    stop_threads = True
                    while not fr.fetch_q.empty():
                        try:
                            fr.fetch_q.get(False)
                        except Empty: continue
                        fr.fetch_q.task_done()
                fr.join()
            utils.echo_msg('ran fetch module {} on region {} ({}/{})...\
            '.format(fetch_mod, regions.region_format(this_region, 'str'), rn+1, len(these_regions)))
            if want_proc: vdatumfun.vdatum_clean_result()

### End
