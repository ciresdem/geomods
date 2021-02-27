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
## current fetch modules: dc, nos, mb, charts, usace, srtm_cgiar, srtm_plus, tnm, gmrt, emodnet, mar_grav, osm
##
## Possible BAG errors with GDAL >= 3
##
### Code:

import os
import sys
import time

import requests
import ftplib
import urllib

import lxml.html as lh
import lxml.etree
import json
import copy

import ogr
import gdal

import threading
try:
    import Queue as queue
except: import queue as queue

from geomods import utils
from geomods import regions
from geomods import gdalfun
from geomods import xyzfun

_version = '0.6.0'

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## =============================================================================
r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(_version) }

def fetch_queue(q, p, s):
    '''fetch queue `q` of fetch results'''
    
    while True:
        fetch_args = q.get()
        this_region = fetch_args[2]
        this_dt = fetch_args[4].lower()
        fetch_args[2] = None
        if not fetch_args[3]():
            if p is None and s is None:
                if fetch_args[0].split(':')[0] == 'ftp':
                    fetch_ftp_file(*tuple(fetch_args))
                else: fetch_file(*tuple(fetch_args))
            else:
                if s is None:
                    if not os.path.exists(os.path.dirname(fetch_args[1])):
                        try:
                            os.makedirs(os.path.dirname(fetch_args[1]))
                        except: pass
                    o_x_fn = '.'.join(fetch_args[1].split('.')[:-1]) + '.xyz'
                    utils.echo_msg('processing local file: {}'.format(o_x_fn))
                    if not os.path.exists(o_x_fn):
                        with open(o_x_fn, 'w') as out_xyz:
                            p._dump_xyz([fetch_args[0], fetch_args[1], fetch_args[-1]], epsg = 4326, dst_port = out_xyz)
                else: s._dump_xyz([fetch_args[0], fetch_args[1], fetch_args[-1]], epsg = 4326)
                    
        q.task_done()

def fetch_ftp_file(src_url, dst_fn, params = None, callback = None, datatype = None, verbose = False):
    '''fetch an ftp file via urllib'''
    
    status = 0
    f = None
    halt = callback
    
    if verbose: utils.echo_msg('fetching remote ftp file: {}...'.format(src_url))
    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    try:
        f = urllib.request.urlopen(src_url)
    except: status - 1
    
    if f is not None:
        with open(dst_fn, 'wb') as local_file:
             local_file.write(f.read())
    if verbose: utils.echo_msg('fetched remote ftp file: {}.'.format(os.path.basename(src_url)))
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

def fetch_nos_xml(src_url, timeout = 2, read_timeout = 10, verbose = False):
    '''fetch src_url and return it as an XML object'''

    results = lxml.etree.fromstring('<?xml version="1.0"?><!DOCTYPE _[<!ELEMENT _ EMPTY>]><_/>'.encode('utf-8'))
    try:
        req = fetch_req(src_url, timeout = timeout, read_timeout = read_timeout)
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

class fetch_results(threading.Thread):
    '''fetch results gathered from a fetch module.
    results is a list of URLs with data type'''
    
    def __init__(self, results, region, out_dir, proc = None, stream = None, vdc = None, callback = lambda: False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.results = results
        self.region = region
        self._outdir = out_dir
        self.stop_threads = callback
        self.proc = proc
        self.stream = stream
        
    def run(self):
        for _ in range(3):
            t = threading.Thread(target = fetch_queue, args = (self.fetch_q, self.proc, self.stream))
            t.daemon = True
            t.start()

        fq = [[row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2]] for row in self.results]
        [self.fetch_q.put([row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2]]) for row in self.results]
        self.fetch_q.join()
        
## =============================================================================
##
## ISO XML Metadata parsing
##
## =============================================================================
class iso_xml:
    def __init__(self, xml_url, timeout = 2, read_timeout = 10):
        self.url = xml_url
        self.xml_doc = self._fetch(timeout = timeout, read_timeout = read_timeout)
        self.namespaces = {
            'gmd': 'http://www.isotc211.org/2005/gmd', 
            'gmi': 'http://www.isotc211.org/2005/gmi', 
            'gco': 'http://www.isotc211.org/2005/gco',
            'gml': 'http://www.isotc211.org/2005/gml',
            'th': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0',
        }
        
    def _fetch(self, timeout = 2, read_timeout = 10):
        return(fetch_nos_xml(self.url, timeout = timeout, read_timeout = read_timeout))

    def title(self):
        t = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString', namespaces = self.namespaces)
        return(t.text if t is not None else 'Unknown')
        
    def bounds(self, geom = True):
        wl = self.xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        el = self.xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = self.namespaces)
        sl = self.xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = self.namespaces)                            
        nl = self.xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = self.namespaces)                           
        if wl is not None and el is not None and sl is not None and nl is not None:
            region = [float(wl.text), float(el.text), float(sl.text), float(nl.text)]
            if geom: return(gdalfun.gdal_region2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)]))
            else: return(region)
        else: return(None)

    def polygon(self, geom = True):
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces = self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces = self.namespaces)
            [opoly.append([float(x) for x in node.text.split()]) for node in nodes]
            if geom: return(gdalfun.gdal_wkt2geom(gdalfun.gdal_create_polygon(opoly)))
            else: return(opoly)
        else: return(None)
        
    def date(self):
        dt = self.xml_doc.find('.//gmd:date/gco:Date', namespaces = self.namespaces)
        if dt is None:
            dt = self.xml_doc.find('.//gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date', namespaces = self.namespaces)
        return(dt.text[:4] if dt is not None else '0000')

    def xml_date(self):
        mddate = self.xml_doc.find('.//gmd:dateStamp/gco:DateTime', namespaces = self.namespaces)
        return(utils.this_date() if mddate is None else mddate.text)
        
    def reference_system(self):
        ref_s = self.xml_doc.findall('.//gmd:MD_ReferenceSystem', namespaces = self.namespaces)
        if ref_s is None or len(ref_s) == 0: return(None, None)
        h_epsg = ref_s[0].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
        if h_epsg is not None: h_epsg = h_epsg.text.split(':')[-1]
        if len(ref_s) > 1:
            v_epsg = ref_s[1].find('.//gmd:code/gco:CharacterString', namespaces = self.namespaces)
            if v_epsg is not None: v_epsg = v_epsg.text.split(':')[-1]
        else: v_epsg = None
            
        return(h_epsg, v_epsg)

    def abstract(self):
        try:
            abstract = self.xml_doc.find('.//gmd:abstract/gco:CharacterString', namespaces = self.namespaces)
            abstract = '' if abstract is None else abstract.text
        except: abstract = ''
        return(abstract)

    def linkages(self):
        linkage = self.xml_doc.find('.//{*}linkage/{*}URL', namespaces = self.namespaces)
        if linkage is not None: linkage = linkage.text
        return(linkage)
    
    def data_links(self):
        dd = {}
        
        dfs = self.xml_doc.findall('.//gmd:MD_Format/gmd:name/gco:CharacterString', namespaces = self.namespaces)
        dus = self.xml_doc.findall('.//gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces =  self.namespaces)

        if dfs is not None:
            for i,j in enumerate(dfs):
                dd[j.text] = dus[i].text

        return(dd)

## =============================================================================
##
## Fetches Remote Elevation Datalist (FRED)
##
## the reference vector location and related functions
##
## the reference vector has the following fields:
## Name, ID, Date, Metadata, DataLink, IndexLink, DataType, DataSource
##
## =============================================================================
this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

class FRED:
    def __init__(self, verbose = False):
        self._verbose = verbose
        self.fetchdata = os.path.join(this_dir, 'data')
        self.driver = ogr.GetDriverByName('GeoJSON')
        self.fetch_v = 'FRED.geojson'
        if os.path.exists(self.fetch_v):
            self.FREDloc = self.fetch_v
        elif os.path.exists(os.path.join(self.fetchdata, self.fetch_v)):
            self.FREDloc = os.path.join(self.fetchdata, self.fetch_v)
        else: self.FREDloc = self.fetch_v
        if self._verbose: utils.echo_msg('using {}'.format(self.FREDloc))
        self.ds = None
        self.layer = None
        
    def _create_ds(self):
        utils.remove_glob(self.FREDloc)
        self.ds = self.driver.CreateDataSource(self.FREDloc)
        
        self.layer = self.ds.CreateLayer('FRED', None, ogr.wkbMultiPolygon)
        ldfn = self.layer.GetLayerDefn()
        self.layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        self.layer.CreateField(ogr.FieldDefn('MetadataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataDate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('IndexLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataType', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataSource', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('HorizontalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('VerticalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('LastUpdate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Info', ogr.OFTString))

    def _open_ds(self, mode = 0):
        try:
            self.ds = self.driver.Open(self.FREDloc, mode)
        except: self.ds = None
        if self.ds is None:
            self._create_ds()
        self.layer = self.ds.GetLayer()
        
    def _close_ds(self):
        self.layer = None
        self.ds = None

    def _get_fields(self):
        schema = []
        ldefn = self.layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
        return(schema)
        
    def _add_feature(self, survey):
        '''add a survey to the reference vector layer'''

        layer_defn = self.layer.GetLayerDefn()
        feat = ogr.Feature(layer_defn)
        geom = ogr.CreateGeometryFromJson(survey[1])
        feat.SetGeometry(geom)
        [feat.SetField(key, survey[0][key]) for key in survey[0].keys()]
        self.layer.CreateFeature(feat)
        feat = None

    def _add_surveys(self, surveys):
        '''update or create a reference vector using a list of surveys'''

        self._open_ds(1)
        if self.layer is not None:
            [self.layer.SetFeature(ff) for ff in self.layer]
            [self._add_feature(i) for i in surveys]
        
    def _filter(self, region, where = [], layers = []):
        '''Search for data in the reference vector file'''

        _boundsGeom = gdalfun.gdal_region2geom(region)
        _results = []

        if self._verbose: utils.echo_msg('filtering {}...'.format(self.FREDloc))
        self._open_ds()

        f = 0
        
        for layer in layers:
            this_layer = self.layer
            this_layer.SetAttributeFilter("DataSource = '{}'".format(layer))
            [this_layer.SetAttributeFilter('{}'.format(filt)) for filt in where]
            for feat in this_layer:
                geom = feat.GetGeometryRef()
                if geom is not None:
                    if _boundsGeom.Intersects(geom):
                        _results.append({})
                        f_j = json.loads(feat.ExportToJson())
                        for key in f_j['properties'].keys():
                            _results[-1][key] = feat.GetField(key)
                        
            this_layer = None
        self._close_ds()
        if self._verbose: utils.echo_msg('filtered \033[1m{}\033[m data files from FRED'.format(len(_results)))
        return(_results)

## heaps of thanks to https://github.com/fitnr/stateplane
FIPS_TO_EPSG = {
    "0101": "26929", "0102": "26930", "0201": "26948", "0202": "26949",
    "0203": "26950", "0301": "26951", "0302": "26952", "0401": "26941",
    "0402": "26942", "0403": "26943", "0404": "26944", "0405": "26945",
    "0406": "26946", "0501": "26953", "0502": "26954", "0503": "26955",
    "0600": "26956", "0700": "26957", "0901": "26958", "0902": "26959",
    "0903": "26960", "1001": "26966", "1002": "26967", "1101": "26968",
    "1102": "26969", "1103": "26970", "1201": "26971", "1202": "26972",
    "1301": "26973", "1302": "26974", "1401": "26975", "1402": "26976",
    "1501": "26977", "1502": "26978", "1600": "03088", "1601": "2205",
    "1602": "26980", "1701": "26981", "1702": "26982", "1703": "32199",
    "1801": "26983", "1802": "26984", "1900": "26985", "2001": "26986",
    "2002": "26987", "2111": "26988", "2112": "26989", "2113": "26990",
    "2201": "26991", "2202": "26992", "2203": "26993", "2301": "26994",
    "2302": "26995", "2401": "26996", "2402": "26997", "2403": "26998",
    "2500": "32100", "2600": "32104", "2701": "32107", "2702": "32108",
    "2703": "32109", "2800": "32110", "2900": "32111", "3001": "32112",
    "3002": "32113", "3003": "32114", "3101": "32115", "3102": "32116",
    "3103": "32117", "3104": "32118", "3200": "32119", "3301": "32120",
    "3302": "32121", "3401": "32122", "3402": "32123", "3501": "32124",
    "3502": "32125", "3601": "32126", "3602": "32127", "3701": "32128",
    "3702": "32129", "3800": "32130", "3900": "32133", "4001": "32134",
    "4002": "32135", "4100": "32136", "4201": "32137", "4202": "32138",
    "4203": "32139", "4204": "32140", "4205": "32141", "4301": "32142",
    "4302": "32143", "4303": "32144", "4400": "32145", "4501": "32146",
    "4502": "32147", "4601": "32148", "4602": "32149", "4701": "32150",
    "4702": "32151", "4801": "32152", "4802": "32153", "4803": "32154",
    "4901": "32155", "4902": "32156", "4903": "32157", "4904": "32158",
    "5001": "26931", "5002": "26932", "5003": "26933", "5004": "26934",
    "5005": "26935", "5006": "26936", "5007": "26937", "5008": "26938",
    "5009": "26939", "5010": "26940", "5101": "26961", "5102": "26962",
    "5103": "26963", "5104": "26964", "5105": "26965", "5200": "32161"
}
    
## =============================================================================
##
## Fetches Modules
## Each module should have at class with a self._update function that returns a list of survey dictionaries,
## a self._parse_results function that returns a list of urls to remote data files parsed from the reference vector,
## and a self._yield_xyz function to yield the xyz data from the downloaded remote data files.
##
## self._update returns list of surveys
## -- surveys is a list of dictionaries/geometry with field-value pairs:
## [ {'Name': 'fetch_data', etc.} wkt_geom ...]
## self._parse_results returns list of urls
## -- urls is a list of [[url, file-name, file-type] ...]
## self._yield_xyz yields xyz data
##
## =============================================================================
_fetch_modules = {'dc': lambda r, f, c, v: dc(r, f, c, v),
                  'nos': lambda r, f, c, v: nos(r, f, c, v),
                  'charts': lambda r, f, c, v: charts(r, f, c, v),
                  'ncei_thredds': lambda r, f, c, v: ncei_thredds(r, f, c, v),
                  'usace': lambda r, f, c, v: usace(r, f, c, v),
                  'tnm': lambda r, f, c, v: tnm(r, f, c, v),
                  'gmrt': lambda r, f, c, v: gmrt(r, f, c, v),
                  'mb': lambda r, f, c, v: mb(r, f, c, v),
                  'mar_grav': lambda r, f, c, v: mar_grav(r, f, c, v),
                  'srtm_plus': lambda r, f, c, v: srtm_plus(r, f, c, v),
                  'emodnet': lambda r, f, c, v: emodnet(r, f, c, v),
                  'ngs': lambda r, f, c, v: ngs(r, f, c, v),
                  'chs': lambda r, f, c, v: chs(r, f, c, v),
                  'hrdem': lambda r, f, c, v: hrdem(r, f, c, v),
                  }
_fetch_long_desc = lambda x: 'fetches modules:\n% fetches ... <mod>:key=val:key=val...\n\n  ' + '\n  '.join(['{:14}{}\n'.format(key, x[key](None, [], None, False)._desc) for key in x]) + '\n'
_fetch_short_desc = lambda x: ', '.join(['{}'.format(key) for key in x])

## =============================================================================
##
## Digital Coast ('dc')
##
## =============================================================================

class dc:
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading DC fetch module...')
        
        ## ==============================================
        ## Digital Coast urls and directories
        ## ==============================================
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'lidar3_z', 'lidar4_z', 'raster1', 'raster2', 'raster5']
        self._outdir = os.path.join(os.getcwd(), 'dc')
        
        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''NOAA Digital Coast
        Lidar and Raster data from NOAA's Digital Coast

        < dc >'''

    def _update(self):
        '''Update the FRED reference vector after scanning
        the relevant metadata from Digital Coast.'''
        
        if self._verbose: utils.echo_msg('updating Digital Coast fetch module')
        self.FRED._open_ds()
        for ld in self._dc_dirs:
            cols = []
            page = fetch_html(self._dc_htdata_url + ld)
            if page is None: continue
            tables = page.xpath('//table')
            tr = tables[0].xpath('.//tr')
            if len(tr) <= 0: continue
            [cols.append(i.text_content()) for i in tr[0]]
            if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [    ].'.format(len(tr), ld))
            for i in range(1, len(tr)):
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
                if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(tr), ld, int((float(i)/len(tr)) * 100), dc['ID #']))
                self.FRED.layer.SetAttributeFilter("ID = 'DC-{}'".format(dc['ID #']))
                if len(self.FRED.layer) == 0:
                    if 'Metadata' in dc.keys():
                        this_xml = iso_xml(dc['Metadata'])
                        h_epsg, v_epsg = this_xml.reference_system()
                        geom = this_xml.bounds(geom=True)
                        if geom is not None:
                            self._surveys.append([{'Name': dc['Dataset Name'], 
                                                   'ID': 'DC-{}'.format(dc['ID #']), 
                                                   'Date': this_xml.date(), 
                                                   'MetadataLink': dc['Metadata'],
                                                   'MetadataDate': this_xml.xml_date(), 
                                                   'DataLink': dc['https'],
                                                   'IndexLink': dc['Tile Index'],
                                                   'DataType': ld.split("_")[0],
                                                   'DataSource': 'dc',
                                                   'HorizontalDatum': h_epsg,
                                                   'VerticalDatum': v_epsg,
                                                   'LastUpdate': utils.this_date(),
                                                   'Info': this_xml.abstract()}, geom.ExportToJson()])
            if self._verbose: utils.echo_msg('scanning {} surveys in {} [ OK ].'.format(len(tr), ld))
            
        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['dc']))
    
    def _parse_results(self, r):
        for surv in r:
            surv_shp_zip = os.path.basename(surv['IndexLink'])
            if fetch_file(surv['IndexLink'], surv_shp_zip, verbose = self._verbose) == 0:
                v_shps = utils.p_unzip(surv_shp_zip, ['.shp', '.shx', '.dbf', '.prj'])
                v_shp = None
                for v in v_shps:
                    if '.shp' in v: v_shp = v
                try:
                    v_ds = ogr.Open(v_shp)
                except: v_ds = None
                if v_ds is not None:
                    slay1 = v_ds.GetLayer(0)
                    for sf1 in slay1:
                        geom = sf1.GetGeometryRef()
                        if geom.Intersects(self._boundsGeom):
                            tile_url = sf1.GetField('URL').strip()
                            self._data_urls.append([tile_url, '{}/{}'.format(surv['ID'], tile_url.split('/')[-1]), surv['DataType']])
                    v_ds = slay1 = None
                [utils.remove_glob(v) for v in v_shps]
            utils.remove_glob(surv_shp_zip)
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1].lower()
        if src_ext == 'laz' or src_ext == 'las': dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img': dt = 'raster'
        else: dt = None
        if dt == 'lidar':
            if fetch_file(entry[0], src_dc, callback = lambda: False, verbose = False) == 0:
                xyz_dat = utils.yield_cmd('las2txt -verbose -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                '.format(regions.region_format(self.region, 'te'), '2 29', src_dc), verbose = False)
                xyzc = copy.deepcopy(xyzfun._xyz_config)
                xyzc['name'] = src_dc
                xyzc['epsg'] = 4326
                xyzc['warp'] = epsg
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
                srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.gdal_region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_get_epsg(src_ds)))
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose):
                    yield(xyz)
            src_ds = None
        utils.remove_glob(src_dc)

    def _dump_xyz(self, entry, epsg = 4326, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading NOS fetch module...')

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
        self._outdir = os.path.join(os.getcwd(), 'nos')
        self._nos_fmts = ['.xyz.gz', '.bag.gz', '.bag']

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''NOAA NOS Bathymetric Data
        Bathymetry surveys and data (xyz & BAG)

        < nos >'''

    def _update(self):
        '''Crawl the NOS database and update/generate the NOS reference vector.'''

        ## ==============================================
        ## load the reference vector if it exists and
        ## scan the remote nos directories
        ## ==============================================
        self.FRED._open_ds()
        
        for j in self._nos_directories:
            nosdir = j
            xml_catalog = self._nos_xml_url(nosdir)
            page = fetch_html(xml_catalog)
            rows = page.xpath('//a[contains(@href, ".xml")]/@href')
            if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [    ].'.format(len(rows), nosdir))

            ## ==============================================
            ## Parse each survey found in the directory
            ## and append it to the surveys list
            ## ==============================================
            for i, survey in enumerate(rows):
                sid = survey[:-4]
                perc = int((float(i)/len(rows)) * 100)
                if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(rows), nosdir, perc, sid))
                self.FRED.layer.SetAttributeFilter("ID = 'NOS-{}'".format(sid))

                if len(self.FRED.layer) == 0:
                    xml_url = xml_catalog + survey
                    this_xml = iso_xml(xml_url)
                    h_epsg, v_epsg = this_xml.reference_system()
                    this_data = this_xml.data_links()

                    d_links = []
                    d_types = []
                    for key in this_data.keys():
                        if key in ['GEODAS_XYZ', 'BAG', 'GRID_BAG']:
                            d_links.append(this_data[key])
                            d_types.append(key)
                    geom = this_xml.bounds(geom=True)
                    if geom is not None:
                        self._surveys.append([{'Name': this_xml.title(), 
                                               'ID': 'NOS-{}'.format(sid), 
                                               'Date': this_xml.date(), 
                                               'MetadataLink': this_xml.url,
                                               'MetadataDate': this_xml.xml_date(), 
                                               'DataLink': ','.join(list(set(d_links))),
                                               'IndexLink': '',
                                               'DataType': ','.join(list(set(d_types))),
                                               'DataSource': 'nos',
                                               'HorizontalDatum': h_epsg,
                                               'VerticalDatum': v_epsg,
                                               'LastUpdate': utils.this_date(),
                                               'Info': this_xml.abstract()}, geom.ExportToJson()])
            if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [ OK ].\n'.format(len(rows), nosdir))

        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
    
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['nos']))
                
    def _parse_results(self, r):
        '''Search the NOS reference vector and append the results
        to the results list.'''

        for surv in r:
            fldata = surv['DataLink'].split(',')
            [self._data_urls.append([i, i.split('/')[-1], surv['DataType']]) if i != '' else None for i in fldata]
        return(self._data_urls)
            
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, datatype = None, epsg = 4326):
        if xyzc is None: xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_nos = os.path.basename(entry[1])
        dt = None
        if fetch_file(entry[0], src_nos, callback = lambda: False, verbose = False) == 0:
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
                nos_f_r = nos_f

                xyzc['name'] = nos_f_r
                xyzc['delim'] = ','
                xyzc['skip'] = 1
                xyzc['xpos'] = 2
                xyzc['ypos'] = 1
                xyzc['zpos'] = 3
                xyzc['z-scale'] = -1
                xyzc['epsg'] = 4326
                xyzc['warp'] = epsg
                if os.path.exists(nos_f_r):
                    with open(nos_f_r, 'r') as in_n:
                        for xyz in xyzfun.xyz_parse(in_n, xyz_c = xyzc, region = self.region, verbose = self._verbose):
                            yield(xyz)
                    
            elif dt == 'grid_bag':
                src_bag, nos_zips = utils.procs_unzip(src_nos, ['gdal', 'tif', 'img', 'bag', 'asc'])
                nos_f = '{}.tmp'.format(os.path.basename(src_bag).split('.')[0])
                xyzc['name'] = src_bag
                xyzc['delim'] = ' '
                xyzc['skip'] = 0
                xyzc['xpos'] = 0
                xyzc['ypos'] = 1
                xyzc['zpos'] = 2
                xyzc['z-scale'] = 1
                xyzc['warp'] = None
                xyzc['epsg'] = None
                src_ds = gdal.Open(src_bag)
                if src_ds is not None:
                    srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.gdal_region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_get_epsg(src_ds)))
                    with open(nos_f, 'w') as cx:
                        for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg):
                            yield(xyz)
                    src_ds = None
                utils.remove_glob(src_bag)
            utils.remove_glob(nos_f)
        utils.remove_glob(src_nos)
        
    def _dump_xyz(self, entry, epsg = 4326, dst_port = sys.stdout):
         for xyz in self._yield_xyz(entry, epsg = epsg):
             xyzfun.xyz_line(xyz, dst_port, False)

## =============================================================================
##
## Chart Fetch - ENC & RNC
##
## Fetch digital charts from NOAA, including ENC and RNC
##
## =============================================================================
class charts():
    '''Fetch digital chart data from NOAA'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading CHARTS fetch module...')

        ## ==============================================
        ## Chart URLs and directories
        ## ==============================================
        self._enc_data_catalog = 'http://www.charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'http://www.charts.noaa.gov/RNCs/RNCProdCat_19115.xml'
        self._outdir = os.path.join(os.getcwd(), 'charts')
        self._dt_xml = { 'ENC':self._enc_data_catalog,
                         'RNC':self._rnc_data_catalog }

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''NOAA Nautical CHARTS (RNC & ENC)
        Raster and Vector U.S. Nautical Charts

        < charts >'''
        
        self._info = ''

    def _update(self):
        '''Update or create the reference vector file'''
        
        self.FRED._open_ds()
        for dt in self._dt_xml.keys():
            utils.echo_msg('updating {}'.format(dt))

            this_xml = iso_xml(self._dt_xml[dt], timeout = 1000, read_timeout = 2000)
            charts = this_xml.xml_doc.findall('.//{*}has', namespaces = this_xml.namespaces)
            for i, chart in enumerate(charts):
                this_xml.xml_doc = chart
                title = this_xml.title()
                perc = int((float(i)/len(charts)) * 100)
                if self._verbose: utils.echo_msg_inline('scanning {} surveys in {} [{:3}%] - {}.'.format(len(charts), dt, perc, title))
                self.FRED.layer.SetAttributeFilter("ID = 'CHARTS-{}'".format(title))
                if len(self.FRED.layer) == 0:
                    h_epsg, v_epsg = this_xml.reference_system()
                    this_data = this_xml.linkages()
                    geom = this_xml.polygon(geom=True)
                    if geom is not None:
                        self._surveys.append([{'Name': title,
                                               'ID': 'CHARTS-{}'.format(title), 
                                               'Date': this_xml.date(), 
                                               'MetadataLink': this_xml.url,
                                               'MetadataDate': this_xml.xml_date(), 
                                               'DataLink': this_data,
                                               'IndexLink': '',
                                               'DataType': dt,
                                               'DataSource': 'charts',
                                               'HorizontalDatum': h_epsg,
                                               'VerticalDatum': v_epsg,
                                               'LastUpdate': utils.this_date(),
                                               'Info': this_xml.abstract()}, geom.ExportToJson()])
                    
            if self._verbose: utils.echo_msg('scanning {} surveys in {} [ OK ].\n'.format(len(charts), dt))                                    
        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['charts']))

    def _parse_results(self, r):
        '''Search for data in the reference vector file'''

        for surv in r:
            fldata = surv['DataLink'].split(',')
            [self._data_urls.append([i, i.split('/')[-1], surv['DataType']]) if i != '' else None for i in fldata]
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None):
        if xyzc is None: xyzc = xyzfun._xyz_config
        xyzc['z-scale'] = -1
        xyzc['warp'] = epsg
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
                                xyzfun.xyz_line([float(x) for x in xyz], o_xyz, False)

                ds_ogr = layer_s = None
                ch_f_r = dst_xyz

                if os.path.exists(ch_f_r):
                    with open(ch_f_r, 'r') as in_c:
                        for xyz in xyzfun.xyz_parse(in_c, xyz_c = xyzc, verbose = self._verbose):
                            yield(xyz)

                utils.remove_glob(src_ch)
                utils.remove_glob(dst_xyz)
                utils._clean_zips(src_zips)
        utils.remove_glob(src_zip)
        
    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

## =============================================================================
##
## NCEI THREDDS Catalog (CUDEM)
##
## =============================================================================
class ncei_thredds:
    '''Fetch DEMs from NCEI THREDDS Catalog'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading NCEI_THREDDS fetch module...')
        
        ## ==============================================
        ## NOS urls and directories
        ## ==============================================
        self._nt_catalog = "https://www.ngdc.noaa.gov/thredds/demCatalog.xml"
        self._ngdc_url = "https://www.ngdc.noaa.gov"
        self._outdir = os.path.join(os.getcwd(), 'ncei_thredds')
        
        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = ''
        
        self._info = ''

    def _parse_catalog(self, catalog_url):

        ntCatalog = iso_xml(catalog_url)
        ntCatRefs = ntCatalog.xml_doc.findall('.//th:catalogRef', namespaces = ntCatalog.namespaces)
        for ntCatRef in ntCatRefs:
            ntCatHref =ntCatRef.attrib['{http://www.w3.org/1999/xlink}href']
            if ntCatHref[0] == "/":
                ntCatUrl = '{}{}'.format(self._ngdc_url, ntCatHref)
            else: ntCatUrl = '{}/{}'.format(os.path.dirname(catalog_url), ntCatHref)
            self._parse_dataset(ntCatUrl)
            
    def _parse_dataset(self, catalog_url):
        ntCatXml = iso_xml(catalog_url)
        this_ds = ntCatXml.xml_doc.findall('.//th:dataset', namespaces = ntCatXml.namespaces)
        this_ds_services = ntCatXml.xml_doc.findall('.//th:service', namespaces = ntCatXml.namespaces)
        for i, node in enumerate(this_ds):
            this_title = node.attrib['name']
            this_id = node.attrib['ID']

            perc = int((float(i)/len(this_ds)) * 100)
            if self._verbose: utils.echo_msg_inline('scanning {} datasets in {} [{:3}%] - {}.'.format(len(this_ds), this_ds[0].attrib['name'], perc, this_title))
            
            self.FRED.layer.SetAttributeFilter("ID = '{}'".format(this_id))
            if len(self.FRED.layer) == 0:            
                subCatRefs = node.findall('.//th:catalogRef', namespaces = ntCatXml.namespaces)
                if len(subCatRefs) > 0:
                    self._parse_catalog(catalog_url)
                    break
                try:
                    ds_path = node.attrib['urlPath']
                except: continue

                iso_url = False
                wcs_url = False
                http_url = False
                for service in this_ds_services:
                    service_name = service.attrib['name']
                    if service_name == 'iso': iso_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                    if service_name == 'wcs': wcs_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)
                    if service_name == 'http': http_url = '{}{}{}'.format(self._ngdc_url, service.attrib['base'], ds_path)

                this_xml = iso_xml(iso_url)
                title = this_xml.title()
                h_epsg, v_epsg = this_xml.reference_system()
                #this_data = this_xml.linkages()
                geom = this_xml.bounds(geom=True)

                zv = this_xml.xml_doc.findall('.//gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString', namespaces = this_xml.namespaces)
                if zv is not None:
                    for zvs in zv:
                        if zvs.text == 'bathy' or zvs.text == 'Band1' or zvs.text == 'z':
                            zvar = zvs.text
                        else: zvar = 'z'
                #zvar = 'z' if zv is None else zv.text
                if geom is not None:
                    self._surveys.append([{'Name': title,
                                           'ID': this_id, 
                                           'Date': this_xml.date(), 
                                           'MetadataLink': this_xml.url,
                                           'MetadataDate': this_xml.xml_date(), 
                                           'DataLink': http_url,
                                           'IndexLink': wcs_url,
                                           'DataType': zvar,
                                           'DataSource': 'ncei_thredds',
                                           'HorizontalDatum': h_epsg,
                                           'VerticalDatum': v_epsg,
                                           'LastUpdate': utils.this_date(),
                                           'Info': this_xml.abstract()}, geom.ExportToJson()])
                    
        if self._verbose: utils.echo_msg('scanning {} datasets in {} [OK]'.format(len(this_ds), this_ds[0].attrib['name']))
            
    def _update(self):
        self.FRED._open_ds()
        self._parse_catalog(self._nt_catalog)
        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
    
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['ncei_thredds']))

    def _parse_results(self, r):
        '''Search for data in the reference vector file'''

        for surv in r:
            fldata = surv['DataLink'].split(',')
            wcs_url = "{}?request=GetCoverage&version=1.0.0&service=WCS&coverage={}&bbox={}&format=NetCDF3"\
                .format(surv['IndexLink'], surv['DataType'], regions.region_format(self.region, 'bbox'))
            #[self._data_urls.append([i, i.split('/')[-1], surv['DataType']]) if i != '' else None for i in fldata]
            self._data_urls.append([surv['IndexLink'], surv['DataLink'].split('/')[-1], surv['DataType']])
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        try:
            src_ds = gdal.Open(entry[0])
        except Exception as e:
            fetch_file(entry[0], src_dc, callback = lambda: False, verbose = False)
            try:
                src_ds = gdal.Open(src_dc)
            except Exception as e:
                utils.echo_error_msg('could not read raster file: {}, {}'.format(entry[0], e))
                src_ds = None
        except Exception as e:
            utils.echo_error_msg('could not read raster file: {}, {}'.format(entry[0], e))
            src_ds = None

        if src_ds is not None:
            srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.gdal_region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = True):
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_dc)

    def _dump_xyz(self, entry, epsg = 4326, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading Multibeam fetch module...')

        ## ==============================================
        ## NCEI MB URLs and directories 
        ## ==============================================
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_metadata_url = "https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/Multibeam/iso/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._mb_autogrid = "https://www.ngdc.noaa.gov/maps/autogrid/"
        self._outdir = os.path.join(os.getcwd(), 'mb')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._stop = callback

        self._desc = '''NOAA MULTIBEAM survey data
        Multibeam Bathymetric Data from NOAA.
        In addition to deepwater data, the Multibeam Bathymetry Database (MBBDB) includes hydrographic multibeam survey data from NOAA's National Ocean Service (NOS).

        https://www.ngdc.noaa.gov/mgg/bathymetry/multibeam.html

        < mb >'''

        self.info = "NCEI is the U.S. national archive for multibeam bathymetric data and holds more than 9 million nautical miles of ship trackline data recorded from over 2400 cruises and received from sources worldwide. In addition to deepwater data, the Multibeam Bathymetry Database (MBBDB) includes hydrographic multibeam survey data from NOAA's National Ocean Service (NOS)."
        
    def _update(self):

        self.FRED._open_ds()
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('MB-1'))
        if len(self.FRED.layer) == 0:
            _surveys = [[
                {
                    'Name': 'NOAA Multibeam', 
                    'ID': 'MB-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': self._mb_metadata_url,
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._mb_search_url,
                    'IndexLink': '',
                    'DataType': 'multibeam',
                    'DataSource': 'mb',
                    'HorizontalDatum': 4326,
                    'VerticalDatum': 1092,
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-180,180,-90,90]).ExportToJson(),
            ]]
            self.FRED._add_surveys(_surveys)
        self.FRED._close_ds()
        return(_surveys)

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['mb']))
        
    def _parse_results(self, r, processed = False):
        '''Run the MB (multibeam) fetching module.'''
        these_surveys = {}
        these_versions = {}
        for surv in r:
            if self.region is None: return([])        
            _req = fetch_req(surv['DataLink'], params = {'geometry': regions.region_format(self.region, 'bbox')}, timeout = 20)
            print(_req.url)
            if _req is not None:
                survey_list = _req.text.split('\n')[:-1]
                for r in survey_list:
                    dst_pfn = r.split(' ')[0]
                    dst_fn = dst_pfn.split('/')[-1:][0]
                    survey = dst_pfn.split('/')[6]
                    dn = r.split(' ')[0].split('/')[:-1]
                    version = dst_pfn.split('/')[9][-1]
                    data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                    if survey in these_surveys.keys():
                        if version in these_surveys[survey].keys():
                            these_surveys[survey][version].append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])
                        else: these_surveys[survey][version] = [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]
                    else: these_surveys[survey] = {version: [[data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb']]}
                    
        for key in these_surveys.keys():
            if '2' in these_surveys[key].keys():
                for v2 in these_surveys[key]['2']:
                    self._data_urls.append(v2)
                #print(these_surveys[key]['2'])
            else:
                for v1 in these_surveys[key]['1']:
                    self._data_urls.append(v1)
                    #self._data_urls.append(these_surveys[key]['1'])
            #self._data_urls.append([data_url.split(' ')[0], '/'.join([survey, dst_fn]), 'mb'])

        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None):
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_mb = os.path.basename(entry[1])

        if fetch_file(entry[0], src_mb, callback = self._stop, verbose = False) == 0:
            src_xyz = os.path.basename(src_mb).split('.')[0] + '.xyz'
            out, status = utils.run_cmd('mblist -MX20 -OXYZ -I{}  > {}'.format(src_mb, src_xyz), verbose = False)
            xyzc['name'] = src_mb
            xyzc['delim'] = '\t'
            xyzc['z-scale'] = 1
            xyzc['epsg'] = 4326
            xyzc['warp'] = epsg
            mb_r = src_xyz

            with open(mb_r, 'r') as in_m:
                for xyz in xyzfun.xyz_parse(in_m, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
            utils.remove_glob(src_xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_mb))
        utils.remove_glob(src_mb)

    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading USACE fetch module...')

        ## ==============================================
        ## USACE eHydro URLs and directories
        ## ==============================================    
        self._usace_gj_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'
        self._outdir = os.path.join(os.getcwd(), 'usace')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._stop = callback
        
        self._desc = '''USACE bathymetry surveys via eHydro
        Bathymetric Channel surveys from USACE - U.S. only.
        The hydrographic surveys provided by this application are to be used for informational purposes only and should not be used as a navigational aid. 
        Channel conditions can change rapidly and the surveys may or may not be accurate.

        https://navigation.usace.army.mil/Survey/Hydro

        < usace >'''

        self.info = ''

    def _update(self):
        self.FRED._open_ds()
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('USACE-1'))
        if len(self.FRED.layer) == 0:
            _surveys = [[
                {
                    'Name': 'USACE E-Hydro', 
                    'ID': 'USACE-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': '',
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._usace_gs_api_url,
                    'IndexLink': self._usace_gj_url,
                    'DataType': 'usace',
                    'DataSource': 'usace',
                    'HorizontalDatum': 'varies',
                    'VerticalDatum': 'varies',
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-162, -60, 16, 73]).ExportToJson(),
            ]]
            self.FRED._add_surveys(_surveys)
        self.FRED._close_ds()
        return(_surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['usace']))
        
    def _parse_results(self, r, stype = None):
        '''Run the USACE fetching module'''

        for surv in r:
            _data = {
                'geometry': regions.region_format(self.region, 'bbox'),
                'inSR':4326,
                'outSR':4326,
                'f':'pjson',
            }
            self._req = fetch_req(surv['DataLink'], params = _data)
            if self._req is not None:
                survey_list = self._req.json()
                for feature in survey_list['features']:
                    fetch_fn = feature['attributes']['SOURCEDATALOCATION']
                    if stype is not None:
                        if feature['attributes']['SURVEYTYPE'].lower() == stype.lower():
                            self._data_urls.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                    else: self._data_urls.append([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None):
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_zip = os.path.basename(entry[1])
        xyzc['warp'] = epsg
        
        if fetch_file(entry[0], src_zip, callback = self._stop, verbose = False) == 0:

            ## attempt to discover the state plane zone from xml
            src_xmls = utils.p_unzip(src_zip, ['.xml', '.XML'])
            for src_xml in src_xmls:
                this_xml = lxml.etree.parse(src_xml)
                if this_xml is not None:
                    try:
                        w = this_xml.find('.//westbc').text
                        e = this_xml.find('.//eastbc').text
                        n = this_xml.find('.//northbc').text
                        s = this_xml.find('.//southbc').text
                        this_xml_region = [float(w), float(e), float(s), float(n)]
                    except:
                        utils.echo_error_msg('could not determine survey bb from {}'.format(src_xml))
                        this_xml_region = self.region
                    try:
                        prj = this_xml.find('.//gridsysn').text
                        szone = this_xml.find('.//spcszone').text
                        utils.echo_msg('zone: {}'.format(szone))                        
                        xyzc['epsg'] = int(FIPS_TO_EPSG[szone])
                    except:
                        utils.echo_error_msg('could not determine state plane zone from {}'.format(src_xml))
                        xyzc['epsg'] = None
                utils.remove_glob(src_xml)

            ## otherwise attempt to discover state plane zone from vector
            if xyzc['epsg'] is None:
                this_geom = gdalfun.gdal_region2geom(this_xml_region)
                sp_fn = os.path.join(fetchdata, 'stateplane.geojson')
                sp = ogr.Open(sp_fn)
                layer = sp.GetLayer()
                
                for feature in layer:
                    geom = feature.GetGeometryRef()
                    if this_geom.Intersects(geom):
                        xyzc['epsg'] = feature.GetField('EPSG')
                sp = None
                        
            src_usaces = utils.p_unzip(src_zip, ['.XYZ', '.xyz', '.dat'])
            for src_usace in src_usaces:
                if os.path.exists(src_usace):
                    with open(src_usace, 'r') as in_c:
                        xyzc['name'] = src_usace.split('.')[0]
                        for xyz in xyzfun.xyz_parse(in_c, xyz_c = xyzc, verbose = self._verbose):
                            yield(xyz)
                    utils.remove_glob(src_usace)

            #utils._clean_zips(src_zips)
        
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
        utils.remove_glob(src_zip)

    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading The National Map fetch module...')
        
        ## ==============================================
        ## TNM URLs and directories
        ## ==============================================
        self._tnm_api_url = 'http://tnmaccess.nationalmap.gov/api/v1'
        self._tnm_dataset_url = 'https://tnmaccess.nationalmap.gov/api/v1/datasets?'
        self._tnm_product_url = 'https://tnmaccess.nationalmap.gov/api/v1/products?'
        self._tnm_meta_base = 'https://www.sciencebase.gov/catalog/item/'
        self._outdir = os.path.join(os.getcwd(), 'tnm')
        self._tnm_ds = []
        self._tnm_df = []

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = ''

        self.info = ''

    def _datasets(self):
        _req = fetch_req(self._tnm_dataset_url)
        if _req is not None:
            try:
                _datasets = _req.json()
            except Exception as e:
                utils.echo_error_msg('try again, {}'.format(e))
        else: _datasets = None

        return(_datasets)

    def _parse_dsTag(self, dsTag):
        sbDTag = dsTag['sbDatasetTag']
        meta = '{}{}?format=iso'.format(self._tnm_meta_base, dsTag['id'])

        this_xml = iso_xml(meta)
        h_epsg, v_epsg = this_xml.reference_system()
        if h_epsg is None: h_epsg = 'varies'
        if v_epsg is None: v_epsg = 'varies'
        _data = {
            'max': 10000,
            'datasets': sbDTag,
        }
        try:
            url_enc = urllib.urlencode(_data)
        except: url_enc = urllib.parse.urlencode(_data)
        data_link = '{}{}'.format(self._tnm_product_url, url_enc)
        this_geom = this_xml.bounds(geom=True)
        this_date = this_xml.date()
        
        if this_geom is not None:
            return([{'Name': dsTag['title'],
                     'ID': sbDTag, 
                     'Date': this_date,
                     'MetadataLink': meta,
                     'MetadataDate': this_xml.xml_date(), 
                     'DataLink': data_link,
                     'IndexLink': '',
                     'DataSource': 'tnm',
                     'HorizontalDatum': h_epsg,
                     'VerticalDatum': v_epsg,
                     'LastUpdate': utils.this_date(),
                     'Info': this_xml.abstract()}, this_geom.ExportToJson()])
        else: return(None)

    def _parse_ds(self, ds):
        _surveys = []
        try:
            dsTags = ds['tags']
        except: dsTags = []
        if len(dsTags) > 0:
            for tag in dsTags:
                _surveys.append(self._parse_ds(tag))
        else: _surveys.append(self._parse_dsTag(ds))
        return(_surveys)
        
    def _update(self):
        _surveys = []
        _s = {}
        self.FRED._open_ds()
        tnm_ds = self._datasets()
        for i, ds in enumerate(tnm_ds):
            formats = ''
            for f in ds['formats']:
                if 'isDefault' in f.keys():
                    formats = f['value']
                    break
            mapServerUrl = ds['mapServerUrl']
            _s = self._parse_ds(ds)
            for _survey in _s:
                if _survey is not None:
                    #print(_survey)
                    if _survey[0] is not None:
                        if len(_survey) == 1:
                            _survey = _survey[0]
                        _survey[0]['DataType'] = formats
                        _survey[0]['IndexLink'] = mapServerUrl
                        self.FRED.layer.SetAttributeFilter("ID = '{}'".format(ds['sbDatasetTag']))
                        if len(self.FRED.layer) == 0:
                            self._surveys.append(_survey)
                
        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(_surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['tnm']))

    def _parse_results(self, r, f = None, e = None, q = None):
        if self._verbose: utils.echo_msg('filtering TNM dataset results...')
        for surv in r:
    
            _data = {
                'bbox': regions.region_format(self.region, 'bbox'),
            }
            if q is not None: _data['q'] = str(q)
            if f is None:
                _data['prodFormats'] = surv['DataType']
            else: _data['prodFormats'] = ','.join(f)
            if e is None: e = []
            
            _req = fetch_req(surv['DataLink'], params = _data)
            _dataset_results = []

            if _req is not None:
                try:
                    _dataset_results = _req.json()
                except ValueError: utils.echo_error_msg('tnm server error, try again')
                except Exception as e: utils.echo_error_msg('error, {}'.format(e))                

            if len(_dataset_results) > 0:
                for item in _dataset_results['items']:
                    if len(e) > 0:
                        for extent in e:
                            if item['extent'] == extent:
                                f_url = item['downloadURL']
                                if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                                    tnm_ds = 'ned'
                                elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                                    tnm_ds = 'lidar'
                                else: tnm_ds = 'tnm'
                                self._data_urls.append([f_url, f_url.split('/')[-1], tnm_ds])
                    else:
                        f_url = item['downloadURL']
                        if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                            tnm_ds = 'ned'
                        elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                            tnm_ds = 'lidar'
                        else: tnm_ds = 'tnm'
                        self._data_urls.append([f_url, f_url.split('/')[-1], tnm_ds])
        return(self._data_urls)
    
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================                
    def _yield_xyz(self, entry, epsg = None):
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
                    srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.gdal_region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
                    for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose):
                        if xyz[2] != 0:
                            yield(xyz)
                    src_ds = None

                utils.remove_glob(src_tnm)
                utils._clean_zips(src_zips)
        utils.remove_glob(entry[1])

    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)
            
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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading GMRT fetch module...')

        ## ==============================================
        ## GMRT URLs and directories
        ## ==============================================    
        self._gmrt_grid_url = "https://www.gmrt.org:443/services/GridServer?"
        self._gmrt_grid_urls_url = "https://www.gmrt.org:443/services/GridServer/urls?"
        self._gmrt_grid_metadata_url = "https://www.gmrt.org/services/GridServer/metadata?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''The Global Multi-Reosolution Topography Data Synthesis (GMRT) 
        Global Elevation raster dataset via GMRT.

        < gmrt >'''
        
        self.info = 'The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout the global and coastal oceans.'
        
    def _update(self):
        self.FRED._open_ds()
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('GMRT-1'))
        if len(self.FRED.layer) == 0:
            self._surveys = [[
                {
                    'Name': 'GMRT', 
                    'ID': 'GMRT-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': self._gmrt_grid_metadata_url,
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._gmrt_grid_urls_url,
                    'IndexLink': '',
                    'DataType': 'raster',
                    'DataSource': 'gmrt',
                    'HorizontalDatum': 3857,
                    'VerticalDatum': 1092,
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-180,180,-90,90]).ExportToJson(),
            ]]
            self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['gmrt']))
    
    def _parse_results(self, r, fmt = 'geotiff', res = 'max', layer = 'topo'):
        '''Run the GMRT fetching module'''
        if layer != 'topo' and layer != 'topo-mask': layer = 'topo'
        for surv in r:
            _data = {
                'north':self.region[3],
                'west':self.region[0],
                'south':self.region[2],
                'east':self.region[1],
                'mformat':'json',
                'resolution':res,
                'format':fmt,
            }
            #                'layer': layer,
            _req = fetch_req(surv['DataLink'], params = _data, tries = 10, timeout = 2)
            if _req is not None:
                gmrt_urls = _req.json()
                for url in gmrt_urls:
                    opts = {}
                    url_base = url.split('?')[0]
                    for url_opt in url.split('?')[1].split('&'):
                        opt_kp = url_opt.split('=')
                        opts[opt_kp[0]] = opt_kp[1]

                    opts['layer'] = layer
                    try:
                        url_enc = urllib.urlencode(opts)
                    except: url_enc = urllib.parse.urlencode(opts)

                    this_url = '{}?{}'.format(url_base, url_enc)
                    url_region = [float(opts['west']), float(opts['east']), float(opts['south']), float(opts['north'])]
                    outf = 'gmrt_{}_{}.{}'.format(opts['layer'], regions.region_format(url_region, 'fn'), gdalfun.gdal_fext(opts['format']))
                    self._data_urls.append([this_url, outf, 'gmrt'])
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is a an item from self._data_urls
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326):
        src_gmrt = 'gmrt_tmp.tif'
        if fetch_file(entry[0], src_gmrt, callback = lambda: False, verbose = False) == 0:
            try:
                src_ds = gdal.Open(src_gmrt)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose, warp = epsg):
                    yield(xyz)
            src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_gmrt))
        utils.remove_glob(src_gmrt)

    def _dump_xyz(self, src_gmrt, epsg = 4326, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_gmrt, epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading mar_grav fetch module...')

        ## ==============================================
        ## marine gravity URLs and directories
        ## ==============================================
        self._mar_grav_info_url = 'https://topex.ucsd.edu/WWW_html/mar_topo.html'
        self._mar_grav_url = 'https://topex.ucsd.edu/cgi-bin/get_data.cgi'        
        self._outdir = os.path.join(os.getcwd(), 'mar_grav')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''Marine Gravity from Sattelite Altimetry topographic data.
        Elevation data from Scripps Marine Gravity dataset.

        < mar_grav >'''
        
        self.info = ''
        
    def _update(self):

        self.FRED._open_ds()
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('MAR_GRAV-1'))
        if len(self.FRED.layer) == 0:
            self._surveys = [[
                {
                    'Name': 'MAR_GRAV', 
                    'ID': 'MAR_GRAV-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': self._mar_grav_info_url,
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._mar_grav_url,
                    'IndexLink': '',
                    'DataType': 'xyz',
                    'DataSource': 'mar_grav',
                    'HorizontalDatum': 4326,
                    'VerticalDatum': 1092,
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-180,180,-90,90]).ExportToJson(),
            ]]
            self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['mar_grav']))
        
    def _parse_results(self, r):
        '''Run the mar_grav fetching module.'''

        for surv in r:        
            _data = {
                'north':self.region[3],
                'west':self.region[0],
                'south':self.region[2],
                'east':self.region[1],
                'mag':1,
            }

            try:
                url_enc = urllib.urlencode(_data)
            except: url_enc = urllib.parse.urlencode(_data)
            
            url = '{}?{}'.format(surv['DataLink'], url_enc)
            outf = 'mar_grav_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            self._data_urls.append([url, outf, 'mar_grav'])
        return(self._data_urls)
        
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = lambda: False, verbose = self._verbose) == 0:
            xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            xyzc['x-off'] = -360
            xyzc['verbose'] = True
            xyzc['warp'] = epsg
            xyzc['name'] = '<mar_grav data-stream>'
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
    
    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

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
    '''Fetch SRTM15+ data'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading SRTM+ fetch module...')

        ## ==============================================
        ## SRTM URLs and directories
        ## ==============================================
        self._srtm_info_url = 'https://topex.ucsd.edu/WWW_html/srtm15_plus.html'
        self._srtm_url = 'https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi'
        self._outdir = os.path.join(os.getcwd(), 'srtm_plus')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''SRTM15+ elevation data (Scripps)
        Global Bathymetry and Topography at 15 Arc Sec:SRTM15+
    
        < srtm_plus >'''
        
        self.info = ''
        
    def _update(self):
        self.FRED._open_ds()
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('SRTM-1'))
        if len(self.FRED.layer) == 0:
            self._surveys = [[
                {
                    'Name': 'SRTM15+', 
                    'ID': 'SRTM-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': self._srtm_info_url,
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._srtm_url,
                    'IndexLink': '',
                    'DataType': 'xyz',
                    'DataSource': 'srtm_plus',
                    'HorizontalDatum': 4326,
                    'VerticalDatum': 1092,
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-180,180,-90,90]).ExportToJson(),
            ]]
            self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['srtm_plus']))
        
    def _parse_results(self, r):
        '''Run the srtm+ fetching module.'''

        for surv in r:        
            _data = {
                'north':self.region[3],
                'west':self.region[0],
                'south':self.region[2],
                'east':self.region[1],
            }

            try:
                url_enc = urllib.urlencode(_data)
            except: url_enc = urllib.parse.urlencode(_data)
            
            url = '{}?{}'.format(surv['DataLink'], url_enc)
            outf = 'srtm_plus_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            self._data_urls.append([url, outf, 'mar_grav'])
        return(self._data_urls)
        
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = lambda: False, verbose = self._verbose) == 0:
            xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            xyzc['x-off'] = -360
            xyzc['verbose'] = True
            xyzc['warp'] = epsg
            xyzc['name'] = '<mar_grav data-stream>'
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
    
    def _dump_xyz(self, entry, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

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
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading EMODNET fetch module...')

        ## ==============================================
        ## EMODNET URLs and directories
        ## ==============================================    
        self._emodnet_grid_url = 'https://ows.emodnet-bathymetry.eu/wcs?'
        self._emodnet_grid_cap = 'https://ows.emodnet-bathymetry.eu/wms?request=GetCapabilities&service=WMS'
        self._outdir = os.path.join(os.getcwd(), 'emodnet')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''EMODNET Elevation Data
        European Bathymetry/Topographic data from EMODNET

        https://portal.emodnet-bathymetry.eu/help/help.html

        < emodnet >'''
                
        self.info = ''

    def _update(self):

        self.FRED._open_ds()
        xml_cap = iso_xml(self._emodnet_grid_cap).xml_doc.find('{http://www.opengis.net/wms}Capability')
        layer = xml_cap.find('{http://www.opengis.net/wms}Layer')
        layers = layer.findall('{http://www.opengis.net/wms}Layer')
        for l in layers:
            name = l.find('{http://www.opengis.net/wms}Name').text
            if name != 'emodnet:mean': continue
            
            self.FRED.layer.SetAttributeFilter("ID = '{}'".format('{}'.format(name)))
            if len(self.FRED.layer) == 0:
                xml = l.find('{http://www.opengis.net/wms}MetadataURL')
                xml_res = xml.find('{http://www.opengis.net/wms}OnlineResource')
                xml_url = xml_res.get('{http://www.w3.org/1999/xlink}href')
                
                title = l.find('{http://www.opengis.net/wms}Title').text
                abstract = l.find('{http://www.opengis.net/wms}Abstract').text
                bb = l.find('{http://www.opengis.net/wms}EX_GeographicBoundingBox')
                region = [float(x) for x in [bb[0].text, bb[1].text, bb[2].text, bb[3].text]]
                try:
                    dl = urllib.urlencode({'coverage': name, 'Identifier': name})
                except: dl = urllib.parse.urlencode({'coverage': name, 'Identifier': name})
                
                kws = l.find('{http://www.opengis.net/wms}KeywordList')
                dt = 'unknown' if len(kws) == 0 else kws[-1].text
                geom = gdalfun.gdal_region2geom(region)
                if geom is not None:
                    self._surveys.append([{'Name': title, 
                                           'ID': '{}'.format(name), 
                                           'Date': utils.this_year(), 
                                           'MetadataLink': xml_url,
                                           'MetadataDate': utils.this_year(), 
                                           'DataLink': self._emodnet_grid_url,
                                           'IndexLink': '',
                                           'DataType': dt,
                                           'DataSource': 'emodnet',
                                           'HorizontalDatum': 4326,
                                           'VerticalDatum': 1092,
                                           'LastUpdate': utils.this_date(),
                                           'Info': abstract}, geom.ExportToJson()])

        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['emodnet']))
        
    def _parse_results(self, r):
        '''Run the EMODNET fetching module'''

        for surv in r:
            desc_data = {
                'request': 'DescribeCoverage',
                'CoverageID': surv['ID'],
                'version': '2.0.1',
                'service': 'WCS',
                }

            desc_req = fetch_req(surv['DataLink'], params = desc_data)
            desc_results = lxml.etree.fromstring(desc_req.text.encode('utf-8'))
            g_env = desc_results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope')[0]
            hl = [float(x) for x in g_env.find('{http://www.opengis.net/gml/3.2}high').text.split()]

            g_bbox = desc_results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
            lc = [float(x) for x in  g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split()]
            uc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split()]

            ds_region = [lc[1], uc[1], lc[0], uc[0]]
            resx = (uc[1] - lc[1]) / hl[0]
            resy = (uc[0] - lc[0]) / hl[1]

            if regions.regions_intersect_ogr_p(self.region, ds_region):
                data = {
                    'request': 'GetCoverage',
                    'version': '1.0.0',
                    'service': 'WCS',
                    'resx': resx,
                    'resy': resy,
                    'crs': 'EPSG:4326',
                    'format': 'GeoTIFF',
                    'coverage': surv['ID'],
                    'Identifier': surv['ID'],
                    'bbox': regions.region_format(self.region, 'bbox'),
                }

                try:
                    enc_data = urllib.urlencode(data)
                except: enc_data = urllib.parse.urlencode(data)
                emodnet_wcs = '{}{}'.format(surv['DataLink'], enc_data)            
                outf = 'emodnet_{}.tif'.format(regions.region_format(self.region, 'fn'))
                self._data_urls.append([emodnet_wcs, outf, 'emodnet'])

        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None):
        src_emodnet = 'emodnet_tmp.tif'
        if fetch_file(entry[0], src_emodnet, callback = lambda: False, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_emodnet)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
        utils.remove_glob(src_emodnet)

    def _dump_xyz(self, src_emodnet, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_emodnet, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

## =============================================================================
##
## HRDEM Fetch - Canada High Resolution DEM dataset
##
## Fetch Canadian HRDEM data.
## https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6
##
## =============================================================================
class hrdem():
    '''Fetch HRDEM data from Canada (NRCAN)'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading HRDEM fetch module...')

        ## ==============================================
        ## Reference vector and directories
        ## ==============================================
        self._hrdem_footprints_url = 'ftp://ftp.maps.canada.ca/pub/elevation/dem_mne/highresolution_hauteresolution/Datasets_Footprints.zip'
        self._local_ref_vector = 'hrdem.gmt'
        if not os.path.exists(self._local_ref_vector):
            self.FRED = os.path.join(fetchdata, 'hrdem.gmt')
        else: self.FRED = self._local_ref_vector
        self._outdir = os.path.join(os.getcwd(), 'hrdem')
        
        self._has_vector = True if os.path.exists(self.FRED) else False
        if not self._has_vector: self.FRED = 'hrdem.gmt'

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = ''
                
        self.info = ''

    def _update(self):
        self.FRED._open_ds()
        v_zip = os.path.basename(self._hrdem_footprints_url)
        status = fetch_ftp_file(self._hrdem_footprints_url, v_zip, verbose = True)
        v_shps = utils.p_unzip(v_zip, ['.shp', '.shx', '.dbf', '.prj'])
        v_shp = None
        for v in v_shps:
            if '.shp' in v: v_shp = v
        try:
            v_ds = ogr.Open(v_shp)
        except: v_ds = None
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            for f in range(0, fcount):
                feature = layer[f]
                name = feature.GetField('Tile_name')
                perc = (f/fcount) * 100.
                if self._verbose: utils.echo_msg_inline('scanning {} datasets [{:3}%] - {}.'.format(fcount, perc, name))
                self.FRED.layer.SetAttributeFilter("Name = '{}'".format('{}'.format(name)))
                if len(self.FRED.layer) == 0:                
                    data_link = feature.GetField('Ftp_dtm')
                    if data_link is not None:
                        geom = feature.GetGeometryRef()
                        if geom is not None:
                            self._surveys.append([{'Name': name, 
                                                   'ID': feature.GetField('Project'), 
                                                   'Date': utils.this_year(), 
                                                   'MetadataLink': feature.GetField('Meta_dtm'),
                                                   'MetadataDate': utils.this_year(), 
                                                   'DataLink': data_link.replace('http', 'ftp'),
                                                   'IndexLink': self._hrdem_footprints_url,
                                                   'DataType': 'dtm',
                                                   'DataSource': 'hrdem',
                                                   'HorizontalDatum': feature.GetField('Coord_Sys').split(':')[-1],
                                                   'VerticalDatum': 'varies',
                                                   'LastUpdate': utils.this_date(),
                                                   'Info': feature.GetField('Provider')}, geom.ExportToJson()])

            self.FRED._add_surveys(self._surveys)
            if self._verbose: utils.echo_msg_inline('scanning {} datasets [OK]'.format(fcount))
        [utils.remove_glob(v) for v in v_shps]
        utils.remove_glob(v_zip)
        self.FRED._close_ds()
        return(self._surveys)

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['hrdem']))

    def _parse_results(self, r):
        '''Search for data in the reference vector file'''

        for surv in r:
            fldata = surv['DataLink'].split(',')
            [self._data_urls.append([i, i.split('/')[-1], surv['DataType']]) if i != '' else None for i in fldata]
        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326):
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
            srcwin = gdalfun.gdal_srcwin(src_ds, gdalfun.gdal_region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = True):
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_dc)

    def _dump_xyz(self, entry, epsg = 4326, dst_port = sys.stdout):
        for xyz in self._yield_xyz(entry, epsg = epsg):
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
    '''Fetch raster data from CHS'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading CHS fetch module...')

        ## ==============================================
        ## CHS URLs and directories
        ## ==============================================
        #self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._chs_grid_cap = 'https://data.chs-shc.ca/geoserver/wcs?request=GetCapabilities&service=WMS'
        self._outdir = os.path.join(os.getcwd(), 'chs')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = ''''''
                
        self.info = ''

    def _update(self):

        self.FRED._open_ds()
        xml_cap = iso_xml(self._chs_grid_cap).xml_doc.find('{http://www.opengis.net/wms}Capability')
        layer = xml_cap.find('{http://www.opengis.net/wms}Layer')
        layers = layer.findall('{http://www.opengis.net/wms}Layer')
        for l in layers:
            name = l.find('{http://www.opengis.net/wms}Name').text

            self.FRED.layer.SetAttributeFilter("ID = '{}'".format('{}'.format(name)))            
            if len(self.FRED.layer) == 0:
                xml_url = ''
                title = l.find('{http://www.opengis.net/wms}Title').text
                abstract = l.find('{http://www.opengis.net/wms}Abstract').text
                bb = l.find('{http://www.opengis.net/wms}EX_GeographicBoundingBox')
                region = [float(x) for x in [bb[0].text, bb[1].text, bb[2].text, bb[3].text]]
                try:
                    dl = urllib.urlencode({'coverage': name, 'Identifier': name})
                except: dl = urllib.parse.urlencode({'coverage': name, 'Identifier': name})
                
                kws = l.find('{http://www.opengis.net/wms}KeywordList')
                dt = 'unknown' if len(kws) == 0 else kws[-1].text
                if dt == 'CarisTiles': continue
                geom = gdalfun.gdal_region2geom(region)
                if geom is not None:
                    self._surveys.append([{'Name': title, 
                                           'ID': '{}'.format(name), 
                                           'Date': utils.this_year(), 
                                           'MetadataLink': xml_url,
                                           'MetadataDate': utils.this_year(), 
                                           'DataLink': self._chs_url,
                                           'IndexLink': '',
                                           'DataType': dt,
                                           'DataSource': 'chs',
                                           'HorizontalDatum': 4326,
                                           'VerticalDatum': 1092,
                                           'LastUpdate': utils.this_date(),
                                           'Info': abstract}, geom.ExportToJson()])
                    
        self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['chs']))
        
    def _parse_results(self, r):
        '''Run the CHS fetching module'''

        for surv in r:
            desc_data = {
                'request': 'DescribeCoverage',
                'CoverageID': surv['ID'],
                'version': '2.0.1',
                'service': 'WCS',
                }

            desc_req = fetch_req(surv['DataLink'], params = desc_data)
            desc_results = lxml.etree.fromstring(desc_req.text.encode('utf-8'))
            g_env = desc_results.findall('.//{http://www.opengis.net/gml/3.2}GridEnvelope')[0]
            hl = [float(x) for x in g_env.find('{http://www.opengis.net/gml/3.2}high').text.split()]

            g_bbox = desc_results.findall('.//{http://www.opengis.net/gml/3.2}Envelope')[0]
            lc = [float(x) for x in  g_bbox.find('{http://www.opengis.net/gml/3.2}lowerCorner').text.split()]
            uc = [float(x) for x in g_bbox.find('{http://www.opengis.net/gml/3.2}upperCorner').text.split()]

            ds_region = [lc[1], uc[1], lc[0], uc[0]]
            resx = (uc[1] - lc[1]) / hl[0]
            resy = (uc[0] - lc[0]) / hl[1]

            if regions.regions_intersect_ogr_p(self.region, ds_region):
                data = {
                    'request': 'GetCoverage',
                    'version': '1.0.0',
                    'service': 'WCS',
                    'resx': resx,
                    'resy': resy,
                    'crs': 'EPSG:4326',
                    'format': 'GeoTIFF',
                    'coverage': surv['ID'],
                    'Identifier': surv['ID'],
                    'bbox': regions.region_format(self.region, 'bbox'),
                }

                try:
                    enc_data = urllib.urlencode(data)
                except: enc_data = urllib.parse.urlencode(data)
                emodnet_wcs = '{}{}'.format(surv['DataLink'], enc_data)            
                outf = 'chs_{}_{}.tif'.format(''.join('_'.join(surv['ID'].split()).split(':')), regions.region_format(self.region, 'fn'))
                self._data_urls.append([emodnet_wcs, outf, 'chs'])

        return(self._data_urls)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None):
        src_chs = 'chs_tmp.tif'
        if fetch_file(entry[0], src_chs, callback = lambda: False, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_chs)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
        utils.remove_glob(src_emodnet)

    def _dump_xyz(self, src_chs, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_chs, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)
            
## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS - US Only
##
## =============================================================================
class ngs:
    '''Fetch NGS monuments from NOAA'''
    
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading NGS Monument fetch module...')

        ## ==============================================
        ## NGS URLs and directories
        ## ==============================================    
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''NOAA NGS Monument Data
        Monument data from NOAA's Nagional Geodetic Survey (NGS) monument dataset.

        < ngs >'''
                
        self.info = ''

    def _update(self):

        self.FRED._open_ds(1)
        self.FRED.layer.SetAttributeFilter("ID = '{}'".format('NGS-1'))
        if len(self.FRED.layer) == 0:
            self._surveys = [[
                {
                    'Name': 'NGS Monuments', 
                    'ID': 'NGS-1', 
                    'Date': utils.this_year(), 
                    'MetadataLink': '',
                    'MetadataDate': utils.this_year(), 
                    'DataLink': self._ngs_search_url,
                    'IndexLink': '',
                    'DataType': 'raster',
                    'DataSource': 'ngs',
                    'HorizontalDatum': '',
                    'VerticalDatum': '',
                    'LastUpdate': utils.this_date(),
                    'Info': self.info,
                },
                gdalfun.gdal_region2geom([-162, -60, 16, 73]).ExportToJson(),
            ]]
            self.FRED._add_surveys(self._surveys)
        self.FRED._close_ds()
        return(self._surveys)
        
    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['ngs']))
        
    def _parse_results(self, r):
        '''Run the NGS (monuments) fetching module.'''

        for surv in r:
            _data = { 'maxlon':self.region[0],
                      'minlon':self.region[1],
                      'maxlat':self.region[3],
                      'minlat':self.region[2] }

            try:
                ngs_data = urllib.urlencode(_data)
            except: ngs_data = urllib.parse.urlencode(_data)
            ngs_url = '{}{}'.format(surv['DataLink'], ngs_data)
            self._data_urls.append([ngs_url, 'ngs_results_{}.json'.format(regions.region_format(self.region, 'fn')), 'ngs'])
            
        return(self._data_urls)

    def _yield_xyz(self, entry, epsg = None):
        src_ngs = 'ngs_tmp.json'
        r = []
        if fetch_file(entry[0], src_ngs, callback = lambda: False, verbose = self._verbose) == 0:
            with open(src_ngs, 'r') as json_file: r = json.load(json_file)
            for mm in r: yield([mm['lon'], mm['lat'], mm['geoidHt']])
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_ngs))
        utils.remove_glob(src_ngs)

    def _dump_xyz(self, src_emodnet, epsg = None, dst_port = sys.stdout):
        for xyz in self._yield_xyz(src_emodnet, epsg = epsg):
            xyzfun.xyz_line(xyz, dst_port, False)

## =============================================================================
##
## Opens Street Map data (OSM)
##
## Fetch various datasets from OSM/Overpass
## <BETA>
##
## =============================================================================
class osm:
    '''Fetch OSM data'''
    
    def __init__(self, extent = None, where = [], callback = None):
        self._verbose = True
        if self._verbose: utils.echo_msg('loading OSM fetch module...')

        ## ==============================================
        ## OSM URLs and directories
        ## ==============================================    
        self._osm_api = 'https://lz4.overpass-api.de/api/interpreter'
        self._outdir = os.path.join(os.getcwd(), 'osm')
        self.osm_types = {
            'highway': ['LINESTRING'],
            'waterway': ['LINESTRING'],
            'building': ['POLYGON'],
        }

        self.where = where
        self.region = extent
        if self.region is not None:
            self._boundsGeom = gdalfun.gdal_region2geom(self.region)
        self.FRED = FRED(verbose = self._verbose)
        self._data_urls = []
        self._surveys = []
        self._stop = callback
        
        self._desc = '''NOAA NGS Monument Data
        Monument data from NOAA's Nagional Geodetic Survey (NGS) monument dataset.

        < ngs >'''
                
        self.info = ''

    def _update(self):
        pass

    def _filter_results(self):
        return(self.FRED._filter(self.region, self.where, ['ngs']))

    def _parse_results(self):
        return(self._data_urls)
    
    def run(self, osm_type = 'all', proc = False):
        '''Run the OSM fetching module.'''
        
        if self.region is None: return([])

        if osm_type == 'all':
            for key in self.osm_types.keys():
                self._fetch_results(osm_type = key, proc = proc)
        else:
            if osm_type in self.osm_types.keys():
                self._fetch_results(osm_type = osm_type, proc = proc)
                
        return(self._results)
        
    def _fetch_results(self, osm_type = 'highway', proc = False):
        c_bbox = regions.region_format(self.region, 'osm_bbox')
        out_fn = 'osm_{}_{}'.format(osm_type, regions.region_format(self.region, 'fn'))
        
        osm_q = '''
        [out:json];
        (
        way['{}']({});
        );
        out body;
        '''.format(osm_type, c_bbox)

        osm_data = fetch_req(self._osm_api, params = {'data': osm_q}, timeout=3600)
        osm_json = osm_data.json()

        if self._verbose: utils.echo_msg('found {} {} elements'.format(len(osm_json['elements']), osm_type))
        
        if proc:
            if not os.path.exists(self._outdir):
                try:
                    os.makedirs(self._outdir)
                except: pass 

            dst_gmt = open('{}.gmt'.format(os.path.join(self._outdir, out_fn)), 'w')
            dst_gmt.write('# @VGMT1.0 @G{} @Nid\n'.format(self.osm_types[osm_type][0]))
            dst_gmt.write('# @Tdouble\n')
            dst_gmt.write('# @R{}\n'.format(regions.region_format(self.region, 'str')))
            dst_gmt.write('# FEATURE_DATA\n')

        for el_n, el in enumerate(osm_json['elements']):
            utils.echo_msg_inline('fetching osm {} data [{}%]'.format(osm_type, float(el_n) / len(osm_json['elements']) * 100))
            if el['type'] == 'way':

                osm_node_q = '''
                [out:json];
                (
                node(id:{});
                );
                out body;
                '''.format(','.join([str(n) for n in el['nodes']]))

                node_data = fetch_req(self._osm_api, params = {'data': osm_node_q})
                if node_data.status_code == 200:
                    self._results.append([node_data.url, '{}.json'.format(out_fn), 'osm'])
                    
                    if proc:
                        node_json = node_data.json()

                        dst_gmt.write('>\n# @D{}\n'.format(el['id']))

                        for node in el['nodes']:
                            for el_node in node_json['elements']:
                                if el_node['id'] == node:
                                    xyzfun.xyz_line([el_node['lon'], el_node['lat']], dst_gmt)
                                    
        if self._verbose: utils.echo_msg_inline('fetching osm {} data [OK]\n'.format(osm_type))
        
        if proc:
            dst_gmt.close()
            utils.run_cmd('ogr2ogr {}.shp {}.gmt'.format(os.path.join(self._outdir, out_fn), os.path.join(self._outdir, out_fn)), verbose = False)
            
## ==============================================
## fetches processing (datalists fmt:400 - 499)
## ==============================================
def fetch_inf_entry(entry = [], warp = None):
    out_inf = {'minmax': [-180,180,-90,90], 'name': entry[0], 'pts': None, 'wkt': gdalfun.gdal_region2wkt([-180,180,-90,90])}
    if warp is not None:
        out_inf['minmax'] = regions.region_warp(out_inf['minmax'], 4326, warp)
        out_inf['wkt'] = gdalfun.gdal_region2wkt(out_inf['minmax'])
    return(out_inf)

def fetch_yield_entry(entry = ['nos:datatype=xyz'], region = None, warp = None, verbose = False):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''

    region = regions.region_warp(region, src_epsg = warp, dst_epsg = 4326)
    fetch_mod = entry[0].split(':')[0]
    fetch_args = entry[0].split(':')[1:]

    fl = _fetch_modules[fetch_mod](regions.region_buffer(region, 5, pct = True), [], lambda: False, False)
    args_d = utils.args2dict(fetch_args, {})

    r = fl._parse_results(fl._filter_results(), **args_d)
    
    for e in r:
        for xyz in fl._yield_xyz(e, warp):
            yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dump_entry(entry = ['nos:datatype=nos'], dst_port = sys.stdout, region = None, verbose = False):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    
    for xyz in fetch_yield_entry(entry, region, verbose):
        xyz_line(xyz, dst_port, False)
            
## =============================================================================
##
## Run fetches from command-line
##
## =============================================================================
_usage = '''{} [OPTIONS] <module[:parameter=value]* ...>

Fetch geographic elevation data.

General Options:
  -R, --region\t\tSpecifies the desired region to search;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -W, --where\t\trestricted_where: Attribute query (like SQL WHERE)

  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -i, --index\t\tPrint the fetch results.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format in WGS84.
  -d, --dump\t\tDump the XYZ elevation data in WGS84 to stdout.
  -u, --update\t\tUpdate the Fetches Remote Elevation Datalist.

  --help\t\tPrint the usage text
  --version\t\tPrint the version information

Modules (see fetches --modules <module-name> for more info):
  {}

Examples:
 % {} -R -90.75/-88.1/28.7/31.25 nos --where "Date > 2000"
 % {} -R region.shp -p dc nos:datatype=bag charts:datatype=enc
 % {} -R region.shp dc:datatype=lidar -l > dc_lidar.urls
 % {} -R -89.75/-89.5/30.25/30.5 tnm:ds=4:formats=IMG gmrt:res=max:fmt=geotiff

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>
'''.format(os.path.basename(sys.argv[0]),
           _fetch_short_desc(_fetch_modules),
           _fetch_short_desc(_fetch_modules),
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]))

def fetches_cli(argv = sys.argv):
    status = 0
    extent = None
    want_list = False
    want_index = False
    want_proc = False
    want_dump = False
    want_update = False
    stop_threads = False
    these_regions = []
    mod_opts = {}
    verbose = True
    w = []
    
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
        elif arg == '--where' or arg == '-W':
            w.append(argv[i + 1])
            i = i + 1
        elif arg == '--list' or arg == '-l':
            want_list = True
        elif arg == '--index' or arg == '-i':
            want_index = True
        elif arg == '--process' or arg == '-p':
            want_proc = True
        elif arg == '--dump' or arg == '-d':
            want_dump = True
        elif arg == '--update' or arg == '-u':
            want_update = True
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format( _version))
            sys.exit(0)
        elif arg == '--modules' or arg == '-m':
            try:
                if argv[i + 1] in _fetch_modules.keys():
                    sys.stderr.write(_fetch_long_desc({k: _fetch_modules[k] for k in (argv[i + 1],)}))
                else: sys.stderr.write(_fetch_long_desc(_fetch_modules))
            except: sys.stderr.write(_fetch_long_desc(_fetch_modules))
            sys.exit(0)
        else: 
            opts = arg.split(':')
            if opts[0] in _fetch_modules.keys():
                mod_opts[opts[0]] = list(opts[1:])
            else: utils.echo_error_msg('invalid module name `{}`'.format(opts[0]))
        i = i + 1

    ## ==============================================
    ## use all modules if none specified
    ## ==============================================
    if len(mod_opts) == 0:
        for key in _fetch_modules.keys():
            mod_opts[key] = ''
        sys.stderr.write(_usage)
        utils.echo_error_msg('you must select at least one fetch module')
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
    else: these_regions = [[-180, 180, -90, 90]]

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
            fl = _fetch_modules[fetch_mod](regions.region_buffer(this_region, 5, pct = True), w, lambda: stop_threads, verbose)
            args_d = utils.args2dict(args)

            if want_update:
                fl._update()
                continue
            #try:
            ir = fl._filter_results()
            if want_index:
                for result in ir:
                    print(json.dumps(result, indent=4))
                continue
            r = fl._parse_results(ir, **args_d)
            utils.echo_msg('found {} data files.'.format(len(r)))
            if want_list:
                for result in r:
                    print(result[0])
                continue
            #except ValueError as e:
            #    utils.echo_error_msg('something went wrong, {}'.format(e))
            #    sys.exit(-1)
            #except Exception as e:
            #    utils.echo_error_msg('{}'.format(e))
            #    sys.exit(-1)
            
            fr = fetch_results(r, this_region, fl._outdir, fl if want_proc else None, fl if want_dump else None, lambda: stop_threads)
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
                stop_threads = True
                while not fr.fetch_q.empty():
                    try:
                        fr.fetch_q.get(False)
                    except Empty: continue
                    fr.fetch_q.task_done()
            fr.join()
            utils.echo_msg('ran fetch module {} on region {} ({}/{})...\
            '.format(fetch_mod, regions.region_format(this_region, 'str'), rn+1, len(these_regions)))            
### End
