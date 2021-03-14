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
### Code:

import os
import sys
import time
import copy

## ==============================================
## networking/etc
## ==============================================
import requests
import ftplib
import urllib
import lxml.html as lh
import lxml.etree
import json

## ==============================================
## queues and threading
## ==============================================
import threading
try:
    import Queue as queue
except: import queue as queue

## ==============================================
## import gdal/numpy
## ==============================================
import ogr
import gdal

## ==============================================
## import geomods
## ==============================================
from geomods import utils
from geomods import regions
from geomods import gdalfun
from geomods import gmtfun
from geomods import xyzfun
from geomods import mbsfun

_version = '0.6.3'

## =============================================================================
##
## Fetching Functions
##
## Generic fetching and processing functions, etc.
##
## The `fetch_results` class will fetch a list of fetch results [[url, file-name, data-type]...]
## in a queue `fetch_queue` using 3 threads; set `p` and `s` as a fetch module object to processes
## and dump XYZ data from the fetched results, respectively. Use `fetch_file` to fetch single files.
##
## =============================================================================
r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(_version) }

def fetch_queue(q, fo, fg):
    '''fetch queue `q` of fetch results'''
    while True:
        fetch_args = q.get()
        this_region = fetch_args[2]
        this_dt = fetch_args[4].lower()
        fetch_args[2] = None
        if not fetch_args[3]():
            if not fg['proc_p'] and not fg['dump_p']:
                if fetch_args[0].split(':')[0] == 'ftp':
                    fetch_ftp_file(*tuple(fetch_args))
                else: fetch_file(*tuple(fetch_args))
            else:
                if not fg['dump_p']:
                    if not os.path.exists(os.path.dirname(fetch_args[1])):
                        try:
                            os.makedirs(os.path.dirname(fetch_args[1]))
                        except: pass
                    #o_x_fn = '.'.join(fetch_args[1].split('.')[:-1]) + '.xyz'
                    o_x_fn = fetch_args[1] + '.xyz'
                    utils.echo_msg('processing local file: {}'.format(o_x_fn))
                    if not os.path.exists(o_x_fn):
                        with open(o_x_fn, 'w') as out_xyz:
                            fetch_dump_xyz(fetch_args, module = fo, epsg = 4326, z_region = fg['z_region'], inc = fg['inc'], dst_port = out_xyz)
                        if os.stat(o_x_fn).st_size == 0: utils.remove_glob(o_x_fn)
                else: fetch_dump_xyz(fetch_args, module = fo, epsg = 4326)
        q.task_done()

def fetch_ftp_file(src_url, dst_fn, params = None, callback = None, datatype = None, overwrite = False, verbose = False):
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

def fetch_file(src_url, dst_fn, params = None, callback = lambda: False, datatype = None, overwrite = False, verbose = False, timeout = 140, read_timeout = 320):
    '''fetch src_url and save to dst_fn'''
    status = 0
    req = None
    halt = callback

    if verbose: progress = utils._progress('fetching remote file: {}...'.format(os.path.basename(src_url)))
    if not os.path.exists(os.path.dirname(dst_fn)):
        try:
            os.makedirs(os.path.dirname(dst_fn))
        except: pass 
    if not os.path.exists(dst_fn) or overwrite:
        try:
            with requests.get(src_url, stream = True, params = params, headers = r_headers, timeout=(timeout,read_timeout)) as req:
                req_h = req.headers
                if req.status_code == 200:
                    curr_chunk = 0
                    with open(dst_fn, 'wb') as local_file:
                        for chunk in req.iter_content(chunk_size = 8196):
                            if halt(): break
                            if verbose: progress.update()
                            if chunk: local_file.write(chunk)
                else: utils.echo_error_msg('server returned: {}'.format(req.status_code))
        except Exception as e:
            utils.echo_error_msg(e)
            status = -1
    if not os.path.exists(dst_fn) or os.stat(dst_fn).st_size ==  0: status = -1
    if verbose: progress.end(status, 'fetched remote file: {}.'.format(os.path.basename(dst_fn)))
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
    def __init__(self, results, region, out_dir, fetch_obj = None, fg = None, callback = lambda: False):
        threading.Thread.__init__(self)
        self.fetch_q = queue.Queue()
        self.results = results
        self.region = region
        self._outdir = out_dir
        self.stop_threads = callback
        self.fetch_obj = fetch_obj
        self.fg = fg
        
    def run(self):
        for _ in range(3):
            t = threading.Thread(target = fetch_queue, args = (self.fetch_q, self.fetch_obj, self.fg))
            t.daemon = True
            t.start()
        for row in self.results:
            if self.fg['list_p']: print(row[0])
            if self.fg['fetch_p']: self.fetch_q.put([row[0], os.path.join(self._outdir, row[1]), self.region, self.stop_threads, row[2], False, True])
        self.fetch_q.join()
        
## =============================================================================
##
## XML Metadata parsing
##
## =============================================================================
def xml2py(node):
    '''parse an xml file into a python dictionary'''
    texts = {}
    if node is None: return(None)
    for child in list(node):
        child_key = lxml.etree.QName(child).localname
        if 'name' in child.attrib.keys(): child_key = child.attrib['name']
        if '{http://www.w3.org/1999/xlink}href' in child.attrib.keys():
            href = child.attrib['{http://www.w3.org/1999/xlink}href']
        else: href = None
        if child.text is None or child.text.strip() == '':
            if href is not None:
                if child_key in texts.keys():
                    texts[child_key].append(href)
                else: texts[child_key] = [href]
            else:
                if child_key in texts.keys():
                    ck = xml2py(child)
                    texts[child_key][list(ck.keys())[0]].update(ck[list(ck.keys())[0]])
                else: texts[child_key] = xml2py(child)
        else:
            if child_key in texts.keys():
                texts[child_key].append(child.text)
            else: texts[child_key] = [child.text]
    return(texts)

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
            'wms': 'http://www.opengis.net/wms',
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
            if geom: return(regions.region2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)]))
            else: return(region)
        else: return(None)

    def polygon(self, geom = True):
        opoly = []
        polygon = self.xml_doc.find('.//{*}Polygon', namespaces = self.namespaces)
        if polygon is not None:
            nodes = polygon.findall('.//{*}pos', namespaces = self.namespaces)
            [opoly.append([float(x) for x in node.text.split()]) for node in nodes]
            if geom: return(regions.wkt2geom(utils.create_wkt_polygon(opoly)))
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

class WCS:
    def __init__(self, url):
        self.url = url
        self.namespaces = {
            'wms': 'http://www.opengis.net/wms', 'wcs': 'http://www.opengis.net/wcs/2.0',
            'ows': 'http://www.opengis.net/ows/2.0', 'gml': 'http://www.opengis.net/gml/3.2',
            'gmlcov': 'http://www.opengis.net/gmlcov/1.0'}
        self._get_capabilities()
        self._s_version = self._si()['ServiceTypeVersion'][0]

    def _get_capabilities(self):
        _data = {'request': 'GetCapabilities', 'service': 'WCS'}
        c = fetch_req(self.url, params = _data)
        cx = lxml.etree.fromstring(c.text.encode('utf-8'))
        self.service_provider = cx.find('.//ows:ServiceProvider', namespaces = self.namespaces)
        self.service_identification = cx.find('.//ows:ServiceIdentification', namespaces = self.namespaces)
        self.operations_metadata = cx.find('.//ows:OperationsMetadata', namespaces = self.namespaces)
        self.service_metadata = cx.find('.//wcs:ServiceMetadata', namespaces = self.namespaces)
        self.contents = cx.find('.//wcs:Contents', namespaces = self.namespaces)

    def _contents(self):
        c = []
        for coverage in self.contents.xpath('//wcs:CoverageSummary', namespaces = self.namespaces):
            c.append(xml2py(coverage))
        return(c)

    def _om(self):
        return(xml2py(self.operations_metadata))

    def _sp(self):
        return(xml2py(self.service_provider))
    
    def _si(self):
        return(xml2py(self.service_identification))
    
    def fix_coverage_id(self, coverage):
        return(':'.join(coverage.split('__')))

    def unfix_coverage_id(self, coverage):
        return('__'.join(coverage.split(':')))

    def _describe_coverage(self, coverage):
        c_d = {}
        valid = False
        c = self._contents()
        for cc in c:
            if coverage == cc['CoverageId']:
                valid = True
                c_d = cc
                break

        om = self._om()
        url = om['DescribeCoverage']['DCP']['HTTP']['Get'][0]
        _data = {'request': 'DescribeCoverage', 'service': 'WCS',
            'version': self._s_version, 'CoverageID': self.unfix_coverage_id(coverage)}
        d = fetch_req(url, params = _data)
        d_r = lxml.etree.fromstring(d.text.encode('utf-8'))
        cd = d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)
        return(xml2py(d_r.find('.//wcs:CoverageDescription', namespaces = self.namespaces)))

    def _get_coverage_region(self, cov_desc):
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        return([lc[1], uc[1], lc[0], uc[0]])
    
    def _get_coverage_url(self, coverage, region = None):
        dl_coverage = self.fix_coverage_id(coverage)
        cov_desc = self._describe_coverage(coverage)
        fmt = cov_desc["ServiceParameters"]["nativeFormat"][0]        
        hl = [float(x) for x in cov_desc["domainSet"]["RectifiedGrid"]["limits"]["GridEnvelope"]['high'][0].split()]
        uc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["upperCorner"][0].split()]
        lc = [float(x) for x in cov_desc["boundedBy"]["Envelope"]["lowerCorner"][0].split()]
        ds_region = [lc[1], uc[1], lc[0], uc[0]]
        resx = (uc[1] - lc[1]) / hl[0]
        resy = (uc[0] - lc[0]) / hl[1]
        data = {'request': 'GetCoverage', 'version': '1.0.0', 'service': 'WCS',
                'resx': resx, 'resy': resy, 'crs': 'EPSG:4326', 'format': fmt,
                'coverage': coverage, 'Identifier': coverage}
        if region is not None: data['bbox'] = regions.region_format(region, 'bbox')        
        try:
            enc_data = urllib.urlencode(data)
        except: enc_data = urllib.parse.urlencode(data)
        return('{}{}'.format(self.url, enc_data))
    
    def fetch_coverage(coverage, region = None):
        c_url = self._get_coverage_url(coverage, region)
        return(fetch_file(c_url, '{}_{}.tif'.format(coverage, regions.region_format(region, 'fn')), params = data, verbose = True))
    
## =============================================================================
##
## Fetches Remote Elevation Datalist (FRED)
##
## the reference vector location and related functions
##
## =============================================================================
this_dir, this_filename = os.path.split(__file__)
fetchdata = os.path.join(this_dir, 'data')

class FRED:
    def __init__(self, verbose = False, local = False):
        self._verbose = verbose
        self.fetchdata = os.path.join(this_dir, 'data')
        self.driver = ogr.GetDriverByName('GeoJSON')
        self.fetch_v = 'FRED.geojson'
        if local:
            self.FREDloc = self.fetch_v
        elif os.path.exists(self.fetch_v):
            self.FREDloc = self.fetch_v
        elif os.path.exists(os.path.join(self.fetchdata, self.fetch_v)):
            self.FREDloc = os.path.join(self.fetchdata, self.fetch_v)
        else: self.FREDloc = self.fetch_v
        if self._verbose: utils.echo_msg('using {}'.format(self.FREDloc))
        self.ds = None
        self.layer = None
        self.open_p = False
        
        self._fields = ['Name', 'ID', 'Date', 'Agency', 'MetadataLink',
                        'MetadataDate', 'DataLink', 'IndexLink', 'Link',
                        'DataType', 'DataSource', 'Resolution', 'HorizontalDatum',
                        'VerticalDatum', 'LastUpdate', 'Etcetra', 'Info']
        
    def _create_ds(self):
        utils.remove_glob(self.FREDloc)
        self.ds = self.driver.CreateDataSource(self.FREDloc)        
        self.layer = self.ds.CreateLayer('FRED', None, ogr.wkbMultiPolygon)
        ldfn = self.layer.GetLayerDefn()
        self.layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('ID', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Date', ogr.OFTInteger))
        self.layer.CreateField(ogr.FieldDefn('Agency', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('MetadataDate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('IndexLink', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Link', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataType', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('DataSource', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Resolution', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('HorizontalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('VerticalDatum', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('LastUpdate', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Etcetra', ogr.OFTString))
        self.layer.CreateField(ogr.FieldDefn('Info', ogr.OFTString))

    def _add_survey(self, Name = 'fetches', ID = None, Date = None, Agency = None, MetadataLink = None,
                    MetadataDate = None, DataLink = None, IndexLink = None, Link = None, DataType = None,
                    DataSource = None, Resolution = None, HorizontalDatum = None, VerticalDatum = None,
                    Lastupdate = None, Etcetra = None, Info = None, geom = None):
        if not self.open_p: return(None)
        if geom is not None and self.layer is not None:
            self._add_feature([
                {'Name': Name, 'ID': ID, 'Agency': Agency,
                 'Date': Date, 'MetadataLink': MetadataLink,
                 'MetadataDate': MetadataDate, 'DataLink': DataLink,
                 'IndexLink': IndexLink, 'Link': Link, 'DataType': DataType,
                 'DataSource': DataSource, 'Resolution': Resolution,
                 'HorizontalDatum': HorizontalDatum, 'VerticalDatum': VerticalDatum,
                 'LastUpdate': utils.this_date(), 'Etcetra': Etcetra, 'Info': Info},
                geom.ExportToJson()
            ])
            [self.layer.SetFeature(ff) for ff in self.layer]
        else: return(None)
        
    def _open_ds(self, mode = 0):
        if not self.open_p:
            try:
                self.ds = self.driver.Open(self.FREDloc, mode)
            except: self.ds = None
            if self.ds is None:
                self._create_ds()
            self.layer = self.ds.GetLayer()
            self.open_p = True
            return(0)
        else: return(-1)
        
    def _close_ds(self):
        self.layer = self.ds = None
        self.open_p = False
        return(0)
        
    def _get_fields(self):
        if self.open_p:
            schema = []
            ldefn = self.layer.GetLayerDefn()
            for n in range(ldefn.GetFieldCount()):
                fdefn = ldefn.GetFieldDefn(n)
                schema.append(fdefn.name)
            return(schema)
        else: return(-1)
        
    def _add_feature(self, survey):
        '''add a survey to the reference vector layer'''
        if self.open_p:
            layer_defn = self.layer.GetLayerDefn()
            feat = ogr.Feature(layer_defn)
            geom = ogr.CreateGeometryFromJson(survey[1])
            feat.SetGeometry(geom)
            for field in self._fields:
                try:
                    feat.SetField(field, survey[0][field])
                except: feat.SetField(field, -1)
            self.layer.CreateFeature(feat)
            feat = None
            return(0)
        else: return(-1)

    def _edit_feature(self, feature, survey):
        if self.open_p:
            geom = ogr.CreateGeometryFromJson(survey[1])
            feature.SetGeometry(geom)
            for field in self._fields:
                try:
                    feature.SetField(field, survey[0][field])
                except: feature.SetField(field, -1)
            self.layer.SetFeature(feature)
            return(0)
        else: return(-1)

    def _add_surveys(self, surveys):
        '''update or create a reference vector using a list of surveys'''
        if self.open_p:
            if self.layer is not None:
                for survey in surveys:
                    self._add_survey(**survey)
            return(0)
        else: return(-1)

    def _get_region(self, where = [], layers = []):
        out_regions = []
        if self._verbose: _progress = utils._progress('gathering regions from {}...'.format(self.FREDloc))
        for i, layer in enumerate(layers):
            if self._verbose: _progress.update_perc((i, len(layers)))
            this_layer = self.layer
            this_layer.SetAttributeFilter("DataSource = '{}'".format(layer))
            [this_layer.SetAttributeFilter('{}'.format(filt)) for filt in where]
            for feat in this_layer:
                geom = feat.GetGeometryRef()
                wkt = geom.ExportToWkt()
                this_region = ogr.CreateGeometryFromWkt(wkt).GetEnvelope()
                if len(out_regions) > 0:
                    out_regions = regions.regions_merge(out_regions, this_region)
                else: out_regions = this_region
        return(out_regions)

    def _attribute_filter(self, where = []):
        if self.open_p:
            attr_str = []
            [attr_str.append(f) for f in where]
            wf = ' AND '.join(attr_str)
            self.layer.SetAttributeFilter(wf)
            return(0)
        else: return(-1)
        
    def _filter(self, region = None, where = [], layers = []):
        '''Search for data in the reference vector file'''
        _results = []
        if region is not None:
            _boundsGeom = regions.region2geom(region)
        else: _boundsGeom = None

        if self._verbose: _progress = utils._progress('filtering {}...'.format(self.FREDloc))
        if not self.open_p:
            self._open_ds()
            close_p = True
        else: close_p = False

        for i, layer in enumerate(layers):
            if self._verbose: _progress.update_perc((i, len(layers)))
            this_layer = self.layer
            where.append("DataSource = '{}'".format(layer))
            if self._verbose: utils.echo_msg('FRED filter: {}'.format(where))
            self._attribute_filter(where = where)
            for feat in this_layer:
                if _boundsGeom is not None:
                    geom = feat.GetGeometryRef()
                    if geom is not None:
                        if _boundsGeom.Intersects(geom):
                            _results.append({})
                            f_j = json.loads(feat.ExportToJson())
                            for key in f_j['properties'].keys():
                                _results[-1][key] = feat.GetField(key)
                else:
                    _results.append({})
                    f_j = json.loads(feat.ExportToJson())
                    for key in f_j['properties'].keys():
                        _results[-1][key] = feat.GetField(key)
                        
            this_layer = None
        if close_p: self._close_ds()
        if self._verbose: _progress.end(0, 'filtered \033[1m{}\033[m data records from FRED'.format(len(_results)))
        return(_results)

## ==============================================
## lambdas for the FRED using the module object `mod`
## ==============================================
_filter_FRED = lambda mod: mod.FRED._filter(mod.region, mod.where, [mod._name])
_update_FRED = lambda mod, s: mod.FRED._add_surveys(s)
_filter_FRED_index = lambda mod: [utils.echo_msg(json.dumps(f, indent = 2)) for f in _filter_FRED(mod)]

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
## Each module should have at class with a self._update function that updates the FRED reference vector,
## a self._parse_results function that parses the FRED feature record and yields a list of url lists
## e.g. `[url, file-name, datatype]`.
## and a self._yield_xyz function to yield the xyz data from the url list.
##
## self._update updates FRED
## self._parse_results yields [[url, file-name, datatype] ...]
## self._yield_xyz yields xyz data from the results
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

_fetch_long_desc = lambda x: 'fetches modules:\n% fetches ... <mod>:key=val:key=val...\n\n  {}\n'\
    .format('\n  '.join(['\033[1m{:14}\033[0m{}\n\n{}\n\no {}\n\n{}\n'\
                         .format(key, x[key](None, [], None, False)._title,
                                 x[key](None, [], None, False)._info,
                                 '\no '.join(x[key](None, [], None, False)._urls),
                                 x[key](None, [], None, False)._usage)
                         for key in x]))
_fetch_short_desc = lambda x: ' '.join(['{}'.format(key) for key in x])

def fetches_config(name = 'fetches', mods = ['gmrt'], region = None, where = [], fetch_p = True, inc = None,
                   index_p = False, list_p = False, proc_p = False, dump_p = False, update_p = False,
                   z_region = None, verbose = True):
    return({'name': name, 'mods': mods, 'region': region, 'where': where,
            'fetch_p': fetch_p, 'inc': inc, 'index_p': index_p, 'list_p': list_p,
            'proc_p': proc_p, 'dump_p': dump_p, 'update_p': update_p, 'z_region': z_region,
            'verbose': verbose})

## =============================================================================
##
## Digital Coast ('dc')
## fetches-dc
## Raster and Lidar data from NOAAs Digital Coast
##
## FRED holds the dataset level, fetch the index shapefile to parse the individual
## data files...
##
## - check bounds in xml for slr dems
## =============================================================================
class dc:
    def __init__(self, extent = None, where = [], callback = None, verbose = True):
        self._verbose = verbose
        if self._verbose: utils.echo_msg('loading DC fetch module...')
        
        ## ==============================================
        ## Digital Coast urls and directories
        ## ==============================================
        self._dc_url = 'https://coast.noaa.gov'
        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'lidar3_z', 'lidar4_z', 'raster1', 'raster2', 'raster5']
        self._outdir = os.path.join(os.getcwd(), 'dc')
        
        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'dc'
        self._info = '''Lidar and Raster data from NOAA's Digital Coast'''
        self._title = '''NOAA Digital Coast'''
        self._usage = '''< dc >'''
        self._urls = [self._dc_url, self._dc_htdata_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        '''Update the FRED reference vector after scanning
        the relevant metadata from Digital Coast.'''
        self.FRED = FRED(verbose = self._verbose, local = True)
        self.FRED._open_ds(1)
        for ld in self._dc_dirs:
            cols = []
            surveys = []
            page = fetch_html(self._dc_htdata_url + ld)
            if page is None: continue
            tr = page.xpath('//table')[0].xpath('.//tr')
            if len(tr) <= 0: continue
            [cols.append(i.text_content()) for i in tr[0]]
            
            if self._verbose: _prog = utils._progress('scanning {} datasets in {}...'.format(len(tr), ld))
            for i in range(1, len(tr)):
                if self._stop(): break
                if self._verbose: _prog.update_perc((i, len(tr))) #dc['ID #']))
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
                self.FRED._attribute_filter(["ID = '{}'".format(dc['ID #'])])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    if 'Metadata' in dc.keys():
                        this_xml = iso_xml(dc['Metadata'])
                        h_epsg, v_epsg = this_xml.reference_system()
                        geom = this_xml.bounds(geom=True)
                        if geom is not None:
                            if self._verbose: _prog.update_perc((i, len(tr)), msg = '{} ** adding: {} **'.format(_prog.opm, dc['ID #']))
                            surveys.append({'Name': dc['Dataset Name'], 'ID': dc['ID #'], 'Date': this_xml.date(),
                                            'MetadataLink': dc['Metadata'], 'MetadataDate': this_xml.xml_date(),
                                            'DataLink': dc['https'], 'IndexLink': dc['Tile Index'], 'Link': self._dc_url,
                                            'DataType': ld.split('_')[0], 'DataSource': 'dc', 'HorizontalDatum': h_epsg,
                                            'VerticalDatum': v_epsg, 'Info': this_xml.abstract(), 'geom': geom})
            self.FRED._add_surveys(surveys)
            if self._verbose:
                _prog.end(0, 'scanned {} datasets in {}.'.format(len(tr), ld))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), ld))
        self.FRED._close_ds()

    def _parse_results(self):
        for surv in _filter_FRED(self):
            if self._stop(): break
            surv_shp_zip = os.path.basename(surv['IndexLink'])
            if fetch_file(surv['IndexLink'], surv_shp_zip, callback = self._stop, verbose = self._verbose) == 0:
                v_shps = utils.p_unzip(surv_shp_zip, ['.shp', '.shx', '.dbf', '.prj'])
                v_shp = None
                for v in v_shps:
                    if v.split('.')[-1] == '.shp': v_shp = v
                try:
                    v_ds = ogr.Open(v_shp)
                    slay1 = v_ds.GetLayer(0)
                    for sf1 in slay1:
                        if self._stop(): break
                        geom = sf1.GetGeometryRef()
                        if geom.Intersects(regions.region2geom(self.region)):
                            yield([sf1.GetField('URL').strip(), '{}/{}'.format(surv['ID'], tile_url.split('/')[-1]), surv['DataType']])
                    v_ds = slay1 = None
                except: pass
                utils.remove_glob(surv_shp_zip, *v_shps)

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326, z_region = None, inc = None):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1].lower()
        if src_ext == 'laz' or src_ext == 'las': dt = 'lidar'
        elif src_ext == 'tif' or src_ext == 'img': dt = 'raster'
        else: dt = None
        if dt == 'lidar':
            if fetch_file(entry[0], src_dc, callback = self._stop, verbose = self._verbose) == 0:
                xyz_dat = utils.yield_cmd('las2txt -stdout -parse xyz -keep_xy {} -keep_class {} -i {}\
                '.format(regions.region_format(self.region, 'te'), '2 29', src_dc), verbose = False)
                xyzc = copy.deepcopy(xyzfun._xyz_config)
                xyzc['name'] = src_dc
                xyzc['epsg'] = 4326
                xyzc['warp'] = epsg
                if z_region is not None:
                    xyzc['upper_limit'] = z_region[1]
                    xyzc['lower_limit'] = z_region[0]
                for xyz in xyzfun.xyz_parse(xyz_dat, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        elif dt == 'raster':
            try:
                src_ds = gdal.Open(entry[0])
            except Exception as e:
                fetch_file(entry[0], src_dc, callback = self._stop, verbose = self._verbose)
                try:
                    src_ds = gdal.Open(src_dc)
                except Exception as e:
                    utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                    src_ds = None
            except Exception as e:
                utils.echo_error_msg('could not read dc raster file: {}, {}'.format(entry[0], e))
                src_ds = None
            
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, regions.region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_get_epsg(src_ds)))
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = True, z_region = z_region):
                    yield(xyz)
            src_ds = None
        utils.remove_glob(src_dc)

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
        self._nos_url = 'https://www.ngdc.noaa.gov/mgg/bathymetry/hydro.html'
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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'nos'
        self._info = '''Bathymetry surveys and data (xyz & BAG)'''
        self._title = '''NOAA NOS Bathymetric Data'''
        self._usage = '''< nos >'''
        self._urls = [self._nos_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        '''Crawl the NOS database and update/generate the NOS reference vector.'''
        self.FRED._open_ds(1)
        for nosdir in self._nos_directories:
            if self._stop(): break
            surveys = []
            xml_catalog = self._nos_xml_url(nosdir)
            page = fetch_html(xml_catalog)
            rows = page.xpath('//a[contains(@href, ".xml")]/@href')
            if self._verbose: _prog = utils._progress('scanning {} surveys in {}...'.format(len(rows), nosdir))

            for i, survey in enumerate(rows):
                if self._stop(): break
                sid = survey[:-4]
                if self._verbose: _prog.update_perc((i, len(rows)))
                self.FRED._attribute_filter(["ID = '{}'".format(sid)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    this_xml = iso_xml(xml_catalog + survey)
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
                        surveys.append({'Name': this_xml.title(), 'ID': sid, 'Agency': 'NOAA/NOS', 'Date': this_xml.date(),
                                        'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(), 'DataLink': ','.join(list(set(d_links))),
                                        'DataType': ','.join(list(set(d_types))), 'DataSource': 'nos', 'HorizontalDatum': h_epsg,
                                        'VerticalDatum': v_epsg, 'Info': this_xml.abstract(), 'geom': geom})
            if self._verbose:
                _prog.end(0, 'scanned {} surveys in {}.'.format(len(rows), nosdir))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), nosdir))
            self.FRED._add_surveys(surveys)
        self.FRED._close_ds()
    
    def _parse_results(self):
        '''Search the NOS reference vector and append the results
        to the results list.'''
        for surv in _filter_FRED(self):
            yield([self._data_urls.append([i, i.split('/')[-1], surv['DataType']]) if i != '' else None for i in surv['DataLink'].split(',')])
            
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326, z_region = None, inc = None):
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_nos = os.path.basename(entry[1])
        dt = None
        if fetch_file(entry[0], src_nos, callback = self._stop, verbose = self._verbose) == 0:
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
                nos_fns = utils.p_unzip(src_nos, ['xyz', 'dat'])
                xyzc['delim'] = ','
                xyzc['skip'] = 1
                xyzc['xpos'] = 2
                xyzc['ypos'] = 1
                xyzc['zpos'] = 3
                xyzc['z-scale'] = -1
                xyzc['epsg'] = 4326
                xyzc['warp'] = epsg
                if z_region is not None:
                    xyzc['upper_limit'] = z_region[1]
                    xyzc['lower_limit'] = z_region[0]
                for nos_f_r in nos_fns:
                    xyzc['name'] = nos_f_r
                    if os.path.exists(nos_f_r):
                        with open(nos_f_r, 'r') as in_n:
                            for xyz in xyzfun.xyz_parse(in_n, xyz_c = xyzc, region = self.region, verbose = self._verbose):
                                yield(xyz)
                utils.remove_glob(*nos_fns)
                #[utils.remove_glob(x) for x in nos_fns]

            elif dt == 'grid_bag':
                src_bags = utils.p_unzip(src_nos, ['gdal', 'tif', 'img', 'bag', 'asc'])

                nos_f = '{}.tmp'.format(os.path.basename(src_bag).split('.')[0])
                for src_bag in src_bags:
                    try:
                        src_ds = gdal.Open(src_bag)
                        srcwin = gdalfun.gdal_srcwin(src_ds, regions.region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_get_epsg(src_ds)))
                        with open(nos_f, 'w') as cx:
                            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, z_region = z_region):
                                yield(xyz)
                        src_ds = None
                    except: pass
                utils.remove_glob(*src_bags)
        utils.remove_glob(src_nos)
        
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
        self._charts_url = 'http://www.charts.noaa.gov/'
        self._enc_data_catalog = 'http://www.charts.noaa.gov/ENCs/ENCProdCat_19115.xml'
        self._rnc_data_catalog = 'http://www.charts.noaa.gov/RNCs/RNCProdCat_19115.xml'
        self._outdir = os.path.join(os.getcwd(), 'charts')
        self._dt_xml = { 'ENC':self._enc_data_catalog,
                         'RNC':self._rnc_data_catalog }

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'charts'
        self._info = '''Raster and Vector U.S. Nautical Charts'''
        self._title = '''NOAA Nautical CHARTS (RNC & ENC)'''
        self._usage = '''< charts >'''
        self._urls = [self._enc_data_catalog, self._rnc_data_catalog]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()

    def _update(self):
        '''Update or create the reference vector file'''
        self.FRED._open_ds(1)
        for dt in self._dt_xml.keys():
            surveys = []
            this_xml = iso_xml(self._dt_xml[dt], timeout = 1000, read_timeout = 2000)
            charts = this_xml.xml_doc.findall('.//{*}has', namespaces = this_xml.namespaces)
            _prog = utils._progress('scanning {} surveys in {}.'.format(len(charts), dt))
            for i, chart in enumerate(charts):
                this_xml.xml_doc = chart
                title = this_xml.title()
                if self._verbose: _prog.update_perc((i, len(charts)))
                self.FRED._attribute_filter(["ID = '{}'".format(title)])
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    h_epsg, v_epsg = this_xml.reference_system()
                    this_data = this_xml.linkages()
                    geom = this_xml.polygon(geom=True)
                    if geom is not None:
                        surveys.append({'Name': title, 'ID': title, 'Agency': 'NOAA', 'Date': this_xml.date(),
                                        'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(),
                                        'DataLink': this_data, 'Link': self._charts_url, 'DataType': dt,
                                        'DataSource': 'charts', 'HorizontalDatum': h_epsg, 'VerticalDatum': v_epsg,
                                        'Info': this_xml.abstract, 'geom': geom})
            self.FRED._add_surveys(surveys)
            if self._verbose:
                _prog.end(0, 'scanned {} surveys in {}'.format(len(charts), dt))
                utils.echo_msg('added {} surveys from {}'.format(len(surveys), dt))
        self.FRED._close_ds()
        
    def _parse_results(self):
        '''Search for data in the reference vector file'''
        for surv in _filter_FRED(self):
            for i in surv['DataLink'].split(','):
                yield([i, i.split('/')[-1], surv['DataType']])

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        xyzc = xyzfun._xyz_config
        xyzc['z-scale'] = -1
        xyzc['warp'] = epsg
        if z_region is not None:
            xyzc['upper_limit'] = z_region[1]
            xyzc['lower_limit'] = z_region[0]
        src_zip = os.path.basename(entry[1])
        
        if fetch_file(entry[0], src_zip, callback = self._stop, verbose = self._verbose) == 0:
            if entry[-1].lower() == 'enc':
                src_encs = utils.p_unzip(src_zip, ['000'])
                for src_ch in src_encs:
                    dst_xyz = src_ch.split('.')[0] + '.xyz'
                    try:
                        ds_ogr = ogr.Open(src_ch)
                        layer_s = ds_ogr.GetLayerByName('SOUNDG')
                        if layer_s is not None:
                            with open(dst_xyz, 'w') as o_xyz:
                                for f in layer_s:
                                    g = json.loads(f.GetGeometryRef().ExportToJson())
                                    for xyz in g['coordinates']:
                                        xyzfun.xyz_line([float(x) for x in xyz], o_xyz, False)
                        ds_ogr = layer_s = None
                    except: pass
                    if os.path.exists(dst_xyz):
                        with open(dst_xyz, 'r') as in_c:
                            for xyz in xyzfun.xyz_parse(in_c, xyz_c = xyzc, verbose = self._verbose):
                                yield(xyz)
                utils.remove_glob(dst_xyz, *src_encs)
        utils.remove_glob(src_zip)
        
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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'ncei_thredds'
        self._info = '''NOAA NCEI THREDDS DEM Catalog Access.
Digital Elevation Models around the world at various resolutions and extents.
NCEI builds and distributes high-resolution, coastal digital elevation models (DEMs) that integrate ocean 
bathymetry and land topography supporting NOAA's mission to understand and predict changes in Earth's environment, 
and conserve and manage coastal and marine resources to meet our Nation's economic, social, and environmental needs.
\nDEMs are used for coastal process modeling (tsunami inundation, storm surge, sea-level rise, contaminant dispersal, 
etc.), ecosystems management and habitat research, coastal and marine spatial planning, and hazard mitigation and 
community preparedness.'''
        self._title = '''NCEI THREDDS'''
        self._usage = '''< ncei_thredds >'''
        self._urls = [self._nt_catalog, self._ngdc_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
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
        if self._verbose: _prog = utils._progress('scanning {} datasets in {}...'.format(len(this_ds), this_ds[0].attrib['name']))
        for i, node in enumerate(this_ds):
            surveys = []
            this_title = node.attrib['name']
            this_id = node.attrib['ID']
            if self._verbose: _prog.update_perc((i, len(this_ds)))
            self.FRED._attribute_filter(["ID = '{}'".format(this_id)])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
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

                zv = this_xml.xml_doc.findall('.//gmd:dimension/gmd:MD_Band/gmd:sequenceIdentifier/gco:MemberName/gco:aName/gco:CharacterString', namespaces = this_xml.namespaces)
                if zv is not None:
                    for zvs in zv:
                        if zvs.text == 'bathy' or zvs.text == 'Band1' or zvs.text == 'z':
                            zvar = zvs.text
                        else: zvar = 'z'
                geom = this_xml.bounds(geom=True)
                if geom is not None:
                    surveys.append({'Name': title, 'ID': this_id, 'Agency': 'NOAA', 'Date': this_xml.date(),
                                    'MetadataLink': this_xml.url, 'MetadataDate': this_xml.xml_date(),
                                    'DataLink': http_url, 'IndexLink': wcs_url, 'Link': self._nt_catalog,
                                    'DataType': 'raster', 'DataSource': 'ncei_thredds', 'HorizontalDatum': h_epsg,
                                    'VerticalDatum': v_epsg, 'Etcetra': zvar, 'Info': this_xml.abstract(), 'geom': geom})
        self.FRED._add_surveys(surveys)
        if self._verbose:
            _prog.end(0, 'scanned {} datasets in {}.'.format(len(this_ds), this_ds[0].attrib['name']))
            utils.echo_msg('added {} surveys from {}'.format(len(surveys), this_ds[0].attrib['name']))
        
    def _update(self):
        self.FRED._open_ds(1)
        self._parse_catalog(self._nt_catalog)
        self.FRED._close_ds()
    
    def _parse_results(self):
        '''Search for data in the reference vector file'''
        for surv in _filter_FRED(self):
            #wcs_url = "{}?request=GetCoverage&version=1.0.0&service=WCS&coverage={}&bbox={}&format=NetCDF3"\
            #    .format(surv['IndexLink'], surv['DataType'], regions.region_format(self.region, 'bbox'))
            for d in surv['DataLink'].split(','):
                if d != '':
                    yield[d, d.split('/')[-1], surv['DataType']]

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326, z_region = None, inc = None):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        try:
            src_ds = gdal.Open(entry[0])
        except Exception as e:
            fetch_file(entry[0], src_dc, callback = self._stop, verbose = self._verbose)
            try:
                src_ds = gdal.Open(src_dc)
            except Exception as e:
                utils.echo_error_msg('could not read raster file: {}, {}'.format(entry[0], e))
                src_ds = None
        except Exception as e:
            utils.echo_error_msg('could not read raster file: {}, {}'.format(entry[0], e))
            src_ds = None

        if src_ds is not None:
            srcwin = gdalfun.gdal_srcwin(src_ds, regions.region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose, z_region = z_region):
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_dc)
            
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
        self._mb_html = "https://www.ngdc.noaa.gov/"
        self._outdir = os.path.join(os.getcwd(), 'mb')

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'mb'
        self._info = '''NCEI is the U.S. national archive for multibeam bathymetric data and holds 
more than 9 million nautical miles of ship trackline data recorded from over 2400 cruises and received 
from sources worldwide. In addition to deepwater data, the Multibeam Bathymetry Database (MBBDB) includes 
hydrographic multibeam survey data from NOAA's National Ocean Service (NOS).'''
        self._title = '''NOAA MULTIBEAM survey data'''
        self._usage = '''< mb >'''
        self._urls = [self._mb_data_url, self._mb_metadata_url, self._mb_autogrid]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('MB-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'NOAA Multibeam', ID = 'MB-1', Agency = 'NOAA', Date = utils.this_year(),
                                  MetadataLink = self._mb_metadata_url, MetadataDate = utils.this_year(),
                                  DataLink = self._mb_search_url, DataType = 'multibeam', DataSource = 'mb',
                                  HorizontalDatum = 4326, VerticalDatum = 1092, Info = self._info,
                                  geom = regions.region2geom([-180,180,-90,90]))
        self.FRED._close_ds()
        
    def _parse_results(self, processed = False):
        '''Run the MB (multibeam) fetching module.'''
        these_surveys = {}
        these_versions = {}
        for surv in _filter_FRED(self):
            if self.region is None: break
            _req = fetch_req(surv['DataLink'], params = {'geometry': regions.region_format(self.region, 'bbox')}, timeout = 20)
            if _req is not None and _req.status_code == 200:
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
            else: utils.echo_error_msg('{}'.format(_req.reason))
                    
        for key in these_surveys.keys():
            if '2' in these_surveys[key].keys():
                for v2 in these_surveys[key]['2']:
                    yield(v2)
            else:
                for v1 in these_surveys[key]['1']:
                    yield(v1)
                    
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_mb = os.path.basename(entry[1])

        # split_entry = entry[0].split('/')
        # ship = split_entry[6]
        # survey = split_entry[7]
        # survey_html = '{}ships/{}/{}_mb.html'.format(self._mb_html, ship, survey)
        # surv_ = fetch_html(survey_html)

        # #print(surv_.xpath("//t[@id == 'summary']"))
        
        # for i in surv_.get_element_by_id('summary'):
        #     if i.text == 'Multibeam Bathymetry':
        #         for j in i.find_class('content'):
        #             print(j)
        #         #content = i.find_class('content')
        #         #print(content)
        #         #if len(content) > 0:
        #         #print(content)
        #         #tr = i.xpath('//table')[0].xpath('.//tr')
        #         #if len(tr) > 0:
        #         #    print(tr[0].text)
        
        #     #tr = page.xpath('//table')[0].xpath('.//tr')
        # #if len(tr) <= 0: continue

        # #[cols.append(i.text_content()) for i in tr[0]]

        if fetch_file(entry[0], src_mb, callback = self._stop, verbose = self._verbose) == 0:
            #src_xyz = os.path.basename(src_mb).split('.')[0] + '.xyz'
            src_xyz = os.path.basename(src_mb) + '.xyz'
            
            out, status = utils.run_cmd('mblist -MA -OXYZ -I{}  > {}'.format(src_mb, src_xyz), verbose = False)
            if status != 0:
                if fetch_file('{}.inf'.format(entry[0]), '{}.inf'.format(src_mb), callback = self._stop, verbose = self._verbose) == 0:
                    mb_fmt = mbsfun.mb_inf_data_format('{}.inf'.format(src_mb))
                    utils.remove_glob('{}.inf'.format(entry[0]))
                out, status = utils.run_cmd('mblist -F{} -MA -OXYZ -I{}  > {}'.format(mb_fmt, src_mb, src_xyz), verbose = False)
            if status == 0:
                xyzc['name'] = src_mb
                xyzc['z-scale'] = 1
                xyzc['epsg'] = 4326
                xyzc['warp'] = epsg
                if z_region is not None:
                    xyzc['upper_limit'] = z_region[1]
                    xyzc['lower_limit'] = z_region[0]

                xyzc['delim'] = '\t'
                with open(src_xyz, 'r') as in_m:
                    xyz_func = lambda p: xyzfun.xyz_dump(in_m, xyz_c = xyzc, region = self.region, verbose = False, dst_port = p)
                    if inc is not None:
                        xyzc['delim'] = None
                        out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
                        for xyz in utils.yield_cmd('gmt blockmedian -I{:.10f} {} -r -V'.format(inc, regions.region_format(self.region, 'gmt')), verbose = self._verbose, data_fun = xyz_func):
                            yield([float(x) for x in xyz.split()])

                    else:
                        for xyz in xyzfun.xyz_parse(in_m, xyz_c = xyzc, region = self.region, verbose = self._verbose):
                            yield(xyz)
                utils.remove_glob(src_mb, src_xyz)
            else:
                utils.echo_error_msg('failed to process local file, {} [{}]...'.format(src_mb, entry[0]))
                with open('{}'.format(os.path.join(self._outdir, 'fetch_{}_{}.err'.format(self._name, regions.region_format(self.region, 'fn')))), 'a') as mb_err:
                    mb_err.write('{}\n'.format(','.join([src_mb, entry[0]])))
                os.rename(src_mb, os.path.join(self._outdir, src_mb))
                utils.remove_glob(src_xyz)
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_mb))
        
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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'usace'
        self._info = '''Bathymetric Channel surveys from USACE - U.S. only.
The hydrographic surveys provided by this application are to be used for informational purposes only 
and should not be used as a navigational aid. Channel conditions can change rapidly and the surveys
may or may not be accurate.'''
        self._title = '''USACE bathymetry surveys via eHydro'''
        self._usage = '''< usace >'''
        self._urls = [self._usace_gj_url, self._usace_gs_api_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('USACE-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'USACE E-Hydro', ID = 'USACE-1', Agency = 'USACE', Date = utils.this_year(),
                                  DataLink = self._usace_gs_api_url, IndexLink = self._usace_gj_url, DataType = 'xyz',
                                  DataSource = 'usace', Info = self._info, geom = regions.region2geom([-162,-60,16,73]))
        self.FRED._close_ds()
        
    def _parse_results(self, stype = None):
        '''Run the USACE fetching module'''
        for surv in _filter_FRED(self):
            _data = {'geometry': regions.region_format(self.region, 'bbox'), 'inSR':4326, 'outSR':4326, 'f':'pjson'}
            _req = fetch_req(surv['DataLink'], params = _data)
            if _req is not None and _req.status_code == 200:
                survey_list = self._req.json()
                for feature in survey_list['features']:
                    fetch_fn = feature['attributes']['SOURCEDATALOCATION']
                    if stype is not None:
                        if feature['attributes']['SURVEYTYPE'].lower() == stype.lower():
                            yield([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
                    else: yield([fetch_fn, fetch_fn.split('/')[-1], 'usace'])
            else: utils.echo_error_msg('{}'.format(_req.reason))
            
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        xyzc = copy.deepcopy(xyzfun._xyz_config)
        src_zip = os.path.basename(entry[1])
        xyzc['warp'] = epsg
        
        if fetch_file(entry[0], src_zip, callback = self._stop, verbose = self._verbose) == 0:
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

            if xyzc['epsg'] is None:
                this_geom = regions.region2geom(this_xml_region)
                sp_fn = os.path.join(fetchdata, 'stateplane.geojson')
                try:
                    sp = ogr.Open(sp_fn)
                    layer = sp.GetLayer()
                
                    for feature in layer:
                        geom = feature.GetGeometryRef()
                        if this_geom.Intersects(geom):
                            xyzc['epsg'] = feature.GetField('EPSG')
                    sp = None
                except: pass
                        
            src_usaces = utils.p_unzip(src_zip, ['.XYZ', '.xyz', '.dat'])
            for src_usace in src_usaces:
                if os.path.exists(src_usace):
                    with open(src_usace, 'r') as in_c:
                        xyzc['name'] = src_usace.split('.')[0]
                        if z_region is not None:
                            xyzc['upper_limit'] = z_region[1]
                            xyzc['lower_limit'] = z_region[0]
                        for xyz in xyzfun.xyz_parse(in_c, xyz_c = xyzc, verbose = self._verbose):
                            yield(xyz)
                    utils.remove_glob(src_usace)
                    
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(entry[0]))
        utils.remove_glob(src_zip)

## =============================================================================
##
## The National Map (TNM) - USGS
##
## Fetch elevation data from The National Map
## NED, 3DEP, NHD, Etc.
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
        self._elev_ds = ['National Elevation Dataset (NED) 1 arc-second', 'Digital Elevation Model (DEM) 1 meter',
                         'National Elevation Dataset (NED) 1/3 arc-second', 'National Elevation Dataset (NED) 1/9 arc-second',
                         'National Elevation Dataset (NED) Alaska 2 arc-second', 'Alaska IFSAR 5 meter DEM',
                         'Original Product Resolution (OPR) Digital Elevation Model (DEM)', 'Ifsar Digital Surface Model (DSM)',
                         'Ifsar Orthorectified Radar Image (ORI)', 'Lidar Point Cloud (LPC)',
                         'National Hydrography Dataset Plus High Resolution (NHDPlus HR)', 'National Hydrography Dataset (NHD) Best Resolution',
                         'National Watershed Boundary Dataset (WBD)', 'USDA National Agriculture Imagery Program (NAIP)',
                         'Topobathymetric Lidar DEM', 'Topobathymetric Lidar Point Cloud']
        self._outdir = os.path.join(os.getcwd(), 'tnm')

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'tnm'
        self._info = '''Various datasets from USGS's National Map. The National Map is a 
collaborative effort among the USGS and other Federal, State, and local partners to improve
and deliver topographic information for the Nation.'''
        ## \nDatasets: {}'''.format(self._fred_datasets())
        self._title = '''The National Map (TNM) from USGS'''
        self._usage = '''< tnm:formats=fmt,fmt:extents=ext,ext >'''
        self._urls = [self._tnm_api_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _fred_datasets(self):
        ids = {}
        for r in self.FRED._filter(self.region, self.where, ['tnm']):
            ids[r['ID']] = r['DataType']
        return(ids)

    def _datasets(self, dataset = None):
        _req = fetch_req(self._tnm_dataset_url)
        if _req is not None and _req.status_code == 200:
            try:
                _datasets = _req.json()
                if dataset is not None:
                    for ds in _datasets:
                        tags = ds['tags']
                        if len(tags) > 0:
                            for t in tags:
                                if dataset == t['sbDatasetTag']:
                                    _datasets = t
                                    break
                        else:
                            if dataset == ds['sbDatasetTag']:
                                _datasets = ds
                                break
            except Exception as e:
                utils.echo_error_msg('try again, {}'.format(e))
        else: _datasets = None
        return(_datasets)
    
    def _dataset_tags():
        tnm_ds = self._datasets()
        dsTags = []
        for ds in tnm_ds:
            tags = ds['tags']
            if len(tags) > 0:
                for t in tags:
                    dsTags.append(t['sbDatasetTag'])
            else: dsTags.append(ds['sbDatasetTag'])
        return(dsTags)

    def _update_dataset(self, ds, fmt, geom, h_epsg, v_epsg):
        self.FRED._attribute_filter(["ID = '{}'".format(ds['id'])])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            if 'IMG' in fmt or 'TIFF' in fmt:
                datatype = 'raster'
            elif 'LAS' in fmt or 'LAZ' in fmt:
                datatype = 'lidar'
            else: datatype = 'tnm'
            try:
                url_enc = urllib.urlencode({'datasets': ds['sbDatasetTag']})
            except: url_enc = urllib.parse.urlencode({'datasets': ds['sbDatasetTag']})
            try:
                pubDate = ds['lastPublishedDate']
            except: pubDate = utils.this_year()
            try:
                metadataDate = ds['lastUpdatedDate']
            except: metadataDate = utils.this_year()            
            if geom is not None:
                self.FRED._add_survey(Name = ds['sbDatasetTag'], ID = ds['id'], Agency = 'USGS', Date = pubDate[-4:],
                                      MetadataLink = ds['infoUrl'], MetadataDate = metadataDate,
                                      DataLink = '{}{}'.format(self._tnm_product_url, url_enc), Link = ds['dataGovUrl'], Resolution = ','.join(ds['extents']),
                                      DataType = datatype, DataSource = 'tnm', HorizontalDatum = h_epsg,
                                      VerticalDatum = v_epsg, Etcetra = fmt, Info = ds['refreshCycle'], geom = geom)

    ## ==============================================
    ## update the FRED geojson with each TNM dataset
    ## each dataset will have any number of products, which get parsed for the data-link
    ## in _parse_results().
    ## ==============================================
    def _update(self):
        '''Update FRED with each dataset in TNM'''
        datasets = self._datasets()
        self.FRED._open_ds(1)
        if self._verbose: _prog = utils._progress('scanning {} datasets from TNM...'.format(len(datasets)))
        for i, ds in enumerate(datasets):
            if self._verbose: _prog.update_perc((i, len(datasets)))
            for fmt in ds['formats']:
                if 'isDefault' in fmt.keys():
                    fmt = fmt['value']
                    break
            this_xml = iso_xml('{}{}?format=iso'.format(self._tnm_meta_base, ds['id']))
            geom = this_xml.bounds(geom = True)
            h_epsg, v_epsg = this_xml.reference_system()
            tags = ds['tags']
            if len(tags) > 0:
                for tag in tags:
                    self._update_dataset(tag, fmt, geom, h_epsg, v_epsg)
            else: self._update_dataset(ds, fmt, geom, h_epsg, v_epsg)
        if self._verbose: _prog.end(0, 'scanned {} datasets from TNM'.format(len(datasets)))
        self.FRED._close_ds()

    def _parse_results(self, f = None, e = None, q = None):
        '''parse the tnm results from FRED,
        creates a generator yielding [url, file-name, data-type]'''
        e = e.split(',') if e is not None else None
        f = f.split(',') if f is not None else None
        for surv in _filter_FRED(self):
            print(surv)
            offset = 0
            total = 0
            while True:
                _dataset_results = []
                _data = {'bbox': regions.region_format(self.region, 'bbox'), 'max': 100, 'offset': offset}
                if q is not None: _data['q'] = str(q)
                if f is None:
                    _data['prodFormats'] = surv['Etcetra']
                else: _data['prodFormats'] = ','.join(f)
                if e is None: e = []

                _req = fetch_req(surv['DataLink'], params = _data)
                if _req is not None and _req.status_code == 200:
                    try:
                        _dataset_results = _req.json()
                        total = _dataset_results['total']
                    except ValueError: utils.echo_error_msg('tnm server error, try again')
                    except Exception as e: utils.echo_error_msg('error, {}'.format(e))
                
                if len(_dataset_results) > 0:
                    for item in _dataset_results['items']:
                        if _data['prodFormats'] is None:
                            fmts = []
                        else: fmts = _data['prodFormats'].split(',')
                        f_url = None
                        if len(e) > 0:
                            for extent in e:
                                if item['extent'] == extent:
                                    for fmt in fmts:
                                        if fmt in item['urls'].keys():
                                            f_url = item['urls'][fmt]
                                            break
                                    if f_url is None: f_url = item['downloadURL']
                                    yield([f_url, f_url.split('/')[-1], surv['DataType']])
                        else:
                            for fmt in fmts:
                                if fmt in item['urls'].keys():
                                    f_url = item['urls'][fmt]
                                    break
                            if f_url is None:  f_url = item['downloadURL']
                            yield([f_url, f_url.split('/')[-1], surv['DataType']])
                offset += 100
                if offset >= total: break

    ## ==============================================
    ## _update_prods() and _parse_prods_results() will update FRED with every product as a feature, rather than
    ## the default of each feature being a TNM dataset. _update_prods() takes much longer time to gather the
    ## products for each dataset and recording them in FRED, though the parsing of results is much faster.
    ## For our purposes, we wont be using most of what's available on TNM, so it is a bit of a waste to store
    ## all their datasets, which are already stored online, in FRED. This means user-time for fetches TNM is a
    ## bit slower, however storage costs are minimal and fewer updates may be necesary...
    ## ==============================================                
    def _update_prods(self):
        '''updated FRED with each product file available from TNM'''
        for dsTag in self._elev_ds:
            offset = 0
            utils.echo_msg('processing TNM dataset {}...'.format(dsTag))
            _req = fetch_req(self._tnm_product_url, params = {'max': 1, 'datasets': dsTag})
            try:
                _dsTag_results = _req.json()
            except ValueError: utils.echo_error_msg('tnm server error, try again')
            except Exception as e: utils.echo_error_msg('error, {}'.format(e))
            total = _dsTag_results['total']
            if self._verbose: _prog = utils._progress('gathering {} products from {}...'.format(total, dsTag))
            
            ds = self._datasets(dataset = dsTag)
            this_xml = iso_xml('{}{}?format=iso'.format(self._tnm_meta_base, ds['id']))
            h_epsg, v_epsg = this_xml.reference_system()
            
            while True:
                _data = {'max': 100, 'datasets': dsTag, 'offset': offset}
                _req = fetch_req(self._tnm_product_url, params = _data)
                try:
                    _dsTag_results = _req.json()
                except ValueError: utils.echo_error_msg('tnm server error, try again')
                except Exception as e: utils.echo_error_msg('error, {}'.format(e))
                if self._verbose: _prog.update_perc((offset,total), msg = 'gathering {} products from {}...'.format(total, dsTag))
                
                for i, item in enumerate(_dsTag_results['items']):
                    if self._verbose: _prog.update_perc((i+offset,total), msg = 'gathering {} products from {}...'.format(total, dsTag))
                    try:
                        self.FRED.layer.SetAttributeFilter("ID = '{}'".format(item['sourceId']))
                    except: pass
                    if self.FRED.layer is None or len(self.FRED.layer) == 0:
                        bbox = item['boundingBox']
                        geom = regions.region2geom([bbox['minX'], bbox['maxX'], bbox['minY'], bbox['maxY']])

                        if item['format'] == 'IMG' or item['format'] == 'GeoTIFF':
                            tnm_ds = 'raster'
                        elif item['format'] == 'LAZ' or item['format'] == 'LAS':
                            tnm_ds = 'lidar'
                        else: tnm_ds = 'tnm'

                        if geom is not None:
                            self.FRED._add_survey(Name = item['title'], ID = item['sourceId'], Agency = 'USGS', Date = item['publicationDate'],
                                                  MetadataLink = item['metaUrl'], MetadataDate = item['dateCreated'],
                                                  DataLink = item['downloadURL'], Link = item['sourceOriginId'], Resolution = item['extent'],
                                                  DataType = tnm_ds, DataSource = 'tnm', HorizontalDatum = h_epsg,
                                                  VerticalDatum = v_epsg, Etcetra = dsTag, Info = item['moreInfo'], geom = geom)
                offset += 100
                if total - offset <= 0: break
            if self._verbose: _prog.end(0, 'gathered {} products from {}'.format(total, dsTag))
                           
    def _parse_prods_results(self, r, f = None, e = None, q = None):
        for surv in _filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    yield[d, d.split('/')[-1], surv['DataType']]
        
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is list: [url, file-name, datatype]
    ## ==============================================                
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        '''yield the xyz data from the tnm fetch module'''
        if fetch_file(entry[0], entry[1], callback = self._stop, verbose = self._verbose) == 0:
            datatype = entry[-1]
            if datatype == 'raster':
                src_tnms = utils.p_unzip(entry[1], ['tif', 'img', 'gdal', 'asc', 'bag'])
                for src_tnm in src_tnms:
                    try:
                        src_ds = gdal.Open(src_tnm)
                        if src_ds is not None:
                            srcwin = gdalfun.gdal_srcwin(src_ds, regions.region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
                            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose, z_region = z_region):
                                if xyz[2] != 0:
                                    yield(xyz)
                        src_ds = None
                    except:
                        utils.echo_error_msg('could not read tnm data: {}'.format(src_tnm))
                        src_ds = None
                    utils.remove_glob(src_tnm)
        utils.remove_glob(entry[1])
            
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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'gmrt'
        self._info = '''The Global Multi-Resolution Topography (GMRT) synthesis is a multi-resolutional 
compilation of edited multibeam sonar data collected by scientists and institutions worldwide, that is 
reviewed, processed and gridded by the GMRT Team and merged into a single continuously updated compilation 
of global elevation data. The synthesis began in 1992 as the Ridge Multibeam Synthesis (RMBS), was expanded 
to include multibeam bathymetry data from the Southern Ocean, and now includes bathymetry from throughout 
the global and coastal oceans.'''
        self._title = '''The Global Multi-Reosolution Topography Data Synthesis (GMRT)'''
        self._usage = '''< gmrt:layer=gmrt_layer >
 :layer=[topo/topo-mask]'''
        self._urls = [self._gmrt_grid_url, self._gmrt_grid_metadata_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('GMRT-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'GMRT', ID = 'GMRT-1', Date = utils.this_year(),
                                  MetadataLink = self._gmrt_grid_metadata_url, MetadataDate = utils.this_year(),
                                  DataLink = self._gmrt_grid_urls_url, DataType = 'raster', DataSource = 'gmrt',
                                  HorizontalDatum = 3857, VerticalDatum = 1092, Info = self._info,
                                  geom = regions.region2geom([-180,180,-90,90]))
        self.FRED._close_ds()

    def _parse_results(self, fmt = 'geotiff', res = 'max', layer = 'topo'):
        '''Run the GMRT fetching module'''
        if layer != 'topo' and layer != 'topo-mask': layer = 'topo'
        for surv in _filter_FRED(self):
            _data = {'north':self.region[3], 'west':self.region[0],
                     'south':self.region[2], 'east':self.region[1],
                     'mformat':'json', 'resolution':res, 'format':fmt}
            
            _req = fetch_req(surv['DataLink'], params = _data, tries = 10, timeout = 2)
            if _req is not None and _req.status_code == 200:
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
                    outf = 'gmrt_{}_{}.{}'.format(opts['layer'], regions.region_format(url_region, 'fn'), utils.gdal_fext(opts['format']))
                    yield([this_url, outf, 'gmrt'])

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## `entry` is a an item from self._data_urls
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326, z_region = None, inc = None):
        src_gmrt = 'gmrt_tmp_{}.tif'.format(regions.region_format(self.region, 'fn'))
        if fetch_file(entry[0], src_gmrt, callback = self._stop, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_gmrt)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, verbose = self._verbose, warp = epsg, z_region = z_region):
                    yield(xyz)
            src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_gmrt))
        utils.remove_glob(src_gmrt)

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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'mar_grav'
        self._info = 'Elevation data from Scripps Marine Gravity dataset.'
        self._title = '''Marine Gravity from Sattelite Altimetry topographic data.'''
        self._usage = '''< mar_grav >'''
        self._urls = [self._mar_grav_info_url, self._mar_grav_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('MAR_GRAV-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'MAR_GRAV', ID = 'MAR_GRAV-1', Date = utils.this_year(),
                                  MetadataLink = self._mar_grav_info_url, MetadataDate = utils.this_year(),
                                  DataLink = self._mar_grav_url, DataType = 'xyz', DataSource = 'mar_grav',
                                  HorizontalDatum = 4326, VerticalDatum = 1092, Info = self._info,
                                  geom = regions.region2geom([-180,180,-90,90]))
        self.FRED._close_ds()
        
    def _parse_results(self):
        '''Run the mar_grav fetching module.'''
        for surv in _filter_FRED(self):
            _data = {'north':self.region[3], 'west':self.region[0],
                     'south':self.region[2], 'east':self.region[1],
                     'mag':1}
            try:
                url_enc = urllib.urlencode(_data)
            except: url_enc = urllib.parse.urlencode(_data)
            url = '{}?{}'.format(surv['DataLink'], url_enc)
            outf = 'mar_grav_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            yield([url, outf, 'mar_grav'])
        
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = self._stop, verbose = self._verbose) == 0:
            xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            xyzc['x-off'] = -360
            xyzc['verbose'] = True
            xyzc['warp'] = epsg
            xyzc['name'] = '<mar_grav data-stream>'
            if z_region is not None:
                xyzc['upper_limit'] = z_region[1]
                xyzc['lower_limit'] = z_region[0]
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
    
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
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'srtm_plus'
        self._info = 'Global Bathymetry and Topography at 15 Arc Sec:SRTM15+'
        self._title = '''SRTM15+ elevation data (Scripps).'''
        self._usage = '''< srtm_plus >'''
        self._urls = [self._srtm_info_url, self._srtm_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('SRTM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'SRTM15+', ID = 'SRTM-1', Date = utils.this_year(),
                                  MetadataLink = self._srtm_info_url, MetadataDate = utils.this_year(),
                                  DataLink = self._srtm_url, DataType = 'xyz', DataSource = 'srtm_plus',
                                  HorizontalDatum = 4326, VerticalDatum = 1092, Info = self._info,
                                  geom = regions.region2geom([-180,180,-90,90]))
        self.FRED._close_ds()
        
    def _parse_results(self):
        '''Run the srtm+ fetching module.'''
        for surv in _filter_FRED(self):
            _data = {'north':self.region[3], 'west':self.region[0],
                     'south':self.region[2], 'east':self.region[1]}
            try:
                url_enc = urllib.urlencode(_data)
            except: url_enc = urllib.parse.urlencode(_data)
            url = '{}?{}'.format(surv['DataLink'], url_enc)
            outf = 'srtm_plus_{}.xyz'.format(regions.region_format(self.region, 'fn'))
            yield([url, outf, 'mar_grav'])
        
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        if fetch_file(entry[0], os.path.basename(entry[1]), callback = self._stop, verbose = self._verbose) == 0:
            xyzc = copy.deepcopy(xyzfun._xyz_config)
            xyzc['skip'] = 1
            xyzc['x-off'] = -360
            xyzc['verbose'] = True
            xyzc['warp'] = epsg
            xyzc['name'] = '<mar_grav data-stream>'
            if z_region is not None:
                xyzc['upper_limit'] = z_region[1]
                xyzc['lower_limit'] = z_region[0]
            with open(os.path.basename(entry[1]), 'r') as xyzf:
                for xyz in xyzfun.xyz_parse(xyzf, xyz_c = xyzc, verbose = self._verbose):
                    yield(xyz)
        utils.remove_glob(os.path.basename(entry[1]))
        
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
        self._emodnet_help_url = 'https://portal.emodnet-bathymetry.eu/help/help.html'
        self._outdir = os.path.join(os.getcwd(), 'emodnet')

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'emodnet'
        self._info = 'European Bathymetry/Topographic data from EMODNET'
        self._title = '''EMODNET Elevation Data.'''
        self._usage = '''< emodnet >'''
        self._urls = [self._emodnet_help_url, self._emodnet_grid_url, self._emodnet_grid_cap]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()

    def _update(self):
        self.FRED._open_ds(1)
        emod_wcs = WCS(self._emodnet_grid_url)
        contents = emod_wcs._contents()
        if self._verbose: _prog = utils._progress('Scanning {} WCS coverages from {}...'.format(len(contents), self._emodnet_grid_url))
        for i, layer in enumerate(contents):
            if self._verbose: _prog.update_perc((i, len(contents)))
            self.FRED._attribute_filter(["ID = '{}'".format(layer['CoverageId'][0])])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                d = emod_wcs._describe_coverage(layer['CoverageId'][0])
                if d is not None:
                    ds_region = emod_wcs._get_coverage_region(d)
                    geom = regions.region2geom(ds_region)
                    url = emod_wcs._get_coverage_url(layer['CoverageId'][0], region = ds_region)
                    self.FRED._add_survey(Name = d['name'][0], ID = layer['CoverageId'][0], Date = utils.this_year(), MetadataLink = layer['Metadata'],
                                          MetadataDate = utils.this_year(), DataLink = url, DataType = 'raster',
                                          DataSource = 'emodnet', HorizontalDatum = 4326, VerticalDatum = 1092,
                                          Info = layer['Abstract'], geom = geom)
        if self._verbose: _prog.end(0, 'Scanned {} WCS coverages from {}'.format(len(contents), self._emodnet_grid_url))
        self.FRED._close_ds()
        
    def _parse_results(self):        
        emod_wcs = WCS(self._emodnet_grid_url)
        for surv in _filter_FRED(self):
            d = emod_wcs._describe_coverage(surv['ID'])
            if d is not None:
                ds_region = emod_wcs._get_coverage_region(d)
                if regions.regions_intersect_ogr_p(self.region, ds_region):
                    emod_url = emod_wcs._get_coverage_url(surv['ID'], region = self.region)
                    outf = 'emodnet_{}.tif'.format(regions.region_format(self.region, 'fn'))
                    yield([emod_url, outf, surv['DataType']])

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None, z_region = None, inc = None):
        src_emodnet = 'emodnet_tmp.tif'
        if fetch_file(entry[0], src_emodnet, callback = self._stop, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_emodnet)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose, z_region = z_region):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_emodnet))
        utils.remove_glob(src_emodnet)

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
        self._hrdem_info_url = 'https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995#wb-auto-6'
        self._local_ref_vector = 'hrdem.gmt'
        if not os.path.exists(self._local_ref_vector):
            self.FRED = os.path.join(fetchdata, 'hrdem.gmt')
        else: self.FRED = self._local_ref_vector
        self._outdir = os.path.join(os.getcwd(), 'hrdem')
        
        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name = 'hrdem'
        self._info = '''Collection of lidar-derived DTMs across Canada.'''
        self._title = '''High Resolution DEMs from NCAR'''
        self._usage = '''< hrdem >'''
        self._urls = [self._hrdem_info_url, self._hrdem_footprints_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds()
        v_zip = os.path.basename(self._hrdem_footprints_url)
        status = fetch_ftp_file(self._hrdem_footprints_url, v_zip, verbose = True)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        v_shp = None
        for v in v_shps:
            if '.shp' in v: v_shp = v
        shp_regions = regions.gdal_ogr_regions(v_shp)
        shp_region = []
        for this_region in shp_regions:
            if len(shp_region) > 0:
                shp_region = regions.regions_merge(shp_region, this_region)
            else: shp_region = this_region
        geom = regions.region2geom(shp_region)
        
        self.FRED._attribute_filter(["ID = '{}'".format('HRDEM-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'High-Resolution DEM (Canada)', ID = 'HRDEM-1', Agency = 'NRCAN', Date = utils.this_year(),
                                  MetadataLink = self._hrdem_info_url, MetadataDate = utils.this_year(),
                                  DataLink = self._hrdem_footprints_url, IndexLink = self._hrdem_footprints_url,
                                  DataType = 'raster', DataSource = 'hrdem', Info = 'Canada Only', geom = geom)
        utils.remove_glob(v_zip, *v_shps)
        self.FRED._close_ds()

    def _parse_results(self):
        for surv in _filter_FRED(self):
            status = fetch_ftp_file(surv['IndexLink'], v_zip, verbose = self._verbose)
            v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
            v_shp = None
            for v in v_shps:
                if v.split('.')[-1] == 'shp':
                    v_shp = v
                    break
            try:
                v_ds = ogr.Open(v_shp)
            except:
                v_ds = None
                status = -1
            if v_ds is not None:
                layer = v_ds.GetLayer()
                try:
                    self.FRED.layer.SetAttributeFilter("Name = '{}'".format(name))
                except: pass                
                fcount = layer.GetFeatureCount()
                for f in range(0, fcount):
                    feature = layer[f]
                    if data_link is not None:
                        geom = feature.GetGeometryRef()
                        if geom.Intersects(regions.region2geom(self.region)):
                            data_link = feature.GetField('Ftp_dtm').replace('http', 'ftp')
                            yield([data_link, data_link.split('/')[-1], surv['DataType']])
            utils.remove_glob(v_zip, *v_shps)
                                
    ## ==============================================
    ## _update_all() and _parse_results_all() will update FRED with all the data
    ## from the hrdem footprints, which can take a long time and use a lot of FRED
    ## space, which is mostly unncessary, as the footprints are fairly quick to download
    ## and process on the spot for the most part.
    ## use thse *_all functions to revert back to having all the data in the FRED...
    ## ==============================================
    def _update_all(self):
        self.FRED._open_ds(1)
        v_zip = os.path.basename(self._hrdem_footprints_url)
        status = fetch_ftp_file(self._hrdem_footprints_url, v_zip, verbose = True)
        v_shps = utils.p_unzip(v_zip, ['shp', 'shx', 'dbf', 'prj'])
        v_shp = None
        for v in v_shps:
            if '.shp' in v: v_shp = v
        try:
            v_ds = ogr.Open(v_shp)
        except:
            v_ds = None
            status = -1
        if v_ds is not None:
            layer = v_ds.GetLayer()
            fcount = layer.GetFeatureCount()
            if self._verbose: _prog = utils._progress('scanning {} datasets...'.format(fcount))
            for f in range(0, fcount):
                feature = layer[f]
                name = feature.GetField('Tile_name')
                if self._verbose: _prog.update_perc((f, fcount))
                try:
                    self.FRED.layer.SetAttributeFilter("Name = '{}'".format(name))
                except: pass
                if self.FRED.layer is None or len(self.FRED.layer) == 0:
                    data_link = feature.GetField('Ftp_dtm')
                    if data_link is not None:
                        geom = feature.GetGeometryRef()
                        self.FRED._add_survey(Name = name, ID = feature.GetField('Project'), Agency = 'NRCAN', Date = utils.this_year(),
                                              MetadataLink = feature.GetField('Meta_dtm'), MetadataDate = utils.this_year(),
                                              DataLink = data_link.replace('http', 'ftp'), IndexLink = self._hrdem_footprints_url,
                                              DataType = 'raster', DataSource = 'hrdem', HorizontalDatum = feature.GetField('Coord_Sys').split(':')[-1],
                                              Info = feature.GetField('Provider'), geom = geom)

            if self._verbose: _prog.end('scanned {} datasets.'.format(fcount))
        utils.remove_glob(v_zip, *v_shps)
        self.FRED._close_ds()

    def _parse_results_all(self):
        '''Parse the results of a filtered FRED'''
        for surv in _filter_FRED(self):
            for d in surv['DataLink'].split(','):
                if d != '':
                    yield[d, d.split('/')[-1], surv['DataType']]
                    
    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = 4326, z_region = None, inc = None):
        src_dc = os.path.basename(entry[1])
        src_ext = src_dc.split('.')[-1]
        try:
            src_ds = gdal.Open(entry[0])
        except Exception as e:
            fetch_file(entry[0], src_dc, callback = self._stop, verbose = self._verbose)
            try:
                src_ds = gdal.Open(src_dc)
            except Exception as e:
                utils.echo_error_msg('could not read hrdem raster file: {}, {}'.format(entry[0], e))
                src_ds = None
        except Exception as e:
            utils.echo_error_msg('could not read hrdem raster file: {}, {}'.format(entry[0], e))
            src_ds = None

        if src_ds is not None:
            srcwin = gdalfun.gdal_srcwin(src_ds, regions.region_warp(self.region, s_warp = epsg, t_warp = gdalfun.gdal_getEPSG(src_ds)))
            for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose, z_region = z_region):
                yield(xyz)
        src_ds = None
        utils.remove_glob(src_dc)

## =============================================================================
##
## CHS Fetch
##
## fetch bathymetric soundings from the Canadian Hydrographic Service (CHS) - Canada Only
## https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d
##
## NONNA 10 and NONNA 100
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
        self._chs_url = 'https://data.chs-shc.ca/geoserver/wcs?'
        self._chs_grid_cap = 'https://data.chs-shc.ca/geoserver/wcs?request=GetCapabilities&service=WMS'
        self._chs_info_url = 'https://open.canada.ca/data/en/dataset/d3881c4c-650d-4070-bf9b-1e00aabf0a1d'
        self._chs_api_url = "https://geoportal.gc.ca/arcgis/rest/services/FGP/CHS_NONNA_100/MapServer/0/query?"
        self._outdir = os.path.join(os.getcwd(), 'chs')

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._surveys = []
        self._stop = callback

        self._name = 'chs'
        self._info = '''CHS NONNA 10m and 100m Bathymetric Survey Grids; Non-Navigational gridded bathymetric data based on charts and soundings.'''
        self._title = '''Bathymetric data from CHS'''
        self._usage = '''< chs >'''
        self._urls = [self._chs_info_url, self._chs_url, self._chs_grid_cap]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        chs_wcs = WCS(self._chs_url)
        contents = chs_wcs._contents()
        if self._verbose: _prog = utils._progress('Scanning {} WCS coverages from {}...'.format(len(contents), self._chs_url))
        for i, layer in enumerate(contents):
            if self._verbose: _prog.update_perc((i, len(contents)))
            self.FRED._attribute_filter(["ID = '{}'".format(layer['CoverageId'][0])])
            if self.FRED.layer is None or len(self.FRED.layer) == 0:
                
                d = chs_wcs._describe_coverage(layer['CoverageId'][0])
                if d is not None:
                    ds_region = chs_wcs._get_coverage_region(d)
                    geom = regions.region2geom(ds_region)
                    url = chs_wcs._get_coverage_url(layer['CoverageId'][0], region = ds_region)
                    try:
                        name = d['name'][0]
                    except: name = d['CoverageId'][0]
                    try:
                        meta = layer['Metadata']
                    except: meta = None
                    try:
                        info = layer['Abstract']
                    except: info = None
                    self.FRED._add_survey(Name = name, ID = layer['CoverageId'][0], Date = utils.this_year(), MetadataLink = meta,
                                          MetadataDate = utils.this_year(), DataLink = url, DataType = 'raster',
                                          DataSource = 'chs', HorizontalDatum = 4326, VerticalDatum = 1092,
                                          Info = info, geom = geom)
        if self._verbose: _prog.end(0, 'Scanned {} WCS coverages from {}'.format(len(contents), self._chs_url))
        self.FRED._close_ds()
        
    def _parse_results(self):        
        chs_wcs = WCS(self._chs_url)
        for surv in _filter_FRED(self):
            d = chs_wcs._describe_coverage(surv['ID'])
            if d is not None:
                ds_region = chs_wcs._get_coverage_region(d)
                if regions.regions_intersect_ogr_p(self.region, ds_region):
                    chs_url = chs_wcs._get_coverage_url(surv['ID'], region = self.region)
                    outf = 'chs_{}.tif'.format(regions.region_format(self.region, 'fn'))
                    yield([chs_url, outf, surv['DataType']])

    ## ==============================================
    ## Process results to xyz
    ## yield functions are used in waffles/datalists
    ## as well as for processing incoming fetched data.
    ## ==============================================    
    def _yield_xyz(self, entry, epsg = None, z_region = None):
        src_chs = 'chs_tmp.tif'
        if fetch_file(entry[0], src_chs, callback = self._stop, verbose = self._verbose) == 0:
            try:
                src_ds = gdal.Open(src_chs)
            except: src_ds = None
            if src_ds is not None:
                srcwin = gdalfun.gdal_srcwin(src_ds, self.region)
                for xyz in gdalfun.gdal_parse(src_ds, srcwin = srcwin, warp = epsg, verbose = self._verbose, z_region = z_region):
                    yield(xyz)
                src_ds = None
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_chs))
        utils.remove_glob(src_chs)
            
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
        self._ngs_url = 'http://geodesy.noaa.gov'
        self._ngs_search_url = 'https://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')

        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        ## ==============================================
        ## Metadata, usage, etc.
        ## ==============================================
        self._name = 'ngs'
        self._info = '''Monument data from NOAA's Nagional Geodetic Survey (NGS) monument dataset.'''
        self._title = '''NOAA NGS Monuments'''
        self._usage = '''< ngs >'''
        self._urls = [self._ngs_url, self._ngs_search_url]
        self._update_if_not_in_FRED()
        
    def _update_if_not_in_FRED(self):
        self.FRED._open_ds()
        self.FRED._attribute_filter("DataSource = '{}'".format(self._name))
        if len(self.FRED.layer) == 0:
            self.FRED._close_ds()
            self._update()
        self.FRED._close_ds()
        
    def _update(self):
        self.FRED._open_ds(1)
        self.FRED._attribute_filter(["ID = '{}'".format('NGS-1')])
        if self.FRED.layer is None or len(self.FRED.layer) == 0:
            self.FRED._add_survey(Name = 'NGS Monuments', ID = 'NGS-1', Agency = 'NGS', Date = utils.this_year(),
                                  MetadataDate = utils.this_year(), DataLink = self._ngs_search_url,
                                  DataType = 'raster', DataSource = self._name, Info = self._info,
                                  geom = regions.region2geom([-162, -60, 16, 73]))
        self.FRED._close_ds()
        
    def _parse_results(self):
        '''Parse the NGS (monuments) data.'''
        for surv in _filter_FRED(self):
            _data = { 'maxlon':self.region[0], 'minlon':self.region[1],
                      'maxlat':self.region[3], 'minlat':self.region[2] }
            try:
                ngs_data = urllib.urlencode(_data)
            except: ngs_data = urllib.parse.urlencode(_data)
            ngs_url = '{}{}'.format(surv['DataLink'], ngs_data)
            yield([ngs_url, 'ngs_results_{}.json'.format(regions.region_format(self.region, 'fn')), 'ngs'])

    def _yield_xyz(self, entry, epsg = None, z_region = None):
        src_ngs = 'ngs_tmp.json'
        r = []
        if fetch_file(entry[0], src_ngs, callback = self._stop, verbose = self._verbose) == 0:
            with open(src_ngs, 'r') as json_file: r = json.load(json_file)
            for mm in r: yield([mm['lon'], mm['lat'], mm['geoidHt']])
        else: utils.echo_error_msg('failed to fetch remote file, {}...'.format(src_ngs))
        utils.remove_glob(src_ngs)

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
        self.osm_types = {'highway': ['LINESTRING'], 'waterway': ['LINESTRING'], 'building': ['POLYGON']}
        
        self.where = where
        self.region = extent
        self.FRED = FRED(verbose = self._verbose)
        self._stop = callback

        self._name ='osm'
        self._info = '''Various datasets from Open Street Map Overpass API'''
        self._title = '''Open Street Map (OSM) data'''
        self._usage = '''< osm >'''
        self._urls = [self._osm_api]
 
    def _update(self):
        pass

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
## fetches processing (datalists fmt: -4)
## ==============================================
def fetch_inf_entry(entry = [], warp = None):
    out_inf = {'minmax': [-180,180,-90,90], 'name': entry[0], 'pts': None, 'wkt': regions.region2wkt([-180,180,-90,90])}
    if warp is not None:
        out_inf['minmax'] = regions.region_warp(out_inf['minmax'], 4326, warp)
        out_inf['wkt'] = regions.region2wkt(out_inf['minmax'])
    return(out_inf)

def fetch_yield_entry(entry = ['nos:datatype=xyz'], region = None, warp = None, verbose = False, z_region = None):
    '''yield the xyz data from the fetch module datalist entry

    yields [x, y, z, <w, ...>]'''
    fetch_mod = entry[0].split(':')[0]
    fetch_args = entry[0].split(':')[1:]

    if region is None:
        region = _fetch_modules[fetch_mod](None, [], None, False).FRED._get_region([], [fetch_mod])

    region = regions.region_warp(region, src_epsg = warp, dst_epsg = 4326)
    fl = _fetch_modules[fetch_mod](regions.region_buffer(region, 5, pct = True), [], lambda: False, verbose)
    args_d = utils.args2dict(fetch_args, {})

    for e in fl._parse_results():
        for xyz in fl._yield_xyz(e, epsg = warp, z_region = z_region):
            yield(xyz + [entry[2]] if entry[2] is not None else xyz)

def fetch_dump_entry(entry = ['nos:datatype=nos'], dst_port = sys.stdout, region = None, warp = None, verbose = False, z_region = None):
    '''dump the xyz data from the fetch module datalist entry to dst_port'''
    for xyz in fetch_yield_entry(entry = entry, region = region, warp = warp, verbose = verbose, z_region = z_region):
        xyz_line(xyz, dst_port, False)

def fetch_dump_xyz(parsed_entry, module = None, epsg = 4326, z_region = None, inc = None, dst_port = sys.stdout):
    '''dump the parsed entry to xyz
    use <fetch_module>._parse_results() to generate parsed entry'''
    if module == None: module = _fetch_modules[parsed_entry[4]](None, [], None, False)
    for xyz in module._yield_xyz([parsed_entry[0], parsed_entry[1], parsed_entry[-1]], epsg = epsg, z_region = z_region, inc = inc):
        xyzfun.xyz_line(xyz, dst_port, False)

def fetch_filter_results(region = None, where = None, mods = [], verbose = False):
    '''return the filtered results from FRED'''
    return(FRED(verbose)._filter(region, where, mods))

def fetch_update_FRED(surveys):
    FRED(verbose)._add_surveys(surveys)

## ==============================================
## Run the Fetch Module(s)
## ==============================================
def fetch(fg):
    out_results = []
    stop_threads = False
    if fg is None:
        utils.echo_error_msg('invalid configuration, {}'.format(fg))
        sys.exit(-1)
        
    for fetch_mod in fg['mods'].keys():
        if stop_threads: break
        status = 0
        args = tuple(fg['mods'][fetch_mod])
        
        if fg['region'] is None or fg['region'][0] is None:
            this_region = None
        else: this_region = regions.region_buffer(fg['region'], 5, pct = True)
        fl = _fetch_modules[fetch_mod](this_region, fg['where'], lambda: stop_threads, fg['verbose'])
        args_d = utils.args2dict(args)
        if fg['verbose']: _prog = utils._progress('running FETCHES module < {} > with {} [{}]...'.format(fetch_mod, args, args_d))

        ## ==============================================
        ## Run update in a thread to cleanly close FRED
        ## ============================================== 
        if fg['update_p']:
            t = threading.Thread(target = fl._update, args = ())
            t.daemon = True
            try:
                t.start()
                while True:
                    time.sleep(2)
                    if fg['verbose']: _prog.update(msg='updating FETCHES module < {} >...'.format(fetch_mod))
                    if not t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit):
                utils.echo_error_msg('user breakage...please wait while fetches exits.')
                stop_threads = True
                status = -1
            except Exception as e:
                utils.echo_error_msg(e)
                stop_threads = True
                status = -1
            t.join()
            if fg['verbose']: _prog.end(status, 'updated FETCHES module < {} >.'.format(fetch_mod))
            
        if this_region is not None:
            if fg['index_p']:
                _filter_FRED_index(fl)
                continue
            ## ==============================================    
            ## Fetch the data
            ## fetching will be done in a queue with 3 threads fetching the data at a time.
            ##
            ## fetch_p must be true to fetch the data to _outdir
            ## list_p will list the urls
            ## dump_p will dump the xyz data from the fetch module to stdout
            ## proc_p will output the xyz data to file in _outdir
            ## ==============================================
            fr = fetch_results(fl._parse_results(**args_d), this_region, fl._outdir, fl, fg, lambda: stop_threads)
            fr.daemon = True
            try:
                fr.start()
                while True:
                    time.sleep(2)
                    if fg['verbose']: _prog.update()
                    if not fr.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit):
                utils.echo_error_msg('user breakage...please wait for while fetches exits.')
                stop_threads = True
                status = -1
                while not fr.fetch_q.empty():
                    try:
                        fr.fetch_q.get(False)
                    except Empty: continue
                    fr.fetch_q.task_done()
            except Exception as e:
                stop_threads = True
                status = -1
                utils.echo_error_msg(e)
            fr.join()
        if fg['verbose']: _prog.end(status, 'ran FETCHES module < {} > with {} [{}]...'.format(fetch_mod, args, fg))
        
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
  -E, --increment\tBlockmedian/mean CELL-SIZE in native units or GMT-style increments.
  -Z, --z-region\t\tRestrict data processing to records that fall within the z-region
\t\t\tUse '-' to indicate no bounding range; e.g. -Z-/0 will restrict processing to data
\t\t\trecords whose z value is below zero.
  -F, --fg-config\tA fetches config JSON file. If supplied, will overwrite all other options.
\t\t\tgenerate a fetches_config JSON file using the --config flag.

  -l, --list\t\tReturn a list of fetch URLs in the given region.
  -d, --dump\t\tDump the XYZ elevation data in WGS84 to stdout.
  -p, --process\t\tProcess fetched elevation data to ASCII XYZ format in WGS84.
\t\t\tIf -E or -Z are set, processing will use those switches in data processing.
  -u, --update\t\tUpdate the Fetches Remote Elevation Datalist.
  -i, --index\t\tPrint the fetch FRED results in the given region.

  --help\t\tPrint the usage text
  --config\t\tSave the fetches config JSON
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
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]), 
           os.path.basename(sys.argv[0]),
           os.path.basename(sys.argv[0]))

def fetches_cli(argv = sys.argv):
    '''run fetches from command-line
    generates a fetches_config from the command-line options
    and either outputs or runs the fetches_config
    on each region supplied (multiple regions can be supplied
    by using a vector file as the -R option.)
    See `fetches_cli_usage` for full cli options.'''
    fg = fetches_config()
    fg_user = None
    status = 0
    region = None
    mod_opts = {}
    verbose = True
    want_config = False
    use_local = False
    
    ## ==============================================
    ## process Command-Line
    ## ==============================================
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            region = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            region = str(arg[2:])
        elif arg == '--where' or arg == '-W':
            fg['where'].append(argv[i + 1])
            i = i + 1
        elif arg == '--increment' or arg == '-E':
            fg['inc'] = gmtfun.gmt_inc2inc(argv[i+1])
            i = i + 1
        elif arg[:2] == '-E': fg['inc'] = gmtfun.gmt_inc2inc(arg[2:].split(':')[0])
        elif arg == '--z-range' or arg == '-Z':
            zr = argv[i + 1].split('/')
            if len(zr) > 1:
                fg['z_region'] = [None if x == '-' else float(x) for x in zr]
            i = i + 1
        elif arg[:2] == '-Z':
            zr = arg[2:].split('/')
            if len(zr) > 1:
                fg['z_region'] = [None if x == '-' else float(x) for x in zr]
        elif arg == '--list' or arg == '-l':
            fg['list_p'] = True
            fg['fetch_p'] = False
        elif arg == '--index' or arg == '-i':
            fg['index_p'] = True
            fg['fetch_p'] = False
        elif arg == '--process' or arg == '-p':
            fg['proc_p'] = True
            fg['fetch_p'] = True
        elif arg == '--dump' or arg == '-d':
            fg['dump_p'] = True
            fg['fetch_p'] = True
        elif arg == '--update' or arg == '-u':
            fg['update_p'] = True
        elif arg == '--fg-config' or arg == '-F':
            fg_user = argv[i + 1]
            i += 1
        elif arg[:2] == '-F': fg_user = arg[2:]
        elif arg == '--help' or arg == '-h':
            sys.stderr.write(_usage)
            sys.exit(1)
        elif arg == '--version' or arg == '-v':
            sys.stdout.write('{}\n'.format( _version))
            sys.exit(0)
        elif arg == '--config': want_config = True
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
    ## load the user fg json and run fetches with that.
    ## ==============================================
    if fg_user is not None:
        if os.path.exists(fg_user):
            try:
                with open(fg_user, 'r') as fgj:
                    fg = json.load(fgj)
                    fg = fetches_config(**fg)
                    fetch(fg)
                    sys.exit(0)
            except Exception as e:
                fg = fetches_config(**fg)
                utils.echo_error_msg(e)
                sys.exit(-1)
        else:
            utils.echo_error_msg('specified json file does not exist, {}'.format(fg_user))
            sys.exit(0)
        
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

    fg['mods'] = mod_opts

    ## ==============================================
    ## process input region(s)
    ## ==============================================
    if region is not None:
        try:
            these_regions = [[float(x) for x in region.split('/')]]
        except ValueError: these_regions = regions.gdal_ogr_regions(region)
        if len(these_regions) == 0: status = -1
        for this_region in these_regions:
            if not regions.region_valid_p(this_region): status = -1
        utils.echo_msg('loaded {} region(s)'.format(len(these_regions)))
    else: these_regions = [None]
    if len(these_regions) == 0: utils.echo_error_msg('failed to parse region(s), {}'.format(region))

    ## ==============================================
    ## fetch some data in each of the input regions
    ## ==============================================
    these_fgs = []
    for rn, this_region in enumerate(these_regions):
        tfg = fetches_config(**fg)
        if this_region is None: tfg['fetch_p'] = False
        tfg['region'] = this_region
        tfg['name'] = 'fetches_{}_{}'.format('' if this_region is None else regions.region_format(this_region, 'fn'), utils.this_year())
        tfg = fetches_config(**tfg)
        if want_config:
            if tfg is not None:
                utils.echo_msg(json.dumps(tfg, indent = 2, sort_keys = True))
                with open('{}.json'.format(tfg['name']), 'w') as fg_json:
                    utils.echo_msg('generating fetches config file: {}.json'.format(tfg['name']))
                    fg_json.write(json.dumps(tfg, indent = 2, sort_keys = True))
            else: utils.echo_error_msg('could not parse config.')
        else: fetch(tfg)
### End
