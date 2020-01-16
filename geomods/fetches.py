### fetches.py
##
## Copyright (c) 2012 - 2020 Matthew Love <matthew.love@colorado.edu>
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
import datetime
import math

import requests
import lxml.html as lh
import lxml.etree

import zipfile
import gzip
import csv
import json
import threading

try:
    import numpy as np
except ImportError:
    print 'NumPy must be installed'
    sys.exit(-1)

try:
    import osgeo.ogr as ogr
    import osgeo.gdal as gdal
    has_gdalpy = True
except ImportError:
    try:
        import ogr
        import gdal
        has_gdalpy = True
    except ImportError:
        has_gdalpy = False
        print 'GDAL must be installed'
        sys.exit(-1)

import regions
import datalists
import gdalfun
import utils

_version = '0.1.7'

gdal.PushErrorHandler('CPLQuietErrorHandler')

## =============================================================================
##
## Fetching Functions
##
## Generic fetching functions, etc.
##
## =============================================================================

r_headers = { 'User-Agent': 'GeoMods: Fetches v%s' %(_version) }

namespaces = {'gmd': 'http://www.isotc211.org/2005/gmd', 
              'gmi': 'http://www.isotc211.org/2005/gmi', 
              'gco': 'http://www.isotc211.org/2005/gco',
              'gml': 'http://www.isotc211.org/2005/gml'}

def fetch_file(src_url, dst_fn, params = None):
    '''fetch src_url and save to dst_fn'''

    if not os.path.exists(os.path.dirname(dst_fn)):
        os.makedirs(os.path.dirname(dst_fn))

    req = requests.get(src_url, stream = True, params = params, headers = r_headers)

    if req:
        with open(dst_fn, 'wb') as local_file:
            for chunk in req.iter_content(chunk_size = 50000):
                local_file.write(chunk)
        return(0)
    else: return(-1)

def fetch_req(src_url, params = None, tries = 3):
    '''fetch src_url and return the requests object'''

    if tries <= 0: return(None)
    try:
        return(requests.get(src_url, stream = True, params = params, timeout = 1, headers = r_headers))
    except: return(fetch_req(src_url, params = params, tries = tries - 1))

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
        return(list(csv.reader(req.text.split('\n'), delimiter = ',')))
    else: return(None)

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

    if update:
        ds = ogr.GetDriverByName('GMT').Open(src_vec, 1)
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
    def __init__(self, extent = None, callback = lambda: True):
        '''Fetch elevation data from the Digital Coast'''

        self._dc_htdata_url = 'https://coast.noaa.gov/htdata/'
        self._dc_dav_id = 'https://coast.noaa.gov/dataviewer/#/lidar/search/where:ID='
        self._dc_dirs = ['lidar1_z', 'lidar2_z', 'raster2']

        self._ref_vector = os.path.join(fetchdata, 'dc.gmt')
        self._outdir = os.path.join(os.getcwd(), 'dc')
        self._log_fn = 'fetch_dc_%s.log' %(datetime.datetime.now().strftime('%d%Y'))

        self._status = 0
        self._surveys = []
        self._results = []

        self._index = False
        self._want_proc = True

        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)
        else: self._status = -1

        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self.stop = callback
        self.region = extent

    def _log_survey(self, surv_url):
        with open(self._log_fn, 'a') as local_file:
            local_file.write(surv_url + "\n")
                
    def search_gmt(self, filters=[]):
        '''Search for data in the reference vector file'''

        gmt1 = ogr.Open(self._ref_vector)
        layer = gmt1.GetLayer(0)

        if len(filters) > 0:
            for filt in [filters]:
                layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            if not self.stop():
                geom = feature1.GetGeometryRef()
                if self._boundsGeom.Intersects(geom):
                    surv_url = feature1.GetField('Data')
                    surv_dt = feature1.GetField('Datatype')
                    suh = fetch_html(surv_url)
                    if suh is None: self._status = -1
                    
                    if self._status == 0:
                        if 'lidar' in surv_dt:

                            ## ==============================================
                            ## Lidar data has a minmax.csv file to get extent
                            ## ==============================================

                            scsv = suh.xpath('//a[contains(@href, ".csv")]/@href')[0]
                            dc_csv = fetch_csv(surv_url + scsv)

                            for tile in dc_csv:
                                try:
                                    tb = [float(tile[1]), float(tile[2]), float(tile[3]), float(tile[4])]
                                    tile_geom = bounds2geom(tb)
                                    if tile_geom.Intersects(self._boundsGeom):
                                        self._results.append(os.path.join(surv_url, tile[0]))
                                except: pass

                        elif 'raster' in surv_dt:
                            
                            ## ==============================================
                            ## Raster data has a tileindex shapefile 
                            ## to get extent
                            ## ==============================================

                            sshpz = suh.xpath('//a[contains(@href, ".zip")]/@href')[0]
                            fetch_file(surv_url + sshpz, os.path.join('.', sshpz))

                            try:
                                zip_ref = zipfile.ZipFile(sshpz)
                                zip_ref.extractall('dc_tile_index')
                                zip_ref.close()
                            except BadZipFile:
                                self.stop = lambda: True
                                break

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
                                        self._results.append(tile_url)

                                shp1 = slay1 = None

                            for i in ti:
                                ts = os.remove(os.path.join('dc_tile_index/', i))
                            os.removedirs(os.path.join('.', 'dc_tile_index'))
                            os.remove(os.path.join('.', sshpz))
        if len(self._results) == 0: self._status = -1
        gmt1 = layer = None

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

                            wl = xml_doc.find('.//gmd:westBoundLongitude/gco:Decimal', namespaces = namespaces)
                            el = xml_doc.find('.//gmd:eastBoundLongitude/gco:Decimal', namespaces = namespaces)
                            sl = xml_doc.find('.//gmd:southBoundLatitude/gco:Decimal', namespaces = namespaces)
                            nl = xml_doc.find('.//gmd:northBoundLatitude/gco:Decimal', namespaces = namespaces)
                            if wl is not None:
                                obbox = bounds2geom([float(wl.text), float(el.text), float(sl.text), float(nl.text)])

                                try: 
                                    odate = int(dc['Year'][:4])
                                except: odate = 1900

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

    def proc_data(self, s_dir, s_fn, o_fn, s_t):
        '''Process a fetched file to xyz and add it to its datalist.
        uses `las2txt` from lastools for las/laz files and
        gdal (geomods.gdalfun) for raster files.'''

        status = 0
        xyz_dir = os.path.join(self._outdir, s_dir, 'xyz')
        
        if not os.path.exists(xyz_dir):
            os.makedirs(xyz_dir)

        o_fn_bn = os.path.basename(o_fn).split('.')[0]            
        o_fn_p_xyz = os.path.join(self._outdir, s_dir, '{}.xyz'.format(o_fn_bn))
        o_fn_xyz = os.path.join(xyz_dir, '{}_{}.xyz'.format(o_fn_bn, self.region.fn))
            
        if s_t == 'las' or s_t == 'laz':

            ## ==============================================
            ## Convert to XYZ
            ## ==============================================

            out, status = utils.run_cmd('las2txt -verbose -parse xyz -keep_class 2 29 -i {}'.format(o_fn), False, None)

            if status == 0:
                o_fn_bn = os.path.basename(o_fn).split('.')[0]

                o_fn_txt = os.path.join(self._outdir, s_dir, '{}.txt'.format(o_fn_bn))

                ## ==============================================
                ## Blockmedian the data
                ## ==============================================

                out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR=space', False, None)
                out, status = utils.run_cmd('gmt blockmedian {} -I.1111111s {} -r -V > {}'.format(o_fn_txt, self.region.gmt, o_fn_p_xyz), False, None)

                os.remove(o_fn_txt)

        elif s_t == 'tif' or s_t == 'img':

            ## ==============================================
            ## Convert to XYZ
            ## Raster data should first be transformed to
            ## WGS84 before dumping to xyz
            ## ==============================================

            gdalfun.dump(o_fn, o_fn_p_xyz)

        if status == 0:

            ## ==============================================
            ## Move processed xyz file to xyz directory
            ## ==============================================

            os.rename(o_fn_p_xyz, o_fn_xyz)

            ## ==============================================        
            ## Add xyz file to datalist
            ## ==============================================

            sdatalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(s_dir)))
            sdatalist._append_datafile('{}'.format(os.path.basename(o_fn_xyz)), 168, 1)
            sdatalist._reset()

            ## ==============================================
            ## Generate .inf file
            ## ==============================================

            out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, '{}.datalist'.format(s_dir))), False, None)

            os.remove(o_fn)

    def print_results(self):
        '''print the fetch results as a list of urls suitable 
        for use in `wget`, etc.'''
        
        if len(self._results) == 0:
            self._status = -1
        else:
            for row in self._results:
                print(row)
        
    def fetch_results(self):
        '''Fetch and possibly process the found fetch results.'''

        if len(self._results) == 0:
            self._status = -1
        else:
            for row in self._results:
                if not self.stop():

                    ## ==============================================
                    ## Parse the result row
                    ## ==============================================

                    surv_dir = row.split('/')[-2:][0]
                    if surv_dir == '.':
                        surv_dir = row.split('/')[-3:][0]

                    surv_fn = row.split('/')[-2:][1]
                    outf = os.path.join(self._outdir, surv_dir, surv_fn)

                    ## ==============================================
                    ## Fetch and Log
                    ## ==============================================

                    self._status = fetch_file(row, outf)
                    
                    if self._status == 0 and os.path.exists(outf):

                        ## ==============================================                    
                        ## validate downloaded file...
                        ## ==============================================

                        if open(outf, 'rb').read(4) == 'LASF':
                            self._status = 0
                        elif gdal.Open(outf, gdal.GA_ReadOnly):
                            self._status = 0
                        else: self._status = -1

                        self._log_survey(row)

                        ## ==============================================                    
                        ## Process the data if wanted
                        ## ==============================================

                        if self._want_proc:
                            try:
                                surv_t = surv_fn.split('.')[1]
                            except: surv_t = ''

                            t = threading.Thread(target = self.proc_data, args = (surv_dir, surv_fn, outf, surv_t))
                            t.start()

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
    def __init__(self, extent = None, callback = lambda: True):
        '''Fetch NOS BAG and XYZ sounding data from NOAA'''

        self._nos_xml_url = lambda nd: 'https://data.noaa.gov/waf/NOAA/NESDIS/NGDC/MGG/NOS/%siso_u/xml/' %(nd)
        self._nos_directories = ["B00001-B02000/", "D00001-D02000/", "F00001-F02000/", \
                                 "H00001-H02000/", "H02001-H04000/", "H04001-H06000/", \
                                 "H06001-H08000/", "H08001-H10000/", "H10001-H12000/", \
                                 "H12001-H14000/", "L00001-L02000/", "L02001-L04000/", \
                                 "T00001-T02000/", "W00001-W02000/"]

        self._outdir = os.path.join(os.getcwd(), 'nos')
        self._ref_vector = os.path.join(fetchdata, 'nos.gmt')

        if os.path.exists(self._ref_vector): 
            self._has_vector = True
        else: self._has_vector = False

        self._status = 0
        self._surveys = []
        self._xml_results = []
        self._results = []
        self.stop = callback

        self._want_proc = True

        self.region = extent
        if extent is not None: 
            self._bounds = bounds2geom(extent.region)
        else: self._status = -1

    def _parse_nos_xml(self, xml_url, sid):        
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
        if self._has_vector:
            gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
            layer = gmt1.GetLayer()
        else: layer = []

        print('{:79}'.format(nosdir))

        xml_catalog = self._nos_xml_url(nosdir)
        page = fetch_html(xml_catalog)
        rows = page.xpath('//a[contains(@href, ".xml")]/@href')

        for survey in rows:
            if not self.stop():
                sid = survey[:-4]

                if self._has_vector:
                    layer.SetAttributeFilter('ID = "{}"'.format(sid))

                if len(layer) == 0:
                    xml_url = xml_catalog + survey
                    s_entry = self._parse_nos_xml(xml_url, sid)

                    if s_entry[0]:
                        self._surveys.append(s_entry)
        gmt1 = layer = None

    def search_gmt(self, filters=[]):
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()

        for filt in filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            if not self.stop():
                geom = feature.GetGeometryRef()
                if geom.Intersects(self._bounds):
                    fldata = feature.GetField('Data').split(',')

                    for i in fldata:
                        self._results.append(i)

    def _update(self):
        for j in self._nos_directories:
            if not self.stop():
                self._scan_directory(j)
                update_ref_vector(self._ref_vector, self._surveys, self._has_vector)

                if os.path.exists(self._ref_vector): 
                    self._has_vector = True
                else: self._has_vector = False

                self._surveys = []

    def proc_data(self, s_dir, s_fn, o_fn, s_t):
        '''Process a fetched file to xyz (navd88) and add it to its datalist.'''

        status = 0
        xyz_dir = os.path.join(self._outdir, s_dir, s_t.lower(), 'xyz')
        
        if not os.path.exists(xyz_dir):
            os.makedirs(xyz_dir)

        o_fn_bn = os.path.basename(o_fn).split('.')[0]
        o_fn_tmp = os.path.join(self._outdir, s_dir, '{}_tmp.xyz'.format(o_fn_bn.lower()))
        o_fn_xyz = os.path.join(xyz_dir, '{}_{}.xyz'.format(o_fn_bn, self.region.fn))

        ## ==============================================
        ## NOS XYZ data comes as gzipped CSV
        ## ==============================================
        if s_t == 'GEODAS':

            ## ==============================================            
            ## gunzip the and parse the file
            ## ==============================================

            in_f = gzip.open(os.path.join(s_dir, s_fn), 'rb')
            s = in_f.read()
            s_csv = csv.reader(s.split('\n'), delimiter = ',')
            next(s_csv, None)

            with open(os.path.join(o_fn_tmp), 'w') as out_xyz:
                d_csv = csv.writer(out_xyz, delimiter = ' ')
                
                for row in s_csv:
                    if len(row) > 2:
                        this_xyz = [float(row[2]), float(row[1]), float(row[3]) * -1]
                        d_csv.writerow(this_xyz)
            in_f.close()

            if status == 0:

                ## ==============================================
                ## transform processed xyz file to NAVD88 
                ## using vdatum
                ## ==============================================

                if len(self.this_vd.vdatum_paths) > 0:
                    self.this_vd.ivert = 'mllw'
                    self.this_vd.overt = 'navd88'
                    self.this_vd.ds_dir = os.path.relpath(os.path.join(xyz_dir, 'result'))
                    
                    self.this_vd.run_vdatum(os.path.relpath(o_fn_tmp))
                    
                    os.rename(os.path.join(xyz_dir, 'result', os.path.basename(o_fn_tmp)), o_fn_xyz)
                else: os.rename(o_fn_tmp, o_fn_xyz)
            
                ## ==============================================
                ## Move processed xyz file to xyz directory
                ## ==============================================

                os.remove(o_fn_tmp)
                if os.stat(o_fn_xyz).st_size == 0:
                    os.remove(o_fn_xyz)
                else:

                    ## ==============================================
                    ## Add xyz file to datalist
                    ## ==============================================

                    sdatalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(s_t)))
                    sdatalist._append_datafile('{}'.format(os.path.basename(o_fn_xyz)), 168, 1)
                    sdatalist._reset()

                    ## ==============================================
                    ## Generate .inf file
                    ## ==============================================

                    out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, '{}.datalist'.format(s_t))), False, None)
            os.remove(o_fn)

        ## ==============================================
        ## NOS BAG data comes as gzipped BAG
        ## ==============================================

        elif 'BAG' in s_t:
            
            s_gz = os.path.join(s_dir, s_fn)

            if s_gz.split('.')[-1] == 'gz':
                utils.run_cmd('gunzip {}'.format(os.path.join(s_dir, s_fn)), False, 'gunzip')
                s_bag = '.'.join(s_gz.split('.')[:-1])
            else: s_bag = os.path.join(s_dir, s_fn)

            s_tif = os.path.join(s_dir, s_fn.split('.')[0].lower() + '.tif')
            
            s_xyz = s_gz.split('.')[0] + '.xyz'

            utils.run_cmd('gdalwarp {} {} -t_srs EPSG:4326'.format(s_bag, s_tif), True, 'gdalwarp')

            out_chunks = gdalfun.chunks(s_tif, 1000)
            os.remove(s_tif)

            for i in out_chunks:
                i_xyz = i.split('.')[0] + '.xyz'
                o_xyz = os.path.join(xyz_dir, os.path.basename(i_xyz))

                gdalfun.dump(i, i_xyz)

                if os.stat(i_xyz).st_size == 0:
                    os.remove(i_xyz)
                else:
                    ## ==============================================
                    ## transform processed xyz file to NAVD88 
                    ## using vdatum
                    ## ==============================================

                    if len(self.this_vd.vdatum_paths) > 0:
                        self.this_vd.ivert = 'mllw'
                        self.this_vd.overt = 'navd88'
                        self.this_vd.ds_dir = os.path.relpath(os.path.join(xyz_dir, 'result'))
                        
                        self.this_vd.run_vdatum(os.path.relpath(i_xyz))
                        
                        os.rename(os.path.join(xyz_dir, 'result', os.path.basename(i_xyz)), o_xyz)
                    else: os.rename(i_xyz, o_xyz)
            
                    ## ==============================================
                    ## Move processed xyz file to xyz directory
                    ## ==============================================
                    
                    os.remove(i_xyz)
                    os.remove(i)
                    if os.stat(o_xyz).st_size == 0:
                        os.remove(o_xyz)
                    else:

                        ## ==============================================
                        ## Add xyz file to datalist
                        ## ==============================================
                        
                        sdatalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(s_t)))
                        sdatalist._append_datafile('{}'.format(os.path.basename(o_xyz)), 168, 1)
                        sdatalist._reset()

                        ## ==============================================
                        ## Generate .inf file
                        ## ==============================================

                        out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, '{}.datalist'.format(s_t))), False, None)
                
    def print_results(self):
        for row in self._results:
            if row:
                print(row)

    def fetch_results(self):
        #try:

        if self._want_proc:
            self.this_vd = utils.vdatum()

        for row in self._results:
            if row:
                if not self.stop():

                    ## ==============================================
                    ## Fetch the data
                    ## ==============================================
                    
                    fetch_file(row, os.path.join(self._outdir, os.path.basename(row)))
                    
                    ## ==============================================                    
                    ## Process the data if wanted
                    ## ==============================================

                    if self._want_proc:
                        outf = os.path.join(self._outdir, os.path.basename(row))
                        surv_dir = self._outdir
                        surv_fn = os.path.basename(outf)
                        surv_t = row.split('/')[-2]
                        
                    #t = threading.Thread(target = self.proc_data, args = (surv_dir, surv_fn, outf, surv_t))
                    #t.start()
                    self.proc_data(surv_dir, surv_fn, outf, surv_t)

        #except: self._status = -1

## =============================================================================
##
## Chart Fetch - ENC & RNC
##
## Fetch digital charts from NOAA, including ENC and RNC
##
## =============================================================================

class charts:
    def __init__(self, extent = None, callback = None):
        '''Fetch digital chart data from NOAA'''

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

        self._want_proc = True

        self.region = extent
        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)
        else: self._status = -1

    def _parse_charts_xml(self, update = True):
        '''parse the charts xyz and extract the survey results'''

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

    def search_gmt(self, filters=[]):
        '''Search for data in the reference vector file'''

        ds = ogr.Open(self._ref_vector)
        layer = ds.GetLayer(0)

        ## ==============================================
        ## Filter the reference vector and append
        ## any matching results to _results
        ## ==============================================

        for filt in filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature1 in layer:
            geom = feature1.GetGeometryRef()
            if self._boundsGeom.Intersects(geom):
                self._results.append(feature1.GetField('Data'))

        ds = layer = None

    def _update(self):
        '''Update or create the reference vector file'''

        for dt in self._dt_xml.keys():
            self._checks = dt

            ## ==============================================
            ## parse the Charts XML
            ## ==============================================

            self.chart_xml = fetch_nos_xml(self._dt_xml[self._checks])
            self._parse_charts_xml(self._has_vector)

            ## ==============================================
            ## Update the reference vector
            ## ==============================================

            if len(self._chart_feats) > 0:
                update_ref_vector(self._ref_vector, self._chart_feats, self._has_vector)

            if os.path.exists(self._ref_vector): 
                self._has_vector = True
            else: self._has_vector = False

            self._chart_feats = []

    def proc_data(self, s_dir, s_fn, o_fn, s_t):
        '''Process a fetched file to xyz (navd88) and add it to its datalist.'''

        status = 0
        xyz_dir = os.path.join(self._outdir, s_dir, 'xyz')
        
        if not os.path.exists(xyz_dir):
            os.makedirs(xyz_dir)

        o_fn_bn = os.path.basename(o_fn).split('.')[0]
        o_fn_p_xyz = os.path.join(self._outdir, s_dir, '{}.xyz'.format(o_fn_bn))
        o_fn_tmp = os.path.join(self._outdir, s_dir, '{}_tmp.xyz'.format(o_fn_bn.lower()))
        o_fn_xyz = os.path.join(xyz_dir, '{}_{}.xyz'.format(o_fn_bn, self.region.fn))

        if s_t == 'ENCs':

            ## ==============================================
            ## extract downloaded ZIP
            ## ==============================================

            zip_ref = zipfile.ZipFile(o_fn)
            zip_ref.extractall(os.path.join(s_dir, 'enc'))
            zip_ref.close()

            s_fn_000 = os.path.join(s_dir, 'enc/ENC_ROOT/', o_fn_bn, '{}.000'.format(o_fn_bn))

            ## ==============================================
            ## open the s57 vector and write out the xyz data
            ## ==============================================

            ds_000 = ogr.Open(s_fn_000)
            layer_s = ds_000.GetLayerByName('SOUNDG')
            if layer_s is not None:
                o_xyz = open(o_fn_p_xyz, 'w')
                for f in layer_s:
                    g = json.loads(f.GetGeometryRef().ExportToJson())
                    for xyz in g['coordinates']:
                        xyz_l = '{} {} {}\n'.format(xyz[0], xyz[1], xyz[2]*-1)
                        o_xyz.write(xyz_l)
                o_xyz.close()
            else: status = -1

            if os.path.exists(o_fn_p_xyz):

                ## ==============================================
                ## Extract the data in the specified region
                ## ==============================================

                out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR=space', False, None)
                out, status = utils.run_cmd('gmt gmtselect {} {} -V > {}'.format(o_fn_p_xyz, self.region.gmt, o_fn_tmp), False, 'extracting data using gmt gmtselect')

                if os.stat(o_fn_tmp).st_size == 0:
                    status = -1
            else: status = -1
        else: status = -1

        if status == 0:

            ## ==============================================
            ## transform processed xyz file to NAVD88 
            ## using vdatum
            ## ==============================================

            if len(self.this_vd.vdatum_paths) > 0:
                self.this_vd.ivert = 'mllw'
                self.this_vd.overt = 'navd88'
                self.this_vd.ds_dir = os.path.relpath(os.path.join(xyz_dir, 'result'))

                self.this_vd.run_vdatum(os.path.relpath(o_fn_tmp))

                os.rename(os.path.join(xyz_dir, 'result', os.path.basename(o_fn_tmp)), o_fn_xyz)
            else: os.rename(o_fn_tmp, o_fn_xyz)
            
            ## ==============================================
            ## Move processed xyz file to xyz directory
            ## ==============================================

            os.remove(o_fn_p_xyz)

            if os.stat(o_fn_xyz).st_size == 0:
                os.remove(o_fn_xyz)
            else:

                ## ==============================================
                ## Add xyz file to datalist
                ## ==============================================

                sdatalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(s_t)))
                sdatalist._append_datafile('{}'.format(os.path.basename(o_fn_xyz)), 168, 1)
                sdatalist._reset()

                ## ==============================================
                ## Generate .inf file
                ## ==============================================

                out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, '{}.datalist'.format(s_t))), False, None)

        if os.path.exists(o_fn_tmp):
            os.remove(o_fn_tmp)
        os.remove(o_fn)

    def print_results(self):
        '''print the resulting urls to stdout'''

        for row in self._results:
            print(row)

    def fetch_results(self):
        '''fetch the charts data and optionally process it.
        unprocessed data is in S57 (ENC) and KAPP (RNC) format.
        will only process ENC data to XYZ'''

        if self._want_proc:
            self.this_vd = utils.vdatum()

        try:
            for row in self._results:

                ## ==============================================
                ## Fetch the data
                ## ==============================================

                fetch_file(row, os.path.join(self._outdir, os.path.basename(row)))

                ## ==============================================                    
                ## Process the data if wanted
                ## ==============================================

                if self._want_proc:
                    outf = os.path.join(self._outdir, os.path.basename(row))
                    surv_dir = self._outdir
                    surv_fn = os.path.basename(outf)
                    surv_t = row.split('/')[-2]
                    
                    t = threading.Thread(target = self.proc_data, args = (surv_dir, surv_fn, outf, surv_t))
                    t.start()
        except: self._status = -1

## =============================================================================
##
## The National Map (TNM) - USGS
##
## Fetch elevation data from The National Map
## NED, 3DEP, Etc.
##
## =============================================================================

class tnm:
    def __init__(self, extent = None, callback = None):
        self._tnm_api_url = "http://viewer.nationalmap.gov/tnmaccess/"
        self._tnm_dataset_url = "http://viewer.nationalmap.gov/tnmaccess/api/datasets?"
        self._tnm_product_url = "http://viewer.nationalmap.gov/tnmaccess/api/products?"

        self._status = 0
        self._results = None

        self._outdir = os.path.join(os.getcwd(), 'tnm')
        self._dataset_results = None
        self._ref_vector = None

        self._tnm_ds = [1, 2]
        self._tnm_df = ['IMG']

        if extent is not None: 
            self._results = fetch_req(self._tnm_dataset_url)
            self._datasets = self._results.json()

            sbDTags = []
            for ds in self._tnm_ds:
                sbDTags.append(self._datasets[ds]['sbDatasetTag'])
        
            self.data = { 'datasets':sbDTags,
                          'bbox':extent.bbox }
        
            self._results = fetch_req(self._tnm_product_url, params = self.data)
            self._dataset_results = self._results.json()

    def print_results(self):
        for i in self._dataset_results['items']:
            print i['downloadURL']

    def fetch_results(self):
        try:
            for i in self._dataset_results['items']:
                fetch_file(i['downloadURL'], os.path.join(self._outdir, os.path.basename(i['downloadURL'])))
        except: self._status = -1

## =============================================================================
##
## MB Fetch
##
## Fetch Multibeam bathymetric surveys from NOAA
## MBSystem is required to process the resulting data
##
## =============================================================================

class mb:
    def __init__(self, extent = None, callback = None):
        self._mb_data_url = "https://data.ngdc.noaa.gov/platforms/"
        self._mb_search_url = "https://maps.ngdc.noaa.gov/mapviewer-support/multibeam/files.groovy?"
        self._outdir = os.path.join(os.getcwd(), 'mb')

        self._status = 0
        self._results = None
        self._surveys = []
        self._survey_list = []
        self._ref_vector = None

        if extent is not None:
            self.data = { 'geometry':extent.bbox }
            self._results = fetch_req(self._mb_search_url, params = self.data)
            self._survey_list = self._results.content.split('\n')[:-1]
        
    def parse_results(self, local):
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

    def print_results(self):
        for res in self._survey_list:
            print self._mb_data_url + res.split(' ')[0]

    def fetch_results(self):
        try:
            for r in self._survey_list:
                survey = r.split(' ')[0].split('/')[6]
                dn = r.split(' ')[0].split('/')[:-1]
                dst_dn = os.path.join(self._outdir, *dn)

                if not os.path.exists(dst_dn):
                    os.makedirs(dst_dn)

                data_url = self._mb_data_url + '/'.join(r.split('/')[3:])
                dst_fn = r.split(' ')[0].split('/')[-1:][0]
                
                fetch_file(data_url.split(' ')[0], os.path.join(dst_dn, dst_fn))
        except: self._status = -1

## =============================================================================
##
## USACE Fetch
##
## Fetch USACE bathymetric surveys
##
## =============================================================================

class usace:
    def __init__(self, extent = None, callback = None):

        self._usace_gj_api_url = 'https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson'
        self._usace_gs_api_url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?outFields=*&where=1%3D1'

        self._status = 0
        self._outdir = os.path.join(os.getcwd(), 'usace')
        self._results = []
        self._survey_list = []
        self._ref_vector = None

        if extent is not None:
            self.data = { 'geometry':extent.bbox,
                          'inSR':4326,
                          'f':'pjson' }
 
            self._results = fetch_req(self._usace_gs_api_url, params = self.data)
            self._survey_list = self._results.json()

    def print_results(self):
        for feature in self._survey_list['features']:
            print feature['attributes']['SOURCEDATALOCATION']

    def fetch_results(self):
        try:
            for feature in self._survey_list['features']:
                fetch_file(feature['attributes']['SOURCEDATALOCATION'], \
                           os.path.join(self._outdir, os.path.basename(feature['attributes']['SOURCEDATALOCATION'])))
        except: self._status = -1

## =============================================================================
##
## GMRT Fetch
##
## fetch extracts of the GMRT.
##
## =============================================================================

class gmrt:
    def __init__(self, extent = None, callback = None):
        self._gmrt_grid_url = "https://www.gmrt.org/services/GridServer?"
        self._outdir = os.path.join(os.getcwd(), 'gmrt')
        self.outf = None

        self._status = 0
        self._results = None
        self._ref_vector = None

        self._want_proc = True

        self.region = extent
        if extent is not None: 
            self.data = { 'north':extent.north,
                          'west':extent.west,
                          'south':extent.south,
                          'east':extent.east,
                          'layer':'topo',
                          'resolution':'max',
                          'format':'geotiff' }
        
            self._results = fetch_req(self._gmrt_grid_url, params = self.data)

    def proc_data(self, s_dir, s_fn, o_fn, s_t):

        status = 0
        xyz_dir = os.path.join(self._outdir, s_dir, 'xyz')
        
        if not os.path.exists(xyz_dir):
            os.makedirs(xyz_dir)

        o_fn_bn = os.path.basename(o_fn).split('.')[0]
        o_fn_xyz = os.path.join(xyz_dir, '{}_{}.xyz'.format(o_fn_bn, self.region.fn))

        ## chunk
        ## convert to xyz
        gdalfun.dump(o_fn, o_fn_xyz)

        ## ==============================================
        ## Add xyz file to datalist
        ## ==============================================

        sdatalist = datalists.datalist(os.path.join(xyz_dir, '{}.datalist'.format(s_t)))
        sdatalist._append_datafile('{}'.format(os.path.basename(o_fn_xyz)), 168, 1)
        sdatalist._reset()

        ## ==============================================
        ## Generate .inf file
        ## ==============================================

        out, status = utils.run_cmd('mbdatalist -O -I{}'.format(os.path.join(xyz_dir, '{}.datalist'.format(s_t))), False, None)

    def print_results(self):
        print(self._results.url)

    def fetch_results(self):
        #try:
        outf = os.path.join(self._outdir, self._results.headers['content-disposition'].split('=')[1].strip())
        
        if not os.path.exists(os.path.dirname(outf)):
            os.makedirs(os.path.dirname(outf))

        with open(outf, 'wb') as local_file:
            for chunk in self._results.iter_content(chunk_size = 50000):
                local_file.write(chunk)

        if self._want_proc:
            surv_dir = self._outdir
            surv_fn = os.path.basename(outf)

            t = threading.Thread(target = self.proc_data, args = (surv_dir, surv_fn, outf, 'gmrt'))
            t.start()

        #except: self._status = -1
            
## =============================================================================
##
## SRTM Fetch (cgiar)
##
## Fetch srtm tiles from cgiar.
##
## =============================================================================

class srtm_cgiar:
    def __init__(self, extent = None, callback = None):
        self._srtm_url = 'http://srtm.csi.cgiar.org'
        self._srtm_dl_url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/'

        self._status = 0
        self._results = []
        self._outdir = os.path.join(os.getcwd(), 'srtm')
        self._ref_vector = os.path.join(fetchdata, 'srtm.gmt')

        self._boundsGeom = None
        if extent is not None: 
            self._boundsGeom = bounds2geom(extent.region)

    def search_gmt(self, filters=[]):
        gmt1 = ogr.GetDriverByName('GMT').Open(self._ref_vector, 0)
        layer = gmt1.GetLayer()

        for filt in filters:
            layer.SetAttributeFilter('{}'.format(filt))

        for feature in layer:
            geom = feature.GetGeometryRef()

            if geom.Intersects(self._boundsGeom):
                geo_env = geom.GetEnvelope()
                srtm_lon = int(math.ceil(abs((-180 - geo_env[1]) / 5)))
                srtm_lat = int(math.ceil(abs((60 - geo_env[3]) / 5)))
                self._results.append('srtm_{:2}_{:2}.zip'.format(srtm_lon, srtm_lat))

    def print_results(self):
        for row in self._results:
            print('{}{}'.format(self._srtm_dl_url, row))

    def fetch_results(self):
        try:
            for row in self._results:
                fetch_file('{}{}'.format(self._srtm_dl_url, row), os.path.join(self._outdir, os.path.basename(row)))
        except: self._status = -1

## =============================================================================
##
## National Geodetic Survey (NGS)
##
## Fetch NGS monuments from NGS
##
## =============================================================================

class ngs:
    def __init__(self, extent = None, callback = None):
        self._ngs_search_url = 'http://geodesy.noaa.gov/api/nde/bounds?'
        self._outdir = os.path.join(os.getcwd(), 'ngs')
        self._ref_vector = None

        self._status = 0
        self._results = None

        if extent is not None:
            self.data = { 'maxlon':extent.east,
                          'minlon':extent.west,
                          'maxlat':extent.north,
                          'minlat':extent.south }

            self._results = fetch_req(self._ngs_search_url, params = self.data)
        
    def print_results(self):
        print self._results.url

    def fetch_results(self):
        try:
            r = self._results.json()

            if len(r) > 0:
                if not os.path.exists(self._outdir):
                    os.makedirs(self._outdir)

                dt_now_str = datetime.datetime.now().strftime('%d%Y')
                outfile = open(os.path.join(self._outdir, 'ngs_results_{}.csv'.format(dt_now_str)), 'w')
            
                outcsv = csv.writer(outfile)
                outcsv.writerow(r[0].keys())
                for row in r:
                    outcsv.writerow(row.values())

                outfile.close()
        except: self._status = -1

### End
