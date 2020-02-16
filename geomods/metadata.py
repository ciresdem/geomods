### metadata.py
##
## Copyright (c) 2019 - 2020 CIRES Coastal DEM Team
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
import ogr

import datalists
import gdalfun
import utils

import Queue as queue
import threading
import time
    
## =============================================================================
##
## spatial-metadata module:
##
## Generate spatial metadata from a datalist of xyz elevation data.
## The output is a shapefile with the boundaries of each specific
## datalist as a layer feature.
##
## Specify an input datalist, region and cell-size.
## o_name is the output prefix: `{o_name}_sm.shp`
##
## Uses geomods modules: datalists, gdalfun, utils
##
## =============================================================================

class spatial_metadata:
    '''Generate spatial metadata from a datalist of xyz elevation data.
    The output is a shapefile with the unioned boundaries of each specific
    datalist as a layer feature.'''
    
    def __init__(self, i_datalist, i_region, i_inc = 0.000277777, o_name = None, callback = lambda: False, verbose = False):

        self.dl_q = queue.Queue()
        self.datalist = i_datalist
        self.region = i_region
        self.inc = i_inc
        self.dist_region = self.region.buffer(6 * self.inc)

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.v_fields = ['Name', 'Agency', 'Date', 'Type',
                         'Resolution', 'HDatum', 'VDatum', 'URL']
        self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                         ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name
        
    def run(self, epsg = 4269):
        '''Run the spatial-metadata module.
        specify the output project epsg.'''

        self.spatial_metadata(epsg)
        
    def _gather_from_queue(self):
        '''Gather geometries from a queue of [[datalist, layer], ...].'''

        while True:
            sm_args = self.dl_q.get()
            dl = sm_args[0]
            layer = sm_args[1]
            
            if not self.stop():
                self._gather_from_datalist(dl, layer)
            
            self.dl_q.task_done()

    def _gather_from_datalist(self, dl, layer):
        '''gather geometries from datalist `dl` and append
        results to ogr `layer`. Load the datalist, generate
        a NUM-MSK grid, polygonize said NUM-MSK then union
        the polygon and add it to the output layer.'''

        defn = layer.GetLayerDefn()
        this_datalist = datalists.datalist(dl[0], self.dist_region)
        this_o_name = this_datalist._name
        this_datalist._load_data()

        try:
            o_v_fields = [dl[3], dl[4], dl[5], dl[6], dl[7], dl[8], dl[9], dl[10].strip()]
        except: o_v_fields = [this_datalist._path_dl_name, 'Unknown', 0, 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']

        this_mask = this_datalist.mask(inc = i_inc)
        if os.path.exists(this_mask) and not self.stop():
            pb = utils._progress('gathering geometries from \033[1m{}\033[m...'.format(this_datalist._path_basename))
            utils.remove_glob('{}_poly.*'.format(this_o_name))

            tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_o_name))
            tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_o_name), None, ogr.wkbMultiPolygon)
            tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

            pb1 = utils._progress('polygonizing \033[1m{}\033[m...'.format(this_datalist._path_basename))
            gdalfun.polygonize(this_mask, tmp_layer, verbose = True)
            pb1.opm = 'polygonized \033[1m{}\033[m.'.format(this_datalist._path_basename)
            pb1.end(status)
                        
            if len(tmp_layer) > 1:
                out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn)
                for i, f in enumerate(self.v_fields):
                    out_feat.SetField(f, o_v_fields[i])

                layer.CreateFeature(out_feat)

            tmp_ds = tmp_layer = out_feat = None
            utils.remove_glob('{}_poly.*'.format(this_o_name))
            os.remove(this_mask)

            pb.opm = 'gathered geometries from \033[1m{}\033[m.'.format(this_datalist._path_basename)
            pb.end(self.status)

    def spatial_metadata(self, epsg = 4269):
        '''Geneate spatial metadata from the datalist'''

        if self.inc < 0.0000925:
            print 'warning, increments less than 1/3 arc-second may be slow.'
    
        dst_vec = '{}_sm.shp'.format(self.o_name, self.region.fn)
        dst_layername = os.path.basename(dst_vec).split('.')[0]

        utils.remove_glob('{}_sm.*'.format(self.o_name))

        ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
        if ds is not None:
            try:
                int(epsg)
            except: epsg = 4326

            gdalfun._prj_file('{}.prj'.format(os.path.basename(dst_vec).split('.')[0]), int(epsg))
            layer = ds.CreateLayer('{}'.format(os.path.basename(dst_vec).split('.')[0]), None, ogr.wkbMultiPolygon)

            for i, f in enumerate(self.v_fields):
                layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i]))
            
            for feature in layer:
                layer.SetFeature(feature)
    
            for _ in range(3):
                t = threading.Thread(target = self._gather_from_queue, args = ())
                t.daemon = True
                t.start()
                
            if len(self.datalist.datalists) > 0:
                for dl in self.datalist.datalists:
                    #self._gather_from_datalist(dl, layer)
                    self.dl_q.put([dl, layer])
            else:
                self.dl_q.put([[self.datalist._path, -1, 1], layer])

            self.dl_q.join()

        ds = layer = None
        if not os.path.exists(dst_vec):
            dst_vec = None

        return(dst_vec)
### End
