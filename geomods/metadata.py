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

import waffles
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
    
    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None, o_extend = 6, callback = lambda: False, verbose = False):

        self.dl_q = queue.Queue()
        self.datalist = i_datalist
        self.inc = i_inc
        self.region = i_region
        self.extend = int(o_extend)
        self.dist_region = self.region.buffer(self.extend * self.inc)

        self.stop = callback
        self.verbose = verbose
        self.has_gmt = True
        self.want_queue = True
        self.use_bounds = True

        self.v_fields = ['Name', 'Agency', 'Date', 'Type',
                         'Resolution', 'HDatum', 'VDatum', 'URL']
        if self.use_bounds:
            self.t_fields = ['string', 'string', 'string', 'string',
                             'string', 'string', 'string', 'string']
        else:
            self.t_fields = [ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString,
                             ogr.OFTString, ogr.OFTString, ogr.OFTString, ogr.OFTString]

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name

        self.gmt_vers = utils.config_get_vers('GMT')

        self.datalist._load_datalists()
        
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

        #this_datalist = datalists.datalist(dl[0], self.dist_region, verbose = self.verbose)
        this_datalist = datalists.datalist(dl[0], self.region, verbose = self.verbose)
        this_o_name = this_datalist._name        
        this_datalist._load_data()

        if len(this_datalist.datafiles) > 0:
            try:
                #o_v_fields = [dl[3], dl[4], dl[5], dl[6], dl[7], dl[8], dl[9], dl[10].strip()]
                o_v_fields = dl[3]
                if len(o_v_fields) != 8:
                    o_v_fields = [this_datalist._name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
            except: o_v_fields = [this_datalist._name, 'Unknown', '0', 'xyz_elevation', 'Unknown', 'WGS84', 'NAVD88', 'URL']
            pb = utils._progress('gathering geometries from datalist \033[1m{}\033[m...'.format(this_o_name))

            if self.use_bounds:
                o_v_fields = ['\\"{}\\"'.format(x) if ' ' in x else x for x in o_v_fields]
                utils.run_cmd('bounds -k {}/{} -n "{}" -gg --verbose >> {}.gmt\
                '.format(self.inc, self.dist_region.region_string, '|'.join(o_v_fields), layer), self.verbose, self.verbose, this_datalist._dump_data)
            else:
                defn = layer.GetLayerDefn()
                if self.gmt_vers is None:
                    this_mask = this_datalist.mask(self.inc)
                else:
                    this_dem = waffles.dem(this_datalist, this_datalist.region, i_inc = self.inc, o_fmt = 'GTiff', o_extend = self.extend, verbose = self.verbose)
                    #this_dem.o_fmt = 'GTiff'
                    this_mask = this_dem.run('mask')

                if os.path.exists(this_mask) and not self.stop():
                    ## should use more unique name...crashes when 2 datalists have same name at same time...
                    utils.remove_glob('{}_poly.*'.format(this_o_name))

                    tmp_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('{}_poly.shp'.format(this_o_name))
                    tmp_layer = tmp_ds.CreateLayer('{}_poly'.format(this_o_name), None, ogr.wkbMultiPolygon)
                    tmp_layer.CreateField(ogr.FieldDefn('DN', ogr.OFTInteger))

                    gdalfun.polygonize(this_mask, tmp_layer, verbose = self.verbose)

                    if len(tmp_layer) > 1:
                        out_feat = gdalfun.ogr_mask_union(tmp_layer, 'DN', defn, self.stop, verbose = self.verbose)
                        for i, f in enumerate(self.v_fields):
                            out_feat.SetField(f, o_v_fields[i])

                        layer.CreateFeature(out_feat)

                    tmp_ds = tmp_layer = out_feat = None
                    utils.remove_glob('{}_poly.*'.format(this_o_name))
                    utils.remove_glob('{}*'.format(this_mask[:-3]))

            pb.end(0, 'gathered geometries from datalist \033[1m{}\033[m.'.format(this_o_name))

    def run(self, epsg = 4269):
        '''Run the spatial-metadata module and Geneate spatial metadata from the datalist
        specify the output project epsg.'''

        if self.inc < 0.0000925:
            utils._msg('warning, increments less than 1/3 arc-second may be slow.')
    
        dst_vec = '{}_sm.shp'.format(self.o_name)
        dst_layername = '{}_sm'.format(self.o_name)
        utils.remove_glob('{}.*'.format(dst_layername))

        if self.use_bounds:
            ##utils.run_cmd('bounds -ggg > {}.gmt'.format(dst_layername, self.verbose, self.verbose))
            with open('{}.gmt'.format(dst_layername), 'w') as gmtf:
                gmtf.write('# @VGMT1.0 @GMULTIPOLYGON\n# @N{}\n# @T{}\n# FEATURE_DATA\n'.format('|'.join(self.v_fields), '|'.join(self.t_fields)))
            for dl in self.datalist.datalists:
                self._gather_from_datalist(dl, dst_layername)
        else:
            ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(dst_vec)
            if ds is not None:
                gdalfun._prj_file('{}.prj'.format(dst_layername), epsg)
                layer = ds.CreateLayer('{}'.format(dst_layername), None, ogr.wkbMultiPolygon)

                for i, f in enumerate(self.v_fields):
                    layer.CreateField(ogr.FieldDefn('{}'.format(f), self.t_fields[i]))

                for feature in layer:
                    layer.SetFeature(feature)

                if self.want_queue:
                    for _ in range(3):
                           t = threading.Thread(target = self._gather_from_queue, args = ())
                           t.daemon = True
                           t.start()

                    if len(self.datalist.datalists) > 0:
                        for dl in self.datalist.datalists:
                            self.dl_q.put([dl, layer])
                    else:
                        self.dl_q.put([[self.datalist._path, -1, 1], layer])

                    self.dl_q.join()
                else:
                    for dl in self.datalist.datalists:
                        self._gather_from_datalist(dl, layer)
                
        ds = layer = None
        utils.run_cmd('ogr2ogr {} {}.gmt'.format(dst_vec, dst_layername))
        if not os.path.exists(dst_vec):
            dst_vec = None

        return(dst_vec)
### End
