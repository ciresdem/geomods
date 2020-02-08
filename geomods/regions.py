### regions.py
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

import gdalfun

_version = '0.1.0'

## =============================================================================
##
## Region Class
##
## a region is a geographic bounding-box with 4 corners and the
## region object made from the region_string 'east/west/south/north'
##
## =============================================================================

def regions_intersect_p_depr(region_a, region_b):
    '''Return True if region_a and region_b intersect.'''
    
    region_c = [0, 0, 0, 0]
    if region_a.west > region_b.west: 
        region_c[0] = region_a.west
    else: region_c[0] = region_b.west
    
    if region_a.west < region_b.east: 
        region_c[1] = region_a.east
    else: region_c[1] = region_b.east
    
    if region_a.south > region_b.south: 
        region_c[2] = region_a.south
    else: region_c[2] = region_b.south
    
    if region_a.west < region_b.south: 
        region_c[3] = region_a.north
    else: region_c[3] = region_b.north
    
    return(region('/'.join(map(str, region_c)))._valid)

def regions_intersect_p(region_a, region_b):
    geom_a = gdalfun.bounds2geom(region_a.region)
    geom_b = gdalfun.bounds2geom(region_b.region)

    if geom_a.Intersects(geom_b):
        return(True)
    else: return(False)

class region:
    '''geographic bounding box regtions 'w/e/s/n' '''

    def __init__(self, region_string):
        self.region_string = region_string
        self.region = map(float, region_string.split('/'))
        self._reset()

    def _reset(self):
        self.west = self.region[0]
        self.east = self.region[1]
        self.south = self.region[2]
        self.north = self.region[3]
        self._format_gmt()
        self._format_bbox()
        self._format_fn()
        self._valid = self._valid_p()        

    def _valid_p(self):
        '''validate region'''

        if self.west < self.east and self.south < self.north: return(True)
        else: return(False)

    def _format_gmt(self):
        '''format region to GMT string'''

        self.gmt = '-R' + '/'.join(map(str, self.region))

    def _format_bbox(self):
        '''format region to bbox string'''

        self.bbox = ','.join([str(self.west), str(self.south), str(self.east), str(self.north)])

    def _format_fn(self):
        '''format region to filename string'''

        if self.north < 0: ns = 's'
        else: ns = 'n'
        if self.east < 0: ew = 'w'
        else: ew = 'e'
        self.fn = ('{}{:02d}x{:02d}_{}{:03d}x{:02d}'.format(ns, abs(int(self.north)), abs(int(self.north * 100) % 100), 
                                                            ew, abs(int(self.east)), abs(int(self.east * 100) % 100)))

    def gdal2region(self):
        '''extract the region from a GDAL file.'''

        pass

    def buffer(self, bv, percentage = False):
        '''buffer region'''

        if percentage: bv = self.pct(bv)
        region_b = [self.region[0]-bv, self.region[1] + bv, self.region[2] - bv, self.region[3] + bv]

        return(region("/".join(map(str, region_b))))

    def chunk(self, inc, n_chunk = 10):
        '''chunk the region into n_chunk by n_chunk cell regions, given inc.'''
        import math

        i_chunk = 0
        x_i_chunk = 0
        x_chunk = n_chunk
        o_chunks = []
        region_x_size = math.floor((self.east - self.west) / inc)
        region_y_size = math.floor((self.north - self.south) / inc)

        while True:
            y_chunk = n_chunk

            while True:
                this_x_chunk = x_chunk
                this_y_chunk = y_chunk
                this_x_origin = x_chunk - n_chunk
                this_y_origin = y_chunk - n_chunk
                this_x_size = this_x_chunk - this_x_origin
                this_y_size = this_y_chunk - this_y_origin
                geo_x_o = self.west + this_x_origin * inc
                geo_x_t = geo_x_o + this_x_size * inc
                geo_y_o = self.south + this_y_origin * inc
                geo_y_t = geo_y_o + this_y_size * inc

                this_region = region('{}/{}/{}/{}'.format(geo_x_o, geo_x_t, geo_y_o, geo_y_t))
                o_chunks.append(this_region)

                if y_chunk >= region_y_size:
                    break
                else: 
                    y_chunk += n_chunk
                    i_chunk += 1

            if x_chunk >= region_x_size:
                break
            else:
                x_chunk += n_chunk
                x_i_chunk += 1

        return(o_chunks)

    def pct(self, pctv):
        ewp = (self.east - self.west) * (pctv * .01)
        nsp = (self.north - self.south) * (pctv * .01)

        return((ewp + nsp) / 2)

    def split(self, sv):
        '''split region by split-value `sv`'''
        
        split_regions = []
        
        ew = (self.east - self.west) / sv
        ns = (self.north - self.south) / sv

        ewo = self.west
        ewn = self.west+ew
        nso = self.south
        nsn = self.south+ew

        for i in range(1, sv + 1):

            region_sub = [ewo, ewn, nso, nsn]
            split_regions.append(region('/'.join(map(str, region_sub))))

            ewo = ewn
            ewn += ew
            nso = nsn
            nsn += ns

        return(split_regions)

### End
