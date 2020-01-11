### regions.py
##
## Copyright (c) 2012 - 2019 Matthew Love <matthew.love@colorado.edu>
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

_version = '0.0.2'

## =============================================================================
##
## Region Class
##
## a region is a geographic bounding-box with 4 corners and the
## region object made from the region_string 'east/west/south/north'
##
## =============================================================================
class region:
    def __init__(self, region_string):
        self.region_string = region_string
        self.region = map( float, region_string.split('/'))
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

    ## 'Validate' Region
    def _valid_p(self):
        if self.west < self.east and self.south < self.north: return(True)
        else: return(False)

    ## Format region
    def _format_gmt(self):
        self.gmt = '-R' + '/'.join(map(str, self.region))

    def _format_bbox(self):
        self.bbox = ','.join([str(self.west), str(self.south), str(self.east), str(self.north)])

    def _format_fn(self):
        if self.north < 0: ns = 's'
        else: ns = 'n'
        if self.east < 0: ew = 'w'
        else: ew = 'e'
        self.fn = ('{}{:4}x{:2}_{}{:3}x{:2}\
        '.format(ns, abs(int(self.north)), abs(int(self.north * 100) % 100), 
                 ew, abs(int(self.east)), abs(int(self.east * 100) % 100)))

    ## Process Region
    def buffer(self, bv, percentage = False):
        if percentage: bv = self.pct(bv)

        region_b = [self.region[0]-bv, self.region[1] + bv, self.region[2] - bv, self.region[3] + bv]

        return(region("/".join(map(str, region_b))))

    def pct(self, pctv):
        ewp = (self.east - self.west) * (pctv * .01)
        nsp = (self.north - self.south) * (pctv * .01)

        return((ewp + nsp) / 2)

    ## Split Region
    def split(self, sv):
        
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

    def regions_intersect_p(self, region_a, region_b):
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

### End
