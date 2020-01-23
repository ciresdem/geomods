### datalists.py
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
import regions

_version = '0.0.2'

## =============================================================================
##
## Datalist Class
##
## MBSystem style datalists.
## Recurse through a datalist file and process the results.
##
## a datalist '*.datalits' file should be formatted as in MBSystem:
## ~path ~format ~weight
##
## a format of -1 represents a datalist
## a format of 168 represents XYZ data
##
## each xyz file in a datalist should have an associated '*.inf' file 
## for faster processing
##
## 'inf' files can be generated using 'mbdatalist -O -V -I~datalist.datalist'
##
## if 'iregion' is specified, will only process data that falls within
## the given region
##
## =============================================================================
class datalist:
    '''MBSystem style datalists for elevation data.'''

    def __init__(self, idatalist, iregion = None):
        if not os.path.exists(idatalist): 
            open(idatalist, 'a').close()

        self.region = iregion
        self._path = idatalist
        self._path_dirname = os.path.dirname(self._path)
        self._path_basename = os.path.basename(self._path)
        self._path_dl_name = os.path.basename(self._path).split('.')[0]
        self._reset()

    def _reset(self):
        '''reload the datalist'''

        self.datalist = []
        self.datafiles = []
        self._load()
        self._proc(self.datafiles)
        self._valid = self._valid_p()

    def _valid_p(self):
        '''validate the datalist'''

        if len(self.datafiles) > 0: 
            return(True)
        else: return(False)

    def _read_inf(self, path_i):
        '''Read MBSystem .inf file and extract minmax info'''

        minmax = [0, 0, 0, 0]
        if os.path.exists(path_i):
            self.iob = open(path_i)
            for il in self.iob:
                til = il.split()
                if len( til ) > 1:
                    if til[0] == 'Minimum':
                        if til[1] == 'Longitude:':
                            minmax[0] = til[2]
                            minmax[1] = til[5]
                        elif til[1] == 'Latitude:':
                            minmax[2] = til[2]
                            minmax[3] = til[5]

        try: oregion = regions.region('/'.join(minmax))
        except: oregion = False
        return(oregion)
                     
    def _load(self):
        '''load and process a datalist'''

        with open(self._path, 'r') as fob:
            for dl in fob:
                if dl[0] != '#' and dl[0] != '\n' and dl[0] != '':
                    if len(dl.split(' ')) >= 2:
                        self.datalist.append(dl.split(' '))

    def _proc(self, datafiles = []):
        '''Recurse through the datalist and gather xyz data'''

        for i in self.datalist:
            dpath = i[0]
            dformat = int(i[1])

            if len(i) > 2: 
                dweight = i[2]
            else: dweight = 1

            if dformat == -1:
                d = datalist(dpath, self.region)
                d._proc(datafiles)
            elif dformat == 168:
                usable = True
                path_d = os.path.join(self._path_dirname, dpath)
                if not os.path.exists(path_d): usable = False
                dinf = self._read_inf(path_d + '.inf')
                
                if self.region and dinf:
                    usable = regions.regions_intersect_p(self.region, dinf)

                if usable:
                    datafiles.append(path_d)

    def _append_datafile(self, dfile, dformat, weight):
        '''append xyz data to a datalist'''

        with open(self._path, 'a') as outfile:
            outfile.write('{} {} {}\n'.format(dfile, dformat, weight))

    def _echo_datafiles(self, osep='/n'):
        '''return a string of datafiles in the datalist'''

        return(osep.join(self.datafiles))

    def _cat_port(self, dst_port):
        '''Catenate the xyz data from the datalist'''

        for fn in self.datafiles:
            with open(fn) as infile:
                for line in infile:
                    dst_port.write(line)

    def _caty(self):
        '''Catenate the xyz data from the datalist to a generator
        usage: for line in self._caty(): proc(line)'''

        for fn in self.datafiles:
            with open(fn) as infile:
                for line in infile:
                    yield(line)

    def _join(self, osep='\n'):
        df = []
        for fn in self.datafiles:
            with open(fn) as infile:
                for line in infile:
                    df.append(line)

        return osep.join(df)

### End
