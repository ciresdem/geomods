### datalists.py
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

import os
import regions

_version = '0.0.1'

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
    def __init__(self, idatalist, iregion = None):
        if not os.path.exists(idatalist): 
            open(idatalist, 'a').close()

        self.region = iregion
        self._path = idatalist
        self._path_dirname = os.path.dirname( self._path )
        self._reset()

    ## Reload the datalist
    def _reset(self):
        self.datalist = []
        self.datafiles = []
        self._load()
        self._proc( self.datafiles )
        self._valid = self._valid_p()

    ## 'Validate' the datalist (does it have datafiles?)
    def _valid_p(self):
        if len( self.datafiles ) > 0: return( True )
        else: return( False )

    ## Read MBSystem .inf file and extract minmax info
    def _read_inf(self, path_i):
        minmax = [0, 0, 0, 0]
        if os.path.exists( path_i ):
            self.iob = open( path_i )
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

        try: oregion = region( '/'.join( minmax ))
        except: oregion = False
        return oregion
                     
    ## Load and Process a datalist file.
    def _load(self):
        fob = open( self._path, 'r' )
        for dl in fob:
            if dl[0] != '#' and dl[0] != '':
                self.datalist.append( dl.split( ' ' ))
        fob.close()

    ## Recurse through the datalist and gather xyz data
    def _proc(self, datafiles = []):

        for i in self.datalist:
            dpath = i[0]
            dformat = int( i[1] )

            if len(i) > 2: 
                dweight = i[2]
            else: dweight = 1

            if dformat == -1:
                d = datalist( dpath, self.region )
                d._proc( datafiles )
            elif dformat == 168:
                usable = True
                path_d = os.path.join( self._path_dirname, dpath )
                if not os.path.exists( path_d ): usable = False
                dinf = self._read_inf( path_d + '.inf' )

                if self.region and dinf:
                    usable = regions.regions_intersect_p( self.region, dinf )

                if usable:
                    datafiles.append( path_d )

    ## Add XYZ data to datalist
    def _append_datafile(self, dfile, dformat, weight):
        with open( self._path, 'a' ) as outfile:
            outfile.write( '%s %s %s\n' %( dfile, dformat, weight ))

    ## Process XYZ data
    def _echo_datafiles( self, osep='/n' ):
        return osep.join( self.datafiles )      

    ## Catenate the xyz data from the datalist
    def _cat(self, dst_xyz):
        with open( dst_xyz, 'w' ) as outfile:
            for fn in self.datafiles:
                with open(fn) as infile:
                    for line in infile:
                        outfile.write( line )

### End
