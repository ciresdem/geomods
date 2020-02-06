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

import sys
import os
import regions
import utils

_version = '0.1.0'

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
        self.datalist = []
        self.datafiles = []
        self._load()
        self._valid = self._valid_p()

    def _reset(self):
        '''reload the datalist'''

        self.datalist = []
        self.datafiles = []
        self._load()
        self._proc(self.datafiles)
        self._valid = self._valid_p()

    def _valid_p(self):
        '''validate the datalist'''

        if len(self.datalist) > 0: 
            return(True)
        else: return(False)

    def _read_inf(self, path_i):
        '''Read .inf file and extract minmax info.
        the .inf file can either be an MBSystem style inf file
        or the result of `gmt gmtinfo file.xyz -C`, which is
        a 6 column line with minmax info, etc.'''

        minmax = [0, 0, 0, 0]
        if os.path.exists(path_i):
            iob = open(path_i)
            for il in iob:
                til = il.split()
                if len(til) > 1:
                    ## GMT inf
                    try:
                        minmax[0] = float(til[0])
                        minmax[1] = float(til[1])
                        minmax[2] = float(til[2])
                        minmax[3] = float(til[3])
                    ## mbsystem inf
                    except:
                        if til[0] == 'Minimum':
                            if til[1] == 'Longitude:':
                                minmax[0] = til[2]
                                minmax[1] = til[5]
                            elif til[1] == 'Latitude:':
                                minmax[2] = til[2]
                                minmax[3] = til[5]

        try: 
            o_region = regions.region('/'.join(minmax))
        except: o_region = None

        return(o_region)
                     
    def _load(self):
        '''load and process a datalist'''

        with open(self._path, 'r') as fob:
            for dl in fob:
                if dl[0] != '#' and dl[0] != '\n' and dl[0] != '':
                    dl_cols = dl.split(' ')
                    if len(dl_cols) >= 2:
                        dl_cols = [x.strip() for x in dl_cols]
                        self.datalist.append(dl_cols)

    def _load_data(self):
        self._proc(self.datafiles)

    def _proc(self, datafiles = []):
        '''Recurse through the datalist and gather xyz data'''

        for i in self.datalist:
            dpath = i[0]
            dformat = int(i[1])

            if len(i) > 2: 
                dweight = i[2]
            else: dweight = 1

            ## ==============================================
            ## Datalist format -1
            ## Open as a datalist and _proc
            ## ==============================================

            if dformat == -1:
                d = datalist(dpath, self.region)
                d._proc(datafiles)

            ## ==============================================
            ## XYZ format 168
            ## check if passes tests and append to datafiles if so.
            ## ==============================================

            elif dformat == 168:
                usable = True
                path_d = os.path.join(self._path_dirname, dpath)

                if not os.path.exists(path_d): 
                    usable = False

                dinf_region = self._read_inf(path_d + '.inf')
                
                if self.region is not None and dinf_region is not None:
                    if not regions.regions_intersect_p(self.region, dinf_region):
                        usable = False

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

## =============================================================================
##
## Mainline - run datalists from console.
##
## =============================================================================

_usage = '''{} ({}): Process and generate datalists

usage: {} [ -hvR [ args ] ] datalist ...

Options:
  -R, --region\t\tSpecifies the desired region;

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} my_data.datalist -R -90/-89/30/31

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def main():

    status = 0
    i_datalist = None
    i_region = None
    want_verbose = False

    argv = sys.argv
        
    ## ==============================================
    ## parse command line arguments.
    ## ==============================================

    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            i_region = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-R':
            i_region = str(arg[2:])

        elif arg == '--help' or arg == '-h':
            print(_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_usage)
            sys.exit(0)

        else: 
            i_datalist = arg

        i = i + 1

    if i_datalist is None:
        print (_usage)
        sys.exit(1)

    ## ==============================================
    ## check platform and installed software
    ## ==============================================

    cmd_vers = utils._cmd_check()
    
    ## ==============================================
    ## Load the input datalist
    ## ==============================================
    

    pb = utils._progress('loading datalist...')
    this_datalist = datalist(i_datalist, i_region)
    if not this_datalist._valid: 
        status = -1
    pb.opm = 'loaded datalist <<\033[1m{}\033[m>>'.format(this_datalist._path_basename)
    pb.end(status)

    if status == 0:
        this_datalist._load_data()

        ## ==============================================
        ## print the data to stdout
        ## ==============================================

        this_datalist._cat_port(sys.stdout)

            
if __name__ == '__main__':
    main()


### End
