### vdatum.py
##
## Copyright (c) 2010 - 2020 CIRES Coastal DEM Team
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
import ConfigParser
import utils

_version = '0.1.0'

## =============================================================================
##
## VDatum class for working with NOAA's VDatum program
## Currently only compatible with VDatum > 4.0
##
## =============================================================================

class vdatum:
    '''vdatum object to communicate with NOAA's VDatum'''

    def __init__(self, vdatum_path = None, verbose = False):
        self.verbose = verbose
        self.status = 0
        if vdatum_path is None:
            co = ConfigParser.ConfigParser()
            try:
                co.read(utils.CONFIG_FILE)
                self.vdatum_path = co.get('VDATUM', 'jar')
            except:
                self.vdatum_path = self._find_vdatum()[0]
                if self.status == 0:
                    co.add_section('VDATUM')
                    co.set('VDATUM', 'jar', self.vdatum_path)
            
                    with open(utils.CONFIG_FILE, 'w') as conf:
                        co.write(conf)
                else: self.vdatum_path = None
        else: self.vdatum_path = vdatum_path

        self._version = None
        self._get_version()
        
        self.ivert = 'navd88:m:height'
        self.overt = 'mhw:m:height'
        self.ihorz = 'NAD83_2011'
        self.ohorz = 'NAD83_2011'

        self.region = '3'
        self.ft = 'txt'
        self.fd = 'space'

        self.ds_dir = 'result'

    def _get_version(self):
        if self.vdatum_path is not None:
            out, status = utils.run_cmd('java -jar {} {}'.format(self.vdatum_path, '-'), prog = False)
            for i in out.split('\n'):
                if '- v' in i.strip():
                    self._version = i.strip().split('v')[-1]
                    break

    def _find_vdatum(self):
        '''Find the VDatum executable on the system and 
        return a list of found vdatum.jar paths'''

        results = []
        if self.verbose:
            pb = utils._progress('checking for vdatum.')
        for root, dirs, files in os.walk('/'):
            if 'vdatum.jar' in files:
                results.append(os.path.join(root, 'vdatum.jar'))
        if len(results) <= 0:
            self.status = -1
        if self.verbose:
            pb.opm = '{}..{}'.format(pb.opm, results[0])
            pb.end(self.status)

        return(results)

    def run_vdatum(self, src_fn):
        '''Run vdatum on src_fn which is an XYZ file'''
        
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},0,1,2:{}:{} region:{}\
        '.format(self.ihorz, self.ivert, self.ohorz, self.overt, self.fd, src_fn, self.ds_dir, self.region)
        out, status = utils.run_cmd('java -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, self.verbose)
        #out, status = run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(self.vdatum_path, vdc), self.verbose, True)

        return(status)

### End
