### vdatumfun.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
### Code:

import os

## ==============================================
## import geomods
## ==============================================
from geomods import utils

## =============================================================================
##
## VDatum - vdatumfun.py
## wrapper functions for NOAA's VDatum
##
## Currently only compatible with VDatum >= 4.0
##
## TODO: add all vdatum cli options
## =============================================================================
_vd_config = {
    'jar': None,
    'ivert': 'navd88:m:height',
    'overt': 'mhw:m:height',
    'ihorz': 'NAD83_2011',
    'ohorz': 'NAD83_2011',
    'region': '3',
    'fmt': 'txt',
    'xyzl': '0,1,2',
    'skip': '0',
    'delim': 'space',
    'result_dir': 'result',
    'verbose': False,
}

def vdatum_locate_jar():
    '''Find the VDatum executable on the local system.

    returns a list of found vdatum.jar system paths'''
    
    results = []
    for root, dirs, files in os.walk('/'):
        if 'vdatum.jar' in files:
            results.append(os.path.abspath(os.path.join(root, 'vdatum.jar')))
            break
    if len(results) == 0:
        return(None)
    else: return(results)

def vdatum_get_version(vd_config = _vd_config):
    '''run vdatum and attempt to get it's version
    
    return the vdatum version or None'''
    
    if vd_config['jar'] is None:
        vd_config['jar'] = vdatum_locate_jar()
    if vd_config['jar'] is not None:
        out, status = utils.run_cmd('java -jar {} {}'.format(vd_config['jar'], '-'), verbose = self.verbose)
        for i in out.decode('utf-8').split('\n'):
            if '- v' in i.strip():
                return(i.strip().split('v')[-1])
    return(None)

def vdatum_xyz(xyz, vd_config = _vd_config):
    '''run vdatum on an xyz list [x, y, z]

    returns the transformed xyz list'''
    
    if vd_config['jar'] is None: vd_config['jar'] = vdatum_locate_jar()[0]
    if vd_config['jar'] is not None:
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -pt:{},{},{} region:{}\
        '.format(vd_config['ihorz'], vd_config['ivert'], vd_config['ohorz'], vd_config['overt'], \
                 xyz[0], xyz[1], xyz[2], vd_config['region'])
        out, status = utils.run_cmd('java -Djava.awt.headless=false -jar {} {}'.format(vd_config['jar'], vdc), verbose = False)
        for i in out.split('\n'):
            if 'Height/Z' in i:
                z = float(i.split()[2])
                break
        return([xyz[0], xyz[1], z])
    else: return(xyz)

def vdatum_clean_result(result_f = 'result'):
    '''clean the vdatum 'result' folder'''
    
    utils.remove_glob('{}/*'.format(result_f))
    try:
        os.removedirs(result_f)
    except: pass
    
def run_vdatum(src_fn, vd_config = _vd_config):
    '''run vdatum on src_fn which is an XYZ file
    use vd_config to set vdatum parameters.

    returns [command-output, command-return-code]'''
    
    if vd_config['jar'] is None: vd_config['jar'] = vdatum_locate_jar()[0]
    if vd_config['jar'] is not None:
        vdc = 'ihorz:{} ivert:{} ohorz:{} overt:{} -nodata -file:txt:{},{},skip{}:{}:{} region:{}\
        '.format(vd_config['ihorz'], vd_config['ivert'], vd_config['ohorz'], vd_config['overt'], \
                 vd_config['delim'], vd_config['xyzl'], vd_config['skip'], src_fn, vd_config['result_dir'], vd_config['region'])
        #return(utils.run_cmd('java -jar {} {}'.format(vd_config['jar'], vdc), verbose = True))
        return(utils.run_cmd('java -Djava.awt.headless=true -jar {} {}'.format(vd_config['jar'], vdc), verbose = vd_config['verbose']))
    else: return([], -1)

### End
