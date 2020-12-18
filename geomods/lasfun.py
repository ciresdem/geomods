### lasfun.py
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
### Commentary:
### Code:

import sys

from geomods import utils
from geomods import regions
from geomods import xyzfun

## ==============================================
## las-file processing (datalists fmt:300)
## ==============================================
def las_inf(src_las):
    '''scan an xyz file and find it's min/max values and
    write an associated inf file for the src_xyz file.

    returns region [xmin, xmax, ymin, ymax, zmin, zmax] of the src_xyz file.'''

    minmax = []
    out, status = utils.run_cmd('lasinfo -nc -nv -stdout -i {}'.format(src_las), verbose = False)
    for i in out.split('\n'):
        if 'min x y z' in i:
            xyz_min = [float(y) for y in [x.strip() for x in i.split(':')][1].split()]
        if 'max x y z' in i:
            xyz_max = [float(y) for y in [x.strip() for x in i.split(':')][1].split()]

    minmax = [xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2]]

    with open('{}.inf'.format(src_las), 'w') as inf:
        utils.echo_msg('generating inf file for {}'.format(src_las))
        inf.write('{}\n'.format(' '.join([str(x) for x in minmax])))
    
    return(minmax)
    
def las_inf_entry(entry):
    '''find the region of the xyz datalist entry
    
    returns the region [xmin, xmax, ymin, ymax, zmin, zmax] of the xyz entry'''

    return(las_inf(entry[0]))
    
def las_yield_entry(entry, region = None, verbose = False, z_region = None):
    '''yield the xyz data from the xyz datalist entry

    yields [x, y, z, <w, ...>]'''
    ln = 0
    delim = ' '
    if z_region is not None:
        min_z = None if z_region[0] is None else z_region[0]
        max_z = None if z_region[1] is None else z_region[1]
        #z_region = ['-' if x is None else str(x) for x in z_region]
    else: min_z = max_z = None
    #out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose = False)
    for line in utils.yield_cmd('las2txt -parse xyz -stdout -keep_class 2 29 -i {} {} {} {}\
    '.format(entry[0], '' if region is None else '-keep_xy {}'.format(regions.region_format(region, 'sstr')),\
             '' if min_z is None else '-drop_z_below {}'.format(min_z),\
             '' if max_z is None else '-drop_z_above {}'.format(max_z)), data_fun = None, verbose = False):
        ln += 1
        yield(line.split())
    if verbose: utils.echo_msg('read {} data points from {}'.format(ln, entry[0]))

def las_dump_entry(entry, dst_port = sys.stdout, region = None, verbose = False, z_region = None):
    '''dump the las data from the las datalist entry to dst_port'''
    
    for xyz in las_yield_entry(entry, region, verbose, z_region):
        xyzfun.xyz_line(xyz, dst_port, True, None)        

### End
