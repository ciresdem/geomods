### mbsfun.py
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

## =============================================================================
##
## MB-System Wrapper Functions - mbsfun.py
##
## =============================================================================

def mb_inf(src_xyz, src_fmt = 168):
    '''generate an info (.inf) file from a src_xyz file using MBSystem.'''

    out, status = run_cmd('mbdatalist -O -F{} -I{}'.format(src_fmt, src_xyz), False)
    return(out, status)

def run_mbgrid(datalist, region, inc, dst_name, dist = '10/3', tension = 35, extras = False, verbose = False):
    '''Run the MBSystem command `mbgrid` given a datalist, region and increment.'''
    
    status = 0
    if extras:
        e_switch = '-M'
    else: e_switch = ''

    if len(dist.split('/')) == 1:
        dist = dist + '/2'
    
    mbgrid_cmd = ('mbgrid -I{} {} -E{:.10f}/{:.10f}/degrees! -O{} -A2 -G100 -F1 -N -C{} -S0 -X0.1 -T{} {} > mb_proc.txt \
    '.format(datalist._path, region.gmt, inc, inc, dst_name, dist, tension, e_switch))
    out, status = run_cmd(mbgrid_cmd, verbose = verbose)

    return(out, status)

### End
