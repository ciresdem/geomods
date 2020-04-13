### waffles.py
##
## Copyright (c) 2013 - 2020 CIRES Coastal DEM Team
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

import threading
import time

import regions
import datalists
import gdalfun
import metadata
import utils
        
_version = '0.3.0'

def inc2str_inc(inc):
    '''convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)'''
    
    import fractions

    str_inc = str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', '')

    return(str_inc)

def this_year():
    '''return the current year'''

    import datetime

    return(datetime.datetime.now().strftime('%Y'))

## GMT Wrapper Functions

def grd2gdal(src_grd, dst_fmt = 'GTiff', verbose = False):
    '''Convert the grd file to tif using GMT'''
    
    status = 0
    if os.path.exists(src_grd):
        #dst_gdal = os.path.join(src_grd[:-4], gdalfun._fext(dst_fmt))
        dst_gdal = '{}.{}'.format(os.path.basename(src_grd).split('.')[0], gdalfun._fext(dst_fmt))
        #dst_gdal = '{}.{}'.format(len(os.path.basename(src_grd).split('.')[-1]), gdalfun._fext(dst_fmt))
        print dst_gdal
        grd2gdal_cmd = ('gmt grdconvert {} {}=gd+n-9999:{} -V\
        '.format(src_grd, dst_gdal, dst_fmt))
        out, status = utils.run_cmd(grd2gdal_cmd, verbose, verbose)

        if status != 0:
            dst_gdal = None
    else: dst_gdal = None

    return(dst_gdal)

def grdinfo(src_grd, verbose = False):
    '''Return an info list of `src_grd`'''

    status = 0
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose, True)
    if not os.path.exists(src_grd):
        return([])
    else:
        grdinfo_cmd = ('gmt grdinfo {} -C'.format(src_grd))
        out, status = utils.run_cmd(grdinfo_cmd, verbose, False)
        try:
            os.remove('gmt.conf')
        except: pass
        
        if status !=0:
            return([])
        else: return(out.split())

def gmtinfo(src_xyz, verbose = False):
    '''Return an info list of `src_xyz`'''

    status = 0
    out, status = utils.run_cmd('gmt gmtset IO_COL_SEPARATOR = SPACE', verbose, True)
    if not os.path.exists(src_xyz):
        return([])
    else:
        gmtinfo_cmd = ('gmt grdinfo {} -C'.format(src_xyz))
        out, status = utils.run_cmd(gmtinfo_cmd, verbose, False)
        try:
            os.remove('gmt.conf')
        except: pass
        
        if status !=0:
            return([])
        else: return(out.split())        

def grdcut(src_grd, src_region, dst_grd, verbose = False):
    '''Cut `src_grd` to `src_region` '''

    status = 0
    if os.path.exists(src_grd):
        cut_cmd1 = ('gmt grdcut -V {} -G{} {}'.format(src_grd, dst_grd, src_region.gmt))
        out, status = utils.run_cmd(cut_cmd1, verbose, verbose)

    return(status)

def grd2xyz(src_grd, dst_xyz, region = None, mask = None, verbose = False, want_datalist = False):
    '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region.
    Optionally, generate a datalist and inf file for the resultant xyz data.'''

    status = 0
    if mask:
        grdmask_cmd = ('gmt grdmath -N -V {} {} OR = tmp.grd'.format(src_grd, mask))
        out, status = utils.run_cmd(grdmask_cmd, verbose, verbose)
        if status == 0: 
            src_grd = 'tmp.grd'

    if region and region._valid:
        region_str = region.gmt
    else: region_str = ''

    grd2xyz_cmd = ('gmt grd2xyz -V {} -s {} > {}'.format(src_grd, region_str, dst_xyz))
    out, status = utils.run_cmd(grd2xyz_cmd, verbose, verbose)

    if status == 0:
        if mask:
            if os.path.exists('tmp.grd'):
                os.remove('tmp.grd')

        if want_datalist:
            #xyz_dir = os.path.join(os.getcwd(), 'xyz')
            #if not os.path.exists(xyz_dir):
            #    os.makedirs(xyz_dir)

            s_datalist = datalists.datalist('{}.datalist'.format(dst_xyz.split('.')[0]))
            s_datalist._append_datafile(['{}'.format(os.path.basename(dst_xyz)), 168, 1])
            s_datalist._reset()
        
            out, status = utils.run_cmd('mbdatalist -O -I{}'.format(s_datalist._path), False, True)
        
    return(status)

def slope(src_dem, dst_slp, verbose = False):
    '''Generate a Slope grid from a DEM with GMT'''

    status = 0
    o_b_name = '{}'.format(src_dem.split('.')[0])

    slope_cmd0 = ('gmt grdgradient -V -fg {} -S{}_pslp.grd -D -R{}\
    '.format(src_dem, o_name, src_dem))
    out, status = utils.run_cmd(slope_cmd0, verbose, True)

    if status == 0:
        slope_cmd1 = ('gmt grdmath -V {}_pslp.grd ATAN PI DIV 180 MUL = {}\
        '.format(o_b_name, dst_slp))
        out, status = utils.run_cmd(slope_cmd1, verbose, True)
        
    if os.path.exists('{}_pslp.grd'.format(o_b_name)):
        os.remove('{}_pslp.grd'.format(o_b_name))

    return(status)

def proximity(src_num_msk, dst_prox):
    '''Generate a Proximity grid with GDAL from a NUM-MSK'''

    status = 0
    gdalfun.proximity(src_num_msk, dst_prox)

    if not os.path.exists(dst_prox):
        status = -1

    return(status)

def num_msk(num_grd, dst_msk, verbose = False):
    '''Generate a num-msk from a NUM grid.'''

    status = 0

    num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = {}\
    '.format(num_grd, dst_msk))
    out, status = utils.run_cmd(num_msk_cmd, verbose, True)

    return(status)

def xyz2grd(datalist, region, inc, dst_name, a = 'n', node = 'pixel', verbose = False):
    '''Run the GMT command `xyz2grd` given a datalist, region and increment.'''
    
    status = 0
    if node == 'pixel':
        reg_str = '-r'
    else: reg_str = ''
    
    num_cmd0 = ('gmt xyz2grd -V {} -I{:.10f} -G{} -A{} {}\
    '.format(region.gmt, inc, dst_name, a, reg_str))
    out, status = utils.run_cmd(num_cmd0, verbose, verbose, datalist._dump_data)

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
    out, status = utils.run_cmd(mbgrid_cmd, verbose, verbose)

    return(out, status)

## =============================================================================
##
## DEM module: generate a Digital Elevation Model using a variety of methods
## dem modules include: 'mbgrid', 'surface', 'num', 'mean'
##
## Requires MBSystem, GMT and GDAL 
##
## =============================================================================

class dem:
    '''Generate a Digital Elevation Model using one of the dem modules.
    DEM Modules include, `mbgrid`, `surface`, `num`, `mean`, `bathy`'''

    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None,\
                 o_node = 'pixel', o_fmt = 'GTiff', o_extend = 6, clip_ply = None, callback = lambda: False, verbose = False):
        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        
        self.node = o_node
        self.o_fmt = o_fmt
        self.extend = int(o_extend)

        self.proc_region = self.region.buffer((self.extend * 2) * self.inc)
        self.dist_region = self.region.buffer(self.extend * self.inc)

        #self.datalist.region = self.proc_region
        
        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.o_name = o_name

        # if o_b_name is None:
        #     if o_name is None: 
        #         o_name = self.datalist._path_basename.split('.')[0]
        #     else: o_name = os.path.basename(o_name).split('.')[0]

        #     str_inc = inc2str_inc(self.inc)
        #     self.o_name = '{}{}_{}_{}'.format(o_name, str_inc, self.region.fn, this_year())
        # else: self.o_name = os.path.join(os.path.dirname(o_b_name), os.path.basename(o_b_name).split('.')[0])

        self.dem = '{}.grd'.format(self.o_name)
        self.pb = utils._progress()

        self.gc = utils.check_config(False, self.verbose)

        self.clip_ply = clip_ply
        
    def run(self, dem_mod = 'mbgrid', args = ()):
        '''Run the DEM module `dem_mod` using args `dem_mod_args`.'''

        #if dem_mod != 'mbgrid' and dem_mod != 'conversion-grid' and dem_mod != 'spatial-metadata' and dem_mod != 'uncertainty':
        #    if len(self.datalist.datafiles) == 0:
        #        self.datalist._load_data()
        #    if len(self.datalist.datafiles) == 0:
        #        self.status = -1

        #self.dem = '{}_{}.{}'.format(self.o_name, dem_mod, gdalfun._fext(self.o_fmt))

        args_d = {}
        for arg in args:
            p_arg = arg.split('=')
            args_d[p_arg[0]] = p_arg[1]

        self.o_name = '{}_{}'.format(self.o_name, dem_mod)
            
        if self.status == 0:
            _dem_mods[dem_mod][0](self)(**args_d)

        if self.clip_ply is not None:
            clip_args = {}
            cp = self.clip_ply.split(':')
            clip_args['src_ply'] = cp[0]
            cargs = cp[1:]
            for arg in cargs:
                p_arg = arg.split('=')
                clip_args[p_arg[0]] = p_arg[1]

            self.clip_dem(**clip_args)

        if self.status != 0:
            self.dem = None
        else:
            if self.o_fmt != 'GMT':
                gmt_dem = self.dem
                self.dem = grd2gdal(self.dem, self.o_fmt)
                utils.remove_glob(gmt_dem)
            
        return(self.dem)

    ## ==============================================
    ## run mbgrid on the datalist and generate the DEM
    ## note: mbgrid will cause popen to hang if stdout 
    ## is not cleared...should add output file to send to...
    ## ==============================================

    def clip_dem(self, src_ply = None, invert = False):

        if src_ply is not None:
            if os.path.exists(src_ply):
                if invert:
                    gr_inv = '-i'
                else: gr_inv = ''

                gmt_dem = self.dem
                self.dem = grd2gdal(self.dem)
                utils.remove_glob(gmt_dem)
                
                gi = gdalfun._infos(self.dem)
                print src_ply.split('.')
                gr_cmd = 'gdal_rasterize -burn {} {} -l {} {} {}'.format(gi['ndv'], gr_inv, os.path.basename(src_ply).split('.')[0], src_ply, self.dem)
                print gr_cmd
                utils.run_cmd(gr_cmd, self.verbose, self.verbose)

                gmt_dem = self.dem
                self.dem = grd2gdal(self.dem, 'netCDF')
                utils.remove_glob(gmt_dem)
                
            else: ## use gsshg via GMT
                if invert:
                    ns = '-N0/1/0/1/0'
                else: ns = '-N1/0/1/0/1'

                reg_str = ''
                if self.node == 'pixel':
                    reg_str = '-r'
                
                dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V {} {}\
                '.format(self.inc, self.dist_region.gmt, ns, reg_str))
                out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, self.verbose)

                if self.status == 0:
                    dem_landmask_cmd1 = ('gmt grdmath -V {} tmp_lm.grd MUL 0 NAN = {}_clp.grd\
                    '.format(self.dem, self.o_name))
                    out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, self.verbose)

                    os.rename('{}_clp.grd'.format(self.o_name), self.dem)
                    utils.remove_glob('tmp_lm.grd')
        
    def mbgrid(self, dist = '10/3', tension = 35, use_datalists = False):
        '''Generate a DEM with MBSystem's mbgrid program.'''

        if self.gc['MBGRID'] is None:
            utils._error_msg('MBSystem must be installed to use the MBGRID module')
            self.status = -1

        if self.gc['GMT'] is None:
            utils._error_msg('GMT must be installed to use the MBGRID module')
            self.status = -1
            
        if self.status == 0:

            if use_datalists:
                self.datalist._archive(dirname = '.mb_tmp_datalist')
                self.datalist = self.datalist.archive_datalist
            
            out, self.status = run_mbgrid(self.datalist, self.dist_region, self.inc, self.o_name, dist = dist, tension = tension, verbose = self.verbose)

            if self.status == 0:
                if self.node == 'pixel':
                    out, self.status = utils.run_cmd('gmt grdsample -T {}.grd -Gtmp.grd'.format(self.o_name), self.verbose, self.verbose)
                    if self.status == 0:
                        os.rename('tmp.grd', '{}.grd'.format(self.o_name))

                utils.remove_glob('*.cmd')

            self.dem = '{}.grd'.format(self.o_name)

            if use_datalists:
                import shutil
                shutil.rmtree('.mb_tmp_datalist')

    ## ==============================================
    ## Run GMT surface on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## ==============================================

    def surface(self, tension = .35, relaxation = 1.2, lower_limit = 'd', upper_limit = 'd'):
        '''Generate a DEM with GMT surface'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        dem_surf_cmd = ('gmt blockmean {} -I{:.10f} -V {} | gmt surface -V {} -I{:.10f} -G{}_p.grd -T{} -Z{} -Ll{} -Lu{} {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, tension, relaxation, lower_limit, upper_limit, reg_str))
        out, self.status = utils.run_cmd(dem_surf_cmd, self.verbose, self.verbose, self.datalist._dump_data)

        if self.status == 0:
            dem_cut_cmd = ('gmt grdcut -V {}_p.grd -G{}.grd {}\
            '.format(self.o_name, self.o_name, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)

            utils.remove_glob('{}_p.grd'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)
            
    ## ==============================================
    ## Run GMT triangulate on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## output may have nodata values around the edges.
    ## ==============================================

    def triangulate(self):
        '''Generate a DEM with GMT surface'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        dem_tri_cmd = ('gmt blockmean {} -I{:.10f} -V {} | gmt triangulate {} -I{:.10f} -V -G{}_t.grd {}\
        '.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd(dem_tri_cmd, self.verbose, self.verbose, self.datalist._dump_data)

        if self.status == 0:
            dem_cut_cmd = ('gmt grdcut -V {}_t.grd -G{}.grd {}\
            '.format(self.o_name, self.o_name, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)

            utils.remove_glob('{}_t.grd'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)

    ## ==============================================
    ## Run GMT nearneighbor on the datalist
    ## generate dem at 10x cells to account for
    ## edge effects
    ## output may have nodata values around the edges.
    ## ==============================================

    def nearneighbor(self, radius='6s'):
        '''Generate a DEM with GMT surface'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        #dem_nn_cmd = ('gmt blockmean {} -I{:.10f} -V {} | gmt nearneighbor {} -I{:.10f} -S{} -V -G{}_nn.grd {}\
        #'.format(self.proc_region.gmt, self.inc, reg_str, self.proc_region.gmt, self.inc, radius, self.o_name, reg_str))
        #out, self.status = utils.run_cmd(dem_nn_cmd, self.verbose, self.verbose, self.datalist._dump_data)

        dem_nn_cmd = ('gmt nearneighbor {} -I{:.10f} -S{} -V -G{}_nn.grd {}\
        '.format(self.proc_region.gmt, self.inc, radius, self.o_name, reg_str))
        out, self.status = utils.run_cmd(dem_nn_cmd, self.verbose, self.verbose, self.datalist._dump_data)

        if self.status == 0:
            dem_cut_cmd = ('gmt grdcut -V {}_nn.grd -G{}.grd {}\
            '.format(self.o_name, self.o_name, self.dist_region.gmt))
            out, self.status = utils.run_cmd(dem_cut_cmd, self.verbose, self.verbose)

            utils.remove_glob('{}_nn.grd'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)
            
    ## ==============================================
    ## run GMT xyz2grd on the datalist to generate
    ## an uninterpolated grid
    ## see the -A switch in GMT xyz2grd for mode options
    ## 
    ## `[d|f|l|m|n|r|S|s|u|z]
    ## By  default  we  will calculate mean values if multiple entries
    ## fall on the same node. Use -A to change this behavior, except it
    ## is ignored if -Z is given.
    ## Append f or s to simply keep the first or last data point that
    ## was assigned to each node. Append l or u or d to find the lowest
    ## (minimum) or upper (maximum) value or the difference between the
    ## maximum and miminum value at each node, respectively. Append m
    ## or r or S to compute mean or  RMS  value  or standard deviation at
    ## each node, respectively. Append n to simply count the number of
    ## data points that were assigned to each node (this only requires
    ## two input columns x and y as z is not consulted).
    ## Append z to sum multiple values that belong to the same node.`
    ## ==============================================

    def num(self, mode = 'n'):
        '''Generate a num and num-msk grid with GMT'''

        #self.dem = '{}_n{}.grd'.format(self.o_name, mode)
        out, self.status = xyz2grd(self.datalist, self.dist_region, self.inc, '{}.grd'.format(self.o_name), mode, self.node, verbose = self.verbose)

        self.dem = '{}.grd'.format(self.o_name)
            
    def mask(self, use_gmt = True):
        '''Generate a num and num-msk grid'''

        if self.verbose:
            pb = utils._progress('generating datalist MASK grid...')
        
        if self.gc['GMT'] is None:
            if self.node != 'pixel':
                self.dist_region = self.dist_region.buffer(self.inc * .5)
            self.datalist._load_data()
            self.dem = self.datalist.mask(region = self.dist_region.region, inc = self.inc, o_name = self.o_name, o_fmt = 'netCDF')
        else:
            self.num('n')

            if self.status == 0:
                num_msk_cmd = ('gmt grdmath -V {} 0 MUL 1 ADD 0 AND = tmp.grd\
                '.format(self.dem, self.o_name))
                out, status = utils.run_cmd(num_msk_cmd, self.verbose, self.verbose)

                if self.status == 0:
                    utils.remove_glob(self.dem)
                    os.rename('tmp.grd', '{}.grd'.format(self.o_name))

            self.dem = '{}.grd'.format(self.o_name)

        if self.verbose:
            pb.end(self.status, 'generated datalist MASK grid.')

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Inverse Distance DEM
    ## ==============================================

    def invdst(self, power = 2.0, smoothing = 0.0, radius1 = 0.1, radius2 = 0.1, angle = 0.0, max_points = 0, min_points = 0, nodata = 0.0):
        '''Generate an inverse distance grid with GDAL'''

        if self.node != 'pixel':
            self.dist_region = self.dist_region.buffer(self.inc * .5)
        
        out_size_x = int((self.dist_region.east - self.dist_region.west) / self.inc)
        out_size_y = int((self.dist_region.north - self.dist_region.south) / self.inc)

        self.status = datalists.datalist2csv(self.datalist, self.proc_region, self.inc, self.node, self.o_name, self.verbose)

        gg_cmd = 'gdal_grid -zfield "field_3" -txe {} {} -tye {} {} -outsize {} {} \
        -a invdist:power={}:smoothing={}:radius1={}:radius2={}:angle={}:max_points={}:min_points={}:nodata={} -l {} {}.vrt {}.grd -of netCDF --config GDAL_NUM_THREADS ALL_CPUS\
        '.format(self.dist_region.west, self.dist_region.east, self.dist_region.north, self.dist_region.south, out_size_x, out_size_y,
                 power, smoothing, radius1, radius2, angle, max_points, min_points, nodata, self.o_name, self.o_name, self.o_name)

        out, status = utils.run_cmd(gg_cmd, self.verbose, self.verbose)

        utils.remove_glob('{}.vrt'.format(self.o_name))
        utils.remove_glob('{}.csv'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)
        utils.remove_glob('gmt.conf')

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## an Moving Average DEM
    ## ==============================================
    
    def m_average(self, radius1 = 0.01, radius2 = 0.01, angle = 0.0, min_points = 0, nodata = 0.0):
        '''Generate a moving average grid with GDAL'''

        if self.node != 'pixel':
            self.dist_region = self.dist_region.buffer(self.inc * .5)
            
        out_size_x = int((self.dist_region.east - self.dist_region.west) / self.inc)
        out_size_y = int((self.dist_region.north - self.dist_region.south) / self.inc)
        
        self.status = datalists.datalist2csv(self.datalist, self.proc_region, self.inc, self.node, self.o_name, self.verbose)
        
        gg_cmd = 'gdal_grid -zfield "field_3" -txe {} {} -tye {} {} -outsize {} {}\
        -a average:radius1={}:radius2={}:angle={}:min_points={}:nodata={} -l {} {}.vrt {}.grd -of netCDF --config GDAL_NUM_THREADS ALL_CPUS\
        '.format(self.dist_region.west, self.dist_region.east, self.dist_region.north, self.dist_region.south, out_size_x, out_size_y,
                 radius1, radius2, angle, min_points, nodata, self.o_name, self.o_name, self.o_name)
        out, status = utils.run_cmd(gg_cmd, self.verbose, self.verbose)

        utils.remove_glob('{}.vrt'.format(self.o_name))
        utils.remove_glob('{}.csv'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)
        utils.remove_glob('gmt.conf')

    ## ==============================================
    ## run GDAL gdal_grid on the datalist to generate
    ## a linear triangulation DEM
    ## ==============================================
    
    def linear(self, radius = -1, nodata = 0.0):
        '''Generate a llinear triangulation grid with GDAL'''

        if self.node != 'pixel':
            self.dist_region = self.dist_region.buffer(self.inc * .5)
            
        out_size_x = int((self.dist_region.east - self.dist_region.west) / self.inc)
        out_size_y = int((self.dist_region.north - self.dist_region.south) / self.inc)
        
        self.status = datalists.datalist2csv(self.datalist, self.proc_region, self.inc, self.node, self.o_name, self.verbose)
        
        gg_cmd = 'gdal_grid -zfield "field_3" -txe {} {} -tye {} {} -outsize {} {}\
        -a linear:radius={}:nodata={} -l {} {}.vrt {}.grd -of netCDF --config GDAL_NUM_THREADS ALL_CPUS\
        '.format(self.dist_region.west, self.dist_region.east, self.dist_region.north, self.dist_region.south, out_size_x, out_size_y,
                 radius, nodata, self.o_name, self.o_name, self.o_name)
        out, status = utils.run_cmd(gg_cmd, self.verbose, self.verbose)

        utils.remove_glob('{}.vrt'.format(self.o_name))
        utils.remove_glob('{}.csv'.format(self.o_name))

        self.dem = '{}.grd'.format(self.o_name)
        utils.remove_glob('gmt.conf')
        
    ## ==============================================
    ## Bathy-Surface module.
    ##
    ## with GMT, gererate a 'bathy-surface'.
    ## specify the coastline mask shapefile with 'mask = '
    ## if mask is not specified, will use gsshg.
    ## returns a netcdf masked grid and an XYZ file
    ## of the bathy surface points.
    ## ==============================================

    def bathy(self, mask = 'gsshg'):
        '''Generate a masked surface with GMT surface and a breakline mask polygon.
        Use a coastline polygon mask to generatre a `bathy-surface` DEM.'''

        reg_str = ''
        if self.node == 'pixel':
            reg_str = '-r'

        self.dem = '{}_bs.grd'.format(self.o_name)

        b_inc = self.inc * 3.0
        print b_inc
        if mask == 'gsshg':
            dem_landmask_cmd = ('gmt grdlandmask -Gtmp_lm.grd -I{:.7f} {} -Df+ -V -N1/0 {}\
            '.format(b_inc, self.proc_region.gmt, reg_str))
            out, self.status = utils.run_cmd(dem_landmask_cmd, self.verbose, self.verbose)
            mask = 'tmp_lm.grd'

        dem_surf_cmd = ('gmt blockmean {} -I{:.7f} -V {} | gmt surface -V {} -I{:.7f} -G{}_p.grd -T.35 -Z1.2 -Lu0 {}\
        '.format(self.proc_region.gmt, b_inc, reg_str, self.proc_region.gmt, b_inc, self.o_name, reg_str))
        out, self.status = utils.run_cmd(dem_surf_cmd, self.verbose, True, self.datalist._dump_data)
        
        if self.status == 0:
            dem_landmask_cmd1 = ('gmt grdmath -V {}_p.grd {} MUL 0 NAN = {}_p_bathy.grd\
            '.format(self.o_name, mask, self.o_name))
            out, self.status = utils.run_cmd(dem_landmask_cmd1, self.verbose, self.verbose)

            grdcut('{}_p_bathy.grd'.format(self.o_name), self.dist_region, self.dem)
            utils.remove_glob('{}_p*.grd'.format(self.o_name))

        if self.status == 0:

            xyz_dir = os.path.join(os.getcwd(), 'xyz')
            if not os.path.exists(xyz_dir):
                os.makedirs(xyz_dir)
                
            bathy_xyz = 'xyz/{}_bs.xyz'.format(self.o_name)                
            self.status = grd2xyz(self.dem, bathy_xyz, want_datalist = True)

    ## ==============================================
    ## conversion-grid module:
    ## Generate a VDatum conversion grid
    ##
    ## Requires NOAA's VDatum version >= 4.x
    ## ==============================================

    def vdatum(self, ivert = 'navd88', overt = 'mllw', region = '3'):
        '''generate a vertical datum transformation grid with vdatum'''

        import vdatum

        if self.gc['VDATUM'] is None:
            utils._error_msg('NOAAs VDATUM must be installed to use the VDATUM module.')
            self.status = -1

        if self.status == 0:
            reg_str = ''
            if self.node == 'pixel':
                reg_str = '-r'

            gdalfun.null('empty.tif', self.dist_region.region, 0.00083333, nodata = 0)
            with open('empty.xyz', 'w') as mt_xyz:
                gdalfun.dump('empty.tif', dst_xyz = mt_xyz, dump_nodata = True)

            this_vd = vdatum.vdatum(vdatum_path = self.gc['VDATUM'], verbose = self.verbose)
            this_vd.ivert = ivert
            this_vd.overt = overt

            #dem_name = '{}to{}'.format(this_vd.ivert, this_vd.overt)
            #self.dem = '{}_{}.grd'.format(self.o_name.split('.')[0], dem_name)
            self.dem = '{}.grd'.format(self.o_name)

            this_vd.run_vdatum('empty.xyz')

            if os.path.exists('result/empty.xyz') and os.stat('result/empty.xyz').st_size != 0:
                #out, status = utils.run_cmd('gmt gmtinfo result/empty.xyz -C')
                #empty_infos = out.split()
                empty_infos = gmtinfo('result/empty.xyz', self.verbose)

                if empty_infos[4] > 0:
                    ll_switch = '-Lld'
                else: ll_switch = '-Ll0'

                if empty_infos[5] > 0:
                    lu_switch = '-Lud'
                else: lu_switch = '-Lu0'

                gc = 'gmt blockmean result/empty.xyz -V -I{} {} {}\
                | gmt surface -I{} {} -G{} -V -T0 {} {} {}\
                '.format(self.inc, self.dist_region.gmt, reg_str, 
                         self.inc, self.dist_region.gmt, self.dem, ll_switch, lu_switch, reg_str)
                utils.run_cmd(gc, self.verbose, True)
            else: self.status = -1

            try:
                utils.remove_glob('empty.*')
                utils.remove_glob('result/*')
                os.removedirs('result')
            except: pass
        else: self.dem = None
        #return(self.dem)

    def spatial_metadata(self, epsg = 4269):
        #self.datalist.i_fmt = -1
        #self.datalist._load_data()
        sm = metadata.spatial_metadata(self.datalist, self.region, i_inc = self.inc, o_name = self.o_name, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        sm.want_queue = True
        sm.run(epsg)

    def uncertainty(self, dem_mod = 'mbgrid', dem = None, msk = None, prox = None):
        unc = uncertainty.uncertainty(self.datalist, self.region, i_inc = self.inc, o_name = self.o_name, o_node = self.node, o_extend = self.extend, callback = self.stop, verbose = self.verbose)
        unc.run(dem_mod, dem, msk)
    
## =============================================================================
##
## Console Scripts
##
## Scripts are: 
## waffles
##
## =============================================================================

_dem_mods = {
    'mbgrid': [lambda x: x.mbgrid, 'Weighted SPLINE DEM via mbgrid', ':tension=35:dist=3/10:use_datalists=False'],
    'surface': [lambda x: x.surface, 'SPLINE DEM via GMT surface', ':tension=.35:relaxation=1.2:lower_limit=d:upper_limit=d'],
    'triangulate': [lambda x: x.triangulate, 'TRIANGULATION DEM via GMT triangulate', None],
    'nearneighbor': [lambda x: x.nearneighbor, 'NEARNEIGHBOR DEM via GMT', 'radius=6s'],
    'num': [lambda x: x.num, 'Uninterpolated DEM via GMT xyz2grd', ':mode=n'],
    'mask': [lambda x: x.mask, 'Data MASK grid', 'None'],
    'invdst': [lambda x: x.invdst, 'Inverse Distance DEM via gdal_grid', ':power=2.0:smoothing=0.0:radus1=0.1:radius2:0.1'],
    'average': [lambda x: x.m_average, 'Moving Average DEM via gdal_grid', ':radius1=0.01:radius2=0.01'],
    'linear': [lambda x: x.linear, 'Linear triangulation DEM via gdal_grid', ':radius=-1'],
    'vdatum': [lambda x: x.vdatum, 'VDATUM transformation grid', ':ivert=navd88:overt=mhw:region=3'],
}
#'spatial-metadata': [lambda x: x.spatial_metadata, 'DEM SPATIAL METADATA', ':epsg'],
#'uncertainty': [lambda x: x.uncertainty, 'DEM UNCERTAINTY grid <beta>', ':gridding-module'],

def dem_mod_desc(x):
    out = 'Modules:\n'
    for key in x:
        out += '  {:16}\t{} [{}]\n'.format(key, x[key][-2], x[key][-1])
    return(out)

_waffles_usage = '''{} ({}): Process and generate Digital Elevation Models and derivatives

usage: {} [ -ahprsuvCEFIORX [ args ] ] module[:parameter=value]* ...

{}
Options:
  -R, --region\t\tSpecifies the desired REGION;
\t\t\tThis can either be a GMT-style region ( -R xmin/xmax/ymin/ymax )
\t\t\tor an OGR-compatible vector file with regional polygons. 
\t\t\tIf a vector file is supplied it will search each region found therein.
  -I, --datalist\tThe input DATALIST.
  -E, --increment\tGridding CELL-SIZE in native units or GMT-style increments.
  -F, --format\t\tOutput grid FORMAT. [GTiff]
  -O, --output-name\tOutput naming BASENAME.
  -X, --extend\t\tNumber of cells with which to EXTEND the REGION. [6]
  -C, --clip\t\tCLIP the output to the clip polygon. [clip_ply.shp:invert=False]

  -a, --archive\t\tArchive the data from the DATALIST in the REGION.
  -p, --prefix\t\tSet BASENAME to PREFIX (append inc/region/year/module info to output BASENAME).
  -r, --grid-node\tuse grid-node registration, default is pixel-node
  -s, --spat-metadata\tGenerate associated SPATIAL-METADATA for the REGION.
  -u, --uncertainty\tGenerate an associated DEM UNCERTAINTY grid. <beta>

  --help\t\tPrint the usage text
  --version\t\tPrint the version information
  --verbose\t\tIncrease the verbosity

 Examples:
 % {} -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27 -V surface:tension=.7
 % {} -I input.datalist -E .3333333s -X 2 -R input_tiles_ply.shp -V -r -s -u mbgrid
 % {} -R-82.5/-82.25/26.75/27 -E0.0000925 vdatum:i_vdatum=navd88:o_vdatum=mhw:vd_region=3 -O ncei -p -r

CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\
'''.format( os.path.basename(sys.argv[0]), 
            _version, 
            os.path.basename(sys.argv[0]), 
            dem_mod_desc(_dem_mods),
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]), 
            os.path.basename(sys.argv[0]))

def main():
    status = 0
    i_region = None
    i_datalist = None
    i_inc = 0.0000925925
    these_regions = []
    stop_threads = False
    want_verbose = False
    want_sm = False
    want_unc = False
    want_archive = False
    want_prefix = False
    mod_opts = {}
    o_bn = None
    o_fmt = 'GTiff'
    node_reg = 'pixel'
    o_extend = 6
    clip_ply = None

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

        elif arg == '--datalist' or arg == '-I':
            i_datalist = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-I':
            i_datalist = str(arg[2:])

        elif arg == '--increment' or arg == '-E':
            i_inc = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-E':
            i_inc = arg[2:]

        elif arg == '--format' or arg == '-F':
            o_fmt = argv[i + 1]
            i = i + 1
        elif arg[:2] == '-F':
            o_fmt = arg[2:]

        elif arg == '--output-name' or arg == '-O':
            o_bn = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-O':
            o_bn = str(arg[2:])
        
        elif arg == '--extend' or arg == '-X':
            o_extend = int(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-X':
            o_extend = int(arg[2:])

        elif arg == '--clip' or arg == '-C':
            clip_ply = str(argv[i + 1])
            i = i + 1
        elif arg[:2] == '-C':
            clip_ply = str(arg[2:])

        elif arg == '--prefix' or arg == '-p':
            want_prefix = True
            
        elif arg == '--grid-node' or arg == '-r':
            node_reg = 'grid'

        elif arg == '--archive' or arg == '-a':
            want_archive = True

        elif arg == '--spat-metadata' or arg == '-s':
            want_sm = True

        elif arg == '--uncertainty' or arg == '-u':
            want_unc = True

        elif arg == '--help' or arg == '-h':
            print(_waffles_usage)
            sys.exit(1)

        elif arg == '--version' or arg == '-v':
            print('{}, version {}\n{}'.format(os.path.basename(sys.argv[0]), _version, utils._license))
            sys.exit(1)

        elif arg == '--verbose' or arg == '-V':
            want_verbose = True

        elif arg[0] == '-':
            print(_waffles_usage)
            sys.exit(0)

        else: 
            opts = arg.split(':')
            if opts[0] in _dem_mods.keys():
                mod_opts[opts[0]] = list(opts[1:])
            else: utils._progress().err_msg('invalid module name `{}`'.format(opts[0]))

        i = i + 1

    if i_region is None and i_datalist is None:
        utils._error_msg('a region and/or a datalist must be specified')
        print(_waffles_usage)
        sys.exit(1)

    try:
        #i_inc = float(i_inc)
        i_inc = regions._inc(i_inc)
    except:
        utils._error_msg('the increment value should be a number')
        print(_waffles_usage)
        sys.exit(1)

    for key in mod_opts.keys():
        mod_opts[key] = [None if x == '' else x for x in mod_opts[key]]

    ## this crashes a lot; re-think!
    gc = utils.check_config(False, want_verbose)
    
    ## ==============================================
    ## process input region(s) and loop
    ## ==============================================

    if want_verbose: pb = utils._progress('loading region(s)...')
    if i_region is None:
        these_regions = [None]
    else:
        try: 
            these_regions = [regions.region(i_region)]
        except: these_regions = [regions.region('/'.join(map(str, x))) for x in gdalfun._ogr_extents(i_region)]
    
    if len(these_regions) == 0:
        these_regions = [None]

    for this_region in these_regions:
        if this_region is not None:
            if not this_region._valid:
                status = -1
    if want_verbose: pb.end(status, 'loaded \033[1m{}\033[m region(s).'.format(len(these_regions)))
            
    if status == -1:
        utils._error_msg('failed to load region(s)')
        print(_waffles_usage)
        sys.exit(1)

    for rn, this_region in enumerate(these_regions):

        status = 0
        
        ## ==============================================
        ## Load the input datalist
        ## ==============================================

        if stop_threads: break

        if i_datalist is not None:
            this_datalist = datalists.datalist(i_datalist, this_region.buffer((o_extend * 2) * i_inc), verbose = want_verbose)
            #this_datalist = datalists.datalist(i_datalist, this_region, verbose = want_verbose)
            if not this_datalist._valid_p():
                utils._error_msg('invalid datalist')
                status = -1
        else: this_datalist = None

        if status != 0: break

        if this_region is None and this_datalist is not None:
            utils._msg('no region supplied, gathering from datalist...')
            this_datalist.gather_region()
            this_region = this_datalist.region

        if this_region is None and this_datalist is None:
            utils._error_msg('no region or datalist; aborting...')
            status = -1
            break
        
        if o_bn is None:
            if this_datalist is None:
                o_name = 'waffles'
            else: o_name = this_datalist._path_basename.split('.')[0]
        #else: o_name = o_bn.split('.')[0]
        else: o_name = os.path.join(os.path.dirname(o_bn), os.path.basename(o_bn).split('.')[0])
        
        if want_prefix:
            str_inc = inc2str_inc(i_inc)
            o_name = '{}{}_{}_{}'.format(o_name, str_inc, this_region.fn, this_year())
        
        if this_datalist is not None:
            if want_sm:
                if i_inc < 0.0000925925:
                    tmp_inc = i_inc * 3
                    o_extend = 2
                else: tmp_inc = i_inc

                if node_reg != 'pixel':
                    this_region = this_region.buffer(tmp_inc * .5)
                
                sm = metadata.spatial_metadata(this_datalist, this_region, i_inc = tmp_inc, o_name = o_name,\
                                               o_extend = o_extend, callback = lambda: stop_threads, verbose = want_verbose)
                sm.want_queue = True
                s = threading.Thread(target = sm.run, args = ())
                s.start()

            if want_archive:
                a = threading.Thread(target = this_datalist._archive, args = ())
                arc_datalist = a.start()
                
        for dem_mod in mod_opts.keys():

            ## ==============================================
            ## Run the DEM module
            ## ==============================================
            
            ## most modules need a datalist and a region and an increment to be valid, except vdatum doesn't need a datalist...
            if this_datalist is None:
                if dem_mod.lower() != 'vdatum':
                    utils._error_msg('a valid datalist is needed to run module: {}...'.format(dem_mod))
                    break
                else: pass
            
            args = tuple(mod_opts[dem_mod])

            utils._msg('Module: {}'.format(dem_mod))
            utils._msg('Module Options: {}'.format(args))
            utils._msg('Datalist: {}'.format(this_datalist._path_basename))
            utils._msg('Region: {}'.format(this_datalist.region.gmt))
            utils._msg('Increment: {}'.format(i_inc))
            utils._msg('Output: {}'.format(o_name))
            
            pb = utils._progress('running geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string))
            dl = dem(this_datalist, this_region, i_inc = i_inc, o_name = o_name, o_node = node_reg,\
                     o_fmt = o_fmt, o_extend = o_extend, clip_ply = clip_ply, callback = lambda: stop_threads, verbose = want_verbose)
            t = threading.Thread(target = dl.run, args = (dem_mod, args))
            
            try:
                t.start()
                while True:
                    time.sleep(1)
                    pb.update()
                    if not t.is_alive():
                        break
            except (KeyboardInterrupt, SystemExit): 
                dl.status = -1
                stop_threads = True

            t.join()
            pb.end(dl.status, 'ran geomods dem module \033[1m{}\033[m on region ({}/{}): \033[1m{}\033[m...\
            '.format(dem_mod.upper(), rn + 1, len(these_regions), this_region.region_string))

            if want_unc and this_datalist is not None and dem_mod.lower() != 'vdatum':
                import uncertainty
                
                pb = utils._progress('running geomods \033[1mUNCERTAINTY\033[m module on region ({}/{}): \033[1m{}\033[m...\
                '.format(rn + 1, len(these_regions), this_region.region_string))
                dp = None
                cnt = 0

                while True:
                    cnt += 1
                    unc = uncertainty.uncertainty(this_datalist, this_region, i_inc = i_inc, o_name = o_name, o_node = node_reg, o_extend = o_extend, callback = lambda: stop_threads, verbose = want_verbose)
                    ## todo: account for possible dem args for uncertainty...
                    dp = unc.run(dem_mod, args, dl.dem)
                    if len(dp) > 0 or cnt == 10:
                        break
                pb.end(0, 'ran geomods \033[1mUNCERTAINTY\033[m module on region ({}/{}): \033[1m{}\033[m...\
                '.format(rn + 1, len(these_regions), this_region.region_string))

        if this_datalist is not None:
            if want_sm: s.join()
            if want_archive: a.join()

if __name__ == '__main__':
    main()

### END
