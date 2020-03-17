### uncertainty.py
##
## Copyright (c) 2019 - 2020 CIRES Coastal DEM Team
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

import random
import numpy as np

import regions
import datalists
import gdalfun
import utils
from waffles import *

from scipy import optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

## =============================================================================
##
## uncertainty module:
##
## =============================================================================

def err2coeff(my_data, coeff_guess = [0, 0.1, 0.2], dst_name = None):
    '''data is 2 col file with `err dist`'''

    #try: 
    #    my_data = np.loadtxt(data, delimiter=' ')
    #except: sys.exit(2)

    error=my_data[:,0]
    distance=my_data[:,1]

    max_int_dist = np.max(distance)
    nbins = 10

    #coeff_guess=[0, 0.1, 0.2]
    n, _ = np.histogram(distance, bins = nbins)

    # want at least 2 values in each bin?
    while 0 or 1 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins = nbins)

    serror, _ = np.histogram(distance, bins = nbins, weights = error)
    serror2, _ = np.histogram(distance, bins = nbins, weights = error**2)

    mean = serror / n
    std = np.sqrt(serror2 / n - mean * mean)

    ydata = np.insert(std, 0, 0)
    
    bins_orig=(_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)

    fitfunc = lambda p, x: p[0] + p[1] * (abs(x) ** abs(p[2]))
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)

    # fit
    
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    #plt.show()

    if dst_name is None:
        out_png = 'unc_best_fit.png'
    else: out_png = dst_name + 'bf.png'
    plt.savefig(out_png)   # save the figure to file
    plt.close()

    #scatter

    plt.scatter(distance, error)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')

    if dst_name is None:
        out_png = 'unc_scatter.png'
    else: out_png = dst_name + 'scatter.png'
    plt.savefig(out_png)
    plt.close()

    return(out)

def hav_dst(pnt0, pnt1):

    import math
    
    x0=float(pnt0[0])
    y0=float(pnt0[1])
    x1=float(pnt1[0])
    y1=float(pnt1[1])
    
    rad = 637100 # m
    
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    
    a = math.sin(dx / 2) * math.sin(dx / 2) + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy / 2) * math.sin(dy / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = rad * c
    return(d)

def gmtselect_split(o_xyz, sub_region, sub_bn, verbose = False):
    '''split an xyz file into an inner and outer region.'''

    status = 0
    out_inner = None
    out_outer = None

    gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
    out, status = utils.run_cmd(gmt_s_inner, verbose, False)

    if status == 0:
        out_inner = '{}_inner.xyz'.format(sub_bn)

    gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
    out, status = utils.run_cmd(gmt_s_outer, verbose, False)

    if status == 0:
        out_outer = '{}_outer.xyz'.format(sub_bn)

    return([out_inner, out_outer])

class uncertainty:

    def __init__(self, i_datalist, i_region, i_inc = 0.000277777, o_name = None, callback = lambda: False, verbose = False):

        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.proc_region = self.region.buffer(10 * self.inc)
        self.dist_region = self.region.buffer(6 * self.inc)
        self.node = 'pixel'

        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.tw = utils._progress()
        
        self.dem = { 
            'dem': None,
            'num': None,
            'msk': None,
            'prox': None,
            'int-unc': None,
        }

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name

        self.region_info = {}
        self.sub_zones = {}
        
    def run(self, dem_mod = 'mbgrid', dem = None, msk = None):
        
        dp = None
        
        if dem is not None: self.dem['dem'] = dem
        if msk is not None: self.dem['msk'] = msk
        self.datalist._load_data()

        self.set_or_make_dem(dem_mod)
        if self.verbose: self.tw.msg(self.dem)

        ## ==============================================
        ## Calculate the percentage of filled cells
        ## and proximity percentiles.
        ## ==============================================
        
        self.region_info = self.region_analysis()
        if self.verbose:
            for x in self.region_info.keys():
                self.tw.msg('region: {}: {}'.format(x, self.region_info[x]))
        
        utils._progress('running \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m'.format(dem_mod))

        sims = 4
        chnk_lvl = 4
                    
        ## ==============================================
        ## chunk region into sub regions
        ## ==============================================
        
        utils._progress('chunking region into sub-regions...')
        chnk_inc = int(chnk_lvl * self.region_info[self.o_name][4])
        self.sub_regions = self.region.chunk(self.inc, chnk_inc)
        if self.verbose: print([x.region for x in self.sub_regions])
        self.tw.end(self.status, 'chunked region into {} sub-regions.'.format(len(self.sub_regions)))

        utils._progress('analyzing {} sub-regions...'.format(len(self.sub_regions)))
            
        self.sub_zones = self.tile_analysis()
        if self.verbose:
            for x in self.sub_zones.keys():
                self.tw.msg('Sub-region {}: {}'.format(x, self.sub_zones[x]))
                
        utils._progress('analyzed {} sub-regions.'.format(len(self.sub_regions)))
        
        for sim in range(0, sims):
            utils._progress('performing simultion {} out of {}...'.format(sim + 1, sims))
            #utils._progress('processing chunk level {}'.format(c))
            self.status = 0
            
            i_dp = self.interpolation(dem_mod)
            
            if dp is None:
                dp = i_dp
            else:
                try:
                    dp = np.concatenate((dp, i_dp), axis = 0)
                except: dp = None
                
            if len(dp) == 0: self.status = -1
            if self.verbose: self.tw.msg('gathered {} error points'.format(len(dp)))
            self.tw.end(self.status, 'performed simulation {} out of {}.'.format(sim + 1, sims))

        self.tw.end(self.status, 'ran \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m'.format(dem_mod))
        
        if self.status == 0:

            if self.verbose: self.tw.msg('gathered {} error points'.format(len(dp)))
            np.savetxt('test.err', dp, '%f', ' ')
            ec = self.err_plot(dp)
            
            ## ==============================================
            ## apply error coefficient to full proximity grid
            ## ==============================================

            utils._progress('applying coefficient to proximity grid')
            ## USE numpy instead
            math_cmd = 'gmt grdmath {} ABS {} POW {} MUL {} ADD = {}_dst_unc.tif=gd+n-9999:GTiff\
            '.format(self.dem['prox'], ec[2], ec[1], 0, self.o_name)
            utils.run_cmd(math_cmd, self.verbose, self.verbose)
            self.tw.end(0, 'applying coefficient to proximity grid')

    def set_or_make_dem(self, dem_mod = 'mbgrid'):
        '''check if dem dict contains dems, otherwise generate them...'''

        utils._progress('checking for DEMs...')

        if self.dem['dem'] is None:
            utils._progress('generating DEM using {}...'.format(dem_mod))
            self.dem['dem'] = dem(self.datalist, self.region, str(self.inc)).run(dem_mod)
            self.tw.end(self.status, 'generated DEM using {}.'.format(dem_mod))
        self.tw.msg('using DEM {}'.format(self.dem['dem']))

        # if self.dem['num'] is None:
        #     utils._progress('generating NUM grid...')
        #     self.dem['num'] = dem(self.datalist, self.region, str(self.inc)).run('num')
        #     tw.end(self.status, 'generated NUM grid.')
        # tw.msg('using NUM grid {}'.format(self.dem['num']))
            
        if self.dem['msk'] is None:
            utils._progress('generating MASK grid...')
            self.dem['msk'] = self.datalist.mask(region = self.proc_region.region, inc = self.inc, o_name = self.o_name)
            #msk_dem = dem(self.datalist, self.region, str(self.inc))
            #msk_dem.o_fmt = 'GTiff'
            #self.dem['msk'] = msk_dem.run('mask')
            self.tw.end(self.status, 'generated MASK grid.')
        self.tw.msg('using MASK {}'.format(self.dem['msk']))
            
        if self.dem['prox'] is None:
            utils._progress('generating proximity grid...')
            self.dem['prox']  = '{}_prox.tif'.format(self.dem['msk'].split('.')[0]) 
            proximity(self.dem['msk'], self.dem['prox'])
            self.tw.end(self.status, 'generated proximity grid.')
        self.tw.msg('using PROXIMITY grid {}'.format(self.dem['prox']))
            
        self.tw.end(self.status, 'checked for DEMs.')

    def split_sample(self, sub_regions, ss_samp, dem_mod = 'mbgrid', s_dp = None):
        
        if self.stop() or self.status !=0:
            return(s_dp)
        
        for n,sub_region in enumerate(sub_regions):
            if self.verbose: self.tw.msg('processing sub-region {}'.format(sub_region))
            
            this_region = regions.region('/'.join(map(str, sub_region[0])))
            o_xyz = '{}_{}.xyz'.format(self.o_name, n)

            if self.verbose:
                self.tw.msg('initial sampling density: {}'.format(sub_region[3]))
                self.tw.msg('desired sampling density: {}'.format(ss_samp))

            if sub_region[3] < ss_samp:
                ss_samp = None
            
            with open(o_xyz, 'w') as o_fh:
                gdalfun.dump(self.dem['dem'], o_fh, False, gdalfun._srcwin(self.dem['dem'], this_region.buffer(10*self.inc).region), self.dem['msk'])

            if os.stat(o_xyz).st_size == 0:
                self.tw.err_msg('no data in sub-region...')
                self.status = -1
            else:
                ## ==============================================
                ## split sub region XYZ data into outer and
                ## inner (outer is a data buffer)
                ## ==============================================

                s_inner, s_outer = gmtselect_split(o_xyz, this_region, 'sub_{}'.format(n))

                if os.stat(s_inner).st_size != 0:
                    sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                else: sub_xyz = []
                ss_len = len(sub_xyz)

                ## Skip if ss_len is smaller than 5 or so
                #if os.stat(s_outer).st_size != 0:
                #    sub_o_xyz = np.loadtxt(s_outer, ndmin=2, delimiter = ' ')
                #print(len(sub_o_xyz))
                
                if ss_samp is not None:
                    #sx_cnt = int(ss_len * (ss_samp / 100)) + 1
                    sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
                else: sx_cnt = 1
                sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                if self.verbose:
                    self.tw.msg('withholding {} out of {} points for error sampling'.format(ss_len - sx_cnt, ss_len))

                np.random.shuffle(sub_xyz)
                np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                sub_datalist =  datalists.datalist('sub_{}.datalist'.format(n), this_region)
                sub_datalist._append_datafile([s_outer, 168, 1])
                sub_datalist._append_datafile([sub_xyz_head, 168, 1])
                sub_datalist._load_data()

                sub_surf = dem(sub_datalist, this_region, str(self.inc), verbose = self.verbose)
                #sub_surf.o_fmt = 'GTiff'
                sub_dem = sub_surf.run(dem_mod)

                sub_msk = sub_datalist.mask(region = this_region.buffer(10*self.inc).region, inc = self.inc)
                sub_prox = '{}_prox.tif'.format(sub_msk.split('.')[0])
                self.status = proximity(sub_msk, sub_prox)

                sub_xyd = gdalfun.query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                sub_dp = gdalfun.query(sub_xyd, sub_prox, 'zg')

                sys.exit()
                os.remove(sub_xyz_head)
                os.remove(sub_xyz_head + '.inf')
                os.remove(sub_datalist._path)

                if s_dp is None:
                    s_dp = sub_dp
                else:
                    try:
                        s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                    except:
                        self.tw.err_msg(ss_len)
                        self.tw.err_msg(sx_cnt)
                        self.tw.err_msg(s_dp)
                        self.tw.err_msg(sub_dp)
                        s_dp = None

            os.remove(o_xyz)
            utils.remove_glob('sub_{}*'.format(n))
            
        return(s_dp)
    
    def tile_analysis(self):

        sub_count = 0
        sub_zones = {}
        
        for sub_region in self.sub_regions:
            sub_count += 1

            ## cut mask grid to sub region
            gdalfun.cut(self.dem['msk'], gdalfun._srcwin(self.dem['msk'], sub_region.region), 'tmp_msk.tif')
            gdalfun.cut(self.dem['dem'], gdalfun._srcwin(self.dem['dem'], sub_region.region), 'tmp_dem.tif')

            s_gc = gdalfun._infos('tmp_msk.tif')
            s_g_max = float(s_gc['nx'] * s_gc['ny'])
            s_sum = gdalfun.sum('tmp_msk.tif')
            s_perc = (s_sum / s_g_max) * 100

            os.remove('tmp_msk.tif')
            
            s_dc = gdalfun._infos('tmp_dem.tif', True)
            s_min = s_dc['zmin']
            s_max = s_dc['zmax']

            os.remove('tmp_dem.tif')
            
            if s_max < 0:
                zone = 'Bathy'
            elif s_min > 0:
                zone = 'Topo'
            else: zone = 'BathyTopo'

            sub_zones[sub_count] = [sub_region.region, s_g_max, s_sum, s_perc, s_min, s_max, zone]
            
        return(sub_zones)

    def region_analysis(self):

        region_info = {}
        
        num_sum = gdalfun.sum(self.dem['msk'])
        gc = gdalfun._infos(self.dem['msk'])
        g_max = float(gc['nx'] * gc['ny'])
        num_perc = (num_sum / g_max) * 100.
        prox_perc_95 = gdalfun.percentile(self.dem['prox'], 95)

        region_info[self.o_name] = [self.region.region, g_max, num_sum, num_perc, prox_perc_95]

        return(region_info)

    def err_plot(self, dp):
        #try:
        dp = dp[dp[:,1]<self.region_info[self.o_name][4],:]
        dp = dp[dp[:,1]>0,:]
        ec = err2coeff(dp)
        self.tw.msg('error coefficient: {}'.format(ec))
        #except:
        #    self.tw.err_msg('cannot allocate memory')

        return(ec)
    
    def interpolation(self, dem_mod = 'mbgrid'):
        '''calculate the interpolation uncertainty.'''

        this_region = self.region
        dp = None

        s_dens = np.array([self.sub_zones[x][3] for x in self.sub_zones.keys()])
        s_5perc = np.percentile(s_dens, 5)
        self.tw.msg('Sampling density for split-sample is: {}'.format(s_5perc))

        ## ==============================================
        ## BATHY
        ## ==============================================
        
        bathy_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Bathy']
        if len(bathy_tiles) > 0:
            bathy_dens = np.array([x[3] for x in bathy_tiles])
            bathy_50perc = np.percentile(bathy_dens, 50)
        else: bathy_50perc = 0.0
        if self.verbose: self.tw.msg('Minimum sampling for Bathy tiles: {}'.format(bathy_50perc))

        bathy_trainers = []
        for i in bathy_tiles:
            if i[3] >= bathy_50perc:
                bathy_trainers.append(i)
                
        bathy_trainers = np.array(bathy_trainers)
                
        #if self.verbose:
        self.tw.msg('possible BATHY training zones: {}'.format(len(bathy_trainers)))

        np.random.shuffle(bathy_trainers)
        #bathy_trainers_h = bathy_trainers[:int(len(bathy_trainers)*.25) + 1]
        bathy_trainers_h = bathy_trainers[:25]

        bathy_t_dens = np.array([x[3] for x in bathy_trainers_h])
        if len(bathy_t_dens) > 0:
            bathy_t_50perc = np.percentile(bathy_t_dens, 50)
        else: bathy_t_50perc = 0

        # if bathy_t_50perc > 0.0:
        #     n_loops = int(((100 - bathy_t_50perc) / self.region_info[self.o_name][3]) + 1)
        # else: n_loops = 0

        utils._progress('calculating split-sample for {} semi-random \033[1mBATHY\033[m zone training tiles'.format(len(bathy_trainers_h)))
        bathy_dp = self.split_sample(bathy_trainers_h, s_5perc, dem_mod, dp)
        if bathy_dp is not None:
            status = 0
        else: status = -1

        self.tw.end(status, 'calculated split-sample for {} semi-random \033[1mBATHY\033[m zone training tiles'.format(len(bathy_trainers_h)))

        ## ==============================================
        ## BATHY TOPO
        ## ==============================================
        
        bathy_topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'BathyTopo']
        if len(bathy_topo_tiles) > 0:
            bathy_topo_dens = np.array([x[3] for x in bathy_topo_tiles])
            bathy_topo_50perc = np.percentile(bathy_topo_dens, 50)
        else: bathy_topo_50perc = 0.0
        if self.verbose: self.tw.msg('Minimum sampling for BathyTopo tiles: {}'.format(bathy_topo_50perc))
        
        bathy_topo_trainers = []
        for i in bathy_topo_tiles:
            if i[3] >= bathy_topo_50perc:
                bathy_topo_trainers.append(i)
                
        bathy_topo_trainers = np.array(bathy_topo_trainers)
                
        #if self.verbose:
        self.tw.msg('possible BATHY-TOPO training zones: {}'.format(len(bathy_topo_trainers)))

        np.random.shuffle(bathy_topo_trainers)
        #bathy_topo_trainers_h = bathy_topo_trainers[:int(len(bathy_topo_trainers)*.25) + 1]
        bathy_topo_trainers_h = bathy_topo_trainers[:25]

        bathy_topo_t_dens = np.array([x[3] for x in bathy_topo_trainers_h])
        if len(bathy_topo_t_dens) > 0:
            bathy_topo_t_50perc = np.percentile(bathy_topo_t_dens, 50)
        else: bathy_topo_t_50perc = 0

        # if bathy_topo_t_50perc > 0.0:
        #     n_loops = int(((100 - bathy_topo_t_50perc) / self.region_info[self.o_name][3]) + 1)
        # else: n_loops = 0

        utils._progress('calculating split-sample for {} semi-random \033[1mBATHY-TOPO\033[m zone training tiles'.format(len(bathy_topo_trainers_h)))
        bathy_topo_dp = self.split_sample(bathy_topo_trainers_h, s_5perc, dem_mod, bathy_dp)
        if bathy_topo_dp is not None:
            status = 0
        else: status = -1

        self.tw.end(status, 'calculated split-sample for {} random \033[1mBATHY-TOPO\033[m zone training tiles'.format(len(bathy_topo_trainers_h)))

        ## ==============================================
        ## TOPO
        ## ==============================================
        
        topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Topo']
        if len(topo_tiles) > 0:
            topo_dens = np.array([x[3] for x in topo_tiles])
            topo_50perc = np.percentile(topo_dens, 50)
        else: topo_50perc = 0.0
        if self.verbose: self.tw.msg('Minimum sampling for Topo tiles: {}'.format(topo_50perc))

        topo_trainers = []
        for i in topo_tiles:
            if i[3] >= topo_50perc:
                topo_trainers.append(i)

        topo_trainers = np.array(topo_trainers)
        
        #if self.verbose:
        self.tw.msg('possible TOPO training zones: {}'.format(len(topo_trainers)))

        np.random.shuffle(topo_trainers)
        #topo_trainers_h = topo_trainers[:int(len(topo_trainers)*.25) + 1]
        topo_trainers_h = topo_trainers[:25]

        topo_t_dens = np.array([x[3] for x in topo_trainers_h])
        if len(topo_t_dens) > 0:
            topo_t_50perc = np.percentile(topo_t_dens, 50)
        else: topo_t_50perc = 0

        # if topo_t_50perc > 0.0:
        #     n_loops = int(((100 - topo_t_50perc) / self.region_info[self.o_name][3]) + 1)
        # else: n_loops = 0

        utils._progress('calculating split-sample for {} semi-random \033[1mTOPO\033[m zone training tiles'.format(len(topo_trainers_h)))
        dp = self.split_sample(topo_trainers_h, s_5perc, dem_mod, bathy_topo_dp)
        if dp is not None:
            status = 0
        else: status = -1

        self.tw.end(status, 'calculated split-sample for {} semi-random \033[1mTOPO\033[m zone training tiles'.format(len(topo_trainers_h)))

        return(dp)            

### END
