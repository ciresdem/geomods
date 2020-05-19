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

## =============================================================================
##
## uncertainty module - uncertainties.py
##
## datasource and interpolation uncertainty.
##
## =============================================================================

def err2coeff(my_data, coeff_guess = [0, 0.1, 0.2], dst_name = None):
    '''data is 2 col file with `err dist`'''

    from scipy import optimize
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText
    except:
        print('you need to install matplotlib to run uncertainty plots...')
    
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
    else: out_png = '{}_bf.png'.format(dst_name)
    plt.savefig(out_png)   # save the figure to file
    plt.close()

    #scatter

    plt.scatter(distance, error)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')

    if dst_name is None:
        out_png = 'unc_scatter.png'
    else: out_png = '{}_scatter.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

    return(out)

class uncertainty:

    def __init__(self, i_datalist, i_region, i_inc = 0.0000925925, o_name = None, o_node = 'pixel', o_extend = 6, callback = lambda: False, verbose = False):

        self.datalist = i_datalist
        self.region = i_region
        self.inc = float(i_inc)
        self.node = o_node

        self.extend = int(o_extend)

        self.proc_region = self.region.buffer((self.extend * 2) * self.inc)
        self.dist_region = self.region.buffer(self.extend * self.inc)
        
        self.status = 0
        self.stop = callback
        self.verbose = verbose

        self.dem_mod = 'mbgrid'
        
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

        self.zones = ['bathy', 'bathy-topo', 'topo']
        self.region_info = {}
        self.sub_zones = {}
        self.trainers = []
        
    def run(self, dem_mod = 'mbgrid', dem = None, msk = None):
        
        i_dp = None
        opts = dem_mod.split(':')
        self.dem_mod = opts[0]
        self.mod_args = list(opts[1:])
        
        if dem is not None: self.dem['dem'] = dem
        if msk is not None: self.dem['msk'] = msk
                
        ## s_dp = self.source()
        ## v_dp = self.vdatum()
        i_dp = self.interpolation()

        ## self.combine(s_dp, v_dp, i_dp)
        
        return(i_dp)
        
    def set_or_make_dem(self):
        '''check if dem dict contains dems, otherwise generate them...'''

        if self.dem['dem'] is None:
            this_dem = waffles(self.datalist, self.region, str(self.inc), o_b_name = self.o_name)
            this_dem.o_fmt = 'GTiff'
            self.dem['dem'] = this_dem.run(self.dem_mod, self.mod_args)
            echo_msg('generated DEM {} using {}.'.format(self.dem['dem'], self.dem_mod))
        else: echo_msg('found DEM {}'.format(self.dem['dem']))

        if self.dem['msk'] is None:
            self.dem['msk'] = self.datalist.mask(region = self.dist_region.region, inc = self.inc, o_name = self.o_name)
            msk_dem = waffles(self.datalist, self.region, str(self.inc))
            msk_dem.o_fmt = 'GTiff'
            self.dem['msk'] = msk_dem.run('mask')
        else: echo_msg('found Data MASK {}'.format(self.dem['msk']))

        if self.dem['prox'] is None:
            self.dem['prox']  = '{}_prox.tif'.format(self.dem['msk'].split('.')[0]) 
            gdalfun.gdal_proximity(self.dem['msk'], self.dem['prox'])
            echo_msg('generated PROXIMITY grid {}.'.format(self.dem['prox']))
        else: echo_msg('found PROXIMITY grid {}.'.format(self.dem['prox']))
            
    def err_plot(self, dp, d_max):
        '''plot a numpy array of 'err dist' values and return the error coefficient.'''
        
        #dp = dp[dp[:,1]<self.region_info[self.o_name][4],:]
        dp = dp[dp[:,1] < d_max,:]
        dp = dp[dp[:,1] > 0,:]
        ec = err2coeff(dp, dst_name = self.o_name)
        echo_msg('error coefficient: {}'.format(ec))

        return(ec)

    def region_analysis(self):
        '''Analyze the input self.region and return infos.'''
        
        region_info = {}
        
        num_sum = gdalfun.gdal_sum(self.dem['msk'])
        gc = gdalfun._infos(self.dem['msk'])
        g_max = float(gc['nx'] * gc['ny'])
        num_perc = (num_sum / g_max) * 100.
        prox_perc_95 = gdalfun.gdal_percentile(self.dem['prox'], 95)

        region_info[self.o_name] = [self.region.region, g_max, num_sum, num_perc, prox_perc_95]

        return(region_info)
    
    def tile_analysis(self):
        '''Anaylize the chunked regions and return infos about them.'''
        
        sub_count = 0
        sub_zones = {}
        
        for sc, sub_region in enumerate(self.sub_regions):

            gdalfun.gdal_cut(self.dem['msk'], gdalfun._srcwin(self.dem['msk'], sub_region.region), 'tmp_msk.tif')
            gdalfun.gdal_cut(self.dem['dem'], gdalfun._srcwin(self.dem['dem'], sub_region.region), 'tmp_dem.tif')

            s_gc = gdalfun._infos('tmp_msk.tif')
            s_g_max = float(s_gc['nx'] * s_gc['ny'])
            s_sum = gdalfun.gdal_sum('tmp_msk.tif')
            s_perc = (s_sum / s_g_max) * 100
            
            s_dc = gdalfun._infos('tmp_dem.tif', True)
            
            if s_dc['zmax'] < 0:
                zone = 'Bathy'
            elif s_dc['zmin'] > 0:
                zone = 'Topo'
            else: zone = 'BathyTopo'

            sub_zones[sc+1] = [sub_region.region, s_g_max, s_sum, s_perc, s_dc['zmin'], s_dc['zmax'], zone]

            remove_glob('tmp_msk.tif')
            remove_glob('tmp_dem.tif')
            
        return(sub_zones)

    def zone_analysis(self):
        '''Analyze uncertainty zones and select training tiles.'''
        
        trainers = []
        
        bathy_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Bathy']
        bathy_topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'BathyTopo']
        topo_tiles = [self.sub_zones[x] for x in self.sub_zones.keys() if self.sub_zones[x][6] == 'Topo']

        for z, tile_set in enumerate([bathy_tiles, bathy_topo_tiles, topo_tiles]):
            if len(tile_set) > 0:
                t_dens = np.array([x[3] for x in tile_set])
                t_50perc = np.percentile(t_dens, 50)
            else: t_50perc = 0.0
            if self.verbose: echo_msg('Minimum sampling for {} tiles: {}'.format(self.zones[z].upper(), t_50perc))

            t_trainers = [x for x in tile_set if x[3] > t_50perc]
            echo_msg('possible {} training zones: {}'.format(self.zones[z].upper(), len(t_trainers)))
            trainers.append(t_trainers)
                
        return(trainers)

    def zone_sort(self):
        '''sort training tiles by distance'''

        train_sorted = []
        
        for z, train in enumerate(self.trainers):
            train_d = []
            
            np.random.shuffle(train)

            while True:
                if len(train) == 0: break
                
                this_center = region(train[0][0]).center()
                train_d.append(train[0])
                train = train[1:]
                if len(train) == 0: break
                
                dsts = [hav_dst(this_center, region(x[0]).center()) for x in train]
                min_dst = np.percentile(dsts, 50)
                d_t = lambda t: hav_dst(this_center, region(t[0]).center()) > min_dst
                
                np.random.shuffle(train)
                train.sort(reverse=True, key=d_t)
                
            if self.verbose: echo_msg(' '.join([region(x[0]).gmt for x in train_d[:25]]))
            train_sorted.append(train_d)
            
        return(train_sorted)
        
    def split_sample(self, sub_regions, ss_samp, s_dp = None):
        '''perform split-sample analysis on the training tiles and return a list of `error distance` values'''
        
        if self.stop() or self.status !=0: return(s_dp)
        
        for n,sub_region in enumerate(sub_regions):
            #if self.verbose:
            echo_msg('processing sub-region ({}) {}'.format(n, sub_region))
            
            this_region = region(sub_region[0])
            o_xyz = '{}_{}.xyz'.format(self.o_name, n)

            if self.verbose:
                echo_msg('initial sampling density: {}'.format(sub_region[3]))
                echo_msg('desired sampling density: {}'.format(ss_samp))

            if sub_region[3] < ss_samp: ss_samp = None
            
            with open(o_xyz, 'w') as o_fh:
                gdalfun.gdal_dump(self.dem['dem'], o_fh, ' ', None, False, gdalfun._srcwin(self.dem['dem'], this_region.buffer(20*self.inc).region), self.dem['msk'])

            if os.stat(o_xyz).st_size == 0:
                echo_msg('no data in sub-region...')
                self.status = -1
            else:
                ## use gdal instead
                s_inner, s_outer = gmtselect_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = self.verbose)

                if os.stat(s_inner).st_size != 0:
                    sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                else: sub_xyz = []
                ss_len = len(sub_xyz)

                if ss_samp is not None:
                    sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
                else: sx_cnt = 1
                sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                if self.verbose: echo_msg('withholding {} out of {} points for error sampling'.format(ss_len - sx_cnt, ss_len))

                np.random.shuffle(sub_xyz)
                np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                sub_datalist = datalist('sub_{}.datalist'.format(n), this_region, verbose = self.verbose)
                sub_datalist._append_entry([s_outer, 168, 1])
                sub_datalist._append_entry([sub_xyz_head, 168, 1])

                sub_surf = waffles(sub_datalist, this_region, i_inc = str(self.inc), o_name = 'sub_{}'.format(n), verbose = self.verbose)
                sub_dem = sub_surf.run(self.dem_mod, self.mod_args)

                if sub_dem is not None:
                    sub_msk = sub_datalist.mask(region = this_region.buffer(10*self.inc).region, inc = self.inc)
                    sub_prox = '{}_prox.tif'.format(sub_msk.split('.')[0])
                    gdalfun.gdal_proximity(sub_msk, sub_prox)
                    gdalfun.gdal_infos(sub_prox)
                    sub_xyd = gdalfun.gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                    sub_dp = gdalfun.gdal_query(sub_xyd, sub_prox, 'zg')
                else: sub_dp = None

                remove_glob(sub_xyz_head)
                remove_glob(sub_datalist._path)
                sub_xyz = None

                if s_dp is None:
                    s_dp = sub_dp
                else:
                    try:
                        s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                    except:
                        if self.verbose: echo_error_msg('found no error points...')
                        pass

            remove_glob(o_xyz)
            remove_glob('sub_{}*'.format(n))
            
        return(s_dp)
            
    def interpolation(self):
        '''calculate the interpolation uncertainty.'''

        dp = None
        sims = 10
        sim_loops = 1
        
        self.set_or_make_dem()
        if self.verbose: echo_msg(self.dem)

        self.region_info = self.region_analysis()
        if self.verbose:
            for x in self.region_info.keys():
                echo_msg('region: {}: {}'.format(x, self.region_info[x]))

        ## ==============================================
        ## chunk region into sub regions
        ## ==============================================
        chnk_lvl = 4        
        echo_msg('chunking region into sub-regions using chunk level {}...'.format(chnk_lvl))
        chnk_inc = int(chnk_lvl * self.region_info[self.o_name][4])
        self.sub_regions = self.region.chunk(self.inc, chnk_inc)
        #if self.verbose: print([x.region for x in self.sub_regions])
        echo_msg('chunked region into {} sub-regions.'.format(len(self.sub_regions)))

        echo_msg('analyzing {} sub-regions...'.format(len(self.sub_regions)))
        self.sub_zones = self.tile_analysis()
        if self.verbose:
            for x in self.sub_zones.keys():
                echo_msg('Sub-region {}: {}'.format(x, self.sub_zones[x]))
        
        echo_msg('running \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m...'.format(self.dem_mod))
                    
        s_dens = np.array([self.sub_zones[x][3] for x in self.sub_zones.keys()])
        s_5perc = np.percentile(s_dens, 5)
        s_dens = None
        echo_msg('Sampling density for region is: {:.12f}'.format(s_5perc))
                
        self.trainers = self.zone_analysis()
        trains = self.zone_sort()
        echo_msg('sorted training tiles.')
        ## generate shapefile from tain_h
        echo_msg('analyzed {} sub-regions.'.format(len(self.sub_regions)))
        
        for sim in range(0, sims):
            echo_msg('performing INTERPOLATION UNCERTAINTY simulation {} out of {}...'.format(sim + 1, sims))
            self.status = 0
            
            for s_l in range(0, sim_loops):
                for z, train in enumerate(trains):
                    train_h = train[:25]
                    dp = self.split_sample(train_h, s_5perc, dp)
            if len(dp) == 0:
                self.status = -1
                #break
            echo_msg('performed INTERPOLATION UNCERTAINTY simulation {} out of {}; {} error points accumulated.'.format(sim + 1, sims, len(dp)))

        echo_msg('ran \033[1mINTERPOLATION\033[m uncertainty module using \033[1m{}\033[m.'.format(self.dem_mod))
        
        if self.status == 0:

            if self.verbose: echo_msg('gathered {} error points'.format(len(dp)))
            np.savetxt('{}.err'.format(self.o_name), dp, '%f', ' ')
            ec = self.err_plot(dp[:50000000], self.region_info[self.o_name][4])
            
            ## ==============================================
            ## apply error coefficient to full proximity grid
            ## ==============================================

            echo_msg('applying coefficient to proximity grid')
            ## USE numpy instead
            
            math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_dst_unc.tif=gd+n-9999:GTiff\
            '.format(self.dem['prox'], ec[2], ec[1], 0, self.o_name)
            utils.run_cmd(math_cmd, self.verbose, self.verbose)
            echo_msg('applyed coefficient to proximity grid')
            
        return(dp)

    def source(self):
        '''source data uncertainty'''
        
        pass
    
### END
