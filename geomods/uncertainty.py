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
## uncertainty module:
##
## =============================================================================

def err2coeff(my_data):
    '''data is 2 col file with `err dist`'''

    from scipy import optimize
    import matplotlib.pyplot as plt
    #try: 
    #    my_data = np.loadtxt(data, delimiter=' ')
    #except: sys.exit(2)

    error=my_data[:,0]
    distance=my_data[:,1]
        
    max_int_dist = np.max(distance)
    nbins = 10

    coeff_guess=[0, 0.1, 0.2]
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

    fitfunc = lambda p, x: p[0] + p[1] * (abs(x) ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)

    # fit
    
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    #plt.show()

    out_png = 'unc_best_fit.png'
    plt.savefig(out_png)   # save the figure to file
    plt.close()

    #scatter

    plt.scatter(distance, error)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')

    out_png = 'unc_scatter.png'
    plt.savefig(out_png)
    plt.close()

    return(out)

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

        self.dem = { 
            'dem': None,
            'num': None,
            'num-grd': None,
            'num-msk': None,
            'prox': None,
            'int-unc': None,
        }

        if o_name is None:
            self.o_name = self.datalist._name
        else: self.o_name = o_name

    def run(self, dem_mod = 'mbgrid', dem = None, num = None, num_msk = None, prox = None):

        if dem is not None:
            self.dem['dem'] = dem

        if num is not None:
            self.dem['num'] = num

        if num_msk is not None:
            self.dem['num-msk'] = num_msk

        if prox is not None:
            self.dem['prox'] = prox

        self.datalist._load_data()
        self.interpolation_uncertainty(dem_mod)

    def set_or_make_dem(self, dem_mod = 'mbgrid'):
        '''check if dem dict contains dems, otherwise generate them...'''

        tw = utils._progress('checking for DEMs...')

        if self.dem['dem'] is None:
            tw.err_msg('generating dem...')
            dems = dem(self.datalist, self.region, str(self.inc)).run(dem_mod)
            for key in dems.keys():
                self.dem[key] = dems[key]

        if self.dem['num'] is None:
            tw.err_msg('generating NUM grid...')
            dems = dem(self.datalist, self.region, str(self.inc)).run('num')
            for key in dems.keys():
                self.dem[key] = dems[key]

        if self.dem['num-msk'] is None:
            tw.err_msg('generating NUM mask...')
            self.dem['num-grd-msk'] = '{}_msk.grd'.format(self.dem['num'].split('.')[0]) 
            num_msk(self.dem['num'], self.dem['num-grd-msk'])
            self.dem['num-msk'] = grd2tif(self.dem['num-grd-msk'])

        if self.dem['prox'] is None:
            tw.err_msg('generating proximity grid...')
            self.dem['prox']  = '{}_prox.tif'.format(self.dem['num'].split('.')[0]) 
            proximity(self.dem['num-msk'], self.dem['prox'])

        tw.opm = 'checked for DEMs.'
        tw.end(self.status)        

    def gmtselect_split(self, o_xyz, sub_region, sub_bn):
        '''split an xyz file into an inner and outer region.'''

        out_inner = None
        out_outer = None

        gmt_s_inner = 'gmt gmtselect -V {} {} > {}_inner.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
        out, self.status = utils.run_cmd(gmt_s_inner, self.verbose, False)
        
        if self.status == 0:
            out_inner = '{}_inner.xyz'.format(sub_bn)

        gmt_s_outer = 'gmt gmtselect -V {} {} -Ir > {}_outer.xyz'.format(o_xyz, sub_region.gmt, sub_bn)
        out, self.status = utils.run_cmd(gmt_s_outer, self.verbose, False)

        if self.status == 0:
            out_outer = '{}_outer.xyz'.format(sub_bn)

        return([out_inner, out_outer])

    def interpolation_uncertainty(self, dem_mod = 'mbgrid'):
        '''calculate the interpolation uncertainty.'''

        loop_count = 1
        sub_count = 0        
        this_region = self.region
        dp = None

        self.set_or_make_dem(dem_mod)
        print self.dem

        ## ==============================================
        ## Calculate the percentage of filled cells
        ## and proximity percentiles.
        ## ==============================================
        tw = utils._progress('')

        num_sum = gdalfun.sum(self.dem['num-msk'])
        gi = grdinfo(self.dem['num-msk'])
        g_max = int(gi[9]) * int(gi[10])

        num_perc = (num_sum / g_max) * 100
        prox_perc_95 = gdalfun.percentile(self.dem['prox'], 75)
        tw.err_msg('waffles: proximity 95th perc: {}'.format(prox_perc_95))

        sub_regions = this_region.chunk(self.inc, 500)
        tw.err_msg('waffles: chunking into {} regions.'.format(len(sub_regions)))
        #tw.err_msg('waffles: processing {} regions {} times.'.format(len(sub_regions), loop_count))

        for sub_region in sub_regions:
            if self.stop(): break

            self.status = 0
            sub_count += 1
            o_xyz = '{}.xyz'.format(self.o_name)

            tw = utils._progress('processing sub region \033[1m{}\033[m...'.format(sub_count))

            self.status = grd2xyz(self.dem['dem'], o_xyz, region = sub_region.buffer(10*self.inc), mask = self.dem['num-msk'])
            if os.stat(o_xyz).st_size == 0:
                tw.err_msg('waffles: error, no data in sub-region...')
                self.status = -1
            else:
                s_inner, s_outer = self.gmtselect_split(o_xyz, sub_region, 'sub_{}'.format(sub_count))

                if os.stat(s_inner).st_size != 0:
                    sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                else: self.status = -1

                if self.status == 0 and not self.stop():
                    self.status = grdcut(self.dem['num-msk'], sub_region, 'tmp_sub.grd')

                    gi = grdinfo('tmp_sub.grd')
                    num_max = int(gi[9]) * int(gi[10])

                    tw.err_msg('waffles: total cells in subgrid is {}'.format(num_max))
                    num_sub_sum = gdalfun.sum(grd2tif('tmp_sub.grd'))
                    tw.err_msg('waffles: total filled cells in subgrid is {}'.format(num_sub_sum))

                    try:
                        os.remove('tmp_sub.grd')
                        os.remove('tmp_sub.tif')
                    except: pass

                    num_sub_perc = (num_sub_sum / num_max) * 100
                    tw.err_msg('waffles: {}% of cells have data'.format(num_sub_perc))
                    s_size = 100 - num_sub_perc

                    if s_size >= 100:
                        n_loops = 0
                    else: n_loops = int(((s_size / num_perc) * 2) + 1)

                    tw.err_msg('waffles: loops for this subregion is {}'.format(n_loops))

                    ## ==============================================
                    ## Split Sample n_loops times
                    ## ==============================================

                    for i in range(0, n_loops):
                        if self.stop(): break

                        np.random.shuffle(sub_xyz)
                        sx_len = len(sub_xyz)
                        sx_len_pct = int(sx_len * (num_sub_perc / 100))

                        if sx_len_pct == 0:
                            break

                        tw.err_msg('waffles: extracting {} random data points out of {}'.format(sx_len_pct, sx_len))
                        sub_xyz_head = 'sub_{}_head.xyz'.format(sub_count)

                        np.savetxt(sub_xyz_head, sub_xyz[:sx_len_pct], '%f', ' ')

                        sub_datalist =  datalists.datalist('sub_{}.datalist'.format(sub_count), sub_region)
                        sub_datalist._append_datafile(s_outer, 168, 1)
                        sub_datalist._append_datafile(sub_xyz_head, 168, 1)
                        sub_datalist._reset()

                        pb = utils._progress('generating sub-region {} random-sample grid...{}/{}'.format(sub_count, i, n_loops))

                        sub_surf = dem(sub_datalist, sub_region, str(self.inc))
                        sub_surf.run(dem_mod)
                        sub_surf.run('num')
                        sub_dems = sub_surf.dem

                        pb.opm = 'generated sub-region {} random-sample grid.{}/{}'.format(sub_count, i, n_loops)
                        pb.end(sub_surf.status)

                        if sub_dems['dem'] is not None:

                            sub_prox = '{}_prox.tif'.format(sub_dems['num'].split(',')[0])
                            self.status = proximity(sub_dems['num-msk'], sub_prox)
                            
                            sub_xyd = gdalfun.query(sub_xyz[sx_len_pct:], sub_dems['dem'], 'xyd')
                            sub_dp = gdalfun.query(sub_xyd, sub_prox, 'zg')

                            if len(sub_dp) != 0:
                                if dp is None:
                                    dp = sub_dp
                                else: dp = np.concatenate((dp, sub_dp), axis = 0)
                        os.remove(sub_xyz_head)
                        os.remove(sub_datalist._path)

                utils.remove_glob('sub_{}*'.format(sub_count))

            tw.opm = 'processed sub region \033[1m{}\033[m.'.format(sub_count)
            tw.end(self.status)

        ## ==============================================
        ## calculate the error coefficients and plot results
        ## ==============================================

        if not self.stop():
            dp = dp[dp[:,1]<prox_perc_95,:]
            dp = dp[dp[:,1]>0,:]
            ec = err2coeff(dp)
            print('error coefficient: {}'.format(ec))
            np.savetxt('test.err', dp, '%f', ' ')

        ## ==============================================
        ## apply error coefficient to full proximity grid
        ## ==============================================
