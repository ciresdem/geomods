
import geomods
from geomods.waffles import *
from scipy import optimize

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText
except:
    print('you need to install matplotlib to run uncertainty plots...')
    
## =============================================================================
##
## waffles uncertainty module - uncertainties.py
##
## datasource and interpolation uncertainty.
##
## =============================================================================
_unc_config = {
    'wg': waffles_config(),
    'dem': None,
    'msk': None,
    'prox': None,
    'zones': ['bathy', 'bathy-topo', 'topo'],
}

def err_fit_plot(xdata, ydata, out, fitfunc, dst_name = 'unc'):
    '''plot a best fit plot'''
    plt.plot(xdata, ydata, 'o')
    plt.plot(xdata, fitfunc(out, xdata), '-')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    out_png = '{}_bf.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

def err_scatter_plot(error_arr, dist_arr, dst_name = 'unc'):
    '''plot a scatter plot'''
    plt.scatter(dist_arr, error_arr)
    #plt.title('Scatter')
    plt.xlabel('distance')
    plt.ylabel('error (m)')
    out_png = '{}_scatter.png'.format(dst_name)
    plt.savefig(out_png)
    plt.close()

def err2coeff(err_arr, coeff_guess = [0, 0.1, 0.2], dst_name = 'unc'):
    '''calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist`'''
    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    max_int_dist = np.max(distance)
    nbins = 10
    n, _ = np.histogram(distance, bins = nbins)
    # want at least 2 values in each bin?
    while 0 or 1 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins = nbins)

    serror, _ = np.histogram(distance, bins = nbins, weights = error)
    serror2, _ = np.histogram(distance, bins = nbins, weights = error**2)
    mean2 = (serror / n)**2
    std = np.sqrt(serror2 / n - mean2)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)
    fitfunc = lambda p, x: p[0] + p[1] * (abs(x) ** abs(p[2]))
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args = (xdata, ydata), full_output = True)
    err_fit_plot(xdata, ydata, out, fitfunc, dst_name)
    err_scatter_plot(error, distance)
    return(out)

def err_plot(err_arr, d_max):
    '''plot a numpy array of 'err dist' values and return the error coefficient.'''    
    err_arr = err_arr[err_arr[:,1] < d_max,:]
    err_arr = err_arr[err_arr[:,1] > 0,:]
    return(err2coeff(err_arr))

# def unc_split_sample(uc = _unc_config, ss_samp = 50, s_dp = None):
#     '''perform split-sample analysis on the training tiles and return a list of `error distance` values'''
#     o_xyz = '{}_{}.xyz'.format(uc['name'], n)

#     ## ==============================================
#     ## extract the xyz data for the region from the DEM
#     ## ==============================================
#     with open(o_xyz, 'w') as o_fh:
#         gdal_parse(uc['dem'], dst_xyz = o_fh, srcwin = gdal_srcwin(uc['dem'], region_buffer(this_region, (20 * uc['wg']['inc'])), mask = uc['msk']))
        
#     if os.stat(o_xyz).st_size != 0:
#         ## ==============================================
#         ## split the xyz data to inner/outer; outer is
#         ## the data buffer, inner will be randomly sampled
#         ## ==============================================
#         s_inner, s_outer = gmt_select_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = True)
#         if os.stat(s_inner).st_size != 0:
#             sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
#         else: sub_xyz = []
#         ss_len = len(sub_xyz)
#         if ss_samp is not None:
#             sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
#         else: sx_cnt = 1
#         sub_xyz_head = 'sub_{}_head.xyz'.format(n)
#         echo_msg('withholding {} out of {} points for error sampling'.format(ss_len - sx_cnt, ss_len))
#         np.random.shuffle(sub_xyz)
#         np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

#         ## ==============================================
#         ## generate the random-sample DEM
#         ## ==============================================
#         wc = waffles_config()
#         wc['name'] = 'sub_{}'.format(n)
#         wc['datalist'] = datalist_master(['{}:168:1',format(s_outer), '{}:168:1'.format(sub_xyz_head)])
#         wc['region'] = this_region
#         wc['inc'] = uc['inc']
#         wc['mod'] = uc['mod']
#         wc['mod_args'] = uc['mod_args']
#         sub_dem = waffles_run(wc)

#         ## ==============================================
#         ## generate the random-sample data MASK and PROX
#         ## ==============================================        
#         if sub_dem is not None:
#             wc['mod'] = 'num'
#             wc['mod_args'] = ('mode=k')
#             wc['name'] = 'sub_{}_msk'.format(n)
#             wc['region'] = region_buffer(this_region, uc['inc'] * 10)
#             sub_msk = waffles_run(wc)
#             sub_prox = '{}_prox.tif'.format(wc['name'])
#             gdal_proximity(sub_msk, sub_prox)
#             #gdal_infos(sub_prox)
            
#             ## ==============================================
#             ## Calculate the random-sample errors
#             ## ==============================================
#             sub_xyd = gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
#             sub_dp = gdal_query(sub_xyd, sub_prox, 'zg')
#         else: sub_dp = None

#         remove_glob(sub_xyz_head)
#         sub_xyz = None

#         if s_dp is None:
#             s_dp = sub_dp
#         else:
#             try:
#                 s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
#             except:
#                 echo_error_msg('found no error points...')
#                 pass
#     else: echo_error_msg('no data in sub-region...')
    
#     remove_glob(o_xyz)
#     remove_glob('sub_{}*'.format(n))
#     return(s_dp)

def waffles_interpolation_uncertainty(uc = _unc_config):
    '''calculate the interpolation uncertainty.'''
    dp = None
    sims = 20
    sim_loops = 2
    chnk_lvl = 2
    echo_msg('running INTERPOLATION uncertainty module using {}...'.format(uc['wg']['mod']))

    ## ==============================================
    ## region analysis
    ## ==============================================
    region_info = {}
    num_sum, g_max, num_perc = gdal_mask_analysis(mask = uc['msk'])
    prox_perc_95 = gdal_percentile(uc['prox'], 95)
    region_info[uc['wg']['name']] = [uc['wg']['region'], g_max, num_sum, num_perc, prox_perc_95] 
    for x in region_info.keys():
        echo_msg('region: {}: {}'.format(x, region_info[x]))

    ## ==============================================
    ## chunk region into sub regions
    ## ==============================================
    echo_msg('chunking region into sub-regions using chunk level {}...'.format(chnk_lvl))
    chnk_inc = int(region_info[uc['wg']['name']][4] / chnk_lvl)
    print(chnk_inc)
    sub_regions = region_chunk(uc['wg']['region'], uc['wg']['inc'], chnk_inc)
    echo_msg('chunked region into {} sub-regions.'.format(len(sub_regions)))

    ## ==============================================
    ## sub-region analysis
    ## ==============================================
    echo_msg('analyzing {} sub-regions...'.format(len(sub_regions)))
    sub_zones = {}    
    for sc, sub_region in enumerate(sub_regions):
        gdal_cut(uc['msk'], sub_region, 'tmp_msk.tif')
        gdal_cut(uc['dem'], sub_region, 'tmp_dem.tif')
        s_sum, s_g_max, s_perc = gdal_mask_analysis('tmp_msk.tif')
        s_dc = gdal_infos('tmp_dem.tif', True)
        zone = 'Bathy' if s_dc['zmax'] < 0 else 'Topo' if s_dc['zmin'] > 0 else 'BathyTopo'
        sub_zones[sc + 1] = [sub_region, s_g_max, s_sum, s_perc, s_dc['zmin'], s_dc['zmax'], zone]
        remove_glob('tmp_*.tif')
    for x in sub_zones.keys():
        echo_msg('Sub-region {}: {}'.format(x, sub_zones[x]))
        
    s_dens = np.array([sub_zones[x][3] for x in sub_zones.keys()])
    s_5perc = np.percentile(s_dens, 5)
    s_dens = None
    echo_msg('Sampling density for region is: {:.12f}'.format(s_5perc))

    ## ==============================================
    ## zone analysis / generate training regions
    ## ==============================================
    trainers = []
    bathy_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'Bathy']
    bathy_topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'BathyTopo']
    topo_tiles = [sub_zones[x] for x in sub_zones.keys() if sub_zones[x][6] == 'Topo']

    for z, tile_set in enumerate([bathy_tiles, bathy_topo_tiles, topo_tiles]):
        if len(tile_set) > 0:
            t_dens = np.array([x[3] for x in tile_set])
            t_50perc = np.percentile(t_dens, 50)
        else: t_50perc = 0.0
        echo_msg('Minimum sampling for {} tiles: {}'.format(uc['zones'][z].upper(), t_50perc))
        t_trainers = [x for x in tile_set if x[3] > t_50perc]
        echo_msg('possible {} training zones: {}'.format(uc['zones'][z].upper(), len(t_trainers)))
        trainers.append(t_trainers)
    trains = regions_sort(trainers)
    echo_msg('sorted training tiles.')
    echo_msg('analyzed {} sub-regions.'.format(len(sub_regions)))

    ## ==============================================
    ## split-sample simulations and error calculations
    ## ==============================================
    s_dp = None
    for sim in range(0, sims):
        echo_msg('performing SPLIT-SAMPLE simulation {} out of {}...'.format(sim + 1, sims))
        status = 0
        for s_l in range(0, sim_loops):
            for z, train in enumerate(trains):
                train_h = train[:25]
                ss_samp = s_5perc
                
                ## ==============================================
                ## perform split-sample analysis on each training region.
                ## ==============================================
                for n, sub_region in enumerate(train_h):
                    echo_msg('processing sub-region ({}) {}'.format(n, sub_region))
                    this_region = sub_region[0]
                    echo_msg('initial/desired sampling density: {}/{}'.format(sub_region[3], ss_samp))
                    if sub_region[3] < ss_samp: ss_samp = None
                    
                    o_xyz = '{}_{}.xyz'.format(uc['wg']['name'], n)

                    ## ==============================================
                    ## extract the xyz data for the region from the DEM
                    ## ==============================================
                    ds = gdal.Open(uc['dem'])
                    with open(o_xyz, 'w') as o_fh:
                        gdal_parse(ds, dst_xyz = o_fh, srcwin = gdal_srcwin(ds, region_buffer(this_region, (20 * uc['wg']['inc']))), mask = uc['msk'])
                    ds = None

                    if os.stat(o_xyz).st_size != 0:

                        ## ==============================================
                        ## split the xyz data to inner/outer; outer is
                        ## the data buffer, inner will be randomly sampled
                        ## ==============================================
                        s_inner, s_outer = gmt_select_split(o_xyz, this_region, 'sub_{}'.format(n), verbose = uc['wg']['verbose'])
                        if os.stat(s_inner).st_size != 0:
                            sub_xyz = np.loadtxt(s_inner, ndmin=2, delimiter = ' ')
                        else: sub_xyz = []
                        ss_len = len(sub_xyz)
                        if ss_samp is not None:
                            sx_cnt = int(sub_region[1] * (ss_samp / 100)) + 1
                        else: sx_cnt = 1
                        sub_xyz_head = 'sub_{}_head.xyz'.format(n)
                        echo_msg('withholding {} out of {} points for error sampling'.format(ss_len - sx_cnt, ss_len))
                        np.random.shuffle(sub_xyz)
                        np.savetxt(sub_xyz_head, sub_xyz[:sx_cnt], '%f', ' ')

                        ## ==============================================
                        ## generate the random-sample DEM
                        ## ==============================================
                        wc = waffles_config()
                        wc['name'] = 'sub_{}'.format(n)
                        wc['datalists'] = [s_outer, sub_xyz_head]
                        wc['region'] = this_region
                        wc['inc'] = uc['wg']['inc']
                        wc['mod'] = uc['wg']['mod']
                        wc['verbose'] = True
                        wc['mod_args'] = uc['wg']['mod_args']
                        sub_dem = waffles_run(wc)
                        print(sub_dem)
                        ## ==============================================
                        ## generate the random-sample data MASK and PROX
                        ## ==============================================        
                        if sub_dem is not None:
                            wc['mod'] = 'num'
                            wc['mod_args'] = ['mode=k']
                            wc['name'] = 'sub_{}_msk'.format(n)
                            wc['region'] = region_buffer(this_region, uc['wg']['inc'] * 10)
                            sub_msk = waffles_run(wc)
                            sub_prox = '{}_prox.tif'.format(wc['name'])
                            gdal_proximity(sub_msk, sub_prox)
                            #print(gdal_infos(sub_prox))
                            #print(gdal_infos(sub_dem))

                            ## ==============================================
                            ## Calculate the random-sample errors
                            ## ==============================================
                            sub_xyd = gdal_query(sub_xyz[sx_cnt:], sub_dem, 'xyd')
                            sub_dp = gdal_query(sub_xyd, sub_prox, 'zg')
                        else: sub_dp = None
                        #print(sub_dp)
                        remove_glob(sub_xyz_head)
                        sub_xyz = None

                        if s_dp is None:
                            s_dp = sub_dp
                        else:
                            #try:
                            s_dp = np.concatenate((s_dp, sub_dp), axis = 0)
                            #except:
                            #echo_error_msg('found no error points...')
                            #    pass
                    else: echo_error_msg('no data in sub-region...')    
                    remove_glob(o_xyz)
                    remove_glob('sub_{}*'.format(n))
                dp = s_dp    
        echo_msg('performed SPLIT-SAMPLE simulation {} out of {}; {} error points accumulated.'.format(sim + 1, sims, len(dp)))
    echo_msg('ran INTERPOLATION uncertainty module using {}.'.format(uc['wg']['mod']))

    if len(dp) > 0:
        ## ==============================================
        ## save err dist data files
        ## ==============================================
        echo_msg('gathered {} error points'.format(len(dp)))
        np.savetxt('{}.err'.format(uc['wg']['name']), dp, '%f', ' ')
        ec = err_plot(dp[:50000000], region_info[uc['wg']['name']][4])

        ## ==============================================
        ## apply error coefficient to full proximity grid
        ## ==============================================
        echo_msg('applying coefficient to proximity grid')
        ## USE numpy instead
        math_cmd = 'gmt grdmath {} 0 AND ABS {} POW {} MUL {} ADD = {}_dst_unc.tif=gd+n-9999:GTiff\
        '.format(uc['prox'], ec[2], ec[1], 0, uc['wg']['name'])
        run_cmd(math_cmd, verbose = uc['wg']['verbose'])
        echo_msg('applyed coefficient to proximity grid')

    return(dp)

if __name__ == '__main__':

    wg = waffles_config()
    wg['name'] = 'test'
    wg['datalists'] = [sys.argv[1]]
    wg['region'] = [float(x) for x in sys.argv[2].split('/')]
    wg['inc'] = gmt_inc2inc(sys.argv[3])
    wg['mod'] = 'surface'
    #dem = waffles_run(wg)
    dem = '{}.tif'.format(wg['name'])

    wg['name'] = 'test_num'
    wg['mod'] = 'num'
    wg['mod_args'] = ('mode=k',)
    msk = '{}.tif'.format('test_num')
    #msk = waffles_run(wg)
    prox = '{}_prox.tif'.format(wg['name'])
    gdal_proximity(msk, prox)

    wg['mod'] = 'surface'
    wg['mod_args'] = ()
    
    uc = _unc_config
    uc['wg'] = wg
    uc['dem'] = dem
    uc['msk'] = msk
    uc['prox'] = prox
    waffles_interpolation_uncertainty(uc)
    #dp = np.loadtxt('test.err')
    #ec = err_plot(dp, 100)
    #print(ec)
### End
