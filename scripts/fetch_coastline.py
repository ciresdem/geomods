#!/usr/bin/env python
### fetch_coastline.py
##
## Copyright (c) 2019 - 2021 CIRES Coastal DEM Team
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
##
## fetch and process a coastline
##
### Code:

import os
import sys
import glob
from geomods import waffles
from geomods import fetches

cs_region = [-161.000278, -152.999722, 17.999722, 23.000278]
want_landsat = True
want_gshhg = False

waffles.run_cmd('regions {} > region_buff.gmt'.format(waffles.region_format(cs_region, 'gmt')), verbose = True)
waffles.run_cmd('gdal_rasterize -tr 0.0000925 0.0000925 -te {} -burn 0 -ot Int16 -co COMPRESS=DEFLATE region_buff.gmt nhd.tif'.format(waffles.region_format(cs_region, 'te')), verbose = True)
waffles.run_cmd('gdal_rasterize -tr 0.0000925 0.0000925 -te {} -burn 1 -ot Int16 -co COMPRESS=DEFLATE region_buff.gmt landsat.tif'.format(waffles.region_format(cs_region, 'te')), verbose = True)
waffles.run_cmd('gdal_rasterize -tr 0.0000925 0.0000925 -te {} -burn 1 -ot Int16 -co COMPRESS=DEFLATE region_buff.gmt gshhg.tif'.format(waffles.region_format(cs_region, 'te')), verbose = True)

### GSHHG (GMT)
if want_gshhg:
    waffles.run_cmd('gmt grdmath {} -I0.0000925 gshhg.tif LDISTG -Df -V = test.tif=gd:GTiff'.format(waffles.region_format(cs_region, 'gmt')), verbose = True)
    #waffles.run_cmd('gmt grdlandmask {} -Gtest.tif=gd:GTiff -I0.0000925 -Df -V'.format(waffles.region_format(cs_region, 'gmt')), verbose = True)
    sys.exit()

### LANDSAT
if want_landsat:
    ls_wrld = 'https://rmgsc.cr.usgs.gov/outgoing/ecosystems/Global/WorldEcologicalLandUnits2015data.zip'
    ls_wrld_zip = 'WorldEcologicalLandUnits2015data.zip'

    waffles.echo_msg('fetching landsat world ecological land units')
    if fetches.fetch_file(ls_wrld, ls_wrld_zip, overwrite = True) == 0:
        waffles.echo_msg('unzipping landsat zip file')
        src_ls, ls_zips = waffles.procs_unzip(ls_wrld_zip, ['tif'])
        waffles.gdal_set_nodata(src_ls, -2147483647)

    waffles.echo_msg('cutting landsat to user region of {}'.format(cs_region))
    waffles.gdal_cut(src_ls, cs_region, 'tmp_ls.tif')

    waffles.echo_msg('masking landsat raster')
    waffles.gdal_mask('tmp_ls.tif', 'tmp_ls_msk.tif')
    waffles.run_cmd('gdal_polygonize.py tmp_ls_msk.tif ls_coast_ply.shp', verbose = True)
    waffles.run_cmd('ogr2ogr -sql "select ST_Buffer(geometry, 0.01) from ls_coast_ply" -dialect SQLite ls_coast_ply1.shp ls_coast_ply.shp', verbose = True)
    waffles.run_cmd('gdal_rasterize -burn 0 ls_coast_ply1.shp landsat.tif')
else: waffles.run_cmd('gdal_rasterize -burn 0 landsat_all_NA.shp landsat.tif')

### USGS NDB
fl = fetches.fetch_infos['tnm'][0](waffles.region_buffer(cs_region, 5, pct = True), [], None)
r = fl.run(ds = 6, formats = 'FileGDB', extent = 'HU-4 Subregion')
fr = fetches.fetch_results(r, cs_region, fl._outdir, None)
fr.start()
fr.join()

r_shp = []
for result in r:
    try:
        gdb_zip = os.path.join(result[2], result[1])
        gdb_files = waffles.unzip(gdb_zip)
        gdb, gdb_files = waffles.procs_unzip(gdb_zip, ['gdb'])

        gdb_bn = os.path.basename('.'.join(gdb_zip.split('.')[:-1]))
        gdb = gdb_bn + '.gdb'
    
        print(gdb)
        print(gdb_bn)
        waffles.run_cmd('ogr2ogr {}_NHDArea.shp {} NHDArea -overwrite'.format(gdb_bn, gdb), verbose = True)
        r_shp.append('{}_NHDArea.shp'.format(gdb_bn))
        waffles.run_cmd('ogr2ogr {}_NHDPlusBurnWaterBody.shp {} NHDPlusBurnWaterBody -overwrite'.format(gdb_bn, gdb), verbose = True)
        r_shp.append('{}_NHDPlusBurnWaterBody.shp'.format(gdb_bn))
    except: waffles.echo_error_msg('unable to process {}'.format(result))

[waffles.run_cmd('ogr2ogr  -update -append nhdArea_merge.shp {}'.format(shp), verbose = True) for shp in r_shp]
waffles.run_cmd('gdal_rasterize -tr 0.0000925 0.0000925 -te {} -burn 1 -ot Int16 -co COMPRESS=DEFLATE nhdArea_merge.shp nhd_tmp.tif'.format(waffles.region_format(cs_region, 'te')), verbose = True)
waffles.run_cmd('gdal_polygonize.py -8 nhd_tmp.tif nhd_rast.shp', verbose = True)
waffles.run_cmd('ogr2ogr -dialect SQLITE -sql "SELECT * FROM nhd_rast WHERE DN=1 order by ST_AREA(geometry) desc limit 8" nhd_clean.shp nhd_rast.shp', verbose = True)
waffles.run_cmd('gdal_rasterize -tr 0.0000925 0.0000925 -te {} -burn 1 -ot Int16 -co COMPRESS=DEFLATE nhd_clean.shp nhd.tif'.format(waffles.region_format(cs_region, 'te')), verbose = True)

### OTHER

### Combine
waffles.run_cmd('gdal_calc.py -A nhd.tif -B landsat.tif --outfile=combined_coast_sum.tif --calc="A + B" --format=GTiff --overwrite', verbose = True)
waffles.run_cmd('gdal_calc.py -A combined_coast_sum.tif --outfile=combined_coast_rc.tif --calc="1*(A > 0)" --format=GTiff --overwrite', verbose = True)
waffles.run_cmd('gdal_polygonize.py -8 combined_coast_rc.tif combined_coast_all.shp', verbose = True)
waffles.run_cmd('ogr2ogr -dialect SQLITE -sql "SELECT * FROM combined_coast_all WHERE DN=0" combined_coast.shp combined_coast_all.shp')
