#!/usr/bin/env python
### dems.py
##
## Copyright (c) 2013 - 2020 Matthew Love <matthew.love@colorado.edu>
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
## Process DEMs
##
### Code:

import sys
import os
import re
import math
import time
import glob
import subprocess

import random
import numpy as np
import matplotlib.pyplot as plt

try:
    import osgeo.ogr as ogr
    import osgeo.gdal as gdal
    has_gdalpy = True
except ImportError:
    try:
        import ogr
        import gdal
        has_gdalpy = True
    except ImportError:
        has_gdalpy = False

try:
    import arcpy
    has_arcpy = True
except ImportError:
    has_arcpy = False

from regions import region
from datalists import datalist
from gdalfun import *

cudem_version = '0.0.5'

usage = '''dems.py - process DEMs

usage: dems.py [ -hvIER [ args ] ]

 Examples:
 % dems.py -Iinput.datalist -E0.000277777 -R-82.5/-82.25/26.75/27
 % dems.py --datalist input.datalist --increment 0.000277777 --region input_tiles_ply.shp

Report bugs to <matthew dot love at colorado dot edu>
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>'''

## data is 2 col file with 'err dist'
def err2coeff(data):
    try: 
        my_data = np.loadtxt( data, delimiter=' ' )
    except:
        sys.exit( 2 )

    error=my_data[:,0]
    distance=my_data[:,1]
        
    max_int_dist = np.max( distance )
    nbins = 10

    coeff_guess=[0, 0.1, 0.2]
    n, _ = np.histogram( distance, bins = nbins )

    # want at least 2 values in each bin?
    while 0 or 1 in n:
        nbins -= 1
        n, _ = np.histogram( distance, bins = nbins )

    serror, _ = np.histogram( distance, bins = nbins, weights = error )
    serror2, _ = np.histogram( distance, bins = nbins, weights = error**2 )

    mean = serror / n
    std = np.sqrt( serror2 / n - mean * mean )

    ydata = np.insert( std, 0, 0 )
    
    bins_orig=( _[1:] + _[:-1] ) / 2
    xdata = np.insert( bins_orig, 0, 0 )

    fitfunc = lambda p, x: p[0] + p[1] * ( abs( x ) ** p[2] )
    errfunc = lambda p, x, y: y - fitfunc( p, x )
    
    out, cov, infodict, mesg, ier = optimize.leastsq( errfunc, coeff_guess, args = ( xdata, ydata ), full_output = True )
    return out

## =============================================================================
##
## Command execution, et cetra
##
## OS System commands and checks.
## run a command with a progress bar with 'run_cmd'
## check if a command exists on the system with 'cmd_exists'
##
## =============================================================================
cmd_exists = lambda x: any( os.access( os.path.join( path, x ), os.X_OK ) for path in os.environ["PATH"].split( os.pathsep ))

def run_cmd(cmd, prog_message = None):
    'Run a command with or without a progress bar.'

    p = subprocess.Popen( cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
    
    if prog_message:
        tw = prog_bar( prog_message )
        while p.poll() != 0:
            tw._update()
            time.sleep( 2 )

        tw._end( p.returncode )

    out, err = p.communicate()
    return out, p.returncode

## =============================================================================
##
## Progress Bar
##
## with 'prog_message' print a simple progress bar and message.
## use the 'prog_bar.pm' variable to update the message while running.
##
## =============================================================================
class prog_bar:
    def __init__(self, prog_message):
        
        self.tw = 4
        self.count = 0
        self.pc = self.count % self.tw
        self.pm = prog_message
        self.colors = True

        self.cyan = '\033[96m'
        self.green = '\033[92m'
        self.red2 = '\033[91m'
        self.red = '\x1b[31m'
        self.NC2 = '\033[00m'
        self.NC = '\x1b[00m'

        self.spinner = ['\\', '|', '/', '~', '']

        sys.stdout.write('\r[%s] %s' % (" " * self.tw, self.pm))
        sys.stdout.flush()

    def _update(self):
        self.pc = ( self.count % ( self.tw + 1 ))
        sys.stdout.write('\r[%-4s]' %( self.green + ( '*' * self.pc ) + self.red + self.spinner[self.pc] + self.NC + (' ' * (( self.tw-1 ) - self.pc ))))
        sys.stdout.flush()
        self.count += 1

    def _end(self, status):
        if status != 0:
            sys.stdout.write( '\r[%sFAIL%s] %s\n' %( self.red, self.NC, self.pm ))
        else:
            sys.stdout.write( '\r[ %sOK%s ] %s \n' %( self.green, self.NC, self.pm ))

        sys.stdout.flush()

## =============================================================================
##
## GDAL CLASS with GDAL Functions
##
## GDAL functions for processing raster and vector data.
##
## =============================================================================
class gdal_dem():
    def __init__(self):
        self.status = 0

    def con_dec(self, x, dec):
        '''Return a float string with n decimals
        (used for ascii output).'''

        if x is None:
            return

        fstr = "%." + str( dec ) + "f"
        return fstr % x

    # Convert a geographic x,y value to a pixel location of geoTransform
    def geo2pixel(self, geo_x, geo_y, geoTransform ):
        if geoTransform[2] + geoTransform[4] == 0:
            pixel_x = ( geo_x - geoTransform[0] ) / geoTransform[1]
            pixel_y = ( geo_y - geoTransform[3] ) / geoTransform[5]
        else:
            pixel_x, pixel_y = apply_gt( geo_x, geo_y, invert_gt( geoTransform ))
        return int( pixel_x ), int( pixel_y )

    # Convert a pixel location to geographic coordinates given geoTransform
    def pixel2geo( pixel_x, pixel_y, geoTransform ):
        geo_x, geo_y = apply_gt( pixel_x, pixel_y, geoTransform )
        return geo_x, geo_y

    def apply_gt(self, in_x, in_y, geoTransform ):
        out_x = geoTransform[0] + in_x * geoTransform[1] + in_y * geoTransform[2]
        out_y = geoTransform[3] + in_x * geoTransform[4] + in_y * geoTransform[5]
        return out_x, out_y

    def invert_gt(self, geoTransform):
        det = geoTransform[1] * geoTransform[5] - geoTransform[2] * geoTransform[4]
        if abs( det ) < 0.000000000000001:
            return

        invDet = 1.0 / det

        # compute adjoint and divide by determinate
        outGeoTransform = [0, 0, 0, 0, 0, 0]
        outGeoTransform[1] = geoTransform[5] * invDet
        outGeoTransform[4] = -geoTransform[4] * invDet
        outGeoTransform[2] = -geoTransform[2] * invDet
        outGeoTransfrom[5] = geoTransform[1] * invDet
        outGeoTransform[0] = ( geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5] ) * invDet
        outGeoTransform[3] = ( -geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4] ) * invDet
        return outGeoTransform 

    def proximity(self, srcgrd, dstgrd):
        '''Compute a proximity grid via GDAL'''

        dst_ds = None
        prog_func = None
        dst_nodata = -9999

        src_ds = gdal.Open( srcgrd )
    
        if src_ds is None:
            self.status = -1

        srcband = src_ds.GetRasterBand( 1 )

        if dst_ds is None:
            drv = gdal.GetDriverByName( 'GTiff' )
            dst_ds = drv.Create( dstgrd, src_ds.RasterXSize, src_ds.RasterYSize, 1, gdal.GetDataTypeByName( 'Float32' ), [] )

        dst_ds.SetGeoTransform( src_ds.GetGeoTransform() )
        dst_ds.SetProjection( src_ds.GetProjectionRef() )
    
        dstband = dst_ds.GetRasterBand(1)
        dstband.SetNoDataValue( dst_nodata )

        gdal.ComputeProximity( srcband, dstband, ['DISTUNITS=PIXEL'], callback = prog_func )

        dstband = srcband = dst_ds = src_ds = None

    def query(self, srcxyz, srcgrd, out_form):
        '''Query a gdal-compatible grid file with xyz data.'''

        xyzl = []

        # Process the src grid file
        ds = gdal.Open( srcgrd )
        dsgeot = ds.GetGeoTransform()
        dsband = ds.GetRasterBand( 1 )
        dsnodata = dsband.GetNoDataValue()

        dsdt = gdal.GetDataTypeName( dsband.DataType )

        # Load the src grid into a numpy array
        tgrid = dsband.ReadAsArray()

        cellsize = [float( dsgeot[1] ), float( dsgeot[5] )]
        xextent = float( dsgeot[0] )
        yextent = float( dsgeot[3] )
        
        # Process the src xyz data
        for xyz in srcxyz:
            x = xyz[0]
            y = xyz[1]
            try: z = xyz[2]
            except: z = dsnodata

            # Continue if values are reasonable.
            if x > xextent and y < yextent:
                xpos,ypos = self.geo2pixel( x, y, dsgeot )

                # Locate the grid cell and get it's value
                try: g = tgrid[ypos,xpos]
                except: g = dsnodata

                d = c = m = s = dsnodata
                
                if g != dsnodata:
                    d = z - g
                    m = z + g
                    c = self.con_dec( math.fabs( float( d / ( g+0.00000001 ) * 100 )), 2 )
                    s = self.con_dec( math.fabs( d / (z + ( g+0.00000001 ))), 4 )
                    
                d = self.con_dec( d, 4 )

                outs = []
                for i in out_form:
                    outs.append( vars()[i] )
                xyzl.append( np.array( outs, dtype = dsdt ))

        dsband = ds = None
        return np.array( xyzl, dtype = dsdt )

    def toxyz(self, src_gdal, dst_xyz):
        self.status = 0
        band_nums = []
        nodata = ['-9999', 'nan']
        srcwin = None
        delim = ' '
        skip = 1

        if band_nums == []: band_nums = [1]

        srcds = gdal.Open( src_gdal )

        if srcds is None:
            self.status = -1
            return self.status

        bands = []
        for band_num in band_nums: 
            band = srcds.GetRasterBand( band_num )
            if band is None:
                self.status = -1
                return self.status
            bands.append( band )

        gt = srcds.GetGeoTransform()

        if srcwin is None:
            srcwin = (0,0,srcds.RasterXSize,srcds.RasterYSize)

        if dst_xyz is not None:
            dst_fh = open( dst_xyz,'wt' )
        else:
            dst_fh = sys.stdout

        band_format = (("%g" + delim) * len(bands)).rstrip(delim) + '\n'

        if abs( gt[0] ) < 180 and abs( gt[3] ) < 180 \
           and abs( srcds.RasterXSize * gt[1] ) < 180 \
           and abs( srcds.RasterYSize * gt[5] ) < 180:
            format = '%.10g' + delim + '%.10g' + delim + '%s'
        else:
            format = '%.3f' + delim + '%.3f' + delim + '%s'

        for y in range( srcwin[1], srcwin[1] + srcwin[3], skip ):

            data = []
            for band in bands:
                nodata.append( band.GetNoDataValue() )
                band_data = band.ReadAsArray( srcwin[0], y, srcwin[2], 1 )    
                band_data = np.reshape( band_data, ( srcwin[2], ))
                data.append( band_data )

            for x_i in range( 0, srcwin[2], skip ):
                x = x_i + srcwin[0]

                geo_x = gt[0] + (x + 0.5) * gt[1] + (y + 0.5) * gt[2]
                geo_y = gt[3] + (x + 0.5) * gt[4] + (y + 0.5) * gt[5]

                x_i_data = []
                for i in range(len(bands)):
                    x_i_data.append(data[i][x_i])

                band_str = band_format % tuple( x_i_data )

                line = format % ( float( geo_x ),float( geo_y ), band_str )

                if band_str not in nodata:
                    dst_fh.write( line )

## =============================================================================
##
## DEM Class for DEM development et cetra
##
## Generate DEMs with GMT and MBSystem
## idatalist and iregion are datalist and region instances, respectively
## the `dem` dictionary holds the file names of the generated DEM files
## all files are GTiff unless specified as 'grd'
## including 'dem', 'dem-grd', 'num', 'num-grd', 'num-msk', 'prox', 'slope'
##
## =============================================================================
class dem:
    def __init__(self, idatalist, iregion, iinc = "0.000277777", oname = None):
        self.status = 0
        self.dem = {}
        self.inc = float( iinc )
        self.datalist = idatalist
        self.region = iregion
        self.proc_region = self.region.buffer( 10 * self.inc )
        self.dist_region = self.region.buffer( 6 * self.inc )
        self.max_prox = self.max_num = None
        self.oname = oname
        if self.oname is None: 
            self.oname = os.path.basename( self.datalist._path ).split( '.' )[0]

    ## 'Validate' the DEM (check if all files in `dem` exist)
    def _valid_p(self):
        for key in self.dem:
            if not os.path.exists( self.dem[key] ):
                return( False )
        return( True )

    ## Remove keys to DEM filenames if they don't exist
    def _rectify(self):
        for key in self.dem:
            if not os.path.exists( self.dem[key] ):
                del self.dem[key]
                
    ## Set a key/value in the `dem` dictionary
    def _set_dem( self, key, value ):
        self.dem[key] = value

    ## Print out the key/value pairs from the `dem` dictionary
    def _print_dem(self):
        for key in self.dem:
            print( key, self.dem[key] )

    ## Remove file associated with key
    def _dem_remove(self, key):
        try: 
            os.remove( self.dem[key] )
            return 0
        except: return -1

    ## Generic Functions
    def grd2tif(self, src_grd):
        '''Convert the grd file to tif using GMT'''

        if os.path.exists( src_grd ):
            grd2tif_cmd = ( 'gmt grdconvert %s %s.tif=gd+n-9999:GTiff -V' 
                            %( src_grd, os.path.basename( src_grd ).split( '.' )[0] ))
            
            out, self.status = run_cmd( grd2tif_cmd, "converting grd to tif" )
        else: self.status = -1

        return self.status

    def grdinfo(self, src_grd):
        '''Return an info list of `src_grd`'''

        if not os.path.exists( src_grd ):
            self.status = -1
            return self.status
        else:
            grdinfo_cmd = ( 'gmt grdinfo %s -C' %( src_grd ))
            out, self.status = run_cmd( grdinfo_cmd, "snarfing info from grid" )
            return out.split()

    def grdcut(self, src_grd, src_region, dst_grd):
        '''Cut `src_grd` to `src_region` '''

        if os.path.exists( src_grd ):
            self.cut_cmd1 = ( 'gmt grdcut %s -G%s %s -V' 
                              %( src_grd, dst_grd, src_region.gmt ))

    def grd2xyz(self, src_grd, dst_xyz, region = None, mask = None):
        '''Convert `src_grd` to xyz possibly using a nodata mask and/or a region'''

        if mask:
            self.grdmask_cmd = ( 'gmt grdmath %s NOT %s = tmp.grd' %( mask, src_grd ))
            out, self.status = run_cmd( self.grdmask_cmd, "masking DEM" )
            src_grd = 'tmp.grd'

        if region and region._valid:
            self.grd2xyz_cmd = ( 'gmt grd2xyz %s %s -s > %s' %( src_grd, region.gmt, dst_xyz ))
        else:
            self.grd2xyz_cmd = ( 'gmt grd2xyz %s -s > %s' %( igrid, dst_xyz ))
            
        out, self.status = run_cmd(self.grd2xyz_cmd, "converting grd to xyz")

        if mask:
            if os.path.exists('tmp.grd'):
                os.remove('tmp.grd')

        return self.status
        
    ## DEM Et Cetra
    def surface(self):
        '''Generate a DEM with GMT surface'''

        dem_cmd = ( 'cat %s | gmt blockmean %s -I%s -r -V | gmt surface %s -I%s -G%s_p.grd -T.35 -Z1.2 -r -V -Lud -Lld' 
                    %(self.datalist._echo_datafiles( ' ' ), self.proc_region.gmt, self.inc, self.proc_region.gmt, self.inc, self.oname ))
        dem_cmd1 = ( 'gmt grdcut %s_p.grd -G%s.grd %s -V' 
                     %( self.oname, self.oname, self.dist_region.gmt ))

        out, self.status = run_cmd( dem_cmd, 'generating DEM surface' )
        out, self.status = run_cmd( dem_cmd1, 'cutting DEM surface' )

        if os.path.exists( '%s_p.grd' %( self.oname )):
            os.remove( '%s_p.grd' %( self.oname ))

        self.grd2tif( self.dem['grd'] )

        self.dem['dem-grd'] = ( '%s.grd' %( self.oname ))
        self.dem['dem'] = ( '%s.tif' %( self.oname ))

        return self.status

    def num(self):
        '''Generate a num and num-msk grid with GMT'''

        num_cmd0 = ( 'cat %s | gmt xyz2grd %s -I%s -r -V -G%s_num.grd -An'
                     %( self.datalist._echo_datafiles( ' ' ), self.dist_region.gmt, self.inc, self.oname ))
        num_cmd1 = ('gmt grdmath %s_num.grd 0 MUL 1 ADD 0 AND = %s_num_msk.tif=gd+n-9999:GTiff -V' %( self.oname, self.oname ))

        out, self.status = run_cmd( num_cmd0, "generating num grid" )
        out, self.status = run_cmd( num_cmd1, "generating num mask grid" )

        self.grd2tif( '%s_num.grd' %( self.oname ))

        self.dem['num'] = '%s_num.tif' %( self.oname )
        self.dem['num-grd'] = '%s_num.grd' %( self.oname )
        self.dem['num-msk'] = '%s_num_msk.tif' %(self.oname )

        self.max_num = int(self.grdinfo(self.dem['num'])[6])

        return self.status   

    def mbgrid(self, extras = False):
        '''Generate a DEM and num grid with MBSystem'''

        mbgrid_cmd = ( 'mbgrid -I%s %s -E%s/%s/degrees! -O%s -A2 -G3 -F1 -N -C10/3 -S0 -X0.1 -T35 -M' 
                       %( self.datalist._path, self.dist_region.gmt, self.inc, self.inc, self.oname ))

        out, self.status = run_cmd( mbgrid_cmd, 'generating DEM and num grid' )

        self.dem['dem-grd'] = '%s.grd' %( self.oname )
        self.dem['num-grd'] = '%s_num.grd' %( self.oname )

        self.grd2tif( self.dem['dem-grd'] )
        self.dem['dem'] = '%s.tif' %( self.oname )

        self.grd2tif(self.dem['num-grd'])
        self.dem['num'] = '%s_num.tif' %( self.oname )

        try:
            os.remove( '%s.cmd' %( self.dem['dem-grd'] ))
            os.remove( '%s.cmd' %( self.dem['num-grd'] ))
            os.remove( '%s.mb-1' %( self.oname ))
            os.remove( '%s_sd.grd' %( self.oname ))
            os.remove( '%s_sd.grd.cmd' %( self.oname ))
        except: pass

        self.max_num = int( self.grdinfo( self.dem['num-grd'] )[6] )
        
        num_msk_cmd = ( 'gmt grdmath %s 0 MUL 1 ADD 0 AND = %s_num_msk.tif=gd+n-9999:GTiff' %( self.dem['num'], self.oname ))
        out, self.status = run_cmd( num_msk_cmd, 'generating num mask grid' )

        self.dem['num-msk'] = ( '%s_num_msk.tif' %( self.oname ))

        return self.status

    def slope(self):
        '''Generate a Slope grid from the DEM with GMT'''

        try: self.dem['dem-grd']
        except KeyError: 
            self.status = -1
            return self.status

        slope_cmd0 = ( 'gmt grdgradient -fg %s -S%s_pslp.grd -D -R%s -V' 
                       %( self.dem['dem-grd'], self.oname, self.dem['dem-grd'] ))
        slope_cmd1 = ( 'gmt grdmath %s_pslp.grd ATAN PI DIV 180 MUL = %s_slp.tif=gd+n-9999:GTiff'
                       %( self.oname, self.oname ))

        out, self.status = run_cmd( slope_cmd0, "generating DEM gradient" )
        out, self.status = run_cmd( slope_cmd1, "generating DEM slope" )

        if os.path.exists( '%s_pslp.grd' %( self.oname )):
            os.remove( '%s_pslp.grd' %( self.oname ))

        self.dem['slope'] = ( '%s_slp.tif' %( self.oname ))

        return self.status

    def proximity(self):
        '''Generate a Proximity grid with GDAL'''

        try: self.dem['num-msk']
        except KeyError: 
            self.status = -1
            return self.status

        self.dem['prox'] = ( '%s_prox.tif' %( self.oname ))

        gdal_dem().proximity( self.dem['num-msk'], self.dem['prox'] )
        
        self.max_prox = self.grdinfo( self.dem['prox'] )[6]
        print self.max_prox

        return self.status

## =============================================================================
##
## Run cudem.py from command-line
##
## =============================================================================
def main():

    iregion = None
    idatalist = None
    igrid = None
    iinc = '0.000277777'
    status = 0
    these_regions = []
    sub_count = 0
    dem_alg = 'mbgrid'
    dem_algs = ['mbgrid', 'gmt-surface']
    do_unc = False

    argv = sys.argv
        
    ## parse command line arguments.
    i = 1
    while i < len( argv ):
        arg = argv[i]

        if arg == '--region' or arg == '-R':
            iregion = str( argv[i + 1] )
            i = i + 1
        elif arg[:2] == '-R':
            iregion = str( arg[2:] )

        elif arg == '--datalist' or arg == '-I':
            idatalist = str( argv[i + 1] )
            i = i + 1
        elif arg[:2] == '-I':
            idatalist = str( arg[2:] )

        elif arg == '--increment' or arg == '-E':
            iinc = str( argv[i + 1] )
            i = i + 1
        elif arg[:2] == '-E':
            iinc = str( arg[2:] )

        elif arg == '--grid' or arg == '-G':
            igrid = str( argv[i + 1] )
            i = i + 1
        elif arg[:2] == '-G':
            igrid = str( arg[2:] )

        elif arg == '--alogrithm' or arg == '-A':
            dem_alg = str( argv[i + 1] )
            i = i + 1
        elif arg[:2] == '-A':
            dem_alg = str( arg[2:] )

        elif arg == '--uncertainty' or arg == '-u':
            do_unc = True

        elif arg[0] == '-':
            print( usage )
            sys.exit( 0 )

        elif idatalist is None:
            idatalist = arg

        elif iregion is None:
            iregion = arg

        else:
            print( usage )
            sys.exit( 1 )

        i = i + 1

    if iregion is None:
        print( usage )
        sys.exit( 1 )
        #iregion = raw_input( 'input region: ' )

    #if idatalist is None:
        #idatalist = raw_input( 'input datalist: ' )
        #print(usage)
        #sys.exit(1)

    if dem_alg not in dem_algs:
        dem_alg = 'gmt'

    ## check platform
    tw = prog_bar( 'checking platform' )
    platform = sys.platform
    tw.pm = 'checking platform - %s' %( platform )
    tw._end( status )

    ## check for installed software
    tw = prog_bar( 'checking for GMT' )
    if cmd_exists( 'gmt' ): 
        gmt_vers, status = run_cmd( 'gmt --version' )
        tw.pm = 'checking for GMT - %s' %( gmt_vers.rstrip() )
    else: status = -1
    tw._end( status )

    tw = prog_bar( 'checking for MBSystem' )
    if cmd_exists( 'mbgrid' ): 
        mbs_vers, status = run_cmd( 'mbgrid -version' )
        tw.pm = 'checking for MBSystem - %s' %( mbs_vers.split( '\n' )[3].rstrip().split()[2] )
    else: status = -1
    tw._end( status )

    tw = prog_bar( 'checking for GDAL command-line' )
    if cmd_exists( 'gdal-config' ): 
        gdal_vers, status = run_cmd( 'gdal-config --version' )
        tw.pm = 'checking for GDAL command-line - %s' %( gdal_vers.rstrip() )
    else: status = -1
    tw._end( status )

    tw = prog_bar( 'checking for GDAL python bindings' )
    if has_gdalpy: 
        status = 0
        gdal_vers = gdal.__version__
        tw.pm = 'checking for GDAL python bindings - %s' %( gdal_vers )
    else: status = -1
    tw._end( status )

    if platform != 'linux2':
        tw = prog_bar( 'checking for arcpy python bindings' )
        if has_arcpy: 
            status = 0
        else: status = -1
        tw._end( status )

    ## process input region(s)
    tw = prog_bar( "processing region(s)" )
    try: 
        these_regions = [region( iregion )]
    except:
        if os.path.exists( iregion ):
            _poly = ogr.Open( iregion )
            _player = _poly.GetLayer( 0 )
            for pf in _player:
                _pgeom = pf.GetGeometryRef()
                these_regions.append( region( "/".join( map( str, _pgeom.GetEnvelope() ))))

    if len(these_regions) == 0:
        status = -1

    for this_region in these_regions:
        if not this_region._valid: 
            status = -1
    tw.pm = 'processing %s region(s)' %( len( these_regions ))
    tw._end( status )

    if status == -1:
        print( usage )
        sys.exit( 1 )

    for this_region in these_regions:
        
        ## process input datalist
        this_datalist = datalist( idatalist, this_region )
        if not this_datalist._valid:
            break

        ## process DEM
        ## generate DEM and Num grid using full region
        this_surf = dem( this_datalist, this_region, iinc )

        ## MBSystem mbgrid
        if dem_alg == 'mbgrid':
            this_surf.mbgrid()

        ## GMT surface
        elif dem_alg == 'gmt-surface':
            this_surf.surface()
            this_surf.num()

        tw = prog_bar( "validating DEM" )

        this_surf._rectify()
        if this_surf._valid_p():
            if not this_surf.dem['dem'] or not this_surf.dem['num']:
                status = -1
        else:
            status = -1
            
        tw._end( status )

        ## Interpolation Uncertainty
        if do_unc:
        
            tw = prog_bar( "preparing for interpolation uncertainty" )._end( status )
            status = this_surf.proximity()
    
            for sub_region in this_region.split( 4 ):
                sub_count += 1
            
                ## Extract XYZ data for sub-region and randomly sample
                status = this_surf.grd2xyz(this_surf.dem['dem-grd'], '%s.xyz' %( this_surf.oname ), region=sub_region.buffer(10*float(iinc)), mask=this_surf.dem['num-grd'])

                sub_xyz = np.loadtxt( this_surf.dem['xyz'], delimiter = ' ' )
                
                np.random.shuffle( sub_xyz )
                sx_len = len( sub_xyz )
                sx_len_pct = int( sx_len * .5 )

                np.savetxt( 'sub_%d_rest.xyz' %( sub_count ), sub_xyz[sx_len_pct:], '%f', ' ' )
                
                sub_datalist =  datalist( 'sub_%d.datalist' %( sub_count ), sub_region )
                sub_datalist._append_datafile( 'sub_%d_rest.xyz' %( sub_count ), 168, 1 )
                sub_datalist._reset()

                ## Generate the random-sample DEM            
                sub_surf = dem( sub_datalist, sub_region, iinc )
                sub_surf.mbgrid()
            
                status = sub_surf.proximity()
            
                ## Query the random-sample DEM with the withheld data
                sub_xyd = gdal_dem().query( sub_xyz[:sx_len_pct], sub_surf.dem['tif'], 'xyd' )

                ## Query the random-sample proximity grid with the withheld data
                sub_dp = gdal_dem().query( sub_xyd, sub_surf.dem['prox'], 'zg' )
                #print np.max(sub_dp[:,1])
            
                ## Cleanup
                fl = glob.glob('sub_%d*' %(sub_count))
                for f in fl:
                    try: os.remove(f)
                    except: pass
            
if __name__ == '__main__':
    main()

### END
