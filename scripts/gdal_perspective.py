#!/usr/bin/env python

import sys
import os
import math
from osgeo import gdal
import numpy as np

gp_version = '0.0.2'

gg_equat = 111321.543

## Colors
trs = (-11000,-10500,-10000,-9500,-9000,-8500,-8000,-7500,-7000,
        -6500,-6000,-5500,-5000,-4500,-4000,-3500,-3000,-2500,-2000,
        -1500,-1000,-500,-0.001,0,100,200,500,1000,1500,2000,2500,3000,
        3500,4000,4500,5000,5500,6000,6500,7000,7500,8000)

colors = ([10,0,121],[26,0,137],[38,0,152],[27,3,166],[16,6,180],[5,9,193],
          [0,14,203],[0,22,210],[0,30,216],[0,39,223],[12,68,231],[26,102,240],
          [19,117,244],[14,133,249],[21,158,252],[30,178,255],[43,186,255],[55,193,255],
          [65,200,255],[79,210,255],[94,223,255],[138,227,255],[138,227,255],
          [51,102,0],[51,204,102],[187,228,146],[255,220,185],[243,202,137],[230,184,88],
          [217,166,39],[168,154,31],[164,144,25],[162,134,19],[159,123,13],[156,113,7],
          [153,102,0],[162,89,89],[178,118,118],[183,147,147],[194,176,176],[204,204,204],
          [229,229,229],[138,227,255],[51,102,0])
          
elevs = [-20,-19.09090909,-18.18181818,-17.27272727,-16.36363636,-15.45454545,
         -14.54545455,-13.63636364,-12.72727273,-11.81818182,-10.90909091,
         -10,-9.090909091,-8.181818182,-7.272727273,-6.363636364,-5.454545455,
         -4.545454545,-3.636363636,-2.727272727,-1.818181818,-0.909090909,-0.001,
         0,48.06865898,96.13731797,240.3432949,480.6865898,721.0298848,
         961.3731797,1201.716475,1442.05977,1682.403064,1922.746359,2163.089654,
         2403.432949,2643.776244,2884.119539,3124.462834,3364.806129,3605.149424,
         3845.492719]

def scale_el(value, min, max, tr):
    if value > 0 and max > 0:
        return (max * tr) / 8000
    elif value < 0 and min < 0:
        return (min * tr) / -11000
    elif value == 0:
        return 0
    else: return None

def make_cpt(ingrd, grdmm):
    out_file = "%s.cpt" %(ingrd)
    out_cpt = open(out_file, 'w')

    ## NoData
    nd = returnNoData(ingrd)
    nd_string = str(nd) + " 0 0 0\n"
    out_cpt.write(nd_string)

    for i,j in enumerate(elevs):
        elevs[i] = scale_el(elevs[i], grdmm[4], grdmm[5], trs[i])
        if elevs[i] != None:
            outs = elevs[i],colors[i][0],colors[i][1],colors[i][2]
            #print outs
            out_cpt.write(" ".join(map(str, outs)))
            out_cpt.write("\n")

## Get min/max info about `ingrd'
def returnMinMax(ingrd):
    # Process the grid file
    ds = gdal.Open(ingrd)
    comp_geot = ds.GetGeoTransform()

    cellsize = [float(comp_geot[1]), float(comp_geot[5])]
    xmin = float(comp_geot[0])
    ymax = float(comp_geot[3])
    #nodata = ds.GetRasterBand(1).GetNoDataValue()
    #print xmin,ymax

    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    tgrid = band.ReadAsArray()
    print(nodata)
    #tgrid[tgrid==nodata]=float('nan')
    #outarray[np.isnan(outarray)]=ndata
    if np.isnan(nodata):
        tgrid=np.ma.MaskedArray(tgrid, mask=(np.isnan(tgrid)))
    else:
        tgrid=np.ma.MaskedArray(tgrid, mask=(tgrid==nodata))

    #print comp_geot
    xs = ds.RasterXSize
    ys = ds.RasterYSize
    zmin = tgrid.min()
    zmax = tgrid.max()

    xmax = xmin + (cellsize[0] * xs)
    ymin = ymax + (cellsize[1] * ys)
    return [xmin, xmax, ymin, ymax, zmin, zmax, xs, ys]

def returnNoData(ingrd):
    ds = gdal.Open(ingrd)
    return ds.GetRasterBand(1).GetNoDataValue()

def degree2radian(deg):
    return (deg/180.0) * math.pi

def radian2degree(rad):
    return (rad/math.pi) * 180.0

def latlonlen(lat):
    lonlen = math.cos(degree2radian(lat)) * gg_equat
    latlen = 1 * gg_equat
    return [lonlen, latlen]

def pre_perspective(ingrd, grdmm, cpt):
    gd_translate = "gdal_translate -ot UInt16 -of PNG -scale %s %s 0 65535 %s temp.png\n\
gdal_translate -srcwin 1 1 %s %s -of PNG temp.png dem_16bit.png\n\
gdaldem color-relief %s %s temp2.tif\n\
gdal_translate -srcwin 1 1 %s %s temp2.tif rgb.tif\n\
rm temp.* temp2.*\n\
convert -transparent \"rgb(0,0,0)\" rgb.tif rgb.png\n\
convert -size 10000 dem_16bit.png dem_16bit_10000.png\n\
convert -transparent \"rgb(0,0,0)\" -size 10000 rgb.png rgb_10000.png\n" %(grdmm[4], grdmm[5], ingrd, 
                                                                          grdmm[6]-2, grdmm[7]-2, ingrd, 
                                                                          cpt, grdmm[6]-2, grdmm[7]-2)
    return gd_translate

def povray_template(ingrd, grdmm, options):
    lllen = latlonlen(grdmm[2])
    povr_out = "%s.pov" %(ingrd)
    povr_file = open(povr_out, 'w')
    povr_file.write("// DEM\n\
\n\
//global_settings { assumed_gamma 2.2 }\n\
      \n\
#include \"colors.inc\"\n\
     \n\
#declare Bi = 2;\n\
\n\
//\n\
// Custom parameters start here\n\
//\n\
#declare rgb_image = \"rgb.png\"\n\
#declare dem_image = \"dem_16bit.png\"\n\
\n\
#declare xdim = %s;  //number of pixels in X-direction\n\
#declare ydim = %s;  //number of pixels in y-direction\n\
#declare max_y = %s; //maximum latitude extent\n\
#declare min_y = %s; //minimum latitude extent\n\
#declare min_z = %s; //minimum elevation\n\
#declare max_z = %s; //maximum elevation\n\
\n\
// Obtained from http://www.csgnetwork.com/degreelenllavcalc.html  \n\
#declare deg_lat_len = %s; //length of a degree of latitude in meters  \n\
#declare deg_lon_len = %s; //length of a degree of longitude in meters\n\
\n\
// Position of camera             \n\
#declare cam_azimuth = %s;\n\
#declare cam_elevation = %s;\n\
#declare cam_distance = %s; \n\
#declare cam_view_angle = %s;\n\
                     \n\
// Position of the \"sun\"  \n\
#declare light_azimuth = cam_azimuth+90;\n\
#declare light_elevation = %s;\n\
#declare light_distance = %s;             \n\
\n\
#declare vertical_exaggeration = %s;\n\
//\n\
// Custom parameters end here\n\
//\n\
\n\
#declare lon_scale = deg_lon_len / deg_lat_len;\n\
#declare z_scale = (100 * (max_z - min_z)) / (deg_lat_len * (max_y - min_y));\n\
\n\
#declare cam_x = cam_distance * cos(radians(cam_elevation)) * sin(radians(cam_azimuth));\n\
#declare cam_y = cam_distance * sin(radians(cam_elevation));\n\
#declare cam_z = cam_distance * cos(radians(cam_elevation)) * cos(radians(cam_azimuth));\n\
#declare light_x = light_distance * cos(radians(light_elevation)) * sin(radians(light_azimuth));\n\
#declare light_y = light_distance * sin(radians(light_elevation));\n\
#declare light_z = light_distance * cos(radians(light_elevation)) * cos(radians(light_azimuth));\n\
\n\
#declare Texture0 = // Color elevation image (24-bit RGB PNG)\n\
texture {\n\
  pigment{\n\
    image_map { \n\
      png rgb_image map_type 0 once interpolate Bi \n\
    } \n\
  } \n\
  finish { ambient 0.4 diffuse 0.8 } \n\
  rotate x*90  \n\
}\n\
\n\
\n\
height_field { // Unsigned 16-bit PNG DEM\n\
   png dem_image \n\
   smooth\n\
   clipped_by {box { <0, 0, 0>, <0.999, 1, 0.999> } }\n\
   texture { Texture0 }\n\
   translate <-0.5, 0, -0.5>\n\
   scale <100*lon_scale*xdim/ydim,\n\
          vertical_exaggeration*z_scale,  //Vertical exaggeration\n\
          100>\n\
} \n\
\n\
\n\
camera {\n\
   angle cam_view_angle\n\
   location <cam_x, cam_y, cam_z>\n\
   look_at <0, 0, 0> \n\
}\n\
\n\
light_source { <light_x, light_y, light_z> color White shadowless parallel }\n\
        \n\
background { White } \n\
" %(grdmm[6],grdmm[7],grdmm[3],
    grdmm[2],grdmm[4],grdmm[5],
    lllen[0],lllen[1],options[0],
    options[1],options[2],options[3],
    options[4],options[5],options[6]))
    povr_file.close()
#    print('\ngdal_perspective v.%s' %(gp_version))
def Usage():
    print('')
    print('gdal_perspective.py srcfile [-width n] [-height n]')
    print('                    [-azimuth n] [-elevation n]')
    print('                    [-distance n] [-view-angle n]')
    print('                    [-light-elevation n] [-light-distance n]')
    print('                    [-vertical-exaggeration n] [-cpt n] [-q]')
    print('')
    sys.exit( 1 )

if __name__ == "__main__":

    ingrd=None
    pov_wid=1000
    pov_height=800
    pov_azimuth=130
    pov_elevation=27
    pov_distance=235
    pov_view_angle=35
    pov_light_elevation=30
    pov_light_distance=10000
    pov_vertical_exag=2
    quiet=False
    user_cpt=None

    gdal.AllRegister()
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv is None:
        sys.exit(0)
        
    # Parse command line arguments.
    i = 1
    while i < len(argv):
        arg = argv[i]

        if arg == '-width':
            pov_wid = sys.argv[i+1]
            i = i + 1

        elif arg == '-height':
            pov_height = sys.argv[i+1]
            i = i + 1

        elif arg == '-azimuth':
            pov_azimuth = sys.argv[i+1]
            i = i + 1

        elif arg == '-view-angle':
            pov_view_angle = sys.argv[i+1]
            i = i + 1
            
        elif arg == '-elevation':
            pov_elevation = sys.argv[i+1]
            i = i + 1

        elif arg == '-distance':
            pov_distance = sys.argv[i+1]
            i = i + 1

        elif arg == '-light-elevation':
            pov_distance = sys.argv[i+1]
            i = i + 1

        elif arg == '-light-distance':
            pov_distance = sys.argv[i+1]
            i = i + 1

        elif arg == '-vertical-exaggeration':
            pov_vertical_exag = sys.argv[i+1]
            i = i + 1

        elif arg == '-cpt':
            user_cpt = sys.argv[i+1]
            i = i + 1

        elif arg == '-q':
            quiet = True

        elif arg[0] == '-':
            Usage()

        elif ingrd is None:
            ingrd = arg

        else:
            Usage()

        i = i + 1

    if ingrd == None:
        Usage()
        sys.exit(0)

    options=[pov_azimuth,pov_elevation,pov_distance,pov_view_angle,pov_light_elevation,pov_light_distance,pov_vertical_exag]
    if not quiet:
        print("gdal_perspective: Obtaining min/max values from grid")
    grdmm = returnMinMax(ingrd)
    if not quiet:
        print("gdal_perspective: Generating CPT file")
    if user_cpt is None:
        make_cpt(ingrd, grdmm)
        pp_sh = pre_perspective(ingrd,grdmm,"%s.cpt" %(ingrd))
    else:
        pp_sh = pre_perspective(ingrd,grdmm,user_cpt)
    povray_template(ingrd, grdmm, options)
    povr_sh = "povray %s.pov +W%s +H%s" %(ingrd,pov_wid,pov_height)
    if not quiet:
        print("gdal_perspective: processing grid file")
    os.system(pp_sh)
    #print pp_sh
    if not quiet:
        print("gdal_perspective: processing pov-ray file")
    os.system(povr_sh)
    #print povr_sh
    # Cleanup
    #os.system("rm rgb.* dem_*")
#--END

