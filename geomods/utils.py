### utils.py
##
## Copyright (c) 2010 - 2021 CIRES Coastal DEM Team
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
## nothing in utils shall depend on other geomods modules...
##
### Code:
import os
import sys
import time
import subprocess
import glob
import math
import zipfile
import gzip
import datetime

## ==============================================
## import gdal/numpy
## ==============================================
import ogr
import osr
import gdal
import numpy as np

## ==============================================
##
## General utility functions - utils.py
##
## ==============================================
_version = '0.4.1'

def inc2str_inc(inc):
    """convert a WGS84 geographic increment to a str_inc (e.g. 0.0000925 ==> `13`)

    Args:
      inc (float): a gridding increment

    Returns:
      str: a str representation of float(inc)
    """
    
    import fractions
    return(str(fractions.Fraction(str(inc * 3600)).limit_denominator(10)).replace('/', ''))

def this_date():
    """get current data

    Returns:
      str: the current date
    """
    
    return(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))

def this_year():
    """get the current year
    
    Returns
      str: the current year
    """
    
    return(datetime.datetime.now().strftime('%Y'))

def rm_f(f_str):
    """os.remove f_str, pass if error

    Args:
      f_str (str): a pathname string

    Returns:
      int: 0
    """
    
    try:
        if os.path.exists(f_str):
            os.remove(f_str)
    except: pass
    return(0)
        
def remove_glob(*args):
    """glob `glob_str` and os.remove results, pass if error
    
    Args:
      *args (str): any number of pathname or dirname strings

    Returns:
      int: 0
    """
    
    for glob_str in args:
        try:
            globs = glob.glob(glob_str)
            for g in globs:
                if os.path.isdir(g):
                    remove_globs('{}/*'.format(g))
                    os.removedirs(g)
                else: os.remove(g)
        except: pass
    return(0)

def base_name(instr, extension):
    """get the filename basename

    Args:
      instr (str): the pathname string
      extension (str): the pathname extension

    Returns:
      str: the basename of a file-string
    """
    
    return(instr[:-len(extension)])

def args2dict(args, dict_args={}):
    """convert list of arg strings to dict.
    
    Args:
      args (list): a list of ['key=val'] pairs
      dict_args (dict): a dict to append to

    Returns:
      dict: a dictionary of the key/values
    """
    
    for arg in args:
        p_arg = arg.split('=')
        dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else True if p_arg[1].lower() == 'true' else None if p_arg[1].lower() == 'none' else p_arg[1]
    return(dict_args)

def int_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(int(val))
    except: return(or_val)

def float_or(val, or_val = None):
    """return val if val is integer

    Args:
      val (?): input value to test
      or_val (?): value to return if val is not an int

    Returns:
      ?: val as int otherwise returns or_val
    """
    
    try:
        return(float(val))
    except: return(or_val)

def gdal_fext(src_drv_name):
    """find the common file extention given a GDAL driver name
    older versions of gdal can't do this, so fallback to known standards.

    Args:
      src_drv_name (str): a source GDAL driver name

    Returns:
      list: a list of known file extentions or None
    """
    
    fexts = None
    try:
        drv = gdal.GetDriverByName(src_drv_name)
        if drv.GetMetadataItem(gdal.DCAP_RASTER): fexts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        if fexts is not None: return(fexts.split()[0])
        else: return(None)
    except:
        if src_drv_name.lower() == 'gtiff': fext = 'tif'
        elif src_drv_name == 'HFA': fext = 'img'
        elif src_drv_name == 'GMT': fext = 'grd'
        elif src_drv_name.lower() == 'netcdf': fext = 'nc'
        else: fext = 'gdal'
        return(fext)
    
## ==============================================
##
## Archives (zip/gzip/etc.)
##
## ==============================================
def _clean_zips(zip_files):
    """remove all files\directories in `zip_files`
    
    Args:
      zip_files (list): a list of files pathname strings.

    Returns:
      int: 0
    """
    
    for i in zip_files:
        if os.path.isfile(i):
            os.remove(i)
            zip_files = [x for x in zip_files if x != i]
    if len(zip_files) > 0:
        for i in zip_files:
            if os.path.isdir(i):
                try:
                    os.removedirs(i)
                except: pass
    return(0)

def unzip(zip_file):
    """unzip (extract) `zip_file`

    Args:
      zip_file (str): a zip file pathname string

    Returns:
      list: a list of extracted file names.
    """
    
    zip_ref = zipfile.ZipFile(zip_file)
    zip_files = zip_ref.namelist()
    zip_ref.extractall()
    zip_ref.close()
    return(zip_files)

def gunzip(gz_file):
    """gunzip `gz_file`

    Args:
      gz_file (str): a gzip file pathname string.

    Returns:
      str: the extracted file name.
    """
    
    if os.path.exists(gz_file):
        gz_split = gz_file.split('.')[:-1]
        guz_file = '{}.{}'.format(gz_split[0], gz_split[1])
        with gzip.open(gz_file, 'rb') as in_gz, \
             open(guz_file, 'wb') as f:
            while True:
                block = in_gz.read(65536)
                if not block:
                    break
                else: f.write(block)
    else:
        echo_error_msg('{} does not exist'.format(gz_file))
        guz_file = None
    return(guz_file)

def p_unzip(src_file, exts=None):
    """unzip/gunzip src_file based on `exts`
    
    Args:
      src_file (str): a zip/gzip filename string
      exts (list): a list of extensions to extract

    Returns:
      list: a list of the extracted files
    """
    
    src_procs = []
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for ext in exts:
                for zf in zfs:
                    if ext == zf.split('.')[-1]:
                        #if ext in zf:
                        src_procs.append(os.path.basename(zf))
                        with open(os.path.basename(zf), 'wb') as f:
                            f.write(z.read(zf))
    elif src_file.split('.')[-1] == 'gz':
        tmp_proc = gunzip(src_file)
        if tmp_proc is not None:
            for ext in exts:
                if ext == tmp_proc.split('.')[-1]:
                    src_procs.append(os.path.basename(tmp_proc))
                    break
                else: remove_glob(tmp_proc)
    else:
        for ext in exts:
            if ext == src_file.split('.')[-1]:
                src_procs.append(src_file)
                break
    return(src_procs)
    
def procs_unzip(src_file, exts):
    """unzip/gunzip src_file based on `exts`

    this function is depreciated, use p_unzip instead.

    Args:
      src_file (str): a zip/gzip filename string
      exts (list): a list of extensions to extract

    Returns:
      list: a list of the extracted files
    """
    
    zips = []
    src_proc = None
    if src_file.split('.')[-1].lower() == 'zip':
        with zipfile.ZipFile(src_file) as z:
            zfs = z.namelist()
            for ext in exts:
                for zf in zfs:
                    if ext in zf:
                        #src_proc = os.path.join(os.path.dirname(src_file), zf)
                        src_proc = os.path.basename(zf)
                        with open(src_proc, 'wb') as f:
                            f.write(z.read(zf))
                        break
    elif src_file.split('.')[-1] == 'gz':
        tmp_proc = gunzip(src_file)
        if tmp_proc is not None:
            for ext in exts:
                if ext in tmp_proc:
                    src_proc = os.path.basename(tmp_proc)
                    os.rename(tmp_proc, src_proc)
                    break
    else:
        for ext in exts:
            if ext in src_file:
                src_proc = src_file
                break
    return([src_proc, zips])

## ==============================================
##
## spatial and raster functions
##
## ==============================================
def euc_dst(pnt0, pnt1):
    """return the distance between pnt0 and pnt1,
    using the euclidean formula.

    `pnts` are geographic and result is in meters.

    Args:
      pnt0 (list): an xyz data list
      pnt1 (list): an xyz data list

    Returns:
      float: the distance beteween pnt0 and pnt1
    """
    
    rad_m = 637100
    distance = math.sqrt(sum([(a-b) ** 2 for a, b in zip(pnt0, pnt1)]))
    return(rad_m * distance)
    
def hav_dst(pnt0, pnt1):
    """return the distance between pnt0 and pnt1,
    using the haversine formula.

    `pnts` are geographic and result is in meters.

    Args:
      pnt0 (list): an xyz data list
      pnt1 (list): an xyz data list

    Returns:
      float: the distance beteween pnt0 and pnt1
    """
    
    x0 = float(pnt0[0])
    y0 = float(pnt0[1])
    x1 = float(pnt1[0])
    y1 = float(pnt1[1])
    rad_m = 637100
    dx = math.radians(x1 - x0)
    dy = math.radians(y1 - y0)
    a = math.sin(dx / 2) * math.sin(dx / 2) + math.cos(math.radians(x0)) * math.cos(math.radians(x1)) * math.sin(dy / 2) * math.sin(dy / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return(rad_m * c)
    
def _geo2pixel(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    if geoTransform[2] + geoTransform[4] == 0:
        pixel_x = ((geo_x-geoTransform[0]) / geoTransform[1])#w + .5
        pixel_y = ((geo_y-geoTransform[3]) / geoTransform[5])# + .5
    else: pixel_x, pixel_y = _apply_gt(geo_x, geo_y, _invert_gt(geoTransform))
    return(int(pixel_x), int(pixel_y))

def _geo2pixel_affine(geo_x, geo_y, geoTransform):
    """convert a geographic x,y value to a pixel location of geoTransform

    note: use _geo2pixel instead

    Args:
      geo_x (float): geographic x coordinate
      geo_y (float): geographic y coordinate
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a list of the pixel values [pixel-x, pixel-y]
    """
    
    import affine
    forward_transform = affine.Affine.from_gdal(*geoTransform)
    reverse_transform = ~forward_transform
    pixel_x, pixel_y = reverse_transform * (geo_x, geo_y)
    pixel_x, pixel_y = int(pixel_x + 0.5), int(pixel_y + 0.5)
    return(pixel_x, pixel_y)

def _pixel2geo(pixel_x, pixel_y, geoTransform):
    """convert a pixel location to geographic coordinates given geoTransform

    Args:
      pixel_x (int): the x pixel value
      pixel_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster
    
    Returns:
      list: [geographic-x, geographic-y]
    """
    
    geo_x, geo_y = _apply_gt(pixel_x, pixel_y, geoTransform)
    return(geo_x, geo_y)

def _apply_gt(in_x, in_y, geoTransform):
    """apply geotransform to in_x,in_y
    
    Args:
      in_x (int): the x pixel value
      in_y (int): the y pixel value
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: [geographic-x, geographic-y]
    """
    
    out_x = geoTransform[0] + (int(in_x + 0.5)*geoTransform[1]) + (int(in_y + 0.5)*geoTransform[2])
    out_y = geoTransform[3] + (int(in_x + 0.5)*geoTransform[4]) + (int(in_y + 0.5)*geoTransform[5])

    return(out_x, out_y)

def _invert_gt(geoTransform):
    """invert the geotransform
    
    Args:
      geoTransform (list): a geo-transform list describing a raster

    Returns:
      list: a geo-transform list describing a raster
    """
    
    det = (geoTransform[1]*geoTransform[5]) - (geoTransform[2]*geoTransform[4])
    if abs(det) < 0.000000000000001: return
    invDet = 1.0 / det
    outGeoTransform = [0, 0, 0, 0, 0, 0]
    outGeoTransform[1] = geoTransform[5] * invDet
    outGeoTransform[4] = -geoTransform[4] * invDet
    outGeoTransform[2] = -geoTransform[2] * invDet
    outGeoTransfrom[5] = geoTransform[1] * invDet
    outGeoTransform[0] = (geoTransform[2] * geoTransform[3] - geoTransform[0] * geoTransform[5]) * invDet
    outGeoTransform[3] = (-geoTransform[1] * geoTransform[3] + geoTransform[0] * geoTransform[4]) * invDet
    return(outGeoTransform)

def geoms_intersect_p(geom_a, geom_b):
    """check if OGR geometries `geom_a` and `geom_b` intersect

    Args:
      geom_a (ogr-geom): an ogr geometry
      geom_b (ogr-geom): an ogr geometry

    Returns:
      bool: True for intersection
    """
    
    if geom_a is not None and geom_b is not None:
        if geom_a.Intersects(geom_b):
            return(True)
        else: return(False)
    else: return(True)

def create_wkt_polygon(coords, xpos=1, ypos=0):
    """convert coords to Wkt

    Args:
      coords (list): x/y geographic coords
      xpos (int): the position of the x value in coords
      ypos (int): the position of the y value in corrds

    Returns:
      wkt: polygon as wkt
    """
    
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for coord in coords: ring.AddPoint(coord[xpos], coord[ypos])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    poly_wkt = poly.ExportToWkt()
    poly = None
    return(poly_wkt)

def wkt2geom(wkt):
    """transform a wkt to an ogr geometry
    
    Args:
      wkt (wkt): a wkt geometry

    Returns:
      ogr-geom: the ogr geometry
    """
    
    return(ogr.CreateGeometryFromWkt(wkt))

def sr_wkt(epsg, esri = False):
    """convert an epsg code to wkt

    Returns:
      (str): wkt or None
    """
    
    try:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(epsg))
        if esri: sr.MorphToESRI()
        return(sr.ExportToWkt())
    except: return(None)

## ==============================================
##
## error plots and calculations
##
## ==============================================
def err_fit_plot(xdata, ydata, out, fitfunc, dst_name='unc', xa='distance'):
    """plot a best fit plot with matplotlib
    
    Args:
      xdata (list): list of x-axis data
      ydata (list): list of y-axis data

    """
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        plt.plot(xdata, ydata, 'o')
        plt.plot(xdata, fitfunc(out, xdata), '-')
        plt.xlabel(xa)
        plt.ylabel('error (m)')
        out_png = '{}_bf.png'.format(dst_name)
        plt.savefig(out_png)
        plt.close()
        
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

def err_scatter_plot(error_arr, dist_arr, dst_name='unc', xa='distance'):
    """plot a scatter plot with matplotlib
    
    Args:
      error_arr (array): an array of errors
      dist_arr (array): an array of distances

    """
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText

        plt.scatter(dist_arr, error_arr)
        #plt.title('Scatter')
        plt.xlabel(xa)
        plt.ylabel('error (m)')
        out_png = '{}_scatter.png'.format(dst_name)
        plt.savefig(out_png)
        plt.close()
        
    except: echo_error_msg('you need to install matplotlib to run uncertainty plots...')

def err2coeff(err_arr, coeff_guess=[0, 0.1, 0.2], dst_name='unc', xa='distance'):
    """calculate and plot the error coefficient given err_arr which is 
    a 2 col array with `err dist
    
    Args:
      error_arr (array): an array of errors and distances

    Returns:
      list: [coefficient-list]
    """
    
    from scipy import optimize
    error = err_arr[:,0]
    distance = err_arr[:,1]
    
    max_int_dist = np.max(distance)
    nbins = 10
    n, _ = np.histogram(distance, bins = nbins)
    while 0 in n:
        nbins -= 1
        n, _ = np.histogram(distance, bins=nbins)
    serror, _ = np.histogram(distance, bins=nbins, weights=error)
    serror2, _ = np.histogram(distance, bins=nbins, weights=error**2)
    mean = serror / n
    std = np.sqrt(serror2 / n - mean * mean)
    ydata = np.insert(std, 0, 0)
    bins_orig=(_[1:] + _[:-1]) / 2
    xdata = np.insert(bins_orig, 0, 0)
    xdata[xdata - 0 < 0.0001] = 0.0001
    fitfunc = lambda p, x: p[0] + p[1] * (x ** p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x)
    out, cov, infodict, mesg, ier = optimize.leastsq(errfunc, coeff_guess, args=(xdata, ydata), full_output=True)
    try:
        err_fit_plot(xdata, ydata, out, fitfunc, dst_name, xa)
        err_scatter_plot(error, distance, dst_name, xa)
    except: echo_error_msg('unable to generate error plots, please check configs.')
    return(out)

## ==============================================
##
## system cmd verification and configs.
##
## run_cmd - run a system command as a subprocess
## and return the return code and output
##
## yield_cmd - run a system command as a subprocess
## and yield the output
##
## ==============================================
cmd_exists = lambda x: any(os.access(os.path.join(path, x), os.X_OK) for path in os.environ['PATH'].split(os.pathsep))

def run_cmd(cmd, data_fun=None, verbose=False):
    """Run a system command while optionally passing data.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    Args:
      cmd (str): a system command to run
      data_fun (lambda): a lambda function taking an output port as arg.
      verbose (bool): increase verbosity

    Returns:
      list: [command-output, command-return-code]
    """
    
    if verbose: _prog = _progress('running cmd: {}...'.format(cmd.rstrip()))
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    if verbose:
        p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True)
    else: p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)

    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()
    
    while p.poll() is None:
        if verbose:
            time.sleep(2)
            _prog.update()

    out = p.stdout.read()
    if not verbose: p.stderr.close()
    p.stdout.close()
    if verbose: _prog.end(p.returncode, 'ran cmd: {}... and returned {}.'.format(cmd.rstrip(), p.returncode))
    return(out, p.returncode)

def yield_cmd(cmd, data_fun=None, verbose=False):
    """Run a system command while optionally passing data.

    `data_fun` should be a function to write to a file-port:
    >> data_fun = lambda p: datalist_dump(wg, dst_port = p, ...)

    Args:
      cmd (str): a system command to run
      data_fun (lambda): a lambda function taking an output port as arg.
      verbose (bool): increase verbosity

    Yields:
      str: each line of output from the cmd
    """
    
    if verbose: echo_msg('running cmd: {}...'.format(cmd.rstrip()))    
    if data_fun is not None:
        pipe_stdin = subprocess.PIPE
    else: pipe_stdin = None
    p = subprocess.Popen(cmd, shell=True, stdin=pipe_stdin, stdout=subprocess.PIPE, close_fds=True)
    
    if data_fun is not None:
        if verbose: echo_msg('piping data to cmd subprocess...')
        data_fun(p.stdin)
        p.stdin.close()

    while p.poll() is None:
        line = p.stdout.readline().decode('utf-8')
        if not line: break
        else: yield(line)
    line = p.stdout.read().decode('utf-8')
    p.stdout.close()
    if verbose: echo_msg('ran cmd: {} and returned {}.'.format(cmd.rstrip(), p.returncode))

def cmd_check(cmd_str, cmd_vers_str):
    """check system for availability of 'cmd_str' 

    Args:
      cmd_str (str): a command string to run
      cmd_vers_str (str): a command that returns a version

    Returns:
      str: the commands version or None
    """
    
    if cmd_exists(cmd_str): 
        cmd_vers, status = run_cmd('{}'.format(cmd_vers_str))
        return(cmd_vers.rstrip())
    else: return("0")

def config_check(chk_vdatum=False, verbose=False):
    """check for needed waffles external software.

    waffles external software: gdal, gmt, mbsystem
    also checks python version and host OS and 
    records waffles version

    Args:
      chk_vdatum (bool): check for vdatum
      verbose (bool): increase verbosity

    Returns:
      dict: a dictionary of gathered results.
    """
    
    _waff_co = {}
    py_vers = str(sys.version_info[0]),
    host_os = sys.platform
    _waff_co['platform'] = host_os
    _waff_co['python'] = py_vers[0]
    ae = '.exe' if host_os == 'win32' else ''

    #if chk_vdatum: _waff_co['VDATUM'] = vdatum(verbose=verbose).vdatum_path
    _waff_co['GDAL'] = cmd_check('gdal_grid{}'.format(ae), 'gdal-config --version').decode()
    _waff_co['GMT'] = cmd_check('gmt{}'.format(ae), 'gmt --version').decode()
    _waff_co['MBGRID'] = cmd_check('mbgrid{}'.format(ae), 'mbgrid -version 2>&1 | grep Version').decode()
    _waff_co['LASTOOLS'] = cmd_check('las2txt{}'.format(ae), 'las2txt -version 2>&1 | awk \'{print $5}\'').decode()
    _waff_co['GEOMODS'] = str(_version)
    return(_waff_co)
    
## ==============================================
##
## stderr messaging and simple progress indicator
##
## ==============================================
def echo_warning_msg2(msg, prefix='waffles'):
    """echo warning msg to stderr using `prefix`

    >> echo_warning_msg2('message', 'test')
    test: warning, message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1mwarining\033[m, {}\n'.format(prefix, msg))

def echo_error_msg2(msg, prefix='waffles'):
    """echo error msg to stderr using `prefix`

    >> echo_error_msg2('message', 'test')
    test: error, message

    Args:
      msg (str): a message
      prefix (str): a prefix for the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: \033[31m\033[1merror\033[m, {}\n'.format(prefix, msg))

def echo_msg2(msg, prefix='waffles', nl=True):
    """echo `msg` to stderr using `prefix`

    >> echo_msg2('message', 'test')
    test: message
    
    Args:
      msg (str): a message
      prefix (str): a prefix for the message
      nl (bool): append a newline to the message
    """
    
    sys.stderr.flush()
    sys.stderr.write('\x1b[2K\r')
    sys.stderr.write('{}: {}{}'.format(prefix, msg, '\n' if nl else ''))
    sys.stderr.flush()

## ==============================================
## echo message `m` to sys.stderr using
## auto-generated prefix
## lambda runs: echo_msg2(m, prefix = os.path.basename(sys.argv[0]))
## ==============================================
echo_msg = lambda m: echo_msg2(m, prefix=os.path.basename(sys.argv[0]))
echo_msg_inline = lambda m: echo_msg2(m, prefix=os.path.basename(sys.argv[0]), nl = False)

## ==============================================
## echo error message `m` to sys.stderr using
## auto-generated prefix
## ==============================================
echo_error_msg = lambda m: echo_error_msg2(m, prefix=os.path.basename(sys.argv[0]))
echo_warning_msg = lambda m: echo_warning_msg2(m, prefix=os.path.basename(sys.argv[0]))

class _progress:
    """geomods minimal progress indicator"""

    def __init__(self, message=None):
        self.tw = 7
        self.count = 0
        self.pc = self.count % self.tw
        self.opm = message
        self.add_one = lambda x: x + 1
        self.sub_one = lambda x: x - 1
        self.spin_way = self.add_one
        self.spinner = ['*     ', '**    ', '***   ', ' ***  ', '  *** ', '   ***', '    **', '     *']
        
        self.perc = lambda p: ((p[0]/p[1]) * 100.)
        
        if self.opm is not None:
            self._clear_stderr()
            sys.stderr.write('\r {}  {:40}\n'.format(" " * (self.tw - 1), self.opm))
        
    def _switch_way(self):
        self.spin_way = self.sub_one if self.spin_way == self.add_one else self.add_one

    def _clear_stderr(self, slen = 79):
        sys.stderr.write('\x1b[2K\r')
        sys.stderr.flush()

    def update_perc(self, p, msg=None):
        if len(p) == 2:
            self._clear_stderr()
            sys.stderr.write('\r[\033[36m{:^5.2f}%\033[m] {:40}\r'.format(self.perc(p), msg if msg is not None else self.opm))
        else: self.update()
        
    def update(self, msg=None):

        self.pc = (self.count % self.tw)
        self.sc = (self.count % (self.tw + 1))
            
        self._clear_stderr()
        sys.stderr.write('\r[\033[36m{:6}\033[m] {:40}\r'.format(self.spinner[self.sc], msg if msg is not None else self.opm))
        
        if self.count == self.tw: self.spin_way = self.sub_one
        if self.count == 0: self.spin_way = self.add_one
        self.count = self.spin_way(self.count)

    def end(self, status, end_msg=None):
        self._clear_stderr()
        if end_msg is None: end_msg = self.opm
        if status != 0:
            sys.stderr.write('\r[\033[31m\033[1m{:^6}\033[m] {:40}\n'.format('fail', end_msg))
        else: sys.stderr.write('\r[\033[32m\033[1m{:^6}\033[m] {:40}\n'.format('ok', end_msg))

### End
