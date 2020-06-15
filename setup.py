### setup.py
##
## Copyright (c) 2020 CIRES Coastal DEM Team
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
from setuptools import setup

setup(
    name = 'geomods',
    version = '0.2.0',
    description = 'Modules and scripts for utilizing geographic data Digital Elevation Models',
    license = 'MIT',
    author = 'CIRES Coastal DEM Team',
    url = 'http://ciresgroups.colorado.edu/coastalDEM',
    packages = ['geomods'],  #same as name
    package_data = {'geomods': ['data/*.gmt']},
    install_requires = [
        'GDAL',
        'numpy',
        'requests',
        'lxml',
    ], 
    entry_points = {
        'console_scripts': [
            'waffles = geomods.waffles:waffles_cli',
            'datalists = geomods.waffles:datalists_cli',
            'fetches = geomods.fetches:main',
        ],
    },
    scripts = [
        'scripts/gdal_chunk.py',
        'scripts/gdal_crop.py',
        'scripts/gdal_clip.py',
        'scripts/gdal_null.py',
        'scripts/gdal_mask.py',
        'scripts/gdal_split.py',
        'scripts/vdatum_cmd.py',
        'scripts/smooth_dem_bathy.py',
        'scripts/xyz_clip.py',
        'scripts/xyz2shp.py',
        'scripts/xyz_histogram.py',
        'scripts/x360.py',
	'scripts/gdal_findreplace.py',
        'scripts/spatial-meta.sh',
        'scripts/clip_xyz.sh',
        'scripts/coastline2xyz.sh',
    ],
    #python_requires = '>=2.7, <3',
    #project_urls = {
    #    'Source': 'https://github.com/ciresdem/geomods',
    #},
)
