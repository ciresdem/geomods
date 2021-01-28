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
import setuptools

with open('README', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'geomods',
    version = '0.4.5',
    description = 'Modules and scripts for utilizing geographic data Digital Elevation Models',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    license = 'MIT',
    author = 'CIRES Coastal DEM Team',
    url = 'http://ciresgroups.colorado.edu/coastalDEM',
    packages = setuptools.find_packages(),#['geomods'],  #same as name
    package_data = {'geomods': ['data/*.gmt']},
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI APPROVED :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires = [
        'GDAL',
        'numpy',
        'scipy',
        'requests',
        'lxml',
        'matplotlib',
    ], 
    entry_points = {
        'console_scripts': [
            'waffles = geomods.waffles:waffles_cli',
            'datalists = geomods.datalists:datalists_cli',
            'fetches = geomods.fetches:fetches_cli',
        ],
    },
    scripts = [
        'scripts/gdal_chunk.py',
        'scripts/gdal_crop.py',
        'scripts/gdal_clip.py',
        'scripts/gdal_null.py',
        'scripts/gdal_mask.py',
        'scripts/gdal_split.py',
        'scripts/gdal_query.py',
        'scripts/gdal_gquery.py',
        'scripts/gdal_perspective.py',
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
        'scripts/create_datalist.sh',
        'scripts/bag2tif2chunks2xyz.sh',
        'scripts/gdal2xyzchunks.py',
        'scripts/error_distance_plots.py',
        'scripts/fetch_coastline.py',
    ],
    #python_requires = '>=2.7, <3',
    #project_urls = {
    #    'Source': 'https://github.com/ciresdem/geomods',
    #},
)
