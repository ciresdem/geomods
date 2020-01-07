### setup.py
##
## Copyright (c) 2020 Matthew Love <matthew.love@colorado.edu>
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
    name='geomods',
    version='1.0',
    description='A useful module',
    license='MIT',
    author='Matthew Love',
    author_email='matthew.love@colorado.edu',
    url='http://ciresgroups.colorado.edu/coastalDEM',
    packages=['geomods'],  #same as name
    package_data={'geomods': ['data/*.gmt']},
    #install_requires=['numpy', 'gdal', 'ogr'], #external packages as dependencies
    scripts=[
        'scripts/fetch.py',
        'scripts/gdal_crop.py' ]
)
