#!/bin/sh
## ------------------------------------------------------------
### coastline2xyz.sh  
## Copyright (c) 2019 - 2020 Matthew Love <matthew.love@colorado.edu>
## This file is liscensed under the GPL v.2 or later and 
## is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details. 
## <http://www.gnu.org/licenses/> 
##--------------------------------------------------------------
##
### Code:

c2x_version="0.0.1"
usage="usage: ocastline2xyz.sh [options]
\n\n
Options:\n
-I\tinput coastline shapefile\n
-Z\toutput z value (default is 0)\n
\n
coastline2xyz.sh v.${c2x_version}\n"

## Getopts
while getopts ":I:Z:" options; do
    case $options in
	I ) in_coast=$OPTARG;;
	Z ) out_zed=$OPTARG;;
	h ) echo -e $usage;;
	\? ) echo -e $usage
	exit 1;;
	* ) echo -e $usage
	    exit 1;;
    esac
done

if [ ! $in_coast ] ; then
    echo -e $usage
    exit 1
fi

mkdir .c2x_tmp
cp $(basename $in_coast .shp).* .c2x_tmp
cd .c2x_tmp
rm $(basename $in_coast .shp).dbf
ogr2ogr tmp.csv $in_coast -f CSV -lco GEOMETRY=AS_WKT -overwrite
if [ ! $out_zed ] ; then
    cat tmp.csv | sed -e '1,1d' | tr ',' '\n' | sed 's/[A-Za-z"()]*//g' | tr ' ' ',' | sed 's/^,//' | awk -F, '{print $1,$2,$3}' > ../$(basename $in_coast .shp).xyz
else
    cat tmp.csv | sed -e '1,1d' | tr ',' '\n' | sed 's/[A-Za-z"()]*//g' | tr ' ' ',' | sed 's/^,//' | awk -v zed="$out_zed" -F, '{print $1,$2,zed}' > ../$(basename $in_coast .shp).xyz
fi
cd ..
rm -rf .c2x_tmp
### End
