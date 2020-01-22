#!/usr/bin/python
#
#  Description: Generate a histogram
#
# 	$Id: xyz_histogram.py,v 1.2 2011/03/07 21:02:15 mlove Exp mlove $	
#--

#--
import sys
import os
import re
import numpy as np
from math import *
#--

def checkColumnNumType(fileob, delim):
	
    colNumType = []
    cols = fileob.readline().split(delim)
    for i in xrange(0, len(cols)):
        try: 
            float(cols[i])
            try:
                int(cols[i])
                colgrp = 2
            except:
                colgrp = 1
        except:
            colgrp = 0
        colNumType.append(colgrp)
    #fileob.seek(0)
    return colNumType

# =============================================================================XYZ-File Class

class xyzFile:

    # =============================================================================
    def __init__(self, infile, delim=",", pt_form="xyz"):
        self.infile = infile
        self.delim = delim
        if infile == sys.stdin:
            self.fileob = infile
        else:
            self.fileob = open(infile, 'r')
        self.fiet = 0

        self.bnames = tuple(re.findall(r'(\D)', pt_form))

        colnum = checkColumnNumType(self.fileob, self.delim)
        if len(self.bnames) > len(colnum):
            self.bnames = self.bnames[:-(len(self.bnames) - len(colnum))]
            
        self.btypes = []
        for i in colnum:
            if i == 1:
                self.btypes.append('f8')
            elif i == 0:
                self.btypes.append('S16')
            elif i == 2:
                self.btypes.append('i4')
        
        self.btypes = tuple(self.btypes)
        if len(self.btypes) > len(self.bnames):
            self.btypes = self.btypes[:-(len(self.btypes) - len(self.bnames))]

    def xyzClose(self):
        self.fileob.close()

    def getXyzFileOb(self):
        return self.fileob

    def getX(self):
        xpos = list(self.bnames).index('x')
        if self.fiet == 0:
            self.xyzreadline()
        elif self.fiet is None:
            return None
        return float(self.fiet[xpos])
    
    def getY(self):
        ypos = list(self.bnames).index('y')
        if self.fiet == 0:
            self.xyzreadline()
        elif self.fiet is None:
            return None
        return float(self.fiet[ypos])

    def getValue(self, value):
        vpos = list(self.bnames).index(str(value))
        if self.fiet == 0:
            self.xyzreadline()
        elif self.fiet is None:
            return None
        return float(self.fiet[vpos])

    def xyzReset(self):
        self.fileob.seek(0)

    def xyzSeek(self, val):
        self.fileob.seek(val)

    def xyzreadline(self):
        field = self.fileob.readline()
        try:
            self.fiet = field.strip().split(self.delim)
        except:
            print >> sys.stderr, "xyzlib: Failed to parse line."
        if self.fiet == ['']:
            return None
        else:
            return self.fiet

    def asciiRow(self):
        '''Return the first row of an ascii file for
        inspection.'''
        
        fi = open(self.infile, 'r')
        field = fi.readline()
        fi.close()
        return field.strip().split(self.delim)
    
    def getPtPositions(self, innames):
        '''return a list of the position of points in an ascii file.'''
        bn_positions = []
        for i in innames:
            bn_positions.append([j for j,x in enumerate(innames) if x == i])

        for i in range(0, len(bn_positions)):
            if bn_positions[i] == []:
                bn_positions[i] = 99
            else:
                bn_positions[i] = bn_positions[i][0]
				
        return bn_positions

    def readAsArray(self):
        '''Read an xyz ascii files point records into a
        NumPy Rec-Array. '''
        
        bn_positions = self.getPtPositions(self.bnames)

        cast = np.cast
        data = [[] for dummy in xrange(len(self.btypes))]
        fi = open(self.infile, 'r')
        j=0
        for line in fi:
            field = line.strip().split(self.delim)
            fields = []

            for i in bn_positions:
                if i == 99:
                    fields.append(0)
                else:
                    fields.append(field[i])

            for i, number in enumerate(fields):
                data[i].append(number)
            print >> sys.stderr, "xyzlib: " + str(j), "\r",
            j += 1
        fi.close()
        for i in xrange(len(self.btypes)):
           data[i] = cast[self.btypes[i]](data[i])
           
        return np.rec.array(data, {'names': self.bnames, 'formats': self.btypes})
            
# =============================================================================END


#--
def conDec(x, dec):
    '''Return a float string with n decimals
    (used for ascii output).'''
    
    fstr = "%." + str(dec) + "f"
    return fstr % x
#--

#--
def my_median(xarray):
    xarray.sort()
    m_num = xarray.size / 2
    if m_num % 2 == 0:
        return xarray[m_num]
    else:
        return np.mean(xarray[m_num-1:m_num+1])
#--

#--
def quartile(xarray):
    xarray.sort()
    lower,upper = np.array_split(xarray, 2)
    lower = lower[:-1]
    return np.median(lower),np.median(upper)
#--

#--
def kurtosis(xarray, hmean):
    mean_dev = xarray - float(hmean) #np.mean(xarray)
    #md4 = mean_dev**4
    md4st = mean_dev/np.std(mean_dev)
    md4m = np.mean(md4st**4)
    return md4m-3
#--

#--
def skewness(xarray, hmean):
    mean_dev = xarray - float(hmean) #np.mean(xarray)
    print mean_dev.shape
    md3 = mean_dev**3
    md3sum = np.sum(md3)
    return md3sum / ((md3.size-1) * np.std(xarray)**3)
#--

#--
def rmseme(xarray):
    dz2 = xarray**2
    hmean2 = np.mean(dz2)
    rmse = sqrt(float(hmean2))
    return rmse
#--

#--
def genHistogram(infile, outfile, delim=" ", p_form='xyz', vcol='z',
                 quickly=False, gf="/usr/share/texmf-dist/fonts/truetype/public/libertine/fxlr.ttf",
                 use_stats=False, gloc=False, bcolor="purple", verbose=False):
    if gloc is True:
        grphval = 0.05
    else:
        grphval = 0.65

    input = xyzFile(infile, delim, p_form)

    np_zs = np.array([])
    this_rec = input.xyzreadline()
    j = 0

    if quickly:
        try:
            xyzf = input.readAsArray()
            np_zs = xyzf[vcol]
        except:
            print >> sys.stderr, "xyz-histogram: Error\n" 
            quickly = None
            
    if quickly is None:
        while this_rec is not None:
            j += 1
            if verbose:
                print >> sys.stderr, "xyz-histogram: " + str(j), "\r",
            this_val = input.getValue(vcol)
            np_zs = np.append(np_zs, this_val)
            this_rec = input.xyzreadline()
    
    hbin,hist = np.histogram(np_zs,100,new=True)

    hmax = conDec(np_zs.max(), 5)
    hmin = conDec(np_zs.min(), 5)
    
    tmax = hmax
    tmin = hmin

    hist_file = "hist_" + str(infile) + ".txt"
    tmp_file = open(hist_file, "w")
    for i,j in enumerate(hbin):
        tmp_val = str(hist[i]) + " " + str(j) + "\n"
        tmp_file.write( tmp_val )
    tmp_file.close()

    tmp_gp = open("_tmp_gp.gnu", "w")

    ## Generate the GnuPlot script
    if use_stats is True:
        
        htot = np_zs.size
        hmean = conDec(np.mean(np_zs), 5)
        hstd = conDec(np.std(np_zs), 5)
        hmed = conDec(np.median(np_zs), 5)
        hskew = conDec(skewness(np_zs, hmean), 5)
        hkurt = conDec(kurtosis(np_zs, hmean), 5)
        hlowerq,hupperq = quartile(np_zs)
        rmse = rmseme(np_zs)

        gplot_file = "#!/usr/bin/gnuplot \
        \nset label 1 at graph " + str(grphval) + ", graph 0.95 \nset label 1 \"Total: " + str(htot) + "\" \
        \nset label 2 at graph " + str(grphval) + ", graph 0.9 \nset label 2 \"Minimum: " + hmin + "\" \
        \nset label 3 at graph " + str(grphval) + ", graph 0.85 \nset label 3 \"Maximum: " + hmax + "\" \
        \nset label 4 at graph " + str(grphval) + ", graph 0.8 \nset label 4 \"Mean: " + hmean + "\" \
        \nset label 5 at graph " + str(grphval) + ", graph 0.75 \nset label 5 \"Median: " + hmed + "\" \
        \nset label 6 at graph " + str(grphval) + ", graph 0.7 \nset label 6 \"Std. Dev.: " + hstd + "\" \
        \nset label 7 at graph " + str(grphval) + ", graph 0.65 \nset label 7 \"Skewness: " + hskew + "\" \
        \nset label 8 at graph " + str(grphval) + ", graph 0.6 \nset label 8 \"Kurtosis: " + hkurt + "\" \
        \nset label 9 at graph " + str(grphval) + ", graph 0.55 \nset label 9 \"1st Quartile: " + conDec(hlowerq, 5) + "\" \
        \nset label 10 at graph " + str(grphval) + ", graph 0.5 \nset label 10 \"3rd Quartile: " + conDec(hupperq, 5) + "\" \
        \nset label 11 at graph " + str(grphval) + ", graph 0.45 \nset label 11 \"RMSE: " + conDec(rmse, 5) + "\" \
        \nset xrange[" + str(float(tmin)-0.5) + ":" + str(float(tmax)+0.5) + "] \
        \nset ylabel \"Frequency\" \nset term png font \"" + gf + "\" 10\
        \nset style histogram cluster gap 1 \
        \nset style fill solid 0.7 border -1 \
        \nset style line 5 lt rgb \"" + bcolor + "\" \
        \nset output \"" + outfile + "\" \nplot \"" + hist_file + "\" w boxes ls 5 notitle"

        
        gplot_file = "#!/usr/bin/gnuplot \
        \nset label 1 at graph " + str(grphval) + ", graph 0.95 \nset label 1 \"Total: " + str(htot) + "\" \
        \nset label 2 at graph " + str(grphval) + ", graph 0.9 \nset label 2 \"Minimum: " + hmin + "\" \
        \nset label 3 at graph " + str(grphval) + ", graph 0.85 \nset label 3 \"Maximum: " + hmax + "\" \
        \nset label 4 at graph " + str(grphval) + ", graph 0.8 \nset label 4 \"Mean: " + hmean + "\" \
        \nset label 5 at graph " + str(grphval) + ", graph 0.75 \nset label 5 \"Median: " + hmed + "\" \
        \nset label 6 at graph " + str(grphval) + ", graph 0.7 \nset label 6 \"Std. Dev.: " + hstd + "\" \
        \nset label 7 at graph " + str(grphval) + ", graph 0.65 \nset label 7 \"Skewness: " + hskew + "\" \
        \nset label 8 at graph " + str(grphval) + ", graph 0.6 \nset label 8 \"Kurtosis: " + hkurt + "\" \
        \nset label 9 at graph " + str(grphval) + ", graph 0.55 \nset label 9 \"1st Quartile: " + conDec(hlowerq, 5) + "\" \
        \nset label 10 at graph " + str(grphval) + ", graph 0.5 \nset label 10 \"3rd Quartile: " + conDec(hupperq, 5) + "\" \
        \nset label 11 at graph " + str(grphval) + ", graph 0.45 \nset label 11 \"RMSE: " + conDec(rmse, 5) + "\" \
        \nset xrange[" + str(float(tmin)-0.5) + ":" + str(float(tmax)+0.5) + "] \
        \nset ylabel \"Frequency\" \nset term png font \"" + gf + "\" 10\
        \nset style histogram cluster gap 1 \
        \nset style fill solid 0.7 border -1 \
        \nset style line 5 lt rgb \"" + bcolor + "\" \
        \nset output \"" + outfile + "\" \nplot \"" + hist_file + "\" w boxes ls 5 notitle"
        # boxes ls 5 notitle"
    else:
        gplot_file = "#!/usr/bin/gnuplot \
        \nset xrange[" + str(float(tmin)-0.5) + ":" + str(float(tmax)+0.5) + "] \
        \nset ylabel \"Frequency\" \nset term png font \"" + gf + "\" 10\
        \nset style histogram cluster \
        \nset style fill solid border -1 \
        \nset style line 5 lt rgb \"" + bcolor + "\" \
        \nset output \"" + outfile + "\" \nplot \"" + hist_file + "\" w boxes ls 5 notitle"
        # boxes ls 5 notitle"

    tmp_gp.write(gplot_file)
    tmp_gp.close()
    import socket
    if socket.gethostname() == "indian":
        gplot_cmd = "/nfs/mgg_apps/MMT/bin/gnuplot _tmp_gp.gnu"
    else:
        gplot_cmd = "gnuplot _tmp_gp.gnu"
    os.system(gplot_cmd)
#--

def Usage():
    print('xyz_histogram.py [-delimiter char] [-p_form string] [-column string] [-font path]')
    print('                 [-bar_color color] [-statistics] [-stat_location]')
    print('                 [-quickly] [-verbose] output [input]')
    print('')
    print('-delimiter       The delimiter separating the input ascii columns')
    print('-p_form          The format of the input ascii data as a string: e.g. xyzgd')
    print('-column          The column of the input ascii data to graph, should be a character')
    print('                 corresponding to one of the characters mentioned in p_form.')
    print('-font            The path to the font to use, default uses liberation sans regular')
    print('-bar_color       The color of the histogram bars (can be hex/html/color name).')
    print('-statistics      If set, will print statistics on the histogram.')
    print('-stat_location   If set, will move the statistics from the right side of the graph')
    print('                 to the left side.')
    print('-verbose         Increase verbosity')
    print('-quickly         Uses xyzlib.py to load the input ascii data, speeding up performance')
    print('                 while using more memory.')
    print('output           The output png image of the plotted histogram.')
    print('input            The input ascii table data.')
    print('')

#--
#
# Mainline
#
#--
if __name__ == "__main__":

    delim=" "
    zloc="z"
    p_form="xyz"
    gf="/usr/share/fonts/liberation/LiberationSans-Regular.ttf"
    bcolor="purple"
    ust=False
    gloc=False
    npd=None
    verbose=False
    input=None
    output=None

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg == '-help':
            Usage()
            sys.exit(0)

        if arg == '-delimiter':
            delim = str(sys.argv[i+1])
            i = i + 1

        elif arg == '-p_form':
            p_form = str(sys.argv[i+1])
            i = i + 1

        elif arg == '-column':
            zloc = str(sys.argv[i+1])
            i = i + 1

        elif arg == '-font':
            gf = str(sys.argv[i+1])
            i = i + 1

        elif arg == '-bar_color':
            bcolor = sys.argv[i+1]
            i = i + 1
            
        elif arg == '-statistics':
            ust = True

        elif arg == '-stat_location':
            gloc = True            

        elif arg == '-quickly':
            npd = True

        elif arg == '-verbose':
            verbose = True

        elif output is None:
            output = arg

        elif input is None:
            input = arg

        elif arg[0] == '-':
            Usage()

        else:
            Usage()

        i = i + 1    

    infile,outfile = input,output

    if infile is None:
        infile = sys.stdin

    if outfile is None:
        Usage()
        sys.exit(0)

    genHistogram(infile, outfile, delim, p_form, zloc, npd, gf, ust, gloc, bcolor, verbose)
#--END
