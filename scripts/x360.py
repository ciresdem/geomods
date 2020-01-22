#!/usr/bin/env python2

import sys
import os
import re
try:
    import numpy as np
			
except ImportError:
    sys.exit('''Apparently, you do not have numpy installed correctly on your system.  Please install to proceed.''')

# =============================================================================
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
        self.fileob.seek(0)
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

if __name__ == '__main__':
        
        which_way=None
        
        try:
                infile=sys.argv[1]
                outfile=sys.argv[2]
        except:
                print "you must enter an in and outfile; x360.py input output"
                sys.exit(0)
        if len(sys.argv) > 3:
                which_way=sys.argv[3]

        inxyz=xyzFile(infile, " ", "xyz")
        outxyz=open(outfile,'w')
        this_line=inxyz.xyzreadline()
        inxyz.xyzReset()

        while this_line is not None:
                x,y,z=inxyz.getX(),inxyz.getY(),inxyz.getValue('z')
                #print x,y,z
                if which_way is None:
                        x2 = x-360
                else:
                        x2 = x+360
                #print x2
                outline = " ".join([str(x2),str(y),str(z)])
                outxyz.write(outline)
                outxyz.write("\n")
                this_line=inxyz.xyzreadline()	

        inxyz.xyzClose()
        outxyz.close()


