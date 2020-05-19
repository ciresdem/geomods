### lasf.py
##
## Copyright (c) 2007 - 2020 CIRES Coastal DEM Team
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

import os
import sys
import struct
import time
import re
import numpy as np

_version = '0.1.5'

def las_file_p(inlas):
    '''Check quickly whether file is an LAS
    file.  The first four(4) bytes	of an
    LAS file should have the characters "LASF"'''

    return open(inlas, 'rb').read(4) == "LASF"

def myBaseName(inlas, extension):
    '''Return the basename of a file-string'''

    return inlas[:-len(extension)]


class LasHeader:

    def __init__(self, inlas="None"):
        self.inlas = inlas
        try:
            if os.path.exists(self.inlas):
                self.nostream = False
            else: self.nostream = True
        except:
            self.nostream = True

    ## ==============================================
    ## Header structs
    ## The format of these header structures is:
    ## (description, bytes, struct-type, array-size)
    ## 
    ## header_struct defenition originally written by Matthew Perry
    ## under the MIT License. (http://code.google.com/p/pylas/source/browse/trunk/pylas.py)
    ## ==============================================
    
    header_struct = (('filesig',4,'c',4),
                     ('filesourceid',2,'H',1),
                     ('reserved',2,'H',1),
                     ('guid1',4,'L',1),
                     ('guid2',2,'H',1),
                     ('guid3',2,'H',1),
                     ('guid4',8,'B',8),
                     ('vermajor',1,'B',1),
                     ('verminor',1,'B',1),
                     ('sysid',32,'c',32),
                     ('gensoftware',32,'c',32),
                     ('fileday',2,'H',1),
                     ('fileyear',2,'H',1),
                     ('headersize',2,'H',1),
                     ('offset',4,'L',1),
                     ('numvlrecords',4,'L',1),
                     ('pointformat',1,'B',1),
                     ('pointreclen',2,'H',1),
                     ('numptrecords',4,'L',1),
                     ('numptbyreturn',20,'L',5),
                     ('xscale',8,'d',1),
                     ('yscale',8,'d',1),
                     ('zscale',8,'d',1),
                     ('xoffset',8,'d',1),
                     ('yoffset',8,'d',1),
                     ('zoffset',8,'d',1),
                     ('xmax',8,'d',1),
                     ('xmin',8,'d',1),
                     ('ymax',8,'d',1),
                     ('ymin',8,'d',1),
                     ('zmax',8,'d',1),
                     ('zmin',8,'d',1))

    vlheader_struct = (('reserved',2,'H',1),
                       ('userid',16,'c',16),
                       ('recordid',2,'H',1),
                       ('reclen',2,'H',1),
                       ('desc',32,'c',32))

    ## ==============================================
    ## parseHeader function originally written by Matthew Perry
    ## under the MIT License. (http://code.google.com/p/pylas/source/browse/trunk/pylas.py)
    ## ==============================================

    def parseHeader(self):
        '''Parse the LAS Header information into a
        python dictionary.  See the header_struct
        defenition above for the names to use to
        call certain values.'''

        if self.nostream:
            return "You cannot parse what you do not have, try createHeader first..."

        ds = open(self.inlas, 'rb')
        lheader = {}
        for i in self.header_struct:
            try:
                if i[2] == 'c':
                    value = ds.read(i[1])
                elif i[3] > 1:
                    value = struct.unpack("=" + str(i[3]) + i[2], ds.read(i[1]))
                else:
                    value = struct.unpack("=" + i[2], ds.read(i[1]))[0]
            except:
                value = ds.read(i[1])
            lheader[i[0]] = value
        ds.close()
        return lheader

    def dumpHeader(self):
        '''Dump the LAS file public header to the screen'''
        lash = self.parseHeader()
        for i in self.header_struct:
            print i[0] + ":" + "\t", lash[i[0]]

    ## ==============================================
    ## Variable Length Header(s)
    ## ==============================================

    def getSeek(self):
        '''Retreive and return the propper seek
        position, which is just the value of the
        size of the header in question. - initial
        seek position, that is.'''
		
        return self.parseHeader()['headersize']

    def getVlNum(self):
        '''Return the number of variable length
        records present in LAS file.'''
		
        return self.parseHeader()['numvlrecords']

    def parseVlr(self, seeknum):
        '''Parse a variable length header from an
        LAS file.  The existence of any variable
        length headers should be mentioned in the
        main header of the LAS file.'''

        if self.nostream:
            return "You cannot parse what you do not have, try createHeader first..."
        
        ds = open(self.inlas,'rb')
        ds.seek(seeknum)
        vlhead = {}
        for i in self.vlheader_struct:
            if i[2] == 'c':
                value = ds.read(i[1])
            elif i[2] == 'H':
                value = struct.unpack("=" + i[2] , ds.read(i[1]) )[0]
            vlhead[i[0]] = value
        ds.close()
        return vlhead
    
    def getVlRecs(self):
        '''Parse variable length record from an
        LAS file, based on the information in
        the LAS variable length header.'''
		
        seeknum = self.getSeek()
        vlnum = self.getVlNum()
        vlrecs = []
        for i in range(0, vlnum):
            vlr = self.parseVlr(seeknum)
            rid = vlr['recordid']
            if len(str(rid)) < 4:
                rid = vlr['reserved']
            seeknum = seeknum + 54 + vlr['reclen']
            recs = (rid, vlr['reclen'], seeknum)
            vlrecs.append(recs)
        return vlrecs

    def readVlRec(self, seeknum, recid, reclen):
        ds = open(self.inlas, 'rb')
        # Las Spec Records
        if int(recid) == 0:
            '''classification lookup'''
        elif int(recid) == 2:
            '''histogram'''
        elif int(recid) == 3:
            ds.seek(seeknum+54)
            return ds.read(reclen)
        # Projection Info
        elif int(recid) == 34735:
            '''do something'''
        elif int(recid) == 34736:
            '''do something'''
        elif int(recid) == 34737:
            ds.seek(seeknum+54)
            return ds.read(reclen)
            
    def getScale(self):
        '''Retreive and return the x, y and z scale
        information from the LAS header.'''

        h = self.parseHeader()
        return h['xscale'], h['yscale'], h['zscale']

    def getOffsets(self):
        '''Return the x,y and z offsets from the header'''
        
        h = self.parseHeader()
        return h['xoffset'], h['yoffset'], h['zoffset']

    def getReturnGrp(self, instream):
        '''Return the return group information.
        This is used in header creation.'''

        firstret = 0
        secondret = 0
        thirdret = 0
        fourthret = 0
        fifthret = 0
        for row in instream:
            return_grp = row['r']
            return_num = return_grp & 7
            if return_num == 1 or return_num == 0:
                firstret += 1
            elif return_num == 2:
                secondret +=1
            elif return_num == 3:
                thirdret += 1
            elif return_num == 4:
                fourthret += 1
            elif return_num == 5:
                fifthret += 1
				
        return firstret,secondret,thirdret,fourthret,fifthret

    def getPtStats(self, returns=False):
        '''Return some statistics about the point-
        cloud, for header creation.'''

        if self.nostream:
            instream = self.inlas
        else:
            instream = LasFile(self.inlas).readAsArray()
        maxx = instream['x'].max()
        minx = instream['x'].min()
        maxy = instream['y'].max()
        miny = instream['y'].min()
        maxz = instream['z'].max()
        minz = instream['z'].min()
        tot = len(instream)
        if returns:
            firstret,secondret,thirdret,fourthret,fifthret = self.getReturnGrp(instream)
        else:
            firstret,secondret,thirdret,fourthret,fifthret = 0,0,0,0,0
        return maxx,maxy,maxz,minx,miny,minz,tot,\
               firstret,secondret,thirdret,fourthret,fifthret

    def getPtStatsAsList(self, returns=False, reclen=10):
        '''Return some statistics about the point-
        cloud, for header creation.'''

        if self.nostream:
            instream = self.inlas
        else:
            instream = LasFile(inlas).readAsList()
        maxx = 0
        maxy = 0
        minx = 0
        miny = 0
        maxz = 0
        minz = 0
        tot = 0

        for i in range(0, len(instream)):
            for j in range(0, len(instream[i]), reclen):
                if instream[i][j] > maxx:
                    maxx = instream[i][j]
                elif instream[i][j] < minx:
                    minx = instream[i][j]
                if instream[i][j+1] > maxy:
                    maxy = instream[i][j+1]
                elif instream[i][j+1] < miny:
                    miny = instream[i][j+1]
                if instream[i][j+2] > maxz:
                    maxz = instream[i][j+2]
                elif instream[i][j+2] < minz:
                    minz = instream[i][j+2]
                tot += 1

        if returns:
            firstret,secondret,thirdret,fourthret,fifthret = self.getReturnGrp(instream)
        else:
            firstret,secondret,thirdret,fourthret,fifthret = 0,0,0,0,0
        return maxx,maxy,maxz,minx,miny,minz,tot,\
               firstret,secondret,thirdret,fourthret,fifthret

    def createHeader(self, returns=False, scale=[0.0001,0.0001,0.0001], offsets=[0,0,0], reclen=10):
        '''Create the las header, which will output
        as version 1.1 whether the input LAS file
        is 1.0 or 1.1, other versions are not yet
        supported (v1.2 & v1.3dev).  If returns equals
        False, then the return numbers wont be calculated
        for the header.'''

        lheader = {}
        fday = time.localtime()[7]
        fyear = time.localtime()[0]
        lstats = self.getPtStatsAsList(returns)
        h_size = 227
        h_offset = 229

        if reclen == 10:
            ptformat = 1
            ptlen = 28
        else:
            ptformat = 0
            ptlen = 28

        values = ["LASF", 0, 0, 0, 0, 0, (0, 0, 0, 0, 0, 0, 0, 0), 1, 1, \
                  "Python " + str(sys.version[0:5]), "ML-LAS " + _version, \
                  fday, fyear, h_size, h_offset, 0, ptformat, ptlen, lstats[6], \
                  (lstats[7], lstats[8], lstats[9], lstats[10], lstats[11]), \
                  scale[0], scale[1], scale[2], offsets[0], offsets[1], offsets[2], \
                  lstats[0] * scale[0], lstats[3] * scale[0], lstats[1] * scale[1], \
                  lstats[4] * scale[1], lstats[2] * scale[2], lstats[5] * scale[2]]

        j = 0
        for i in values:
            lheader[self.header_struct[j][0]] = i
            j += 1

        return lheader

    def createNpHeader(self, returns=False, scale=[0.0001,0.0001,0.0001], offsets=[0,0,0]):
        '''Create an LAS File header as a numpy array'''

        dt = "a4, u2, u2, u4, u2, u2, (8,)u1, u1, u1, a32, a32, u2, \
        u2, u2, u4, u4, u1, u2, u4, (5,)u4, f8, f8, f8, f8, f8, f8, \
        f8, f8, f8, f8, f8, f8"
		
        fday = time.localtime()[7]
        fyear = time.localtime()[0]
        lstats = self.getPtStats()
        h_size = 227
        h_offset = 229
        ptformat = 1
        ptlen = 28
        
        if self.nostream is False:
            offsets = self.getOffsets()
            scale = self.getScale()
        

        values = ("LASF", 0, 0, 0, 0, 0, (0, 0, 0, 0, 0, 0, 0, 0), 1, 1, \
                  "Python " + str(sys.version[0:5]), "ML-LAS " + _version, \
                  fday, fyear, h_size, h_offset, 0, ptformat, ptlen, lstats[6], \
                  (lstats[7], lstats[8], lstats[9], lstats[10], lstats[11]), \
                  scale[0], scale[1], scale[2], offsets[0], offsets[1], offsets[2], \
                  lstats[0] * scale[0], lstats[3] * scale[0], lstats[1] * scale[1], \
                  lstats[4] * scale[1], lstats[2] * scale[2], lstats[5] * scale[2])
        
        return np.array(values, dt)

    def createVlrHeader(self, type="proj", proj_string="NIL"):
        '''Create a variable-length header for
        an LAS file in version 1.1 of the standard.'''
		
        vheader = {}
        j = 0
        if type == "proj":
            rec_id = 34735
            user_id = "LASF_Projection"
            reclen = 0
            desc = proj_string
        elif type == "class":
            rec_id = 0
            user_id = "LASF_Spec"
            desc = "LASF Classifications"
            reclen = 0
            values = [0, user_id, rec_id, reclen, desc]
        for i in values:
            vheader[vlheader_struct[j][0]] = i
            j += 1
        return vheader

    def createVlrRecord(self):
        '''do something'''

class LasFile:

    def __init__(self, inlas):

        if not isLas:
            sys.exit("This does not appear to be a proper LAS file. Please check your source file.")
        
        self.inlas = inlas
        self.hc = LasHeader(inlas)
        self.h = self.hc.parseHeader()
        
        self.bsize = 1000
        self.bnames = ('x','y','z','i','r','c','s','u','p','g')
        self.btypes = ('i4','i4','i4','u2','u1','u1','u1','u1','u2','f8')

        if self.h['pointformat'] == 1:
            self.reclen = 10
        elif self.h['pointformat'] == 0:
            self.reclen = 9
            self.bnames = self.bnames[:-1]
            self.btypes = self.btypes[:-1]
        else:
            sys.exit("Support for point-format %s is not yet supported...") % (str(h['pointformat']))

    def validPoint(self, pt_rec, use_recs="xyzi", clf=None):
        '''Validate a point record in terms of user
        preference, for use in dumpPoints function.
        use_recs is a string containing any of the
        following letters: "xyziecaupg". clf is the
        classifcation field to use, if set, otherwise
        will validate any class records'''

        point = {}
        outpoint = []

        scale = self.hc.getScale()
        offsets = self.hc.getOffsets()
        pt_form = tuple(re.findall(r'(\D)', use_recs))

        for form in pt_form:
            if form == 'x':
                outpoint.append(pt_rec[0] * scale[0] + offsets[0])
            elif form == 'y':
                outpoint.append(pt_rec[1] * scale[1] + offsets[1])
            elif form == 'z':
                outpoint.append(pt_rec[2] * scale[2] + offsets[2])
            else:
                outpoint.append(pt_rec[list(self.bnames).index(form)])

        return outpoint

    def dumpPoints(self, use_recs="xyz", clf=None, delim=","):

        if clf is not None:
            clfs = clf.split(",")
        
        las = open(self.inlas, 'rb')
        las.seek(self.h['offset'])

        BUFF = self.h['pointreclen']

        BUFFERSIZE = BUFF * self.bsize
        
        if BUFF == 28:
            UNPACKFORMAT = 'lllHbBBBHd'
        if BUFF == 20:
            UNPACKFORMAT = 'lllHbBBBH'

        unpackstring = "=" + UNPACKFORMAT

        
        #print >> sys.stderr, BUFF
        #while True:
        for i in range(self.h['numptrecords']):
            thisP = struct.unpack(unpackstring, las.read(BUFF))

            if clf is not None and str(thisP[5]) in clfs:
                print str(delim).join(map(str, self.validPoint(thisP, use_recs, clf)))
            elif clf is None:
                print str(delim).join(map(str, self.validPoint(thisP, use_recs, clf)))

    def readAsArray(self, israw=False):
        '''Read a binary LAS-files point records into a
        NumPy Rec-Array. '''

        lasf = open(self.inlas, 'rb')
        lasf.seek(self.h['offset'])
        scale = self.hc.getScale()
        offsets = self.hc.getOffsets()

        u = np.fromfile(lasf, {'names': self.bnames, 'formats': self.btypes})

        if not israw:
            self.rawbtypes = list(self.btypes)
            self.rawbtypes[0] = 'f8'
            self.rawbtypes[1] = 'f8'
            self.rawbtypes[2] = 'f8'
            self.rawbtypes = tuple(self.rawbtypes)
            uf = np.zeros(len(u), {'names': self.bnames, 'formats': self.rawbtypes})
            for i in self.bnames:
                if i == 'x':
                    uf[i] = u[i] * scale[0] + offsets[0]
                elif i == 'y':
                    uf[i] = u[i] * scale[1] + offsets[1]
                elif i == 'z':
                    uf[i] = u[i] * scale[2] + offsets[2]
                elif i == 'r':
                    uf[i] = u[i] & 7
                else:
                    uf[i] = u[i]
            u = uf
            
        lasf.close()
        return u     

    def readAsList(self):
        '''Load an LAS file into a python list.  This is
        the old default function, now depreciated, though
        it`s abilities are still in working order.'''
		
        las = open(self.inlas, 'rb')
        u = []
        las.seek(self.h['offset'])
        if self.h['pointformat'] == 1:
            BUFFERSIZE = 28 * self.bsize
            UNPACKFORMAT = 'lllHbBBBHd'
        elif self.h['pointformat'] == 0:
            BUFFERSIZE = 20 * self.bsize
            UNPACKFORMAT = 'lllHbBBBH'
        pointlen = self.h['pointreclen']
        while True:
            try:
                data = las.read(BUFFERSIZE)
            except:
                break
            if not data:
                break
            unpackstring = "=" + UNPACKFORMAT * (len(data) / pointlen)
            u.append(struct.unpack(unpackstring, data))
        las.close()
        return u

    def readAsSubSet(self, col="c", clf=2):
        '''Return a subset of the las point array, based
        on the field (col) and the field-value (clf).'''

        clfs = clf.split(",")

        lasp = self.readAsArray()

        if len(clfs) > 1:
            classindex = np.empty(0, 'i4')
            for clf in clfs:
                classindex_tmp = np.where(lasp[col] == int(clf))
                classindex = np.append(classindex, classindex_tmp)
            classindex.sort()
        else:
            classindex = np.where(lasp[col] == int(clf))[0]
        
        newArray = np.zeros(len(classindex), {'names': self.bnames, 'formats': self.rawbtypes})
        for num,item in enumerate(classindex):
            newArray[num] = lasp[item]
        return newArray

    def writeFile(self, inarray, outfile):
        '''Write an LAS file to disk'''
        
        if isnp:
            LasWrite(inarray, outfile).writeLasFile()
        else:
            LasWrite(inarray, outfile).writeLasFileFromList()

class LasWrite:

    def __init__(self, inarray, outfile):

        self.inarray = inarray
        try:
            self.inarray.shape
            self.isnp = True
        except:
            self.isnp = False
        if self.isnp:
            self.inarrayRowlen = len(self.inarray[0])
        else:
            self.inarrayRowlen = 10
        self.outlas = open(outfile, 'wb')
        self.lh = LasHeader(self.inarray)

    def writeHeader(self, xyzscale=[0.01,0.01,0.01], xyzoffsets=[0,0,0], reclen=10):
        h = LasHeader(self.inarray).createHeader(False, xyzscale, xyzoffsets, reclen)
        values = []
        for i in self.lh.header_struct:
            try:
                charfinder = h[i[0]][0] + 1
                for j in h[i[0]]:
                    values.append(j)
            except:
                values.append(h[i[0]])
        FORMATSTRING = '= 4s H H L H H 8B B B 32s \
        32s H H H L L B H L 5L d d d d d d d d d d d d'
        
        s = struct.Struct(FORMATSTRING)
        packed_data = s.pack(*values)
        self.outlas.write(packed_data)
        self.outlas.write(struct.pack('2x'))

    def writePointRecords(self, nstream):
        PACKFORMAT = 'lllHbBBBHd'
        if self.isnp:
            nstream = self.inarray.tolist()
        for i in nstream:
            FORMATSTRING =  "=" + (PACKFORMAT * (len(i) / self.inarrayRowlen))
            s = struct.Struct(FORMATSTRING)
            packed_records = s.pack(*i)
            self.outlas.write(packed_records)

    def writeLasFile(self, xyzscale=[0.01,0.01,0.01], xyzoffsets=[0,0,0]):
        '''Write an LAS binary file based on the
        NumPy Rec-Array of the LAS point cloud '''

        harray = self.lh.createNpHeader(False,xyzscale,xyzoffsets)
        harray.tofile(self.outlas)
        self.outlas.write(struct.pack('2x'))
        self.inarray.tofile(self.outlas)
        self.outlas.close()
		
    def writeLasFileFromList(self):
        
        self.writeHeader()
        self.writePointRecords(self.inarray)
        self.outlas.close()

    def lhead2nphead(self, lhead):
        '''Convert a header dictionary to a numpy array. '''

        dt = "a4, u2, u2, u4, u2, u2, (8,)u1, u1, u1, a32, a32, u2, \
        u2, u2, u4, u4, u1, u2, u4, (5,)u4, f8, f8, f8, f8, f8, f8, \
        f8, f8, f8, f8, f8, f8"
        nh = []
        for i in header_struct:
            nh.append(lhead[i[0]])
        return np.array(tuple(nh), dt)

### End
