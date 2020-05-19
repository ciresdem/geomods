#!/usr/bin/env python
##
import sys
import waffles

if __name__ == '__main__':
    '''adjust the x value in an xyz file by 360'''
    which_way = None
    try:
        infile = sys.argv[1]
        outfile = sys.argv[2]
    except:
        print "you must enter an in and outfile; x360.py input output"
        sys.exit(0)
        
    if len(sys.argv) > 3: which_way = sys.argv[3]
    with open(infile, 'r') as iob:
        with open(outfile, 'w') as oob:
            for xyz in waffles.xyz_parse(iob):
                if which_way is None: x = xyz[0] -= 360
                else: x = xyz[0] += 360
                waffles.xyz_line([x, xyz[1], xyz[2]], oob)

### End

