#!/usr/bin/env python
__description__ = \
"""
Parse the output from diehdral-recorder.tcl and make it pretty and R-readable.
"""
__author__ = "Michael J. Harms"
__usage__ = "parse-dihedral.py diehdral-output"
__date__ = "110413"

import sys

class WaterError(Exception):
    """
    General error class for this module.
    """
    
    pass

def parseWater(input_file):
    """
    Go through a vmd output file and grab the specified dihedrals.
    """ 

    f = open(input_file,'r')
   
    out = ["%12s%12s%12s\n" % (" ","frame","inter_hoh")] 
    for l in f.readlines():
        if not l.startswith("->"):
            continue

        entry = l.split("|")

        out.append("%12i%12i" % (int(entry[1]),int(entry[1])))
        num_internal = len(entry[5].split())
        out.append("%12i\n" % num_internal)

    f.close()

    return out

def main(argv=None):
    """
    Main function!
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n" % __usage__
        raise WaterError(err)

    out = parseWater(input_file)

    print "".join(out)


if __name__ == "__main__":
    main()
