#!/usr/bin/env python
__description__ = \
"""
Parse the output from diehdral-recorder.tcl and make it pretty and R-readable.
"""
__author__ = "Michael J. Harms"
__usage__ = "parse-dihedral.py diehdral-output"
__date__ = "110413"

import sys

class DihedralError(Exception):
    """
    General error class for this module.
    """
    
    pass

def parseDihedral(input_file,place_in_bins=True):
    """
    Go through a vmd output file and grab the specified dihedrals.
    """ 

    f = open(input_file,'r')
   
    out = [] 
    for l in f.readlines():
        if l.startswith("-*"):
            column_names = l.split("|")[1:]

        if not l.startswith("->"):
            continue

        entry = l.split("|")

        out.append("%12i%12i" % (int(entry[1]),int(entry[1])))

        if place_in_bins:

            for e in entry[2:]:
                dihed = float(e) % 360
                if dihed >= 0 and dihed < 120:
                    out.append("%12i" % 0)
                elif dihed >= 120 and dihed < 240:
                    out.append("%12i" % 1)
                else:
                    out.append("%12i" % 2)
        else:
            out.extend(["%12.3f" % (float(e) % 360) for e in entry[2:]])
        out.append("\n")

    f.close()

    header = "%12s%12s" % (" ","frame")
    header = header + "".join(["%12s" % c for c in column_names]) + "\n"

    out.insert(0,header)

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
        raise DihedralError(err)

    out = parseDihedral(input_file,place_in_bins=False)

    print "".join(out)


if __name__ == "__main__":
    main()
