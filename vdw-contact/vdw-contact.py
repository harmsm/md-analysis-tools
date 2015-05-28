#!/usr/bin/env python
__description__ =\
"""
Parses the output from count-contacts.tcl and makes it pretty and R-readable.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "count-contacts.py vmd_output_file"

import sys

class CountContactsError(Exception):
    """
    General error class for this module.
    """
    
    pass

 
def parseVMDLine(line):
    """
    """

    # Parse the line
    split_line = line.split("|")
  
    # Frame number
    out = [int(split_line[1])]
    out.extend([int(v) for v in split_line[2:]])

    return tuple(out)


def main(argv=None):
    """
    """

    # Parse command line
    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
    except IndexError:
        err = __usage__
        raise CountContactsError(err)

    # Read the input file
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    lines = [l for l in lines if l.startswith("->")]

    # Parse each line in the file
    out = []
    for l in lines:
        out.append(("%10i%10i%10i%10i%10i%10i\n") % parseVMDLine(l))

    # Add line numbers and header
    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%10s%10s%10s%10s%10s%10s\n" % 
               (" ","frame","vdw200","vdw225","vdw250","vdw275","vdw300"))

    print "".join(out)

 
if __name__ == "__main__":
    
    main()     
