#!/usr/bin/env python
__description__ =\
"""
Parses the output from specific-packing.tcl and makes it pretty and R-readable.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "specific-packing.py vmd_output_file"

import sys

class SpecificPackingError(Exception):
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
    out.append(float(split_line[2]))
    out.append(float(split_line[3]))
    
    for x in split_line[4:]:
        cont = x.split("{")[-1]
        cont = cont.strip()[1:]
        
        num_atoms = len(cont.split())
        out.append(num_atoms)

        #if num_atoms > 0:
        #    out.append(1)
        #else:
        #    out.append(0)


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
        raise SpecificPackingError(err)

    # Read the input file
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    lines = [l for l in lines if l.startswith("->")]

    # Parse each line in the file
    out = []
    for l in lines:
        out.append(("%10i%10.3f%10.3f%10i%10i%10i%10i%10i%10i\n") % parseVMDLine(l))

    # Add line numbers and header
    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n" % 
               (" ","frame","sel1_sasa","sel2_sasa","c15","c25","c35","c45","c55","c65"))

    print "".join(out)

 
if __name__ == "__main__":
    
    main()     
