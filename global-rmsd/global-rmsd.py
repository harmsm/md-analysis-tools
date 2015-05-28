#!/usr/bin/env python
__description__ =\
"""
Parses the output from calc-rmsd.tcl and makes it pretty and R-readable.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "calc-rsmd.py vmd_output_file"

import sys

HEADER = ["frame","ca_rmsd","lig_rmsd"]

class ParseVmdError(Exception):
    """
    General error class for this module.
    """
    
    pass

 
def parseVMDLine(line,fmt_string="%10.3f"):
    """
    Parse a line of output from VMD.
    """

    # Parse the line
    split_line = line.split("|")
  
    # Frame number 
    out = [int(split_line[1])]
    out.extend([float(v) for v in split_line[2:]])

    fmt = "%s%s\n" % ("%10i",(len(out[1:])*fmt_string))
    
    return fmt % tuple(out)


def main(argv=None):
    """
    Main function.  Parse command line, read file, return output string.  
    """

    if argv == None:
        argv = sys.argv[1:]

    # Parse command line
    try:
        input_file = argv[0]
    except IndexError:
        err = __usage__
        raise ParseVmdError(err)

    # Read the input file
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    # Grab appropriate lines
    lines = [l for l in lines if l.startswith("->")]

    # Parse each line in the file
    out = []
    for l in lines:
        out.append(parseVMDLine(l,"%10.3f"))


    # Add line numbers and header
    header = [" "]
    header.extend(HEADER)
    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%10s%10s%10s\n" % tuple(header))

    return "".join(out)

 
if __name__ == "__main__":
    
    print main()     
