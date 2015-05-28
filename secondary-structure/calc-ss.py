#!/usr/bin/env python
__description__ =\
"""
Parses the output from calc-ss.tcl and makes it pretty and R-readable.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "calc-ss.py vmd_output_file"

import sys

ss_equiv_dict = {"C":1, # coil
                 "T":2, # turn
                 "E":3, # extended
                 "B":3, # isolated bridge
                 "b":3, # isolated bridge
                 "H":4, # alpha helix
                 "G":4, # 310 helix
                 "I":4} # pi helix

class CalcSSError(Exception):
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
    out.extend([ss_equiv_dict[v] for v in split_line[2].split()])

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
        raise CalcSSError(err)

    # Read the input file
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    lines = [l for l in lines if l.startswith("->")]

    # Parse each line in the file
    out = []
    for l in lines:
        out.append(parseVMDLine(l))
    num_residues = len(out[0]) 
    fmt = num_residues*"%6i" 

    # Add line numbers 
    out = [fmt % o for o in out]
    out = ["%6i%s\n" % (i,o) for i, o in enumerate(out)]
    
    # Create header
    header_fmt = "%6s" + num_residues*"%6s" + "\n"
    header = ["r%i" % i for i in range(num_residues-1)]
    header.insert(0,"frame")
    header.insert(0," ")
    header = header_fmt % tuple(header)
    out.insert(0,header)
    
    print "".join(out)

 
if __name__ == "__main__":
    
    main()     
