#!/usr/bin/env python
__description__ =\
"""
Parses the output from calc-distances.tcl and makes it pretty and R-readable.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "calc-distances.py vmd_output_file"

import sys

# MUST MATCH THE RESIDUES OF INTEREST IN CALC_DISTANCES.TCL
single_names = ["H3N","H3M","H3C","H5N","H5M","H5C","H7N","H7C","H10N","H10M",
                "H10C","H12N","H12C","LC3","LC17"]


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
    out.extend([float(v) for v in split_line[2:] if v.strip() != ""])

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
        v = parseVMDLine(l)

        out.append("%10i" % v[0])
        out.append(((len(v)-1)*"%10.3f") % tuple(v[1:]))
        out.append("\n")

    out = "".join(out)
    out = out.split("\n")

    pairs = []
    for i in range(len(single_names)):
        for j in range(i+1,len(single_names)):
            pairs.append("%s-%s" % (single_names[i],single_names[j]))

    # Add line numbers and header
    out = ["%10i%s\n" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%10s" % (" ","frame"))
    out.insert(1,(len(pairs)*"%10s") % tuple(pairs))
    out.insert(2,"\n")

    print "".join(out[:-1])

 
if __name__ == "__main__":
    
    main()     
