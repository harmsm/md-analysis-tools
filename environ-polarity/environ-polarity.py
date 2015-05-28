#!/usr/bin/env python
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "080513"
__usage__ = "outputProcessor.py vmd_output_dir"

import sys, os
from math import sqrt

def mean(some_list):
    """
    Calculate the mean value of a list.
    """

    return float(sum(some_list))/len(some_list)

def stdev(some_list):
    """
    Calculate the standard deviation of a list.
    """

    m = mean(some_list)
    var = mean([(v - m)**2 for v in some_list])
    return sqrt(var)
    

def readFile(vmd_output):
    """
    """

    fmt = "%10.3F%10.3F%10.3F%10.3F%10.3F\n"

    # Read file
    f = open(vmd_output,'r')
    lines = f.readlines()
    f.close()

    # Find lines with actual data
    hash = [l[0:3] for l in lines]
    start_index = hash.index("BEG") + 1
    end_index = hash.index("END")
    lines = lines[start_index:end_index]

    # Find unique distance cutoffs
    columns = [l.split("\t") for l in lines] 

    cutoffs = [float(l[1]) for l in columns]
    cutoffs = dict([(c,[]) for c in cutoffs])
    cutoffs = cutoffs.keys()
    cutoffs.sort()
   
    # Extract mean and standard deviation for this trajectory 
    out = []
    for c in cutoffs:

        protein = [int(l[2]) for l in columns if float(l[1]) == c] 
        prot_mean = mean(protein)
        prot_stdev = stdev(protein)

        solvent = [int(l[3]) for l in columns if float(l[1]) == c] 
        solv_mean = mean(solvent)
        solv_stdev = stdev(solvent)

        out.append(fmt % (c,prot_mean,prot_stdev,solv_mean,solv_stdev))  

    return out


def main():
    """
    Function to call if run from command line.
    """

    # Parse command line.
    try:
        file = sys.argv[1]
    except IndexError:
        print __usage__
        sys.exit()

    out = readFile(file)

    # Add counter to each line of output and add header
    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]   
    out.insert(0,"%10s%10s%10s%10s%10s%10s\n" % (" ","cutoff",
                "prot_mean","prot_sd","solv_mean","solv_sd"))

    print "".join(out)

if __name__ == "__main__":
    main()
