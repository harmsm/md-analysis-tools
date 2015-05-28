#!/usr/bin/env python
__description__ =\
"""
Calculate water residence times and spit out in R-readable format.  Analysis
assumes that same water never leaves and comes back; only time this will cause
a problem is if a long-residence water pops out, then happens to come back and
become a long-residence water again.  Otherwise, short-term associations will
have little effect on results.  Residence times are calculated in frames.
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "water-reisdence-time.py vmd_output_file"

import sys

class WaterResidenceTimeError(Exception):
    """
    General error class for this module.
    """
    
    pass

def waterResidenceTime(input_file,frame_cutoff=0):
    """
    Calculate the residence time for waters observed within a specific cutoff
    of the protein (recorded by water-residence-time.tcl).
    """

    # Read the input file
    f = open(input_file,'r')

    # Parse each line in the file.  For each water seen, update the 
    # water_counter dictionary.
    water_counter = {}
    for i, l in enumerate(f.readlines()):

        if not l.startswith("->"):
            continue

        entry = l.split("|")
        frame = int(entry[1])
        if frame < frame_cutoff:
            continue

        waters = entry[2].split()
        for w in waters:
            try:
                water_counter[w] += 1
            except KeyError:
                water_counter[w] = 1

    for w in waters:
        water_counter[w] = water_counter[w]

    f.close()
  
    all_counts = [(water_counter[k],k) for k in water_counter.keys()]
    all_counts.sort(reverse=True)

    out = ["%10s%12s%12s\n" % (" ","water","residence")]
    for i in range(len(all_counts)):
        out.append("%10i%12s%12i\n" % (i,all_counts[i][1],
                                         all_counts[i][0]))

    return out

def main(argv=None):
    """
    Function to call to run analysis.
    """

    # Parse command line
    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
    except IndexError:
        err = __usage__
        raise WaterResidenceTimeError(err)

    out = waterResidenceTime(input_file)

    print "".join(out)


 
if __name__ == "__main__":
    
    main()     
