#!/usr/bin/env python
"""
combine_gr.py

Generate R-readable output from VMD calculated radial distribution functions.
"""


__author__ = "Michael J. Harms"
__date__ = "070214"
__usage__ = "combine_gr.py inputfile OR directory with input files"

import sys, os
from numpy import mean, std

def parseFile(input_file):
    """
    Take the output from calc_gr.tcl and parse it, returning a list of series
    in the file.
    """

    # Read in file
    f = open(input_file,'r')
    contents = f.readlines()
    f.close()
    
    # Find line that contains data (it will all be on one line, with the series
    # separated by { SERIES }.
    start_index = contents.index("BEGIN\n") + 1
    data = contents[start_index]

    # Create list of series data
    series_list = data.split("{")
    series_list = ["".join([x for x in s if x != "}"]) for s in series_list]
    series_list = [[float(x) for x in s.split()] for s in series_list]

    # Strip blank first line and last "frames processed" entries
    series_list.pop(0)
    series_list.pop(-1)

    return series_list

def combineGr(file_list):
    """
    Combines output from all calc_gr.tcl output files in file_list.
    """

    # Parse all input files in file_list
    # all_series will have data format: 
    #   all_series[trajectory][series_type][radius]   
    all_series = [] 
    for f in file_list:
        all_series.append(parseFile(f))

    # Pull out ranges for trajectories, series, and radii
    trajectories = range(len(all_series))
    series = range(1,len(all_series[0]))
    radii = range(len(all_series[0][0]))

    # Create output file header
    out = []
    out.append(5*"%14s" % (" ","r","g_r","norm_hist","hist"))
    out.append("\n")

    # Append individual trajectory data
    index = 0
    for t in trajectories:
        for r in radii:
            out.append("%14i%14.6E" % (index,all_series[0][0][r]))
            for s in series:
                out.append("%14.6E" % all_series[t][s][r])
            out.append("\n")
            index += 1

    return out


def main():
    """
    Main function to call if run from command line.
    """

    # Pull input off of command line
    try:
        input = sys.argv[1:]
    except IndexError:
        print __usage__
        sys.exit()

    # Deal with all arguments on command line (single files or directories)
    file_list = []
    for single_input in input:
        if os.path.isfile(single_input):
            file_list.append(single_input)
        elif os.path.isdir(single_input):
            tmp_file_list = os.listdir(single_input)
            tmp_file_list = [os.path.join(single_input,f)
                             for f in tmp_file_list]
            tmp_file_list = [f for f in tmp_file_list if os.path.isfile(f)]
            file_list.extend(tmp_file_list)
        else:
            print __usage__
            sys.exit()
    
    # Combine files
    out = combineGr(file_list)

    # Dump output to std out
    print "".join(out)


# Call if run from command line
if __name__ == "__main__":
    main()
