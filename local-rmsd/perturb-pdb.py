#!/usr/bin/env python
__description__ = """
Take a pdb file and an R-formatted file containing the columns x_mean, y_mean,
and z_mean, then shift the coordinates of the residues in the pdb file by x_mean,
y_mean, and z_mean.  It loads the 10 x the total magnitude of the shift into the
b-factor column of the pdb.  The script is kinda dumb: it applies the 
perturbations sequentially to each residue in the pdb fie, so make sure that your
pdb and data file have exactly the same residues.  
"""
__author__ = "Michael J. Harms"
__date__ = "110520  (the day before the end of the world)"
__usage__ = "perturb-pdb.py pdb_file data_file"

import sys
from math import sqrt

class PerturbPdbError(Exception):
    """
    General error class for this module.
    """

    pass

def readDataFile(data_file):
    """
    Read an R-formatted data file and return a list of tuples containing x,y,z
    perturbtations.
    """
        
    fmt_err = """
    Data file does not have correct format!  Data file should be R-formatted
    and have x_mean, y_mean, and z_mean columns!
    """

    f = open(data_file,'r')
    lines = f.readlines()
    f.close()

    # Strip comments and blank lines
    lines = [l for l in lines if (not l.startswith("#")) and l.strip() != ""]

    # Figure out which columns to take
    columns = lines[0].split()
    try:
        x = columns.index("x_mean") + 1
        y = columns.index("y_mean") + 1
        z = columns.index("z_mean") + 1
    except ValueError:
        err = "\n\n%s\n\n" % fmt_err
        raise PerturbPdbError(err)

    # Take data
    out = []
    for l in lines[1:]:
        col = l.split()
        out.append([float(col[x]),float(col[y]),float(col[z])])
   
    return out 


def perturbPdb(pdb_file,data_file):
    """
    Take a pdb file and perturb its coordinates by the values in data file.
    """

    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    data = readDataFile(data_file)
    
    delta = [sqrt((d[0]**2+d[1]**2+d[2]**2)) for d in data]
    delta = [10*d/sum(delta) for d in delta]

    current_residue = None
    out = []
    for line in pdb:
        if line[0:6] == "ATOM  ":
            if current_residue == None:
                current_residue = line[21:26]
                residue_index = 0 
            
            if current_residue != line[21:26]:
                current_residue = line[21:26]
                residue_index += 1

            new_x = float(line[30:38]) + data[residue_index][0]
            new_y = float(line[38:46]) + data[residue_index][1]
            new_z = float(line[46:54]) + data[residue_index][2]
            new_b = delta[residue_index]

            new_line = "%s%8.3f%8.3f%8.3f%s%6.2f%s" % (line[:30],   
                                                       new_x,new_y,new_z,
                                                       line[54:60],
                                                       new_b,line[66:])
            out.append(new_line)
        
        else:
            out.append(line)

    return out


def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]
    
    try:
        pdb_file = argv[0]
        data_file = argv[1]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n" % __usage__
        raise PerturbPdbError(err)

    out = perturbPdb(pdb_file,data_file)

    return "".join(out)
    

if __name__ == "__main__":

    print main()
