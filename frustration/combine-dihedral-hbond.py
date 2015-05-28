#!/usr/bin/env python
__description__ = \
"""
This determines the hydrogen bonds made by a set of important atoms over a set
of trajectories.  It first reads all frames from all trajectories to generate 
a list of all populated hydrogen bonds.  It then goes through the trajectories
again and records, for a given trajectory, which hydrogen bonds are populated.

Three inputs are required: 
    1. output from hbond-recorder.tcl
    2. output from water-recorder.tcl
    3. gro file used for these scripts

It records unique donor-acceptor pairs for atoms specified in 
global_important_atoms.  It also records which class any waters involved in 
hydrogen bonds fall into.  
"""
__author__ = "Michael J. Harms"
__usage__ = "./combind-dihedral-hbond.py [USER-SETTABLE PARAMS AT TOP OF SCRIPT]"
__date__ = "100426"

import sys, os

# -------------------- USER-SETTABLE PARAMETERS ----------------------------- #

global_important_atoms = [("NH1","82ARG"),
                          ("NH2","82ARG"),
                          ("SD" ,"75MET"),
                          ("OE1","41GLN"),
                          ("NE2","41GLN"),
                          ("OE1","41GLU"),
                          ("OE2","41GLU"),
                          ("O3" ,"249N06"),
                          ("O3" ,"249N15"),
                          ("O3" ,"249NPR")]

atom_aliases = {("NH1","82ARG"): ("NH","82ARG"),
                ("NH2","82ARG"): ("NH","82ARG"),
                ("OE1","41GLU"): ("OE","41GLU"),
                ("OE2","41GLU"): ("OE","41GLU"),
                ("OE1","11GLU"): ("OE","11GLU"),
                ("OE2","11GLU"): ("OE","11GLU"),
                ("OE1","13GLU"): ("OE","13GLU"),
                ("OE2","13GLU"): ("OE","13GLU"),
                ("O3" ,"249N06"):("O3","249LIG"),
                ("O3" ,"249N15"):("O3","249LIG"),
                ("O3" ,"249NPR"):("O3","249LIG")}
        
water_type_list = ["ante","lower","upper","inter"]
    
file_list = ["SR2-Q41E-M75L_N0600_run000",
             "SR2-Q41E-M75L_N0600_run001",
             "SR2-Q41E-M75L_N0600_run002",
             "SR2_N0600_run000",
             "SR2_N0600_run001", 
             "SR2_N0600_run002", 
             "SR2_N0600_run003", 
             "SR2_N0600_run004", 
             "SR2-Q41E-M75L_N1500_run000",
             "SR2-Q41E-M75L_N1500_run001",
             "SR2-Q41E-M75L_N1500_run002",
             "SR2_N1500_run000", 
             "SR2_N1500_run001", 
             "SR2_N1500_run002", 
             "SR2_NPR_run000",
             "SR2_NPR_run001",
             "SR2_NPR_run002"] 

gro_location = "../../local-run-output/"
output_dir = "trajectory-data/"


# --------------------------------------------------------------------------- #

class BinConformationsError(Exception):
    """
    General error class for this module.
    """

    pass


def smoothSeries(series,smooth_window=5):
    """
    Smooth a series with a window of fixed size.
    """
    
    # Error check
    if smooth_window <= 0:
        err = "Smoothing window should be an odd integer >= 1!\n"
        raise BinConformationsError(err)

    # Force the smooth window to be an integer
    smooth_window = int(smooth_window)

    # make the smoothing window odd if the user specifies an even number
    if smooth_window % 2 == 0:
        smooth_window += 1

    new_series = []
    half_window = int(smooth_window)/2

    # The first half window is not smoothed
    new_series.extend(series[:half_window])

    # The middle range is smoothed over smooth_window
    indexes = range(half_window,len(series[half_window:-half_window])+half_window)
    for i in indexes:
        new_s = 0
        for k in range(-half_window,half_window+1):
            new_s += series[i+k]

        new_s = float(new_s)/smooth_window
        new_series.append(new_s)

    # The last half-window is not smoothed
    for i in range(len(series)-half_window,len(series)):
        new_series.append(series[i])

    for i, v in enumerate(new_series):
        new_series[i] = int(round(new_series[i],0))

    # A useful test hack for making sure that the smoothing is behaving as expected
    #for i in range(len(series)):
    #    print i, (len(series[0])*"%5i") % tuple(series[i]),
    #    print (len(series[0])*"%5.2f") % tuple(new_series[i])
    #sys.exit()

    return new_series

def parseGroFile(gro_file):
    """
    Parse a gro file, generating a set of dictionaries mapping index to an 
    named atom/residue pair and atom/residue pair to index.
    """

    f = open(gro_file,'r')
    lines = f.readlines()
    f.close()

    atom_to_index = {}
    index_to_atom = []
    for l in lines[2:]:
        index = int(l[15:20]) - 1
        atom = l[11:16].strip()
        residue = l[0:8].strip()

        index_to_atom.append((atom,residue))
        atom_to_index[(atom,residue)] = index
       
    return atom_to_index, index_to_atom
    

def binConformations(input_root,gro_file,specified_keys=None,
                     build_key_list=False):
    """
    """

    important_atoms = global_important_atoms

    # Read hbond file
    f = open("%s_hbond.txt" % input_root,'r')
    hbond_lines = f.readlines()
    f.close()
    hbond_lines = [l for l in hbond_lines if l.startswith("->")] 
   
    # Read water file 
    f = open("%s_water.txt" % input_root,'r')
    water_lines = f.readlines()
    f.close()
    water_lines = [l for l in water_lines if l.startswith("->")] 

    # Parse the GRO file so we can link atom indexes to atom names
    atom_to_index, index_to_atom = parseGroFile(gro_file)

    # Go through the important atoms, grab their indexes, and discard
    # important atoms that do not exist in this file. 
    to_remove = []
    important_indexes = []
    for a in important_atoms:
        try:
            important_indexes.append(atom_to_index[a])
        except KeyError:
            to_remove.append(a)
    important_atoms = [a for a in important_atoms if a not in to_remove]
  
    # Print out the important atoms (a useful way to make sure the script 
    # is working). 
    print important_atoms
 
    # If we're building a key list, just make a list of keys. 
    if build_key_list:
        key_list = []

    # If we already know the keys, make a dictionary using those keys.  
    # Otherwise, remain agnostic about which keys will be used.
    if specified_keys == None:
        obs = {}
    else:
        obs = dict([(k,[0 for j in range(len(hbond_lines))])
                    for k in specified_keys])

    # Walk through the vmd output file and generate a dictionary containing
    # all hydrogen bonds involving important atoms as a function of frame.
    for line_index, l in enumerate(hbond_lines):
        entry = l.split("|")[1:]
        frame = int(entry[0])

        # Create lists of donor and acceptor indexes from the vmd output 
        # line
        if entry[1].startswith(" {"):
            bonds = entry[1].split("}")
            donors = [int(a) for a in bonds[0].strip(" {").split()]
            acceptors = [int (a) for a in bonds[1].strip(" {").split()]
        else:
            bonds = entry[1].split()
            donors = int(bonds[0])
            acceptors = int(bonds[1])
            
        # Go through all donor/acceptor pairs
        for i in range(len(donors)):

            # Only take hydrogen bonds with interesting atoms
            if donors[i] in important_indexes or \
               acceptors[i] in important_indexes:

                # Unique donor/acceptor name
                d = index_to_atom[donors[i]]
                a = index_to_atom[acceptors[i]]
                
                # Skip spurious carbon hydrogen bonds due to slopply h-bond
                # cutoff
                if d[0][0] == "C" or a[0][0] == "C":
                    continue

                # Collapse indistinguishable atoms (Arg NH1/NH2, etc.)
                if d in atom_aliases.keys():
                    d = atom_aliases[d]
                if a in atom_aliases.keys():
                    a = atom_aliases[a]

                # Collapse solvent based on classification in water-recorder 
                if d[1][-3:] == "SOL" or a[1][-3:] == "SOL":

                    # Parse the water line
                    relevant_waters = water_lines[line_index]
                    water_entry = water_lines[line_index].split("|")[1:]
                    water_frame = int(entry[0])
                    if water_frame != frame:
                        err = "Mismatch between frames in water and hbond!\n"
                        raise BinConformationsError(err)

                    water_entry = water_entry[1:]

                    # Now assign a water type based on whether the donor or 
                    # acceptor was assigned a type by water-recorder.tcl
                    water_type = None
                    for j in range(len(water_entry)):
                        waters = [int(w) for w in water_entry[j].split()]
                        if donors[i] in waters:
                            water_type = water_type_list[j]
                            break

                        if acceptors[i] in waters:
                            water_type = water_type_list[j]
                            break

                    # If no type was assigned, this is a bulk water 
                    if water_type == None:
                        water_type = "bulk"

                    # Record the water type
                    if d[1][-3:] == "SOL":
                        d = ("SOL",water_type)  
                    if a[1][-3:] == "SOL":
                        a = ("SOL",water_type)  

                # If we get here, we've officially got an interesting hydrogen
                # bond!  Record the donor/acceptor pair.
                key = (d,a)

                if build_key_list == True:
                    key_list.append(key)
                else:

                    # Record that there is a hydrogen bond between d and a in
                    # this frame.  If no d->a hydrogen bonds have been 
                    # recorded so far, make a list of 0s that shows all frames
                    # beside this one have no hydrogen bond formed.

                    if key not in obs.keys():
                        
                        if specified_keys != None:
                            err = "%r not in specified keys!\n" % key
                            raise BinConformationsError(err)
                        else:
                            obs[key] = [0 for j in range(len(hbond_lines))]

                    obs[key][frame] += 1

    # If we're just building a key list, we're done at this point.  Return the
    # list of keys.
    if build_key_list:
        return key_list
   
    # Smooth the observations 
    for k in obs.keys():
        obs[k] = smoothSeries(obs[k])

    # Make pretty (R-readable) output for each donor/acceptor pair.
    out = ["%20s%20s" % (" ","frame")]
    for k in obs.keys():
        n = "%s.%s_%s.%s" % tuple([k[0][0],k[0][1],k[1][0],k[1][1]])
        out.append("%20s" % n)
    out.append("\n")

    for i in range(len(hbond_lines)):
        out.append("%20i%20i" % (i,i))
        for k in obs.keys():
            out.append("%20i" % obs[k][i])
        out.append("\n")

    return "".join(out)



def main(argv=None):
    """
    Main function for this module.
    """

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Build list of all unique hydrogen bonds observed in all trajectories
    print "Pass 1: building list of all possible hydrogen bonds"
    keys = []
    for f in file_list:
        print f
        gro = os.path.join(gro_location,("%s.gro" % f))
        keys.extend(binConformations(f,gro,build_key_list=True))
    keys = dict([(k,[]) for k in keys]).keys()
    

    # Count those hydrogen bonds over all trajectories and write to files
    print "Pass 2: recording hydrogen bonds for each trajectory"
    for f in file_list:
        print f
        gro = os.path.join(gro_location,("%s.gro" % f))
        out = binConformations(f,gro,specified_keys=keys)

        g = open(os.path.join(output_dir,"%s_hbond-processed.txt" % f),'w')
        g.write(out)
        g.close()


if __name__ == "__main__":
    main()
