#   vmd -e visualize.tcl my_psf my_dcd 
#

# Load a trajectory and visualize it in a way that I like

set ligand "CA ZN"

# Delete default representation
mol delrep 0 top

# Align protein from all frames to protein from the first frame
set sel1 [atomselect top "protein and name CA" frame 0]
set sel2 [atomselect top "protein and name CA" frame 0]
set move_sel [atomselect top "all" frame 0]

set num_frames [molinfo top get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    $sel2 frame $frame_number
    $move_sel frame $frame_number

    set transformation_matrix [measure fit $sel2 $sel1]

    $move_sel move $transformation_matrix
    
}

# Create new representations
mol representation Lines
mol selection "protein"
mol addrep top

mol representation Tube
mol selection "protein"
mol addrep top

mol representation VDW
mol selection "resname ${ligand}"
mol addrep top

#mol selection "resid 17" 
#mol addrep top

#mol selection "solvent within 4 protein"
#mol addrep top
#mol selupdate 3 top on

#mol representation HBonds
#mol selection "(resname ${ligand})"
#mol addrep top
#mol selupdate 4 top on

mol smoothrep top 0 5
mol smoothrep top 1 5
#mol smoothrep top 2 5
#mol smoothrep top 3 5
#mol smoothrep top 4 5

