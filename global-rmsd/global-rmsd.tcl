#   vmd -e visualize.tcl my_psf my_dcd -args ligand
#

# Load a trajectory and visualize it in a way that I like

set ligand "E40 N15 N06 EST NPR"

# Delete default representation
mol delrep 0 top

# Align protein from all frames to protein from the first frame
set sel1 [atomselect top "protein and name CA" frame 0]
set sel2 [atomselect top "protein and name CA" frame 0]

set lig_sel1 [atomselect top "resname $ligand" frame 0]
set lig_sel2 [atomselect top "resname $ligand" frame 0]

set move_sel [atomselect top "all" frame 0]

set num_frames [molinfo top get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    $sel2 frame $frame_number
    $lig_sel2 frame $frame_number
    $move_sel frame $frame_number

    # Align the CA atoms of the protein to the first frame
    set transformation_matrix [measure fit $sel2 $sel1]
    $move_sel move $transformation_matrix

    # Measure the rmsd of the protein and ligand
    puts "-> | $frame_number | [ measure rmsd $sel2 $sel1 ] | [ measure rmsd $lig_sel2 $lig_sel1 ] "

    
}

exit

