# Find the number of atoms within 2-3 A of the ligand. 

set ligand_name "N15 N06 E40 NPR EST"

set vdw_200 [atomselect 0 "protein within 2.00 of resname $ligand_name"]
set vdw_225 [atomselect 0 "protein within 2.25 of resname $ligand_name"]
set vdw_250 [atomselect 0 "protein within 2.50 of resname $ligand_name"]
set vdw_275 [atomselect 0 "protein within 2.75 of resname $ligand_name"]
set vdw_300 [atomselect 0 "protein within 3.00 of resname $ligand_name"]

# Now, go through every frame of the simulation...
set num_frames [molinfo 0 get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    # Update frame on all of these selections
    $vdw_200 frame $frame_number; $vdw_200 update
    $vdw_225 frame $frame_number; $vdw_225 update
    $vdw_250 frame $frame_number; $vdw_250 update
    $vdw_275 frame $frame_number; $vdw_275 update
    $vdw_300 frame $frame_number; $vdw_300 update

    puts "-> | $frame_number | [ $vdw_200 num ] | [ $vdw_225 num ] | [ $vdw_250 num ] | [ $vdw_275 num ] | [ $vdw_300 num ] "

    
}

exit

