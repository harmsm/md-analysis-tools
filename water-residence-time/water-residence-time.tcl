
# Some useful base selections
set ligand_sel "resname EST E40 N15 N06 NPR"
set distance_cutoff 3.5

set water_sel [atomselect 0 "(solvent within $distance_cutoff of protein) and not $ligand_sel"]

set num_frames [molinfo 0 get numframes]
for { set frame_number 0 } { $frame_number < $num_frames } { incr frame_number } {

    $water_sel frame $frame_number
    $water_sel update

    # Record the indexes of nearby waters
    puts "-> | $frame_number | [$water_sel get index] "

}

exit
