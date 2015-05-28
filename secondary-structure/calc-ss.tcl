
set prot [atomselect 0 "protein and name CA" frame 0]

# Now, go through every frame of the simulation...
set num_frames [molinfo 0 get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    animate goto $frame_number
    mol ssrecalc 0

    # Update frame on all of these selections
    #$prot frame $frame_number; $prot update 

    puts "-> | $frame_number | [ $prot get structure ]"

}

exit

