
set aromat1_sel [atomselect 0 "resid 206 and name CG CD1 CE1 CZ CD2 CE2"]
set aromat2_sel [atomselect 0 "resid 209 and name CG CD1 CE1 CZ CD2 CE2"]
set sulfur_sel [atomselect 0 "resid 230 and name SD"]

# Now, go through every frame of the simulation...
set num_frames [molinfo 0 get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    # Update frame on all of these selections
    $aromat1_sel frame $frame_number; $aromat1_sel update
    $aromat2_sel frame $frame_number; $aromat2_sel update
    $sulfur_sel frame $frame_number; $sulfur_sel update

    # Grab position of sulfur and aromatic centroid
    set aromat1_pos [measure center $aromat1_sel]
    set aromat2_pos [measure center $aromat2_sel]
    set sulfur_pos [$sulfur_sel get "x y z"] 

    # Calculate distance between them
    set dx2 [expr ([lindex $aromat1_pos 0] - [lindex [lindex $sulfur_pos 0] 0])**2]
    set dy2 [expr ([lindex $aromat1_pos 1] - [lindex [lindex $sulfur_pos 0] 1])**2]
    set dz2 [expr ([lindex $aromat1_pos 2] - [lindex [lindex $sulfur_pos 0] 2])**2]
    set d1  [expr sqrt($dx2 + $dy2 + $dz2) ]
    
    set dx2 [expr ([lindex $aromat2_pos 0] - [lindex [lindex $sulfur_pos 0] 0])**2]
    set dy2 [expr ([lindex $aromat2_pos 1] - [lindex [lindex $sulfur_pos 0] 1])**2]
    set dz2 [expr ([lindex $aromat2_pos 2] - [lindex [lindex $sulfur_pos 0] 2])**2]
    set d2  [expr sqrt($dx2 + $dy2 + $dz2) ]

    # Write out
    puts "-> | $frame_number | $d1 | $d2"
    
}

exit

