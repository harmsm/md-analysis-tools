
set frame_cutoff 0

# Align CA atoms over the trajectory
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

unset sel1 sel2 move_sel num_frames


# Now calculate the shift in the positions of each CA atom relative to the
# starting state.

set residue_indicies [[atomselect top "protein and name CA"] get index]
set num_residues [ llength $residue_indicies]
set num_frames [molinfo top get numframes]


for { set i 0 } { $i < $num_residues } { incr i } {

    set ref_x 0
    set ref_y 0
    set ref_z 0

    set k [ lindex $residue_indicies $i ]

    set ref_sel [atomselect top "index $k" frame 0]
    set ref_x [$ref_sel get x]
    set ref_y [$ref_sel get y]
    set ref_z [$ref_sel get z]
    
    set this_sel [atomselect top "index $k" frame 0]

    set x_m 0
    set y_m 0
    set z_m 0
    for { set j 0 } { $j < $num_frames } { incr j } {

        if { $j < $frame_cutoff } {
            continue
        }

        $this_sel frame $j
    
        set x($j) [expr [$this_sel get x] - $ref_x]
        set y($j) [expr [$this_sel get y] - $ref_y]
        set z($j) [expr [$this_sel get z] - $ref_z]

        set x_m [expr $x_m + $x($j)]
        set y_m [expr $y_m + $y($j)]
        set z_m [expr $z_m + $z($j)]

    }

    set x_m [expr $x_m/$num_frames]
    set y_m [expr $y_m/$num_frames]
    set z_m [expr $z_m/$num_frames]
   
    set x_s 0 
    set y_s 0
    set z_s 0

    for { set j 0 } { $j < $num_frames } { incr j } {
        
        if { $j < $frame_cutoff } {
            continue
        }

        $this_sel frame $j

        set x_s [expr $x_s + (pow(($x($j) - $x_m),2))]
        set y_s [expr $y_s + (pow(($y($j) - $y_m),2))]
        set z_s [expr $z_s + (pow(($z($j) - $z_m),2))]

        unset x($j) y($j) z($j)

    }

    set x_s [expr sqrt($x_s/($num_frames-$frame_cutoff))]
    set y_s [expr sqrt($y_s/($num_frames-$frame_cutoff))]
    set z_s [expr sqrt($z_s/($num_frames-$frame_cutoff))]

    puts "-> | $i | $x_m | $y_m | $z_m | $x_s | $y_s | $z_s "

    unset ref_sel this_sel k
    unset ref_x ref_y ref_z
    unset x_s y_s z_s
    unset x_m y_m z_m

}

exit
