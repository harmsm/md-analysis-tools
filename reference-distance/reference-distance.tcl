# Calculate the distance between three atoms and spit out as a function of
# trajectory.

# vmd -dispdev text -e calc-distances.tcl -m gro_file xtc_file

# Residues of interest (MUST MATCH THE RESIDUES OF INTEREST IN CALC_DISTANCES.PY)
set residues "27 41 51 60 75 87 110 127 179 195 214 223 234 L1 L2"

# Hideous code that creates a list of distances to calculate
set pair_list [ list ]
for { set i 0 } { $i < [ llength $residues ] } { incr i } {

    set res1 [ lindex $residues $i ]
 
    if {$res1 == "L1"} {
        set atom1 [ [ atomselect top "resname N06 N15 EST E40 and name C3" ] get index ]
    } elseif {$res1 == "L2"} {
        set atom1 [ [ atomselect top "resname N06 N15 EST E40 and name C17" ] get index ]
    } else {
        set atom1 [ [ atomselect top "resid $res1 and name CA" ] get index ]
    }

    for { set j [expr $i + 1] } { $j < [llength $residues] } { incr j } {

        puts "$i $j"

        set res2 [ lindex $residues $j ]
        if {$res2 == "L1"} {
            set atom2 [ [ atomselect top "resname N06 N15 EST E40 and name C3" ] get index ]
        } elseif {$res2 == "L2"} {
            set atom2 [ [ atomselect top "resname N06 N15 EST E40 and name C17" ] get index ]
        } else {
            set atom2 [ [ atomselect top "resid $res2 and name CA" ] get index ]
        }

        lappend pair_list [ list $atom1 $atom2 ]

    }
}

foreach pair $pair_list { 
    puts $pair
}

# Caculate distances for the guys we're interested in and write to standard out
set num_frames [molinfo top get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    puts -nonewline "-> | $frame_number | "
    foreach pair $pair_list {
        puts -nonewline "[measure bond $pair frame $frame_number] |"
    }
    puts ""
    
}

exit

