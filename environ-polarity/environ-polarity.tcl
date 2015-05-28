# Count the number of polar atoms in the vicinity of a user specified selection
# over the course of a trajectory.

set USAGE "environ-polarity.tcl selection-of-interest"

# Parse command line
set sel [lindex $argv]

if { $sel == "" } { 
    puts "Incorrect number of arguments!"
    puts $USAGE
    exit
}

set min_dist 3
set max_dist 6
set self "($sel)"
set prot "((oxygen or nitrogen) and protein)"
set solv "(oxygen and solvent)"

# Go through every frame
set n [molinfo top get numframes]
puts "BEGIN"
for {set i 0} {$i < $n} {incr i} {
  
    # Advance the frame 
    molinfo top set frame $i

    for {set d $min_dist} {$d < [expr $max_dist + 1] } {incr d} {
    
        set polar_prot [atomselect top \
                        "($prot and within $d of ($self)) and not ($self)"]
        set polar_solv [atomselect top \
                        "($solv and within $d of ($self)) and not ($self)"]
     
        set num_polar_prot [llength [$polar_prot get name]] 
        set num_polar_solv [llength [$polar_solv get name]] 
           
        puts "$i\t$d\t$num_polar_prot\t$num_polar_solv"

    }
 
}

# Write to standard out
puts "END"

exit
