# usage:
# vmd -e hbond-recorder.tcl -dispdev text -f GRO XTC -args cutoff angle 
#
# Calculate whether a whole slough of hydrogen bonds are formed or not.  This
# can then be used to bin frames into specific conformational types.
#

# Parse command line
#set cutoff [lindex $argv 0 ]
#set angle [lindex $argv 1 ]

set cutoff 3.5
set angle 45

# Some useful base selections
set ligand_sel "resname EST E40 N15 N06 NPR"

# A-ring components

set hbond_sel [atomselect 0 "protein or (solvent within 5.5 of protein)"]

puts "-* Cutoff: $cutoff A, angle: $angle degrees"

set num_frames [molinfo 0 get numframes]
for { set frame_number 0 } { $frame_number < $num_frames } { incr frame_number } {

    $hbond_sel frame $frame_number
    $hbond_sel update

    # Record
    puts "-> | $frame_number | [measure hbonds $cutoff $angle $hbond_sel] "

}

