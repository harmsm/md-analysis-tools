# usage:
# vmd -e hbond-recorder.tcl -dispdev text -f GRO XTC -args cutoff angle 
#
# Calculate whether a whole slough of hydrogen bonds are formed or not.  This
# can then be used to bin frames into specific conformational types.
#

# Some useful base selections
set ante_water "(name OW) and (within 5 of (resid 249 and name C3)) and (within 10 of ((resid 41 75 72) and name CA))"
set lower_water "(name OW) and (within 5 of (resid 94 and name CA)) and (within 9 of (resid 92 and name CA))"
set upper_water "(name OW) and (within 8 of (resid 41 and name CA)) and (within 8 of (resid 45 and name CA)) and (within 7 of (resid 74 and name CA))"
set inter_water "(name OW) and (within 8 of (resid 45 and name CA)) and (within 8 of (resid 67 and name CA)) and (within 7 of (resid 71 and name CA))"

set ante_water_sel  [atomselect 0 $ante_water]
set lower_water_sel [atomselect 0 $lower_water]
set upper_water_sel [atomselect 0 $upper_water]
set inter_water_sel [atomselect 0 $inter_water]

set num_frames [molinfo 0 get numframes]
for { set frame_number 0 } { $frame_number < $num_frames } { incr frame_number } {

    # Update selections
    $ante_water_sel frame $frame_number;  $ante_water_sel update
    $lower_water_sel frame $frame_number; $lower_water_sel update
    $upper_water_sel frame $frame_number; $upper_water_sel update
    $inter_water_sel frame $frame_number; $inter_water_sel update

    # Record
    puts -nonewline "-> | $frame_number "
    puts -nonewline "| [ $ante_water_sel get index] "
    puts -nonewline "| [$lower_water_sel get index] "
    puts -nonewline "| [$upper_water_sel get index] "
    puts            "| [$inter_water_sel get index] "

}

