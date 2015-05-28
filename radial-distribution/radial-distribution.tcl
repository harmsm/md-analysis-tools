# calc_gr.tcl
# Use VMD to determine the pairwise radial distribution function between two
# selections for an entire MD trajectory.
#   
#   vmd -dispdev text -e calc_gr.tcl my_psf my_dcd -args "sel1" "sel2"
# 
# Note: the quotes are required on the selections to make the argument parsing
# work.  sel1 and sel2 are vmd style selections like "resid 41"
#

set USAGE "vmd -dispdev text -e radial-distribution.tcl gro xtc -args sel1 sel2"

# Parse command line
set obj1 [lindex $argv 0]
set obj2 [lindex $argv 1]

if { $obj1 == "" || $obj2 == "" } { 
    puts "Incorrect number of arguments!"
    puts $USAGE
    exit
}

# A procedure that sets the unit cell to the same value over the entire
# trajectory.
proc set_unitcell {a b c {molid top} {alpha 90.0} {beta 90.0} {gamma 90.0}} {

    # set the unitcell parameters for the whole trajectory
    # <akohlmey 02.07.2003 10:38:17 yello.theochem.ruhr-uni-bochum.de>
    # Copyright (c) 2003 by <Axel.Kohlmeyer@theochem.ruhr-uni-bochum.de>

    # arguments: 
    #  molid = molecule id where the unitcell is added to (default top)

    if {![string compare $molid top]} {
        set molid [molinfo top]
    }

    set n   [molinfo $molid get numframes]

    for {set i 0} {$i < $n} {incr i} {
        molinfo $molid set frame $i
        molinfo $molid set {a b c alpha beta gamma} \
            [list $a $b $c $alpha $beta $gamma]
    }
}

# Find dimensions of the entire protein/water/ion box
set all [atomselect top "all"]
set minmax [measure minmax $all]
foreach {min max} $minmax { break }
foreach {xmin ymin zmin} $min { break }
foreach {xmax ymax zmax} $max { break }

# Use these dimensions to determine the unit cell that will be used in the 
# g(r) calculation.
set set_a [expr $xmax - $xmin]
set set_b [expr $ymax - $ymin]
set set_c [expr $zmax - $zmin]
set_unitcell $set_a $set_b $set_c

set obj1_sel [atomselect top $obj1]
set obj2_sel [atomselect top $obj2]
set gr [measure gofr $obj1_sel $obj2_sel first 0 last -1 step 1 usepbc TRUE]

# Write to standard out
puts "BEGIN"
puts $gr
puts "END"

exit
