# vmd -e ligand-orientation.tcl -dispdev text -m reference_structure -f GRO XTC -args ligand_name
#
# Measure changes in ligand orientation over a trajectory relative to a
# a reference structure. 
#
# Steps: 
# 1. Align the reference ligand to it's principal axes, bringing the reference
#    protein along for the ride.
# 2. Go through the trajectory, aligning the protein at each frame to the 
#    principal-axis-aligned reference protein.
# 3. Measure the transformation matrix between the frame ligand and the 
#    principal-axis-aligned reference ligand.
# 4. Write out transformation matrix
#
# Currently set up for the reference protein to be 1ere chain A and the frame
# protein to be SR2 (or non-indel mutants thereof).  To modify change the
# ref_prot_sel and frame_prot_sel so that the residues match between the 
# reference and target proteins.
#
# This script also assumes the reference ligand is EST.  To change this, change
# ref_ligand_name.  If you want to make the ligand align differently (e.g. to 
# not select B, C, D rings), change ring_sel.

package require Orient
namespace import Orient::orient

set ligand_name "EST N15 N06 E40" 
#[lindex $argv 0]

# Both the reference and target ligands must have the same ring selection to 
# calculate transformation matrix correctly
set ref_ligand_name "EST"

# B, C, D rings
set ring_sel "name C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18"

# Reference vectors

# X (V0->V1)
set V0_sel "name C13"
set V1_sel "name C9"

# Y (V2->V3)
set V2_sel "name C12"
set V3_sel "name C14"

# Z (V4->V5)
set V4_sel "name C13"
set V5_sel "name C18"

# These selections must match one another residue by residue to perform correct
# alignment of frame to reference protein
set ref_prot_sel "(protein and name CA) and ((resid 312 to 416) or (resid 418 to 455) or (resid 460 to 546))"
set frame_prot_sel "(protein and name CA) and ((resid 0 to 23) or (resid 30 to 87) or (resid 89 to 92) or (resid 94 to 214) or (resid 218 to 231))"

# -----------------------------------------------------------------------------
# No user modification should be required past this point
# -----------------------------------------------------------------------------

proc measureOrientation { frame_number ligand V0 V1 V2 V3 V4 V5 } {


    # Measure the principal axes of the frame ligand
    set F [draw principalaxes $ligand]
    set center_mass [measure center $ligand weight mass]
    set V0_coord [ $V0 get { x y z } ]   
    set V1_coord [ $V1 get { x y z } ]   
    set V2_coord [ $V2 get { x y z } ]   
    set V3_coord [ $V3 get { x y z } ]   
    set V4_coord [ $V4 get { x y z } ]   
    set V5_coord [ $V5 get { x y z } ]   
 
    set ref_x [list]
    set ref_y [list]
    set ref_z [list]
    lappend ref_x 0 0 0 
    lappend ref_y 0 0 0 
    lappend ref_z 0 0 0 

    for { set i 0 } { $i < 3 } { incr i } {
        lset ref_x $i [ expr [ lindex $V1_coord 0 $i ] - [ lindex $V0_coord 0 $i ] ]
        lset ref_y $i [ expr [ lindex $V3_coord 0 $i ] - [ lindex $V2_coord 0 $i ] ]
        lset ref_z $i [ expr [ lindex $V5_coord 0 $i ] - [ lindex $V4_coord 0 $i ] ]
    } 

    puts "-> | $frame_number | $center_mass | [lindex $F 0] | [lindex $F 1] | [lindex $F 2] | $ref_x | $ref_y | $ref_z"

    unset F
    unset center_mass 
    unset V0_coord
    unset V1_coord
    unset V2_coord
    unset V3_coord
    unset V4_coord
    unset V5_coord
    unset ref_x
    unset ref_y
    unset ref_z

}


# Create selections for reference and frame ligands
set ref_ligand_sel "resname $ref_ligand_name and $ring_sel"
set frame_ligand_sel "resname $ligand_name and $ring_sel"

puts "*************************************************************************"
puts "Reference ligand: $ref_ligand_sel"
puts "Frame ligand:     $frame_ligand_sel"
puts "*************************************************************************"

# Align the reference ligand protein complex to the principal axes of the
# hormone
set move_sel [atomselect 0 "all"]

# Create selections for reference protein and reference ligand
set ref_prot [atomselect 0 $ref_prot_sel]
set ref_ligand [atomselect 0 $ref_ligand_sel]

# Calculate the center of mass for the reference ligand
set ref_ligand_center [measure center $ref_ligand weight mass]
puts "Reference ligand center of mass: $ref_ligand_center"

# There *must* be a better way to do this.  I hate tcl.  
set move_by [list]
lappend move_by [expr -1*[lindex $ref_ligand_center 0]]
lappend move_by [expr -1*[lindex $ref_ligand_center 1]]
lappend move_by [expr -1*[lindex $ref_ligand_center 2]]
$move_sel moveby $move_by

# Check to make sure it really moved
set ref_ligand_center [measure center $ref_ligand weight mass]
puts "New reference ligand center of mass: $ref_ligand_center"

# Store translation here
set translation [measure center $ref_ligand weight mass]

# Now calculate principal axes and align along x, y, and z
set I [draw principalaxes $ref_ligand]
set A [orient $ref_ligand [lindex $I 2] {0 0 1}]
$move_sel move $A
set I [draw principalaxes $ref_ligand]
set A [orient $ref_ligand [lindex $I 1] {0 1 0}]
$move_sel move $A
set I [draw principalaxes $ref_ligand]

set ref_V0 [atomselect 1 "$ref_ligand_sel and $V0_sel"]
set ref_V1 [atomselect 1 "$ref_ligand_sel and $V1_sel"]
set ref_V2 [atomselect 1 "$ref_ligand_sel and $V2_sel"]
set ref_V3 [atomselect 1 "$ref_ligand_sel and $V3_sel"]
set ref_V4 [atomselect 1 "$ref_ligand_sel and $V4_sel"]
set ref_V5 [atomselect 1 "$ref_ligand_sel and $V5_sel"]

measureOrientation -1 $ref_ligand $ref_V0 $ref_V1 $ref_V2 $ref_V3 $ref_V4 $ref_V5

# Create selections for each frame
set frame_prot [atomselect 1 $frame_prot_sel frame 0]
set frame_ligand [atomselect 1 $frame_ligand_sel frame 0]
set move_sel [atomselect 1 "all" frame 0]
set frame_V0 [atomselect 1 "$frame_ligand_sel and $V0_sel"]
set frame_V1 [atomselect 1 "$frame_ligand_sel and $V1_sel"]
set frame_V2 [atomselect 1 "$frame_ligand_sel and $V2_sel"]
set frame_V3 [atomselect 1 "$frame_ligand_sel and $V3_sel"]
set frame_V4 [atomselect 1 "$frame_ligand_sel and $V4_sel"]
set frame_V5 [atomselect 1 "$frame_ligand_sel and $V5_sel"]

# Now, go through every frame of the simulation...
set num_frames [molinfo 1 get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    # Update frame on all of these selections
    $frame_prot frame $frame_number
    $frame_ligand frame $frame_number
    $move_sel frame $frame_number
    $frame_V0 frame $frame_number
    $frame_V1 frame $frame_number
    $frame_V2 frame $frame_number
    $frame_V3 frame $frame_number
    $frame_V4 frame $frame_number
    $frame_V5 frame $frame_number
    
    # Align the frame protein to the reference protein
    set alignment [measure fit $frame_prot $ref_prot]
    $move_sel move $alignment

    measureOrientation $frame_number $frame_ligand $frame_V0 $frame_V1 $frame_V2 $frame_V3 $frame_V4 $frame_V5

}

exit

