set ligand_sel "resname EST E40 N15 N06 NPR"

# Calculate some important dihedrals 

set r41_chi1 [[atomselect 0 "resid 41 and name N CA CB CG"] get index]
set r41_chi2 [[atomselect 0 "resid 41 and name CA CB CG CD"] get index]
set r41_chi3 [[atomselect 0 "resid 41 and name CB CG CD NE2 OE2"] get index]

set r75_chi1 [[atomselect 0 "resid 75 and name N CA CB CG"] get index]
set r75_chi2 [[atomselect 0 "resid 75 and name CA CB CG SD CD1"] get index]
set r75_chi3 [[atomselect 0 "resid 75 and name CB CG SD CE CD1 CD2"] get index]

set r82_chi1 [[atomselect 0 "resid 82 and name N CA CB CG"] get index]
set r82_chi2 [[atomselect 0 "resid 82 and name CA CB CG CD"] get index]
set r82_chi3 [[atomselect 0 "resid 82 and name CB CG CD NE"] get index]

set lig_chi1 [[atomselect 0 "$ligand_sel and name C4 C3 O3 H3"] get index]
if { [ llength $lig_chi1 ] == 3 } { 
    set lig_chi1 [[atomselect 0 "$ligand_sel and name C5 C4 C3 O3"] get index]
} 

# Write out header
puts -nonewline "-* "
puts -nonewline "| r41_chi1 | r41_chi2 | r41_chi3 "
puts -nonewline "| r75_chi1 | r75_chi2 | r75_chi3 "
puts -nonewline "| r82_chi1 | r82_chi2 | r82_chi3 "
puts            "| lig_chi1"


# Write out dihedrals for each frame
set num_frames [molinfo 0 get numframes]
for { set frame_number 0 } { $frame_number < $num_frames } { incr frame_number } {

    puts -nonewline "-> | $frame_number "
    puts -nonewline "| [measure dihed $r41_chi1 frame $frame_number] "
    puts -nonewline "| [measure dihed $r41_chi2 frame $frame_number] "
    puts -nonewline "| [measure dihed $r41_chi3 frame $frame_number] "
    puts -nonewline "| [measure dihed $r75_chi1 frame $frame_number] "
    puts -nonewline "| [measure dihed $r75_chi2 frame $frame_number] "
    puts -nonewline "| [measure dihed $r75_chi3 frame $frame_number] "
    puts -nonewline "| [measure dihed $r82_chi1 frame $frame_number] "
    puts -nonewline "| [measure dihed $r82_chi2 frame $frame_number] "
    puts -nonewline "| [measure dihed $r82_chi3 frame $frame_number] "
    puts            "| [measure dihed $lig_chi1 frame $frame_number] "

}

