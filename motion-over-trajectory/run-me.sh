#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001"
frame_min=0
frame_max=40
frame_step=10

for x in `echo $file_list`; do
   
    # Make a temporary directory 
    rm -rf tmp 
    mkdir tmp
    cd tmp

    # Extract required frames
    vmd -dispdev text -e ../../../tools/extract_frames.tcl -f ../../../local-run-output/${x}.* -args output $frame_min $frame_max $frame_step

    # Combine all of the output pdb files into a single pdb file
    rm -f ../${x}.pdb
    for y in *.pdb; do
        grep -v CRYST ${y} >> ../${x}.pdb
    done

    # Delete temporary pdb
    cd ..
    rm -rf tmp

done
