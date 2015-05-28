#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001"
for x in `echo $file_list`; do
    echo $x
 
    vmd -dispdev text -e ligand-orientation.tcl -m 1ere_A_1-600.pdb -f ../../local-run-output/${x}.* > ${x}.txt
    python ligand-orientation.py ${x}.txt > ${x}-processed.txt

done

