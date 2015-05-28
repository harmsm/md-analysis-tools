#!/bin/bash 

file_list="SR2_N0600_run000 SR2_N0600_run001"

for x in $file_list; do 

    echo $x

    vmd -dispdev text -e local-rmsd.tcl ../../local-run-output/${x}.* > ${x}.txt 
    python local-rmsd.py ${x}.txt > ${x}-processed.txt

done
