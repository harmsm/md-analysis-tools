#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001"
selection="resid 41"

for x in $file_list; do 
    echo ${x}

    vmd -dispdev text -e hbond-counter.tcl ../../local-run-output/${x}.* -args "$selection" > ${x}.txt 
    python hbond-counter.py ${x}.txt > ${x}-processed.txt

done
