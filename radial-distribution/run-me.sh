#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001"
sel1="resid 41"
sel2="solvent"

for x in `echo $file_list`; do 

    echo ${x}
   
    vmd -dispdev text -e radial-distribution.tcl -f ../../local-run-output/${x}.* -args "$sel1" "$sel2" > ${x}.txt

    python radial-distribution.py ${x}.txt > ${x}-processed.txt

done

