#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001"
for x in `echo $file_list`; do 
    
    echo ${x}

    vmd -dispdev text -e s-aromatic-pi.tcl -f ../../local-run-output/${x}.* > ${x}.txt
    ./s-aromatic-pi.py ${x}.txt > ${x}-processed.txt

done
