#!/bin/bash

file_list="SR2_N0600_run000 SR2_N0600_run001 SR2_N0600_run002 SR2_N1500_run000 SR2_N1500_run001 SR2_N1500_run002 SR2-Q41E-M75L_N0600_run000 SR2-Q41E-M75L_N0600_run001 SR2-Q41E-M75L_N0600_run002 SR2-Q41E-M75L_N1500_run000 SR2-Q41E-M75L_N1500_run001 SR2-Q41E-M75L_N1500_run002"
file_list="SR2_N0600_run002 SR2-Q41E-M75L_N0600_run002"

for x in `echo $file_list`; do 

    echo ${x}
   
    vmd -dispdev text -e water-residence-time.tcl -f ../../local-run-output/${x}.* > ${x}.txt 

    python ./water-residence-time.py ${x}.txt ../../local-run-output/${x}.gro > ${x}-processed.txt

done

