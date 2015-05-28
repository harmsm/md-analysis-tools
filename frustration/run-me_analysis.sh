#!/bin/bash

x=${1}

#vmd -dispdev text -e master-script.tcl -f ../../local-run-output/${x}.* > ${x}.tmp

rm -f xx*

csplit -k ${x}.tmp '/WATER/'; wait $!
mv xx01 ${x}_water.txt
mv xx00 tmp_input

csplit -k tmp_input '/DIHEDRAL/'; wait $!
mv xx00 ${x}_hbond.txt
mv xx01 ${x}_dihed.txt

rm -f tmp_input 
 
python ./parse-dihedral.py ${x}_dihed.txt > ${x}_dihed-raw-processed.txt
python ./parse-water.py ${x}_water.txt > ${x}_water-processed.txt


