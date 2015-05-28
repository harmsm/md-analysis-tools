#!/bin/bash

# Convert all processed dihedral files into gnu-plot readable files.  

for x in *dihed-raw-processed.txt; 
    do awk '{ if ($1 >= 800 && $1 != "frame") { print $0 }} ' ${x} > ${x}.gnu; 
done
