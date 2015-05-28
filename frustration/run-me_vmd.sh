#!/bin/bash

x=${1}

vmd -dispdev text -e master-script.tcl -f ../../local-run-output/${x}.* > ${x}.tmp

