#!/bin/sh

# Go through every line in $1, split the line on spaces, then run a script
# simultaneously on all file_roots listed on the line.  Wait until those
# jobs are finished, then go to the next line.

fname=$1
script_to_run=$2

if [ ! "$fname" ]; then
    echo "Please specify a run input file!"
    exit
fi

if [! "$script_to_run" ]; then
    echo "Please specify a script to run!"
    exit
fi

exec<$fname

while read line; do

    file_list=$line
    for x in `echo $file_list`; do
        echo ${x}
        bash ${script_to_run} ${x} &
    done
    wait $!

done
