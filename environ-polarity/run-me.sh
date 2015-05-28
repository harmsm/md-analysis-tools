
file_list="SR2_N0600_run000 SR2_N0600_run001"
selection="resid 41"

for x in $file_list; do 
    echo $x

    vmd -dispdev text -e environ-polarity.tcl ../../local-run-output/${x}.* -args "$selection" > ${x}.txt

    python environ-polarity.py ${x}.txt > ${x}-processed.txt

done
