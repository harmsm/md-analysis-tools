
set sel1 [atomselect 0 "resid 206 and sidechain"]
set sel2 [atomselect 0 "resname N06 NPR EST E40 N15"]
    
# Now, go through every frame of the simulation...
set num_frames [molinfo 0 get numframes] 
for {set frame_number 0} {$frame_number < $num_frames} {incr frame_number} {

    $sel1 frame $frame_number; $sel1 update
    $sel2 frame $frame_number; $sel2 update

    set sasa_sel1 [measure sasa 1.4 $sel1]
    set sasa_sel2 [measure sasa 1.4 $sel2]
    set cont15 [measure contacts 1.5 $sel1 $sel2] 
    set cont25 [measure contacts 2.5 $sel1 $sel2] 
    set cont35 [measure contacts 3.5 $sel1 $sel2] 
    set cont45 [measure contacts 4.5 $sel1 $sel2] 
    set cont55 [measure contacts 5.5 $sel1 $sel2] 
    set cont65 [measure contacts 6.5 $sel1 $sel2] 

    puts "-> | $frame_number | $sasa_sel1 | $sasa_sel2 | $cont15 | $cont25 | $cont35 | $cont45 | $cont55 | $cont65"
 
    unset sasa_sel1 sasa_sel2 cont15 cont25 cont35 cont45 cont55 cont65
}

exit

