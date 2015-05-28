load SR2_chainA.pdb, SR2

hide eve
color gray50
show ribbon
show sticks, resn STR
color red, (name O* and resn STR)
color blue, (name N* and resn STR)

sele h3, resid 27-51 
sele h5, resid 60-87 
sele h7, resid 110-127
sele h10, resid 179-214 
sele h12, resid 223-238

color pink, h3
color white, h5
color cyan, h7
color limegreen, h10
color orange, h12

show spheres, name CA and (resid 27,41,51,60,75,87,110,127,179,195,214,223,234)
color blue, name CA and (resid 27,60,110,179,223)
color red, name CA and (resid 51,87,127,214,234)
color purple, name CA and (resid 41,75,195)
