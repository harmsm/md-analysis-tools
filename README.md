# md-analysis-tools
scripts for analyzing molecular dynamics trajectories using VMD

## Description
Scripts will run the analysis given a .gro and .xtc/.trr file.  (It's expecting GROMACS output, but could be easily tweaked for other sorts of MD package otuput).  

Most analysis directories (say `local-rmsd`) have an associated tcl file (say `local-rmsd.tcl`)= that is actually passed through vmd, which will spew out a vast quantity of incomprehensible input. `local-rmsd.py` then processes this input and spits it out in a human (and "R") readable format.  The analysis can be run on multiple files using the `run-me.sh` script that is in each directory. Analysis directories sometimes have an associated `.R` file for plotting and/or doing further analysis. 
