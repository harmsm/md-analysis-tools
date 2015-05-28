set xrange [0:360]
set yrange [0:360]
set zrange [0:360]

set xlab "chi1"
set ylab "chi2"
set zlab "chi3"

unset key

splot  'SR2_N0600_run000_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N0600_run001_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N0600_run002_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N0600_run003_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N0600_run004_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N1500_run000_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N1500_run001_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_N1500_run002_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_NPR_run000_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_NPR_run001_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2_NPR_run002_dihed-raw-processed.txt.gnu' using 3:4:5 with points
