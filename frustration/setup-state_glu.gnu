set xrange [0:360]
set yrange [0:360]
set zrange [0:360]

set xlab "chi1"
set ylab "chi2"
set zlab "chi3"

unset key

splot  'SR2-Q41E-M75L_N0600_run000_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2-Q41E-M75L_N0600_run001_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2-Q41E-M75L_N0600_run002_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2-Q41E-M75L_N1500_run000_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2-Q41E-M75L_N1500_run001_dihed-raw-processed.txt.gnu' using 3:4:5 with points
replot 'SR2-Q41E-M75L_N1500_run002_dihed-raw-processed.txt.gnu' using 3:4:5 with points
