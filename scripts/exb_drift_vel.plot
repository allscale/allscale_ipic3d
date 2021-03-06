
set terminal postscript eps color enhanced "Times" 26
set grid xtics ytics
set xlabel "time" font "Times, 26"
set ylabel "velocity" font "Times, 26"
#set xrange [:.67]
set output "| epstopdf --filter > ../imgs/exb_drift_velx.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
set lmargin 9.5

plot "../results/exb_drift.txt" using 2:6 with lines lw 3.0 lc 3 title "ExB Drift"

