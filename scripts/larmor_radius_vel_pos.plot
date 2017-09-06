
set terminal postscript eps color enhanced "Times" 26
set grid xtics ytics
set xlabel "position" font "Times, 26"
set ylabel "velocity" font "Times, 26"
set output "| epstopdf --filter > ../imgs/larmor_radius_pos_vel.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
#set rmargin 1.
#set lmargin 9.5

plot "../results/larmor_radius_simple.txt" using 4:7 with linespoints lw 2.0 lc 3 pt 3 title "Larmor radius"

