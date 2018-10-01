
set terminal postscript eps color enhanced "Times" 26
set grid xtics ytics
set xlabel "position" font "Times, 26"
set ylabel "velocity" font "Times, 26"
set output "| epstopdf --filter > ../imgs/larmor_radius_pos_vel.pdf"

set datafile missing '#'

set for [i=1:5] linetype i dt i

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
#set rmargin 1.
#set lmargin 9.5

# line style
set style line 1 lt 5 lc rgb "red" lw 2

plot "../results/larmor_radius.txt" using 3:6 with line lw 4.0 lc 3 title "theory" , \
     "<(sed -n '1,121p' ../results/larmor_radius.txt)" using 3:6 ls 1 title "numerical"

