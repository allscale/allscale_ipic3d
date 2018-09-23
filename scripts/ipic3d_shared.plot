
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set xlabel "Number of workers" font "Times, 26"
set ylabel "Particles/s" font "Times, 26"
set xrange [.5:33]
#set yrange [0:]
#set logscale xy
set xtics (1, 2, 4, 8, 16, 32)
set output "| epstopdf --filter > ../imgs/180923_ipic3d_shared.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
set lmargin 9.5

plot "ipic3d_shared.dat" using 1:2 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "ipic3d_shared.dat" using 1:3 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "ipic3d_shared.dat" using 1:4 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "ipic3d_shared.dat" using 1:5 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

