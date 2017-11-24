
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
#set xlabel "Number of workers" font "Times, 26"
#set ylabel "Time [secs]" font "Times, 26"
#set xrange [.5:9]
#set yrange [0:]
#set logscale xy
#set xtics (1, 2, 4, 8)
set output "| epstopdf --filter > ../imgs/171124.particles.test.cell.000.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
#set lmargin 9.5

plot "../results/particles.test.cell.000.txt" using 1:2:3
# with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "Ideal scaling" , \

