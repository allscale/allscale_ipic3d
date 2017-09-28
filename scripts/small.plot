
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set xlabel "Number of workers" font "Times, 26"
set ylabel "Time [secs]" font "Times, 26"
#set xrange [.5:9]
#set yrange [0:]
set logscale xy
set xtics (1, 2, 4, 8, 16, 32)
set output "| epstopdf --filter > ../imgs/170926_small_24x24x24_strong_scaling.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
set lmargin 9.5

plot "../results/170926_small_24x24x24_strong_scaling.txt" using 1:(($3)/($1)) with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "Ideal scaling" , \
    "" using 1:($2) with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "AllScale iPIC3D"
#    "" using 1:($4) with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "Original iPIC3D"

