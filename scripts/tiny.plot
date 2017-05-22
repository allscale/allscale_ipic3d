
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set xlabel "Number of cores" font "Times, 26"
set ylabel "Time [secs]" font "Times, 26"
set xrange [0:9]
set yrange [0:]
set output "| epstopdf --filter > ../imgs/170329_tiny_strong_scaling.pdf"

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
set lmargin 6.5

plot "../results/170329_tiny_strong_scaling.txt" using 1:($2) with linespoints lw 3.0 ps 2.0 pt 3 title "AllScale PIC3D" , \
    "" using 1:($4) with points lw 3.0 ps 2.0 pt 4 title "Original iPIC3D"

