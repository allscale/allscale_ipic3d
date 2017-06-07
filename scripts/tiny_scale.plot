
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set xlabel "Number of cores" font "Times, 26"
set ylabel "Time [secs]" font "Times, 26"
set xrange [0:]
set yrange [0:]
set output "| epstopdf --filter > ../imgs/exlu_gops_var5_piv_k420_170105.pdf"

# legend
set key width 0 samplen 1.8
set key top right

# margins
set tmargin .5
set rmargin 1.
set lmargin 6.5

plot "../results/170329_tiny_strong_scaling.txt" using 1:($3/$2) with linespoints lw 3.0 title "AllScale PIC3D"

