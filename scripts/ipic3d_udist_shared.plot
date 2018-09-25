
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set title "Uniform: 1M Particles" font "Times, 26"
set xlabel "Number of workers" font "Times, 26"
set ylabel "Particles/s" font "Times, 26"
set xrange [.5:33]
#set yrange [0:]
#set logscale xy
set xtics (1, 2, 4, 8, 16, 32)
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_1m.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key bottom right

# margins
set tmargin 2
set rmargin 1.
set lmargin 10.5

plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:2 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:8 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:14 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:20 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"


set title "Uniform: 2M Particles" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_2m.pdf"
plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:3 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:9 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:15 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:21 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 4M Particles" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_4m.pdf"
plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:4 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:10 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:16 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:22 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 8M Particles" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_8m.pdf"
plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:5 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:11 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:17 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:23 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 16M Particles" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_16m.pdf"
plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:6 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:12 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:18 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:24 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 32M Particles" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180923_ipic3d_udist_shared_32m.pdf"
plot "../results/single/180923_ipic3d_udist_shared.dat" using 1:7 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:13 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:19 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:25 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

