
set terminal postscript eps color enhanced "Times" 26
set grid noxtics ytics
set title "Uniform: 8M Particles / Node" font "Times, 26"
set xlabel "Number of nodes" font "Times, 26"
set ylabel "Particles/s" font "Times, 26"
#set xrange [.5:33]
#set yrange [0:]
#set logscale xy
set xtics (1, 2, 4, 8, 16, 32, 64, 128, 256)
set output "| epstopdf --filter > ../imgs/180924_ipic3d_udist_dist_8m_reg.pdf"

set datafile missing '#'

# legend
set key width 0 samplen 1.8
set key bottom right

# margins
set tmargin 2
set rmargin 1.6
set lmargin 8.5

plot "../results/multiple/180924_ipic3d_udist.dat" using 1:2 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:6 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:10 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:14 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 16M Particles / Node" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180924_ipic3d_udist_dist_16m_reg.pdf"
plot "../results/multiple/180924_ipic3d_udist.dat" using 1:3 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:7 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:11 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:15 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 32M Particles / Node" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180924_ipic3d_udist_dist_32m_reg.pdf"
plot "../results/multiple/180924_ipic3d_udist.dat" using 1:4 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:8 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:12 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:16 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"

set title "Uniform: 48M Particles / Node" font "Times, 26"
set output "| epstopdf --filter > ../imgs/180924_ipic3d_udist_dist_48m_reg.pdf"
plot "../results/multiple/180924_ipic3d_udist.dat" using 1:5 with linespoints lw 3.0 lc 3 ps 2.0 pt 3 title "MON=0; RES=0" , \
     "" using 1:9 with linespoints lw 3.0 lc 6 ps 2.0 pt 6 title "MON=1; RES=0" , \
     "" using 1:13 with linespoints lw 3.0 lc 4 ps 2.0 pt 4 title "MON=0; RES=1" , \
     "" using 1:17 with linespoints lw 3.0 lc 7 ps 2.0 pt 7 title "MON=1; RES=1"
