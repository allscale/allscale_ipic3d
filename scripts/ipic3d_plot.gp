set title "Simple demo of scatter data conversion to grid data"
#set terminal windows size 1280,1024 title "AllScale ipic3d initial particle distribution"
set terminal X11 size 1280,1024 title "AllScale ipic3d initial particle distribution"
unset hidden3d
set ticslevel 0
set view 60,30
# set autoscale
#set parametric
set style data points
set xlabel "data style point - no dgrid"
set xrange [0:160]
set yrange [0:160]
set zrange [0:160]
set key box
splot "t_begin.txt" pt 7 ps 1 lc rgb "red", \
      "t_end.txt" pt 7 ps 1 lc rgb "green"
pause -1
